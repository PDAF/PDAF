!$Id: compute_rms_smoother.F90 1711 2016-12-20 10:17:58Z lnerger $
!BOP
!
! !ROUTINE: compute_rms_smoother --- Compute RMS errors for smoothed states
!
! !INTERFACE:
SUBROUTINE compute_rms_smoother(step, dim_lag, dim_p, dim_ens, state_p)


! !DESCRIPTION:
! Helper routine for for pre/poststep routine.
! Used in the filters: ETKF/LETKF/ESTKF\LESTKF
! 
! The routine is called by prepoststep_ens_pdaf and prepoststep_etkf_pdaf.
! This is a very simple implementation that just computes the RMS errors
! for all smoothed states and displays the RMS error. (In case of a linear
! model without model noise, the RMS errors of all past smoothed states
! are identical to the true and estimated error of the filtered state.) 
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2012-05 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_memcount, &
       ONLY: memcount
  USE mod_model, &
       ONLY: dim_state, local_dims, dt
  USE mod_parallel, &
       ONLY: mype_filter, npes_filter, COMM_filter, MPI_REAL, &
       MPIerr, MPIstatus
  USE mod_assimilation, &
       ONLY: delt_obs
  USE PDAF_interfaces_module, &    ! Required to get pointer to smoother state array
       ONLY: PDAF_get_smootherens

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
  INTEGER, INTENT(in) :: dim_lag     ! Size of lag for smoother
  INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
  REAL, INTENT(inout) :: state_p(dim_p) ! Array for state vector

! !CALLING SEQUENCE:
! Called by: prepoststep_ens_pdaf  
! Called by: prepoststep_etkf_pdaf  

! *** local variables ***
  INTEGER :: i, j, lag         ! counters
  REAL :: invdim_ens                   ! Inverse ensemble size
  REAL :: invdim_ensm1                 ! Inverse of ensemble size minus 1
  INTEGER, SAVE :: allocflag = 0       ! Flag for memory counting

  REAL, POINTER :: sens_pointer(:,:,:) ! Pointer to smoother ensemble
  INTEGER :: lastlag                   ! Number initialized lags in sens_pointer
  INTEGER :: status                    ! Status flag for PDAF_get_smootherens
  REAL :: rmserror_est                 ! estimated RMS error
  REAL :: rmserror_true                ! true RMS error
  REAL :: rmserror_rel                 ! relative error in estimated error
  REAL, ALLOCATABLE :: variance(:)     ! model state variances
  REAL, ALLOCATABLE :: truevariance(:) ! model state variances
  REAL, ALLOCATABLE :: truefield_p(:)  ! true local model state

  ! Variables for parallelization - local fields
  INTEGER :: offset   ! Row-offset according to domain decomposition
  REAL, ALLOCATABLE :: variance_p(:)     ! local variance
  REAL, ALLOCATABLE :: truevariance_p(:) ! local model state variances


! **********************
! *** INITIALIZATION ***
! **********************

  ! Get access to smoother array
  IF (step > 0 .AND. dim_lag > 0) THEN
     CALL PDAF_get_smootherens(sens_pointer, lastlag, status)
  END IF

  ! Initialize numbers
  invdim_ens    = 1.0 / REAL(dim_ens)  
  invdim_ensm1  = 1.0 / REAL(dim_ens - 1)


! **********************************************
! *** Compute RMS errors for smoothed states ***
! **********************************************

  ALLOCATE(variance_p(dim_p))
  ALLOCATE(variance(dim_state))
  ALLOCATE(truefield_p(dim_p))
  ALLOCATE(truevariance_p(dim_p))
  ALLOCATE(truevariance(dim_state))
  ! count allocated memory
  CALL memcount(3, 'r', 2*dim_state + 3*dim_p)

  ! Initialize screen output
  IF (mype_filter == 0 .AND. lastlag > 0) THEN
     WRITE (*, '(8x, a)') 'Analyze smoothed states'
     WRITE (*, '(12x, a20/, 15x, a3, 5x, a7, 10x, a4, 5x, a24)') &
          'Smoother RMS errors:','lag', 'sampled','true','relative underestimation'
  END IF

  lagloop: DO lag = 1, lastlag
     ! Compute RMS errors for each lag

     rmserror_est  = 0.0
     rmserror_true = 0.0

     ! Compute local mean state 
     state_p = 0.0
     DO i = 1, dim_ens
        DO j = 1, dim_p
           state_p(j) = state_p(j) + sens_pointer(j, i, lag)
        END DO
     END DO
     state_p(:) = invdim_ens * state_p(:)


     ! ************************************
     ! *** Part 1: Estimated RMS errors ***
     ! ************************************

     ! *** Compute local sampled variances ***
     variance_p(:) = 0.0
     DO i = 1, dim_ens
        DO j = 1, dim_p
           variance_p(j) = variance_p(j) &
             + (sens_pointer(j, i, lag) - state_p(j)) &
             * (sens_pointer(j, i, lag) - state_p(j))
        END DO
     END DO
     variance_p(:) = invdim_ensm1 * variance_p(:)

     ! *** Assemble global variance vector on filter PE 0 ***
     PE0_a: IF (mype_filter /= 0) THEN

        ! send sub-fields from PEs /=0
        CALL MPI_send(variance_p(1 : dim_p), dim_p, &
             MPI_REAL,0, mype_filter, COMM_filter, MPIerr)

     ELSE PE0_a
        ! receive and assemble variance field

        ! On PE 0 init variance directly
        variance(1 : dim_p) = variance_p(1 : dim_p)

        ! Receive part of variance field from PEs > 0 into 
        ! correct part of global variance

        offset = 0

        DO i = 2, npes_filter
           ! Increment offset
           offset = offset + local_dims(i - 1)

           ! Receive variance part
           CALL MPI_recv(variance(1 + offset), local_dims(i), &
                MPI_REAL, i - 1, i - 1, COMM_filter, MPIstatus, MPIerr)
        END DO
      
     END IF PE0_a

     ! *** Compute RMS errors according to sampled variances ***
     pe0: IF (mype_filter == 0) THEN
     ! total estimated RMS error
        DO i = 1, dim_state
           rmserror_est = rmserror_est + variance(i)
        ENDDO
        rmserror_est = SQRT(rmserror_est / dim_state)
     END IF pe0


     ! *******************************
     ! *** Part 2: True RMS errors ***
     ! *******************************

     ! Compute true state field
     truefield_p(:) = 1.0
     DO i = 1, ABS(step - lag * delt_obs)
        truefield_p(:) = truefield_p(:) + 1.0 * dt
     END DO

     ! True local variance field
     truevariance_p(:) = 0.0
     DO j = 1, dim_p
        truevariance_p(j) = truevariance_p(j) &
             + (state_p(j) - truefield_p(j)) &
             * (state_p(j)  -truefield_p(j))
     END DO

     ! *** Assemble global variance vector on filter PE 0
     PE0_b: IF (mype_filter /= 0) THEN
      
        ! send sub-fields from PEs /=0
        CALL MPI_send(truevariance_p(1 : dim_p), dim_p, &
             MPI_REAL, 0, mype_filter, COMM_filter, MPIerr)

     ELSE PE0_b
        ! receive and assemble variance field

        ! On PE 0 init variance directly
        truevariance(1 : dim_p) = truevariance_p(1 : dim_p)

        ! Receive part of variance field from PEs > 0 into 
        ! correct part of global variance

        offset = 0

        DO i = 2, npes_filter
           ! Increment offset
           offset = offset + local_dims(i - 1)

           ! Receive variance part
           CALL MPI_recv(truevariance(1 + offset : local_dims(i) + offset), &
                local_dims(i), MPI_REAL, &
                i - 1, i - 1, COMM_filter, MPIstatus, MPIerr)
        END DO
      
     END IF PE0_b

     ! Total true RMS error
     pe0a: IF (mype_filter == 0) THEN
        DO i = 1, dim_state
           rmserror_true = rmserror_true + truevariance(i)
        ENDDO
        rmserror_true = SQRT(rmserror_true / dim_state)
     END IF pe0a

     ! *** Screen IO ***
     IF (mype_filter == 0) THEN
        rmserror_rel = (rmserror_true - rmserror_est) / rmserror_true
       WRITE (*, '(15x, i5, es12.4, 3x, es12.4, 3x, es12.4)') &
            lag, rmserror_est, rmserror_true, rmserror_rel
      IF (lag==lastlag) then
!         WRITE (*, '(a,15x, i5, es12.4, 3x, es12.4, 3x, es12.4)') &
!              'smoother:', lag, rmserror_est, rmserror_true, rmserror_rel
      ENDIF
     END IF

  END DO lagloop


! ********************
! *** finishing up ***
! ********************

  DEALLOCATE(variance, variance_p, truefield_p, truevariance, truevariance_p)

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE compute_rms_smoother
