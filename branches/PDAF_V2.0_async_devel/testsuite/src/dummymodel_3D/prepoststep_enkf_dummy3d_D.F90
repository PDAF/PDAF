!$Id: prepoststep_enkf_dummy3d_D.F90 1549 2015-01-17 09:29:05Z lnerger $
!BOP
!
! !ROUTINE: prepoststep_enkf --- Used-defined Pre/Poststep routine for PDAF
!
! !INTERFACE:
SUBROUTINE prepoststep_enkf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF (EnKF):
! 
! The routine is called for before and after 
! the analysis. Also it is called once at the
! initial time before any forecasts are computed.
! The routine provides full access to the state 
! estimate and the state ensemble to the user.
! Thus, user-controlled pre- and poststep 
! operations can be performed here. For example 
! the forecast and the analysis states and ensemble
! covariance matrix can be analyzed, e.g. by 
! computing the estimated variances. In addition, 
! the estimates can be written to disk. If a user 
! considers to perform adjustments to the 
! estimates (e.g. for balances), this routine is 
! the right place for it.
!
! The routine is called by all filter processes.
!
! For the dummy model with domain-decomposition 
! we compute the estimated and the true estimation 
! errors. These values are then written into a 
! file.
!
! !REVISION HISTORY:
! 2004-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_memcount, &
       ONLY: memcount
  USE mod_parallel, &
       ONLY: mype_filter, npes_filter, COMM_filter, MPI_DOUBLE_PRECISION, &
       MPIerr, MPIstatus
  USE mod_model, &
       ONLY: dims, dim_l, dims_l_all, dt, step_null
  USE mod_assimilation, &
       ONLY: filename

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
     ! (When the routine is called before the analysis -step is provided.)
  INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   ! PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   ! PE-local dimension of observation vector
  REAL, INTENT(inout) :: state_p(dim_p) ! PE-local forecast/analysis mean state
  ! The array 'state' is not generally not initialized in the case of EnKF.
  ! It can be used freely here.
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) ! PE-local state ensemble
  INTEGER, INTENT(in) :: flag        ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_get_state      (as U_prepoststep)
! Called by: PDAF_enkf_update    (as U_prepoststep)
! Calls: memcount
! Calls: MPI_send
! Calls: MPI_recv
!EOP


! *** local variables ***
  INTEGER :: i, j, member           ! Counters
  INTEGER :: dim_var_p              ! Local size of variance field (dims(3)=1 only)
  INTEGER :: offset                 ! Row-offset according to domain decomposition
  REAL :: invdim_ens                ! Inverse of ensemble size
  REAL :: invdim_ensm1              ! Inverse of ensemble size minus 1
  REAL :: rmserror_est              ! Estimated RMS error
  REAL :: rmserror_true             ! True RMS error
  REAL :: rmserror_rel              ! Relative error in estimated error
  INTEGER, SAVE :: allocflag = 0    ! Flag for memory counting
  LOGICAL, SAVE :: firstio = .TRUE. ! File output is peformed for initial time?
  REAL, ALLOCATABLE :: variance(:)       ! Estimated model state variance vector
  REAL, ALLOCATABLE :: truevariance(:)   ! True model state variance vector
  REAL, ALLOCATABLE :: truefield_p(:)    ! PE-local true model state
  REAL, ALLOCATABLE :: variance_p(:)     ! PE-local variance vector
  REAL, ALLOCATABLE :: truevariance_p(:) ! PE-local model state variance vector


! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype_filter == 0) THEN
     IF (step - step_null == 0) THEN
        WRITE (*, '(i7, 3x, a)') step, 'Analyze initial state ensemble'
     ELSE IF (step > 0) THEN
        WRITE (*, '(8x, a)') 'Analyze assimilated state ensemble'
     ELSE IF (step < 0) THEN
        WRITE (*, '(8x, a)') 'Analyze forecasted state ensemble'
     END IF
  END IF

  ! Size of local variable field (only dims(3)=1)
  dim_var_p = dim_l(1) * dim_l(2)

  ! Allocate fields
  ALLOCATE(variance(dims(1) * dims(2)))
  ALLOCATE(variance_p(dim_var_p))
  IF (allocflag == 0) THEN
     ! count allocated memory
     CALL memcount(3, 'r', dims(1) * dims(2) + dim_var_p)
  END IF

  ! Initialize numbers
  rmserror_est  = 0.0
  rmserror_true = 0.0
  rmserror_rel  = 0.0
  invdim_ens    = 1.0 / REAL(dim_ens)
  invdim_ensm1  = 1.0 / REAL(dim_ens - 1)


! **************************
! *** Compute mean state ***
! **************************

  IF (mype_filter == 0) &
       WRITE (*, '(8x, a, i5)') '--- compute ensemble mean'

  ! local 
  state_p = 0.0
  DO member = 1, dim_ens
     DO i = 1, dim_p
        state_p(i) = state_p(i) + ens_p(i, member)
     END DO
  END DO
  state_p(:) = invdim_ens * state_p(:)


! *********************************
! *** Compute sampled variances ***
! *********************************
  variance_p(:) = 0.0
  DO member = 1, dim_ens
     DO j = 1, dim_var_p
        variance_p(j) = variance_p(j) &
             + (ens_p(j, member) - state_p(j)) &
             * (ens_p(j, member) - state_p(j))
     END DO
  END DO
  variance_p(:) = invdim_ensm1 * variance_p(:)

! *** assemble global variance vector on filter PE 0
  PE0_a: IF (mype_filter /= 0) THEN

     ! send sub-fields from PEs /=0
     CALL MPI_send(variance_p(1 : dim_var_p), dim_var_p, &
          MPI_DOUBLE_PRECISION, 0, mype_filter, COMM_filter, MPIerr)

  ELSE PE0_a
     ! receive and assemble variance field

     ! On PE 0 init variance directly
     variance(1 : dim_var_p) = variance_p(1 : dim_var_p)

     ! Receive part of variance field from PEs > 0 into 
     ! correct part of global variance

     offset = 0

     DO i = 2, npes_filter
        ! Increment offset
        offset = offset + dims_l_all(1, i - 1) * dims_l_all(2, i - 1)

        ! Receive variance part
        CALL MPI_recv(variance(1 + offset), dims_l_all(1, i) * dims_l_all(2, i), &
             MPI_DOUBLE_PRECISION, i - 1, i - 1, COMM_filter, MPIstatus, MPIerr)
     END DO
      
  END IF PE0_a

  DEALLOCATE(variance_p)


! ************************************************************
! *** Compute RMS errors according to sampled covar matrix ***
! ************************************************************

  ! *** Only on PE 0 ***
  pe0: IF (mype_filter == 0) THEN
     ! total estimated RMS error
     DO i = 1, dims(1) * dims(2)
        rmserror_est = rmserror_est + variance(i)
     ENDDO
     rmserror_est = SQRT(rmserror_est / (dims(1) * dims(2)))
  END IF pe0

  ! *** Compute true variances
  ALLOCATE(truefield_p(dim_var_p))
  ALLOCATE(truevariance_p(dim_var_p))
  ALLOCATE(truevariance(dims(1) * dims(2)))
  IF (allocflag == 0) THEN
     ! count allocated memory
     CALL memcount(3, 'r', 2 * dim_var_p + dims(1) * dims(2))
     allocflag = 1
  END IF

  truefield_p(:) = 1.0
!   DO i = 1, ABS(step)
!      truefield_p(:) = truefield_p(:) + 1.0 * dt
!   END DO

  truevariance_p(:) = 0.0
  DO j = 1, dim_var_p
     truevariance_p(j) = truevariance_p(j) &
          + (state_p(j) - truefield_p(j)) &
          * (state_p(j) - truefield_p(j))
  END DO

  ! *** assemble global variance vector on filter PE 0
  PE0_b: IF (mype_filter /= 0) THEN

     ! send sub-fields from PEs /=0
     CALL MPI_send(truevariance_p(1 : dim_var_p), dim_var_p, &
          MPI_DOUBLE_PRECISION, 0, mype_filter, COMM_filter, MPIerr)

  ELSE PE0_b
     ! receive and assemble variance field

     ! On PE 0 init variance directly
     truevariance(1 : dim_var_p) = truevariance_p(1 : dim_var_p)

     ! Receive part of variance field from PEs > 0 into 
     ! correct part of global variance

     offset = 0

     DO i = 2, npes_filter
        ! Increment offset
        offset = offset + dims_l_all(1, i - 1) * dims_l_all(2, i - 1)

        ! Receive variance part
        CALL MPI_recv(truevariance(1 + offset &
             : dims_l_all(1, i) * dims_l_all(2, i) + offset), &
             dims_l_all(1, i) * dims_l_all(2, i), MPI_DOUBLE_PRECISION, &
             i - 1, i - 1, COMM_filter, MPIstatus, MPIerr)
     END DO
      
  END IF PE0_b

  pe0a: IF (mype_filter == 0) THEN

     ! total true RMS error
     DO i = 1, dims(1) * dims(2)
        rmserror_true = rmserror_true + truevariance(i)
     ENDDO
     rmserror_true = SQRT(rmserror_true / (dims(1) * dims(2)))

     ! deallocate fields
     DEALLOCATE(truefield_p, truevariance, truevariance_p)

  END IF pe0a


! *****************
! *** Screen IO ***
! *****************

  ! Output RMS errors given by sampled covar matrix
  IF (mype_filter == 0) THEN
     rmserror_rel = (rmserror_true - rmserror_est) / rmserror_true
     WRITE (*, '(12x, a, es12.4)') &
          'RMS error according to sampled variance: ', rmserror_est
     WRITE (*, '(15x, a, es12.4)') &
          'RMS error according to true variance: ', rmserror_true
     WRITE (*, '(14x, a, es12.4)') &
          'Relative underestimation of variances: ', rmserror_rel
  END IF

! *******************
! *** File output ***
! *******************

  IF (mype_filter == 0) THEN
     IF (firstio) THEN
        OPEN(unit = 20, file = filename, status = 'replace')
        firstio = .FALSE.
     ELSE
        OPEN(unit = 20, file = filename, status = 'old', position = 'append')
     END IF

     WRITE (20, *) ABS(step), rmserror_est, rmserror_true
    
     CLOSE(20)
  END IF


! ********************
! *** finishing up ***
! ********************

  DEALLOCATE(variance)

END SUBROUTINE prepoststep_enkf
