!$Id: prepoststep_seek.F90 1590 2015-06-12 11:59:35Z lnerger $
!BOP
!
! !ROUTINE: prepoststep_seek --- Used-defined Pre/Poststep routine for PDAF
!
! !INTERFACE:
SUBROUTINE prepoststep_seek(step, dim_p, dim_eof, dim_eof_p, dim_obs_p, &
     state_p, Uinv, eofV_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF (SEEK):
! 
! The routine is called before the analysis
! and after the re-diagonalization.  Also it 
! is called once at the initial time before 
! any forecasts are computed. 
! The routine provides full access to the state 
! estimate and the state covariance matrix 
! to the user.  Thus, user-controlled pre- and 
! poststep operations can be performed here. 
! For example the forecast and the analysis 
! states and error covariance matrix can be 
! analized, e.g. by computing the estimated 
! variance.  In addition, the estimates can be 
! written to disk.  If a user considers to 
! perform adjustments to the estimates (e.g. 
! for balances), this routine is the right 
! place for it.
!
! The routine is called by all filter processes.
!
! For the dummy model with domain decomposition 
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
       ONLY: dim_state, local_dims, dt, step_null
  USE mod_assimilation, &
       ONLY: incremental, filename

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
     ! (When the routine is called before the analysis -step is provided.)
  INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_eof     ! Number of EOF modes used in SEEK
  INTEGER, INTENT(in) :: dim_eof_p   ! PE-local number of EOF modes
  INTEGER, INTENT(in) :: dim_obs_p   ! PE-local dimension of observation vector
  REAL, INTENT(inout) :: state_p(dim_p)  ! PE-local forecast/analysis state
  ! *** The covariance P is decomposed as P = V U V^T ***
  REAL, INTENT(inout) :: Uinv(dim_eof,dim_eof)   ! Inverse of matrix U
  REAL, INTENT(inout) :: eofV_p(dim_p,dim_eof)   ! PE-local matrix V
  INTEGER, INTENT(in) :: flag        ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: prepoststep_pdaf
! Calls: PDAF_add_increment
! Calls: memcount
! Calls: dgesv (LAPACK)
! Calls: MPI_send
! Calls: MPI_recv
!EOP

! *** local variables ***
  INTEGER :: i, j                      ! counters
  INTEGER, SAVE :: allocflag = 0       ! Flag for memory counting
  REAL, ALLOCATABLE :: variance_p(:)   ! PE-local variance
  REAL, ALLOCATABLE :: variance(:)     ! model state variances
  REAL, ALLOCATABLE :: state_p_tmp(:)  ! temporary state vector
  REAL, ALLOCATABLE :: truevariance_p(:) ! local model state variances
  REAL, ALLOCATABLE :: truevariance(:) ! model state variances
  REAL, ALLOCATABLE :: truefield_p(:)  ! true model state
  REAL, ALLOCATABLE :: temp_Uinv(:,:)  ! temporary storage of Uinv
  REAL, ALLOCATABLE :: eofVt(:,:)      ! Array for tranpose of eofV
  INTEGER, ALLOCATABLE :: ipiv(:)      ! vector of pivot indices for SGESV
  REAL :: rmserror_est                 ! estimated RMS error
  REAL :: rmserror_true                ! true RMS error
  REAL :: rmserror_rel                 ! relative error in estimated error
  LOGICAL, SAVE :: firstio = .TRUE.    ! File output is performed for first time?
  INTEGER :: sgesv_info                ! output flag of SGESV
  REAL, EXTERNAL :: ddot               ! BLAS scalar product

  ! Variables for parallelization - local fields
  INTEGER :: offset   ! Row-offset according to domain decomposition


! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype_filter == 0) THEN
     IF (step - step_null == 0) THEN
        WRITE (*, '(i7, 3x, a)') step, 'Analize initial state ensemble'
     ELSE IF (step > 0) THEN
        WRITE (*, '(8x, a)') 'Analize assimilated state ensemble'
     ELSE IF (step < 0) THEN
        WRITE (*, '(8x, a)') 'Analize forecast state ensemble'
     END IF
  END IF

  ! *** Allocate fields: gather mode matrix on each PE ***
  ALLOCATE(variance(dim_state))
  ALLOCATE(variance_p(dim_p))
  ALLOCATE(temp_Uinv(dim_eof, dim_eof))
  ALLOCATE(eofVt(dim_eof, dim_p))
  ALLOCATE(ipiv(dim_eof))
  IF (allocflag == 0) THEN
     ! count allocated memory
     CALL memcount(3, 'r', 2 * dim_state + dim_eof**2 + dim_eof * dim_state)
     CALL memcount(3, 'i', dim_eof)
  END IF

  ! Initialize numbers
  rmserror_est  = 0.0
  rmserror_true = 0.0
  rmserror_rel  = 0.0


! *********************************
! *** Compute sampled variances ***
! *********************************

  ! save matrix Uinv
  temp_Uinv = Uinv

  ! initialized transposed of eofV
  eofVt = TRANSPOSE(eofV_p)

  ! *** call solver - compute X = U V^T
  CALL dgesv(dim_eof, dim_p, temp_Uinv, dim_eof, ipiv, &
       eofVt, dim_eof, sgesv_info)

  ! *** check if solve was successful
  update: IF (sgesv_info /= 0) THEN
     IF (mype_filter == 0) WRITE (*, '(/2x, a/)') &
          '!!! Problem in solve for state analysis !!!'
  ELSE

     ! compute variances from P = V U V^T = V X
     variance_p(:) = 0.0
     DO i = 1, dim_p
        variance_p(i) = variance_p(i)  &
             + DDOT(dim_eof, eofV_p(i, 1), dim_p, eofVt(1, i), 1)
     END DO

! *** assemble global variance vector on filter PE 0
     PE0_a: IF (mype_filter /= 0) THEN

        ! send sub-fields from PEs /=0
        CALL MPI_send(variance_p(1 : dim_p), dim_p, &
             MPI_DOUBLE_PRECISION, 0, mype_filter, COMM_filter, MPIerr)

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
             MPI_DOUBLE_PRECISION, i - 1, i - 1, COMM_filter, MPIstatus, MPIerr)
        END DO
      
     END IF PE0_a

  END IF update

  DEALLOCATE(variance_p)


! ************************************************************
! *** Compute RMS errors according to sampled covar matrix ***
! ************************************************************

  ! *** Only on PE 0 ***
  pe0: IF (mype_filter == 0) THEN
     ! total estimated RMS error
     DO i = 1, dim_state
        rmserror_est = rmserror_est + variance(i)
     ENDDO
     rmserror_est = SQRT(rmserror_est / dim_state)
  END IF pe0

  ! *** Compute true variances
  ALLOCATE(truefield_p(dim_p))
  ALLOCATE(truevariance_p(dim_p))
  ALLOCATE(truevariance(dim_state))
  ALLOCATE(state_p_tmp(dim_p))
  IF (allocflag == 0) THEN
     ! count allocated memory
     CALL memcount(3, 'r', dim_state + 3 * dim_p)
     allocflag = 1
  END IF

  ! Initialize temporary state vector
  state_p_tmp = state_p

  ! Add analysis state increment if called directly after analysis
  IF (incremental == 1 .AND. step > 0) THEN
     CALL PDAF_add_increment(dim_p, state_p_tmp)
  END IF

  truefield_p(:) = 1.0
  DO i = 1, ABS(step)
     truefield_p(:) = truefield_p(:) + 1.0 * dt
  END DO

  truevariance_p(:) = 0.0
  DO j = 1, dim_p
     truevariance_p(j) = truevariance_p(j) &
         + (state_p_tmp(j) - truefield_p(j)) * (state_p_tmp(j) - truefield_p(j))
  END DO

  ! *** assemble global variance vector on filter PE 0
  PE0_b: IF (mype_filter /= 0) THEN

     ! send sub-fields from PEs /=0
     CALL MPI_send(truevariance_p(1 : dim_p), dim_p, &
          MPI_DOUBLE_PRECISION, 0, mype_filter, COMM_filter, MPIerr)

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
             local_dims(i), MPI_DOUBLE_PRECISION, &
             i - 1, i - 1, COMM_filter, MPIstatus, MPIerr)
     END DO
      
  END IF PE0_b

  pe0a: IF (mype_filter == 0) THEN

     ! total true RMS error
     DO i = 1, dim_state
        rmserror_true = rmserror_true + truevariance(i)
     ENDDO
     rmserror_true = SQRT(rmserror_true / dim_state)

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
     
     WRITE (20,*) ABS(step), rmserror_est, rmserror_true
    
     CLOSE(20)
  END IF


! ********************
! *** finishing up ***
! ********************

  DEALLOCATE(variance, temp_Uinv, eofVt, ipiv)

END SUBROUTINE prepoststep_seek
