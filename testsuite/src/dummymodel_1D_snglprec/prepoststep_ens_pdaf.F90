!$Id: prepoststep_ens_pdaf.F90 1711 2016-12-20 10:17:58Z lnerger $
!BOP
!
! !ROUTINE: prepoststep_ens_pdaf --- Used-defined Pre/Poststep routine for PDAF
!
! !INTERFACE:
SUBROUTINE prepoststep_ens_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: EnKF/SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
! 
! NOTE: For ETKF, subtype=3 (fixed covariance matrix), the routine
! prepoststep_etkf_pdaf has to be used. IT is required because the
! size of Uinv is different in ETKF and SEIK/ESTKF.
!
! The routine is called for all ensemble filters before the analysis
! and after the ensemble transformation. For local filters like LSEIK
! the routine is called before and after the loop over all local
! analysis domains. Also it is called once at the initial time before 
! any forecasts are computed.
! The routine provides full access to the state estimate and the
! state ensemble to the user. Thus, user-controlled pre- and poststep 
! operations can be performed here. For example the forecast and the
! analysis states and ensemble covariance matrix can be analyzed,
! e.g. by computing the estimated variances. In addition, the
! estimates can be written to disk. If a user considers to perform
! adjustments to the estimates (e.g. for balances), this routine is 
! the right place for it.
!
! The routine is called by all filter processes.
!
! For the dummy model with domain decomposition 
! we compute the estimated and the true estimation 
! errors. These values are then written into a 
! file.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_memcount, &
       ONLY: memcount
  USE mod_parallel, &
       ONLY: mype_filter, npes_filter, COMM_filter, MPI_REAL, &
       MPIerr, MPIstatus
  USE mod_model, &
       ONLY: dim_state, local_dims, dt, step_null
  USE mod_assimilation, &
       ONLY: incremental, filename, subtype, covartype, dim_lag

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
     ! (When the routine is called before the analysis -step is provided.)
  INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   ! PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   ! PE-local dimension of observation vector
  REAL, INTENT(inout) :: state_p(dim_p) ! PE-local forecast/analysis state
  ! The array 'state_p' is not generally not initialized in the case of SEIK.
  ! It can be used freely here except for subtype==3, where it must not be changed.
  REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) ! Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      ! PE-local state ensemble
  INTEGER, INTENT(in) :: flag        ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_get_state      (as U_prepoststep)
! Called by: PDAF_seik_update    (as U_prepoststep)
! Called by: PDAF_lseik_update    (as U_prepoststep)
! Calls: PDAF_add_increment
! Calls: PDAF_seik_TtimesA
! Calls: memcount
! Calls: MPI_send
! Calls: MPI_recv
!EOP

! *** local variables ***
  INTEGER :: i, j, member              ! counters
  INTEGER, SAVE :: allocflag = 0       ! Flag for memory counting
  LOGICAL, SAVE :: firstio = .TRUE.    ! File output is peformed for first time?
  LOGICAL, SAVE :: initialstep         ! Whether routine is called at the initial time step
  REAL :: invdim_ens                   ! Inverse ensemble size
  REAL :: invdim_ensm1                 ! Inverse of ensemble size minus 1
  REAL :: rmserror_est                 ! estimated RMS error
  REAL :: rmserror_true                ! true RMS error
  REAL :: rmserror_rel                 ! relative error in estimated error
  REAL, ALLOCATABLE :: variance(:)     ! model state variances
  REAL, ALLOCATABLE :: stateinc_p(:)   ! local temporary vector
  REAL, ALLOCATABLE :: truevariance(:) ! model state variances
  REAL, ALLOCATABLE :: truefield_p(:)  ! true local model state

  ! Variables for parallelization - local fields
  INTEGER :: offset   ! Row-offset according to domain decomposition
  REAL, ALLOCATABLE :: variance_p(:)     ! local variance
  REAL, ALLOCATABLE :: truevariance_p(:) ! local model state variances


! **********************
! *** INITIALIZATION ***
! **********************

  IF (step - step_null == 0) THEN
     IF (mype_filter == 0) &
          WRITE (*, '(i7, 3x, a)') step, 'Analyze initial state ensemble'
     initialstep = .TRUE.
  ELSE IF (step > 0) THEN
     IF (mype_filter == 0) &
          WRITE (*, '(8x, a)') 'Analyze assimilated state ensemble'
     initialstep = .FALSE.
  ELSE IF (step < 0) THEN
     IF (mype_filter == 0) &
          WRITE (*, '(8x, a)') 'Analyze forecasted state ensemble'
     initialstep = .FALSE.
  END IF

  ! Allocate fields
  ALLOCATE(variance(dim_state))
  ALLOCATE(variance_p(dim_p))
  IF (allocflag == 0) THEN
     ! count allocated memory
     CALL memcount(3, 'r', dim_state + dim_p)
     IF (subtype == 3) THEN
        ! Count also memory for special type-3 prepoststep
        CALL memcount(3, 'r', (dim_ens - 1) * (dim_ens - 1) &
             + (dim_ens - 1) * dim_ens &
             + dim_ens * dim_ens + (dim_ens - 1))
     END IF
  END IF

  ! Initialize numbers
  rmserror_est  = 0.0
  rmserror_true = 0.0
  rmserror_rel  = 0.0
  invdim_ens    = 1.0 / REAL(dim_ens)  
  invdim_ensm1  = 1.0 / REAL(dim_ens - 1)


  fsubtype: IF (subtype /= 3) THEN
! **************************************************************
! *** Perform prepoststep for filters with ensemble          ***
! *** transformation. The state and error information is     ***
! *** completely stored in the ensemble.                     ***
! **************************************************************

     ! *** Compute mean state
     IF (mype_filter == 0) &
          WRITE (*, '(8x, a)') '--- compute ensemble mean'

     ! local 
     state_p = 0.0
     DO member = 1, dim_ens
        DO i = 1, dim_p
           state_p(i) = state_p(i) + ens_p(i, member)
        END DO
     END DO
     state_p(:) = invdim_ens * state_p(:)

     ! *** Compute local sampled variances ***
     variance_p(:) = 0.0
     DO member = 1, dim_ens
        DO j = 1, dim_p
           variance_p(j) = variance_p(j) &
             + (ens_p(j, member) - state_p(j)) &
             * (ens_p(j, member) - state_p(j))
        END DO
     END DO
     IF (covartype == 1) THEN
        ! For covariance matrix with factor r^-1 (new SEIK - real ensemble)
        variance_p(:) = invdim_ensm1 * variance_p(:)
     ELSE
        ! For covariance matrix with factor (r+1)^-1 (old SEIK)
        variance_p(:) = invdim_ens * variance_p(:)
     END IF

  ELSE fsubtype
! ***********************************************************************
! *** Perform prepoststep for filters with fixed ensemble (subtype=3) ***
! *** In this case, the ensemble mean is the state estimate, but the  ***
! *** ensemble spread is the initial error estimate. To compute the   ***
! *** estimated analysis error one has to compute the analysis        ***
! *** variances from the ensemble and the the transform matrix        ***
! *** Uinv (for SEEK/SEIK) or Ainv (for ESTKF and ETKF).              ***
! ***                                                                 ***
! *** For local filters, we can in general not compute the analysis   ***
! *** error estimate, because the transform matrix is different for   ***
! *** each local analysis domain and only the transform matrix for    ***
! *** the last analysis domain is stored.                             ***
! ***********************************************************************

     IF (mype_filter == 0) &
          WRITE (*, '(8x, a)') 'Subtype=3 variant to compute variances!'

     ! Compute the variance estimates
     CALL comp_variance_estimate(step, dim_p, dim_ens, dim_ens-1, state_p, &
          Uinv, ens_p, variance_p, initialstep)

  END IF fsubtype


! ******************************************************
! *** Assemble global variance vector on filter PE 0 ***
! ******************************************************
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
  IF (allocflag == 0) THEN
     ! count allocated memory
     CALL memcount(3, 'r', dim_state + 2 * dim_p)
  END IF

  truefield_p(:) = 1.0
  DO i = 1, ABS(step)
     truefield_p(:) = truefield_p(:) + 1.0 * dt
  END DO

  ! Add analysis state increment if called directly after analysis
  ALLOCATE(stateinc_p(dim_p))
  ! count allocated memory
  IF (allocflag == 0) CALL memcount(3, 'r', dim_p)
  stateinc_p = state_p
  IF (incremental == 1 .AND. step > 0) THEN
     CALL PDAF_add_increment(dim_p, stateinc_p)
  END IF

  truevariance_p(:) = 0.0
  DO j = 1, dim_p
     truevariance_p(j) = truevariance_p(j) &
          + (stateinc_p(j) - truefield_p(j)) &
          * (stateinc_p(j)  -truefield_p(j))
  END DO

  DEALLOCATE(stateinc_p)

  ! *** assemble global variance vector on filter PE 0
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

     WRITE (20, *) ABS(step), rmserror_est, rmserror_true

     CLOSE(20)
  END IF


! **********************************************
! *** Compute RMS errors for smoothed states ***
! **********************************************

  IF (dim_lag > 0 .AND. step > 0) THEN
     CALL compute_rms_smoother(step, dim_lag, dim_p, dim_ens, state_p)
  END IF


! ********************
! *** finishing up ***
! ********************

  DEALLOCATE(variance)

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE prepoststep_ens_pdaf
