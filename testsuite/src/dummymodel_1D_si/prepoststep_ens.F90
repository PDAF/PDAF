!$Id: prepoststep_ens.F90 1565 2015-02-28 17:04:41Z lnerger $
!BOP
!
! !ROUTINE: prepoststep_ens --- Used-defined Pre/Poststep routine for PDAF
!
! !INTERFACE:
SUBROUTINE prepoststep_ens(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF (ensemble filters):
! 
! The routine is called for all ensemble filters
! before the analysis and after the ensemble
! transformation. For the local filters (LSEIK and
! LETKF) the routine is called before and after the
! loop  over all local analysis domains. Also it 
! is called once at the initial time before 
! any forecasts are computed.
! The routine provides full access to the state 
! estimate and the state ensemble to the user.
! Thus, user-controlled pre- and poststep 
! operations can be performed here. For example 
! the forecast and the analysis states and ensemble
! covariance matrix can be analized, e.g. by 
! computing the estimated variances. In addition, 
! the estimates can be written to disk. If a user 
! considers to perform adjustments to the 
! estimates (e.g. for balances), this routine is 
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
       ONLY: mype_filter, npes_filter, COMM_filter, MPI_DOUBLE_PRECISION, &
       MPIerr, MPIstatus
  USE mod_model, &
       ONLY: dim_state, local_dims, dt, step_null
  USE mod_assimilation, &
       ONLY: incremental, filename, subtype, covartype

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
  ! It can be used freely here.
  REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) ! Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      ! PE-local state ensemble
  INTEGER, INTENT(in) :: flag        ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: prepoststep_pdaf
! Calls: PDAF_add_increment
! Calls: PDAF_seik_TtimesA
! Calls: memcount
! Calls: dgemm (BLAS)
! Calls: dgesv (LAPACK)
! Calls: MPI_send
! Calls: MPI_recv
!EOP

! *** local variables ***
  INTEGER :: i, j, row, member         ! counters
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

  ! Variables for type-3 prepoststep
  INTEGER :: dgesv_info                ! output flag for DGESV
  REAL :: rdim_ens                     ! dim_ens in real format
  REAL, ALLOCATABLE :: Ttrans(:,:)     ! matrix T^T
  REAL, ALLOCATABLE :: TUT(:,:)        ! temporary matrix TUT^T
  INTEGER, ALLOCATABLE :: ipiv(:)      ! vector of pivot indices for SGESV
  REAL, ALLOCATABLE :: tempUinv(:,:)   ! temporary matrix Uinv
  


! **********************
! *** INITIALIZATION ***
! **********************

  IF (step - step_null == 0) THEN
     IF (mype_filter == 0) &
          WRITE (*, '(i7, 3x, a)') step, 'Analize initial state ensemble'
     initialstep = .TRUE.
  ELSE IF (step > 0) THEN
     IF (mype_filter == 0) &
          WRITE (*, '(8x, a)') 'Analize assimilated state ensemble'
     initialstep = .FALSE.
  ELSE IF (step < 0) THEN
     IF (mype_filter == 0) &
          WRITE (*, '(8x, a)') 'Analize forecasted state ensemble'
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
! *** Perform prepoststep for SEIK with re-inititialization. ***
! *** The state and error information is completely in the   ***
! *** ensemble.                                              ***
! *** Also performed for SEIK without re-init at the initial ***
! *** time.                                                  ***
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
! ******************************************************************
! *** Perform prepoststep for SEIK without re-init (subtype=3)   ***
! *** The state and error information is stored in the           ***
! *** forecast ensemble, the state estimate and the matrix Uinv. ***
! ***                                                            ***    
! *** We compute the variance of the covariance matrix which     ***
! *** is given as                                                ***
! ***                            T    T                          ***
! ***                 P = X T U T (X )                           ***
! ***                  i   i   i    i                            ***
! ******************************************************************

     IF (mype_filter == 0) &
          WRITE (*, '(8x, a)') 'Type-3 variant to compute variances!'

     IF (initialstep) THEN
        ! local 
        state_p = 0.0
        DO member = 1, dim_ens
           DO i = 1, dim_p
              state_p(i) = state_p(i) + ens_p(i, member)
           END DO
        END DO
        state_p(:) = invdim_ens * state_p(:)
     END IF

     ! Allocate fields
     ALLOCATE(tempUinv(dim_ens - 1, dim_ens - 1))
     ALLOCATE(Ttrans(dim_ens - 1, dim_ens))
     ALLOCATE(TUT(dim_ens, dim_ens))
     ALLOCATE(ipiv(dim_ens - 1)) 

     ! Initialize matrix T^T
     DO i = 1, dim_ens - 1
        DO j = 1, dim_ens
           Ttrans(i, j) = -invdim_ens
        END DO
     END DO
     DO i = 1, dim_ens - 1
        Ttrans(i, i) = Ttrans(i, i) + 1.0
     END DO

     IF (step > 0 .AND. (step - step_null /= 0)) THEN
        ! Initialize temporary Uinv (We must not change Uinv here!)
        tempUinv(:, :) = Uinv(:, :)
     ELSE IF (step < 0 .OR. (step - step_null == 0)) THEN
        ! Initialize invariant Uinv (dim_ens T T^T)
        IF (covartype == 1) THEN
           ! For covariance matrix with factor r^-1 (new SEIK - real ensemble)
           rdim_ens = REAL(dim_ens - 1)
        ELSE
           ! For covariance matrix with factor (r+1)^-1 (old SEIK)
           rdim_ens = REAL(dim_ens)
        END IF
        CALL dgemm('n', 't', dim_ens - 1, dim_ens - 1, dim_ens, &
             rdim_ens, Ttrans, dim_ens - 1, Ttrans, dim_ens - 1, &
             0.0, tempUinv, dim_ens - 1)
     END IF

     ! call solver - compute W = U T^T
     CALL dgesv(dim_ens - 1, dim_ens, tempUinv, dim_ens - 1, ipiv, &
          Ttrans, dim_ens - 1, dgesv_info) 

     ! Compute T W = T U T^T using operation in subroutine
     CALL PDAF_seik_TtimesA(dim_ens - 1, dim_ens, Ttrans, TUT)

     ! Compute local sampled variances
     variance_p(:) = 0.0
     DO i = 1, dim_ens
        DO j = 1, dim_ens
           DO row = 1, dim_p
              variance_p(row) = variance_p(row) &
                   + ens_p(row, j) * ens_p(row, i) * TUT(i, j)
           END DO
        END DO
     END DO
    
     DEALLOCATE(Ttrans, TUT)

  END IF fsubtype


! ******************************************************
! *** Assemble global variance vector on filter PE 0 ***
! ******************************************************
  PE0_a: IF (mype_filter /= 0) THEN

     ! send sub-fields from PEs /=0
     CALL MPI_send(variance_p(1 : dim_p), dim_p, &
          MPI_DOUBLE_PRECISION,0, mype_filter, COMM_filter, MPIerr)

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

     WRITE (20, *) ABS(step), rmserror_est, rmserror_true

     CLOSE(20)
  END IF

! ********************
! *** finishing up ***
! ********************

  DEALLOCATE(variance)

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE prepoststep_ens
