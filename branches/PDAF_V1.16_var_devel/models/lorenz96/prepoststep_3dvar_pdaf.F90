!$Id$
!BOP
!
! !ROUTINE: prepoststep_3dvar_pdaf --- Used-defined Pre/Poststep routine for PDAF
!
! !INTERFACE:
SUBROUTINE prepoststep_3dvar_pdaf(step, dim, dim_ens_g, dim_ens, dim_obs, &
     state, Uinv, ens, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF (all filters):
! 
! The routine is called for global filters before
! the analysis and after the ensemble transformation.
! For local filters, the routine is called before
! and after the loop over all local analysis domains.
! Also it is called once at the initial time before 
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
! This variant is for the Lorenz96 model without
! parallelization.
! We compute the estimated errors. These error
! estimated as well as the state estimate are 
! written into a file.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_memcount, &
       ONLY: memcount
  USE mod_model, &
       ONLY: step_null
  USE mod_modeltime, &
       ONLY: time
  USE mod_assimilation, &
       ONLY: subtype, covartype, stepnull_means, dim_lag, Vmat, dim_cvec
  USE output_netcdf_asml, &
       ONLY: write_netcdf_asml

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
     ! (When the routine is called before the analysis -step is provided.)
  INTEGER, INTENT(in) :: dim         ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens_g   ! Global size of state ensemble
  INTEGER, INTENT(in) :: dim_ens     ! Size of ensemble (=global size here)
  INTEGER, INTENT(in) :: dim_obs     ! Dimension of observation vector
  REAL, INTENT(inout) :: state(dim)  ! Forecast/analysis state
  ! The array 'state' is not generally not initialized in the case of SEIK.
  ! It can be used freely here.
  REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) ! Inverse of matrix U
  REAL, INTENT(inout) :: ens(dim, dim_ens)          ! State ensemble
  INTEGER, INTENT(in) :: flag        ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_get_state      (as U_prepoststep)
! Called by: PDAF_seik_update    (as U_prepoststep)
! Called by: PDAF_lseik_update    (as U_prepoststep)
! Calls: PDAF_seik_TtimesA
! Calls: memcount
! Calls: dgemm (BLAS)
! Calls: dgesv (LAPACK)
!EOP

! *** local variables ***
  INTEGER :: i, j, row, member         ! counters
  INTEGER, SAVE :: allocflag = 0       ! Flag for memory counting
  CHARACTER(len=3) :: calltype         ! Whether routine is called at 
             !(ini) initial time, (for) after forecast, or (ana) after analysis
  REAL :: rmse_est                     ! estimated RMS error
  REAL :: rmse_true                    ! true RMS error
  REAL, ALLOCATABLE :: variance(:)     ! model state variances
  INTEGER, SAVE, ALLOCATABLE :: hist_true(:,:) ! Array for rank histogram about true state
  INTEGER, SAVE, ALLOCATABLE :: hist_mean(:,:) ! Array for rank histogram about ensemble mean
  REAL :: fact                         ! Scaling factor
  ! Variables for mean errors from step 0
  REAL, SAVE :: sum_rmse_est_null(2) = 0.0  ! RMS error estimate accumulated over time
  REAL, SAVE :: sum_rmse_true_null(2) = 0.0 ! True RMS error accumulated over time
  INTEGER, SAVE :: nsum_null(2) = 0         ! Length of sums over time
  REAL :: mrmse_est_null = 0.0              ! Time-mean of estimated RMS error
  REAL :: mrmse_true_null = 0.0             ! Time-mean of true RMS error
  ! Variables for sum from step stepnull_means
  REAL, SAVE :: sum_rmse_est_step(2) = 0.0  ! RMS error estimate accumulated over time
  REAL, SAVE :: sum_rmse_true_step(2) = 0.0 ! True RMS error accumulated over time
  INTEGER, SAVE :: nsum_step(2) = 0         ! Length of sums over time
  REAL :: mrmse_est_step = 0.0              ! Time-mean of estimated RMS error
  REAL :: mrmse_true_step = 0.0             ! Time-mean of true RMS error
  REAL :: skewness                          ! Skewness of ensemble
  REAL :: kurtosis                          ! Kurtosis of ensemble
  ! Variables for smoother erros
  REAL, Allocatable :: rmse_s(:)        ! estimated RMS error of smoothed states
  REAL, ALLOCATABLE :: trmse_s(:)       ! true RMS error of smoothed states
  REAL, ALLOCATABLE :: mrmse_s_null(:)  ! Time-mean of estimated smoother RMS error
  REAL, ALLOCATABLE :: mtrmse_s_null(:) ! Time-mean of true smoother RMS error
  REAL, ALLOCATABLE :: mrmse_s_step(:)  ! Time-mean of estimated smootherRMS error
  REAL, ALLOCATABLE :: mtrmse_s_step(:) ! Time-mean of true smoother RMS error

  ! Variables for type-3 prepoststep
  INTEGER :: dgesv_info                ! output flag for DGESV
  INTEGER, SAVE :: allocflag2 = 0      ! Flag for memory counting
  REAL :: rdim_ens                     ! dim_ens in real format
  REAL, ALLOCATABLE :: Ttrans(:,:)     ! matrix T^T
  REAL, ALLOCATABLE :: TUT(:,:)        ! temporary matrix TUT^T
  INTEGER, ALLOCATABLE :: ipiv(:)      ! vector of pivot indices for SGESV
  REAL, ALLOCATABLE :: tempUinv(:,:)   ! temporary matrix Uinv


! **********************
! *** INITIALIZATION ***
! **********************

  IF (step - step_null == 0) THEN
     WRITE (*, '(i7, 3x, a)') step, 'Analize initial state - for 3D-Var'
     calltype = 'ini'
  ELSE IF (step - step_null > 0) THEN
     WRITE (*, '(8x, a)') 'Analize assimilated state - for 3D-Var'
     calltype = 'ana'
  ELSE IF (step - step_null < 0) THEN
     WRITE (*, '(8x, a)') 'Analize forecasted state - for 3D-Var'
     calltype = 'for'
  END IF

  ! Allocate fields
  ALLOCATE(variance(dim))
  IF (allocflag == 0) THEN
     ! count allocated memory
     CALL memcount(3, 'r', dim)
  END IF

  ! Fields for rank histograms
  IF (calltype == 'ini') THEN
     allocate(hist_true(dim_ens+1, 2))
     allocate(hist_mean(dim_ens+1, 2))

     hist_true = 0
     hist_mean = 0
  END IF

  ! Initialize numbers
  rmse_est  = 0.0
  rmse_true = 0.0


! **************************************************************
! *** Perform prepoststep for 3D-Var in which dim_ens=1      ***
! *** The sampled error is here computed from B^(1/2)        ***
! **************************************************************

  ! *** Initialize state estimate
  state(:) = ens(:,1)

  ! *** Compute local sampled variances ***
  variance(:) = 0.0
  DO member = 1, dim_cvec
     DO j = 1, dim
        variance(j) = variance(j) &
             + Vmat(j,member) * Vmat(j,member)
     END DO
  END DO


! ************************************************************
! *** Compute RMS errors according to sampled covar matrix ***
! ************************************************************

  ! total estimated RMS error
  DO i = 1, dim
     rmse_est = rmse_est + variance(i)
  ENDDO
  rmse_est = SQRT(rmse_est / dim)


! *******************************
! *** Compute true RMS errors ***
! *******************************

  CALL compute_truermse(calltype, step, time, dim, state, &
       rmse_true, dim_ens, ens, hist_true, hist_mean, &
       skewness, kurtosis)


! *******************************************
! *** Compute sum of RMS errors over time ***
! *******************************************

  ! Sums from step 0
  IF (calltype == 'for') THEN
     sum_rmse_est_null(1) = sum_rmse_est_null(1) + rmse_est
     sum_rmse_true_null(1) = sum_rmse_true_null(1) + rmse_true
     nsum_null(1) = nsum_null(1) + 1
     mrmse_est_null = sum_rmse_est_null(1) / REAL(nsum_null(1))
     mrmse_true_null = sum_rmse_true_null(1) / REAL(nsum_null(1))
  ELSE IF (calltype == 'ana') THEN
     sum_rmse_est_null(2) = sum_rmse_est_null(2) + rmse_est
     sum_rmse_true_null(2) = sum_rmse_true_null(2) + rmse_true
     nsum_null(2) = nsum_null(2) + 1
     mrmse_est_null = sum_rmse_est_null(2) / REAL(nsum_null(2))
     mrmse_true_null = sum_rmse_true_null(2) / REAL(nsum_null(2))
  END IF

  ! Sums from step stepnull_means
  IF (ABS(step) >= stepnull_means) THEN
     IF (calltype == 'for') THEN
        sum_rmse_est_step(1) = sum_rmse_est_step(1) + rmse_est
        sum_rmse_true_step(1) = sum_rmse_true_step(1) + rmse_true
        nsum_step(1) = nsum_step(1) + 1
        mrmse_est_step = sum_rmse_est_step(1) / REAL(nsum_step(1))
        mrmse_true_step = sum_rmse_true_step(1) / REAL(nsum_step(1))
     ELSE IF (calltype == 'ana') THEN
        sum_rmse_est_step(2) = sum_rmse_est_step(2) + rmse_est
        sum_rmse_true_step(2) = sum_rmse_true_step(2) + rmse_true
        nsum_step(2) = nsum_step(2) + 1
        mrmse_est_step = sum_rmse_est_step(2) / REAL(nsum_step(2))
        mrmse_true_step = sum_rmse_true_step(2) / REAL(nsum_step(2))
     END IF
  END IF


! *********************
! *** Screen output ***
! *********************

  ! Output RMS errors given by sampled covar matrix
  WRITE (*, '(8x,a)') '--- RMS errors'
  WRITE (*, '(12x, a, es12.4, a, i6, a, 2es12.4)') &
       'Estimated: ', rmse_est, ' time mean (from 0,', &
       stepnull_means, '): ', mrmse_est_null, mrmse_est_step
  WRITE (*, '(17x, a, es12.4, a, i6, a, 2es12.4)') &
       'True: ', rmse_true, ' time mean (from 0,', &
       stepnull_means, '): ', mrmse_true_null, mrmse_true_step


! **********************************************
! *** Compute RMS errors for smoothed states ***
! **********************************************

  smoother: IF (dim_lag > 0 .AND. calltype=='ana') THEN

     ALLOCATE(rmse_s(dim_lag))
     ALLOCATE(trmse_s(dim_lag))
     ALLOCATE(mrmse_s_null(dim_lag))
     ALLOCATE(mtrmse_s_null(dim_lag))
     ALLOCATE(mrmse_s_step(dim_lag))
     ALLOCATE(mtrmse_s_step(dim_lag))
     IF (allocflag2 == 0) THEN 
        CALL memcount(3, 'r', 6 * dim_lag)
        allocflag2 = 1
     END IF

     mrmse_s_null = 0.0
     mtrmse_s_null = 0.0
     mrmse_s_step = 0.0
     mtrmse_s_step = 0.0

     CALL compute_rms_smoother(step, dim_lag, dim, dim_ens, state, &
          variance, rmse_s, trmse_s, mrmse_s_null, mtrmse_s_null, &
          mrmse_s_step, mtrmse_s_step)

  ELSE
     ALLOCATE(rmse_s(1))
     ALLOCATE(trmse_s(1))
     ALLOCATE(mrmse_s_null(1))
     ALLOCATE(mtrmse_s_null(1))
     ALLOCATE(mrmse_s_step(1))
     ALLOCATE(mtrmse_s_step(1))
  END IF smoother

  
! *******************
! *** File output ***
! *******************

  CALL write_netcdf_asml(calltype, step, time, dim, state, &
       rmse_est, rmse_true, mrmse_est_null, mrmse_true_null, mrmse_est_step, &
       mrmse_true_step, dim_ens, ens, hist_true, hist_mean, skewness, &
       kurtosis, dim_lag, rmse_s, trmse_s, mrmse_s_null, mtrmse_s_null, &
       mrmse_s_step, mtrmse_s_step)


! ********************
! *** finishing up ***
! ********************

  DEALLOCATE(variance)

  IF (dim_lag > 0 .AND. calltype=='ana') THEN
     DEALLOCATE(rmse_s, trmse_s)
     DEALLOCATE(mrmse_s_null, mtrmse_s_null, mrmse_s_step, mtrmse_s_step)
  END IF

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE prepoststep_3dvar_pdaf
