!>  Interface routine to call initialization of PDAF
!!
!! This routine collects the initialization of variables for PDAF.
!! In addition, the initialization routine PDAF_init is called
!! to perform the internal initialization of PDAF.
!!
!! This variant is for the online mode of PDAF
!! and only for the 3D-Var variants!
!!
!! This routine is generic. However, it assumes a constant observation
!! error (rms_obs_X, etc.). Further, with parallelization the local state
!! dimension dim_state_p is used.
!!
!! __Revision history:__
!! * 2008-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_pdaf()

  USE PDAF                        ! PDAF interface definitions
  USE mod_parallel_model, &       ! Parallelization variables for model
       ONLY: mype_world, mype_model, npes_model, COMM_model, abort_parallel
  USE mod_parallel_pdaf, &        ! Parallelization variables fro assimilation
       ONLY: n_modeltasks, task_id, COMM_filter, COMM_couple, filterpe
  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: dim_state_p, dim_state, screen, filtertype, subtype, dim_ens, &
       type_forget, forget, &
       locweight, cradius, sradius, &
       type_trans, type_sqrt, delt_obs, &
       type_opt, dim_cvec, dim_cvec_ens, mcols_cvec_ens, beta_3dvar, &
       dims_cv_ens_p, off_cv_ens_p, dims_cv_p, off_cv_p, &
       observe_ens, type_obs_init, do_omi_obsstats
  USE mod_model, &                ! Model variables
       ONLY: nx, ny, nx_p
  USE obs_A_pdafomi, &            ! Variables for observation type A
       ONLY: assim_A, rms_obs_A
  USE obs_B_pdafomi, &            ! Variables for observation type B
       ONLY: assim_B, rms_obs_B

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: filter_param_i(5) ! Integer parameter array for filter
  REAL    :: filter_param_r(2) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag
  INTEGER :: i                 ! Counter
  INTEGER :: dim_cvec_ens_p    ! PE-local dimension of ensemble control vector
  INTEGER :: dim_cvec_p        ! PE-local dimension of control vector

! *** External subroutines ***
  EXTERNAL :: init_ens_pdaf            ! Ensemble initialization
  EXTERNAL :: next_observation_pdaf, & ! Provide time step, model time, 
                                       ! and dimension of next observation
       distribute_state_pdaf, &        ! Routine to distribute a state vector to model fields
       prepoststep_ens_pdaf            ! User supplied pre/poststep routine
  EXTERNAL :: init_3dvar_pdaf, &       ! Initialize state and B-matrix for 3D-Var
       prepoststep_3dvar_pdaf          ! User supplied pre/poststep routine
  

! ***************************
! ***   Initialize PDAF   ***
! ***************************

  IF (mype_world == 0) THEN
     WRITE (*,'(/1x,a)') 'INITIALIZE PDAF - ONLINE MODE'
  END IF

  ! *** Define state dimension ***
  dim_state_p = nx_p * ny  ! Local state dimension
  dim_state = nx * ny      ! Global state dimension


! **********************************************************
! ***   CONTROL OF PDAF - used in call to PDAF_init      ***
! **********************************************************

! *** IO options ***
  screen = 2         ! Write screen output (1) for output, (2) add timings

! *** Size of control vector and ensemble size  ***
  dim_ens = n_modeltasks  ! Size of ensemble for all ensemble filters
                     !   We use n_modeltasks here, initialized in init_parallel_pdaf
  dim_cvec = dim_ens ! dimension of control vector (parameterized part)
  mcols_cvec_ens = 1 ! Multiplication factor for ensemble control vector (to simulate localization)

! *** Options for 3D-Var method

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++
  ! +++ For available options see MOD_ASSIMILATION +++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++

  filtertype = 200   ! Type of DA method
  subtype = 0        ! Subtype of 3D-Var

  forget  = 1.0      ! Forgetting factor value for inflation in EnVar and hybrid 

  type_opt = 1       ! Type of minimizer for 3DVar
                     !   (1) LBFGS, (2) CG+, (3) plain CG
                     !   (12) CG+ parallel, (13) plain CG parallel
  beta_3dvar = 0.5   ! Hybrid weight for hybrid 3D-Var


! *********************************************************************
! ***   Settings for analysis steps  - used in call-back routines   ***
! *********************************************************************

! *** Forecast length (time interval between analysis steps) ***
  delt_obs = 2     ! Number of time steps between analysis/assimilation steps

! *** Which observation type to assimilate
  assim_A = .true.
  assim_B = .false.

! *** specifications for observations ***
  rms_obs_A = 0.5    ! Observation error standard deviation for observation A
  rms_obs_B = 0.5    ! Observation error standard deviation for observation B

! *** Localization settings
  locweight = 0     ! Type of localizating weighting
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRADIUS
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
  cradius = 0.0     ! Cut-off radius for observation domain in local filters
  sradius = cradius ! Support radius for 5th-order polynomial
                    ! or radius for 1/e for exponential weighting


! ***********************************
! *** Some optional functionality ***
! ***********************************

! *** Parse command line options   ***
! *** This is optional, but useful ***

  call init_pdaf_parse()

  ! Set size of control vector for ensemble 3D-Var
  ! Using mcols_cvec_ens simulates the effect when localization would be applied
  dim_cvec_ens = dim_ens * mcols_cvec_ens


! *** Initial Screen output ***
! *** This is optional      ***

  IF (mype_world == 0) call init_pdaf_info()


! **************************************************
! *** Initialize decomposition of control vector ***
! **************************************************

  ! Parameterized part

  IF (filtertype==200 .AND. subtype>0 .AND. (type_opt==12 .OR. type_opt==13)) THEN

     ! split control vector
     ALLOCATE (dims_cv_ens_p(npes_model))
     ALLOCATE (off_cv_ens_p(npes_model))

     dims_cv_ens_p = FLOOR(REAL(dim_cvec_ens) / REAL(npes_model))
     DO i = 1, (dim_cvec_ens - npes_model * dims_cv_ens_p(1))
        dims_cv_ens_p(i) = dims_cv_ens_p(i) + 1
     END DO

     off_cv_ens_p(1) = 0
     DO i = 2, npes_model
        off_cv_ens_p(i) = off_cv_ens_p(i-1) + dims_cv_ens_p(i-1)
     END DO

     IF (mype_world == 0) THEN
        WRITE (*, '(/2x, a, i3, a)') &
             '-- Decomposition of control vector over', npes_model, ' PEs'
        DO i = 1, npes_model
           WRITE (*, '(5x, a, i3, a, i3, a, 2i5)') &
                'task ', task_id, ' PE(model) ', i-1, &
                ' dims_cv_ens_p, off_cv_ens_p: ', dims_cv_ens_p(i), off_cv_ens_p(i)
        END DO
     END IF

     ! Set dimension of control vector for my PE-local domain
     dim_cvec_ens_p = dims_cv_ens_p(mype_model + 1)
  ELSE
     dim_cvec_ens_p = dim_cvec_ens
  END IF

  ! Ensemble part of control vector

  IF (filtertype==200 .AND. (subtype==0 .OR. subtype==3 .OR. subtype==4) &
       .AND. (type_opt==12 .OR. type_opt==13)) THEN

     ! split control vector
     ALLOCATE (dims_cv_p(npes_model))
     ALLOCATE (off_cv_p(npes_model))

     dims_cv_p = FLOOR(REAL(dim_cvec) / REAL(npes_model))
     DO i = 1, (dim_cvec - npes_model * dims_cv_p(1))
        dims_cv_p(i) = dims_cv_p(i) + 1
     END DO

     off_cv_p(1) = 0
     DO i = 2, npes_model
        off_cv_p(i) = off_cv_p(i-1) + dims_cv_p(i-1)
     END DO

     IF (mype_world == 0) THEN
        WRITE (*, '(/2x, a, i3, a)') &
             '-- Decomposition of control vector over', npes_model, ' PEs'
        DO i = 1, npes_model
           WRITE (*, '(5x, a, i3, a, i3, a, 2i5)') &
                'task ', task_id, ' PE(model) ', i-1, &
                ' dims_cv_p, off_cv_p: ', dims_cv_p(i), off_cv_p(i)
        END DO
     END IF

     ! Set dimension of control vector for my PE-local domain
     dim_cvec_p = dims_cv_p(mype_model + 1)
  ELSE
     dim_cvec_p = dim_cvec
  END IF


! *****************************************************
! *** Call PDAF initialization routine on all PEs.  ***
! ***                                               ***
! *** Here, only the calls for 3D-Var schemes are   ***
! *** implemented.                                  ***
! ***                                               ***
! *** For all methods, first the arrays of integer  ***
! *** and real number parameters are initialized.   ***
! *** Subsequently, PDAF_init is called.            ***
! ***                                               ***
! *** As for the ensemble filters, one can after    ***
! *** call to PDAF_init insert calls to             ***
! *** PDAF_set_iparam and PDAF_set_rparam to set    ***
! *** further optional parameters for PDAF.         ***
! *****************************************************

  ! *** 3D-Var ***

  filter_param_i(1) = dim_state_p    ! State dimension
  filter_param_i(2) = dim_ens        ! Size of ensemble
  filter_param_i(3) = type_opt       ! Choose type of optimizer
  filter_param_i(4) = dim_cvec_p     ! Dimension of control vector (parameterized part)
  filter_param_i(5) = dim_cvec_ens_p ! Dimension of control vector (ensemble part)
  filter_param_r(1) = forget         ! Forgetting factor
  filter_param_r(2) = beta_3dvar     ! Hybrid weight for hybrid 3D-Var

  IF (subtype==0) THEN
     ! parameterized 3D-Var
     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 5,&
          filter_param_r, 1, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_3dvar_pdaf, &
          screen, status_pdaf)
  ELSE
     ! Ensemble or hybrid 3D-Var
     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 5,&
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_pdaf, &
          screen, status_pdaf)
  END IF


! *** Check whether initialization of PDAF was successful ***
  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in initialization of PDAF - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF


! **********************************
! *** Prepare ensemble forecasts ***
! **********************************

  IF (.NOT. (filtertype==200 .AND. subtype==0)) THEN
     ! For 3D ensemble Var and hybrid Var
     CALL PDAF_init_forecast(next_observation_pdaf, distribute_state_pdaf, &
          prepoststep_ens_pdaf, status_pdaf)
  ELSE
     ! For parameterized 3D-Var
     CALL PDAF_init_forecast(next_observation_pdaf, distribute_state_pdaf, &
          prepoststep_3dvar_pdaf, status_pdaf)
  END IF

END SUBROUTINE init_pdaf
