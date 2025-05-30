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
  USE mod_parallel_pdaf, &        ! Parallelization variables
       ONLY: mype_world, n_modeltasks, task_id, &
       COMM_model, COMM_filter, COMM_couple, filterpe, abort_parallel
  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: dim_state_p, screen, filtertype, subtype, dim_ens, &
       type_iau, steps_iau, forget, locweight, cradius, sradius, delt_obs, &
       type_opt, dim_cvec, dim_cvec_ens, mcols_cvec_ens, beta_3dvar, &
       solver_iparam1, solver_iparam2, solver_rparam1, solver_rparam2
  USE obs_OBSTYPE_pdafomi, &      ! Variables for observation OBSTYPE
       ONLY: assim_OBSTYPE, rms_obs_OBSTYPE

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: filter_param_i(7) ! Integer parameter array for filter
  REAL    :: filter_param_r(4) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag

! *** External subroutines ***
  EXTERNAL :: init_ens_pdaf            ! Ensemble initialization
  EXTERNAL :: next_observation_pdaf, & ! Provide time step, model time, 
                                       ! and dimension of next observation
       distribute_state_pdaf, &        ! Routine to distribute a state vector to model fields
       prepoststep_pdaf                ! User supplied pre/poststep routine
  EXTERNAL :: init_3dvar_pdaf, &       ! Initialize state and B-matrix for 3D-Var
       prepoststep_3dvar_pdaf          ! User supplied pre/poststep routine
  

! ***************************
! ***   Initialize PDAF   ***
! ***************************

  IF (mype_world == 0) THEN
     WRITE (*,'(/1x,a)') 'INITIALIZE PDAF - ONLINE MODE'
  END IF

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE init_pdaf.F90: Initialize state dimension here!'

  ! *** Define state dimension ***
!  dim_state = ?
  dim_state_p = 10  ! + Dummy value to be able to compile


! **********************************************************
! ***   CONTROL OF PDAF - used in call to PDAF_init      ***
! **********************************************************

! *** IO options ***
  screen = 2         ! Write screen output (1) for output, (2) add timings

! *** Size of control vector and ensemble size  ***
  dim_ens = 9         ! Size of ensemble for ensemble and hybrid Var
  dim_cvec = dim_ens  ! dimension of control vector (parameterized part)
  mcols_cvec_ens = 1  ! Multiplication factor for ensemble control vector (to simulate localization)

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

  ! Set parameters for solver; one could also try to use the defaults
  IF (type_opt==1) THEN
     ! Solver: LBFGS
     solver_iparam1 = 5      ! Number of corrections used in limited memory matrix; 3<=m<=20
     solver_iparam2 = 0      ! -Not used-
     solver_rparam1 = 1.0e-5 ! Parameter 'pgtol'; limit for stopping iterations
     solver_rparam2 = 1.0e+7 ! Parameter 'factr'; tolerance in termination test
  ELSEIF (type_opt==2 .OR. type_opt==12) THEN
     ! Solver: CG+
     solver_iparam1 = 2      ! Parameter 'method'; (1) Fletcher-Reeves, (2) Polak-Ribiere, (3) positive Polak-Ribiere
     solver_iparam2 = 1      ! Parameter 'irest'; (0) no restarts; n>0 restart every n steps
     solver_rparam1 = 1.0e-5 ! Convergence parameter 'eps'
     solver_rparam2 = 0.0    ! -Not used-
  ELSEIF (type_opt==3 .OR. type_opt==13) THEN
     ! Solver: CG+
     solver_iparam1 = 200    ! Maximum number of iterations
     solver_iparam2 = 0      ! -Not used-
     solver_rparam1 = 1.0e-7 ! Convergence parameter 'eps'
     solver_rparam2 = 0.0    ! -Not used-
  END IF

  ! Incremental updating (IAU)
  type_iau = 0       ! Type of incremental updating
  steps_iau = 1      ! Number of time steps over which IAU is applied


! *********************************************************************
! ***   Settings for analysis steps  - used in call-back routines   ***
! *********************************************************************

! *** Forecast length (time interval between analysis steps) ***
  delt_obs = 2      ! This should be set according to the data availability

! *** Which observation type to assimilate
  assim_OBSTYPE = .true.

! *** specifications for observations ***
  rms_obs_OBSTYPE = 0.5    ! Observation error standard deviation

! *** Localization settings
  locweight = 0     ! Type of localizating weighting
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRADIUS
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
  cradius = 2.0     ! Cut-off radius for observation domain in local filters
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


! *****************************************************
! *** Call PDAF initialization routine on all PEs.  ***
! ***                                               ***
! *** Here, only the calls for 3D-Var schemes are   ***
! *** implemented.                                  ***
! ***                                               ***
! *** For all methods, first the arrays of integer  ***
! *** and real number parameters are initialized.   ***
! *** Subsequently, PDAF_init is called.            ***
! *****************************************************

  ! *** 3D-Var ***

  filter_param_i(1) = dim_state_p    ! State dimension
  filter_param_i(2) = dim_ens        ! Size of ensemble
  filter_param_i(3) = type_opt       ! Choose type of optimizer
  filter_param_i(4) = dim_cvec       ! Dimension of control vector (parameterized part)
  filter_param_i(5) = dim_cvec_ens   ! Dimension of control vector (ensemble part)
  filter_param_i(6) = solver_iparam1 ! Parameter setting for solver
  filter_param_i(7) = solver_iparam2 ! Parameter setting for solver
  filter_param_r(1) = forget         ! Forgetting factor
  filter_param_r(2) = beta_3dvar     ! Hybrid weight for hybrid 3D-Var
  filter_param_i(3) = solver_rparam1 ! Parameter setting for solver
  filter_param_i(4) = solver_rparam2 ! Parameter setting for solver

  IF (subtype==0) THEN
     ! parameterized 3D-Var
     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 7,&
          filter_param_r, 4, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_3dvar_pdaf, &
          screen, status_pdaf)
  ELSE
     ! Ensemble or hybrid 3D-Var
     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 7,&
          filter_param_r, 4, &
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


! **********************
! *** Initialize IAU ***
! **********************
  
  CALL PDAF_iau_init(type_iau, steps_iau, status_pdaf)


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
