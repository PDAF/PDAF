!>  Interface routine to call initialization of PDAF
!!
!! This routine collects the initialization of variables for PDAF.
!! In addition, the initialization routine PDAF_init is called
!! to perform the internal initialization of PDAF.
!!
!! This variant is for the offline mode of PDAF
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

  USE PDAF, &                     ! Interface definitions to PDAF core routines
       ONLY: PDAF_init, PDAF_set_offline_mode
  USE mod_parallel_pdaf, &        ! Parallelization variables
       ONLY: mype_world, n_modeltasks, task_id, &
       COMM_model, COMM_filter, COMM_couple, filterpe, abort_parallel
  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: dim_state_p, screen, filtertype, subtype, dim_ens, &
       incremental, type_forget, forget, &
       locweight, cradius, sradius, &
       type_trans, type_sqrt, &
       type_opt, dim_cvec, dim_cvec_ens, mcols_cvec_ens, &
       beta_3dvar
  USE obs_A_pdafomi, &            ! Variables for observation type A
       ONLY: assim_A, rms_obs_A
  USE obs_B_pdafomi, &            ! Variables for observation type B
       ONLY: assim_B, rms_obs_B
  USE obs_C_pdafomi, &            ! Variables for observation type C
       ONLY: assim_C, rms_obs_C

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: filter_param_i(5) ! Integer parameter array for filter
  REAL    :: filter_param_r(2) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag

! *** External subroutines ***
  EXTERNAL :: init_ens_offline     ! Ensemble initialization
  EXTERNAL :: init_3dvar_offline   ! Initialize state and B-matrix for 3D-Var
  

! ***************************
! ***   Initialize PDAF   ***
! ***************************

  IF (mype_world == 0) THEN
     WRITE (*,'(/1x,a)') 'INITIALIZE PDAF - OFFLINE MODE'
  END IF


! **********************************************************
! ***   CONTROL OF PDAF - used in call to PDAF_init      ***
! **********************************************************

! *** IO options ***
  screen = 2         ! Write screen output (1) for output, (2) add timings

! *** Size of control vector and ensemble size  ***
  dim_ens = 1         ! Size of ensemble for ensemble and hybrid Var
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


! *********************************************************************
! ***   Settings for analysis steps  - used in call-back routines   ***
! *********************************************************************

! *** Which observation type to assimilate
  assim_A = .true.
  assim_B = .false.
  assim_C = .false.

! *** specifications for observations ***
  rms_obs_A = 0.5    ! Observation error standard deviation for observation A
  rms_obs_B = 0.5    ! Observation error standard deviation for observation B
  rms_obs_C = 0.5    ! Observation error standard deviation for observation C

! *** Localization settings
  locweight = 0     ! Type of localizating weighting
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRADIUS
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
  cradius = 0       ! Cut-off radius in grid points for observation domain in local filters
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
! *** Here, only the call for 3D-Var schemes is     ***
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
  filter_param_r(1) = forget         ! Forgetting factor
  filter_param_r(2) = beta_3dvar     ! Hybrid weight for hybrid 3D-Var

  IF (subtype==0) THEN
     ! parameterized 3D-Var
     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 5,&
          filter_param_r, 1, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_3dvar_offline, &
          screen, status_pdaf)
  ELSE
     ! Ensemble or hybrid 3D-Var
     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 5,&
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_offline, &
          screen, status_pdaf)
  END IF


! *** Check whether initialization of PDAF was successful ***
  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in initialization of PDAF - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF


! *************************************
! *** Activate offline mode of PDAF ***
! *************************************

  CALL PDAF_set_offline_mode(screen)

END SUBROUTINE init_pdaf
