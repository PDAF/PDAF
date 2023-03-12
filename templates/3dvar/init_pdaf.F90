!$Id: init_pdaf.F90 906 2021-12-01 17:26:32Z lnerger $
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
!! error (rms_obs_A, etc.). Further, with parallelization the local state
!! dimension dim_state_p is used.
!!
!! __Revision history:__
!! * 2008-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_pdaf()

  USE pdaf_interfaces_module, &   ! Interface definitions to PDAF core routines
       ONLY: PDAF_init, PDAF_get_state
  USE mod_parallel_pdaf, &        ! Parallelization variables
       ONLY: mype_world, n_modeltasks, task_id, &
       COMM_model, COMM_filter, COMM_couple, filterpe, abort_parallel
  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: dim_state_p, screen, filtertype, subtype, dim_ens, &
       incremental, forget, locweight, local_range, srange, &
       filename, delt_obs, &
       type_opt, dim_cvec, dim_cvec_ens, mcols_cvec_ens, beta_3dvar
  USE obs_OBSTYPE_pdafomi, &      ! Variables for observation OBSTYPE
       ONLY: assim_OBSTYPE, rms_obs_OBSTYPE
!   USE mod_model, &                ! Model variables
!        ONLY: nx, ny

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: filter_param_i(7) ! Integer parameter array for filter
  REAL    :: filter_param_r(3) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag
  INTEGER :: doexit, steps     ! Not used in this implementation
  REAL    :: timenow           ! Not used in this implementation

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

  WRITE (*,*) 'TEMPLATE init_pdaf.F90: Initialize state dimension here!'

  ! *** Define state dimension ***
!  dim_state = ?
  dim_state_p = 10


! **********************************************************
! ***   CONTROL OF PDAF - used in call to PDAF_init      ***
! **********************************************************

! *** IO options ***
  screen      = 2  ! Write screen output (1) for output, (2) add timings

! *** Filter specific variables
  filtertype = 200  ! Type of filter
                    !   (200) 3D-Var schemes
  dim_ens = 9       ! Size of ensemble for all ensemble/hybrid 3D-Var
  subtype = 0       ! subtype of 3D-Var: 
                    !   (0) parameterized 3D-Var
                    !   (1) 3D Ensemble Var using LESTKF for ensemble update
                    !   (4) 3D Ensemble Var using ESTKF for ensemble update
                    !   (6) hybrid 3D-Var using LESTKF for ensemble update
                    !   (7) hybrid 3D-Var using ESTKF for ensemble update
  forget  = 1.0     ! Forgetting factor
  incremental = 0   ! (1) to perform incremental updating (only in SEIK/LSEIK!)
  type_opt = 1      ! Type of minimizer for 3DVar
                    !   (1) LBFGS, (2) CG+, (3) plain CG
                    !   (12) CG+ parallel, (13) plain CG parallel
  dim_cvec = dim_ens  ! dimension of control vector (parameterized part)
  mcols_cvec_ens = 1  ! Multiplication factor for ensemble control vector (to simulate localization)
  beta_3dvar = 0.5  ! Hybrid weight for hybrid 3D-Var


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
                    !   (1) exponentially decreasing with SRANGE
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
  local_range = 2.0     ! Range in grid points for observation domain in local filters
  srange = local_range  ! Support range for 5th-order polynomial
                    ! or range for 1/e for exponential weighting

! *** File names
  filename = 'output.dat'


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
     CALL PDAF_get_state(steps, timenow, doexit, next_observation_pdaf, &
          distribute_state_pdaf, prepoststep_ens_pdaf, status_pdaf)
  ELSE
     ! For parameterized 3D-Var
     CALL PDAF_get_state(steps, timenow, doexit, next_observation_pdaf, &
          distribute_state_pdaf, prepoststep_3dvar_pdaf, status_pdaf)
  END IF

END SUBROUTINE init_pdaf
