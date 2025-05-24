!>  Interface routine to call initialization of PDAF
!!
!! This routine collects the initialization of variables for PDAF.
!! In addition, the initialization routine PDAF_init is called
!! to perform the internal initialization of PDAF.
!!
!! This variant is for the online mode of PDAF.
!!
!! This routine is generic. However, it assumes a constant observation
!! error (rms_obs). Further, with parallelization the local state
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
       COMM_model, COMM_filter, COMM_couple, filterpe, abort_parallel, &
       mype_filter
  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: dim_state_p, dim_state, screen, filtertype, subtype, dim_ens, &
       delt_obs, type_iau, steps_iau, &
       type_forget, forget, &
       rank_ana_enkf, locweight, cradius, sradius, &
       type_trans, type_sqrt, &
       observe_ens, type_obs_init, do_omi_obsstats, &
       ensgroup
  USE mod_model, &                ! Model variables
       ONLY: nx, ny, nx_p
  USE obs_A_pdafomi, &            ! Variables for observation type A
       ONLY: assim_A, rms_obs_A
  USE obs_B_pdafomi, &            ! Variables for observation type B
       ONLY: assim_B, rms_obs_B

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: filter_param_i(2) ! Integer parameter array for filter
  REAL    :: filter_param_r(1) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag
  REAL    :: lim_coords(2,2)   ! limiting coordinates of process sub-domain
  INTEGER :: i, off_nx         ! Counters

! *** External subroutines ***
  EXTERNAL :: init_ens_pdaf            ! Ensemble initialization
  EXTERNAL :: next_observation_pdaf, & ! Provide time step, model time, 
                                       ! and dimension of next observation
       distribute_state_pdaf, &        ! Routine to distribute a state vector to model fields
       prepoststep_pdaf                ! User supplied pre/poststep routine
  

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

! *** Ensemble size ***
  dim_ens = n_modeltasks  ! Size of ensemble for all ensemble filters
                     !   We use n_modeltasks here, initialized in init_parallel_pdaf

! *** Options for filter method

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++
  ! +++ For available options see MOD_ASSIMILATION +++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++

  filtertype = 6     ! Type of filter
  subtype = 0        ! Subtype of filter

  forget  = 1.0      ! Forgetting factor value for inflation
  type_forget = 0    ! Type of forgetting factor

  type_trans = 0     ! Type of ensemble transformation (deterministic or random)
  type_sqrt = 0      ! SEIK/LSEIK/ESTKF/LESTKF: Type of transform matrix square-root


! *********************************************************************
! ***   Settings for analysis steps  - used in call-back routines   ***
! *********************************************************************

! *** Type of initial ensemble ***
  ensgroup = 1       ! (1) for ensemble from true state; (2) rotated ensemble by 90 degrees

! *** Forecast length (time interval between analysis steps) ***
  delt_obs = 2       ! This should be set according to the data availability

! *** Which observation type to assimilate
  assim_A = .true.
  assim_B = .false.

! *** specifications for observations ***
  rms_obs_A = 0.5    ! Observation error standard deviation for observation A
  rms_obs_B = 0.5    ! Observation error standard deviation for observation B

! *** Localization settings
  locweight = 0     ! Type of localizating weighting
  cradius = 0.0     ! Cut-off radius in grid points for observation domain in local filters
  sradius = cradius ! Support radius for 5th-order polynomial
                    ! or radius for 1/e for exponential weighting


! ***********************************
! *** Some optional functionality ***
! ***********************************

! *** Parse command line options   ***
! *** This is optional, but useful ***

  call init_pdaf_parse()

! *** Initial Screen output ***
! *** This is optional      ***

  IF (mype_world == 0) call init_pdaf_info()


! *****************************************************
! *** Call PDAF initialization routine on all PEs.  ***
! ***                                               ***
! *** Here, the full selection of filters is        ***
! *** implemented. In a real implementation, one    ***
! *** reduces this to selected filters.             ***
! ***                                               ***
! *** For all filters, PDAF_init is first called    ***
! *** specifying only the required parameters.      ***
! *** Further settings are done afterwards using    ***
! *** calls to PDAF_set_iparam & PDAF_set_rparam.   ***
! *****************************************************

  ! *** Here we specify only the required integer and real parameters
  ! *** Other parameters are set using calls to PDAF_set_iparam/PDAF_set_rparam
  filter_param_i(1) = dim_state_p ! State dimension
  filter_param_i(2) = dim_ens     ! Size of ensemble
  filter_param_r(1) = forget      ! Forgetting factor

  CALL PDAF_init(filtertype, subtype, 0, &
       filter_param_i, 2,&
       filter_param_r, 1, &
       COMM_model, COMM_filter, COMM_couple, &
       task_id, n_modeltasks, filterpe, init_ens_pdaf, &
       screen, status_pdaf)

  ! *** Additional parameter specifications ***
  ! *** -- These are all optional --        ***

  ! Generic settings
  CALL PDAF_set_iparam(5, type_forget, status_pdaf)      ! Type of forgetting factor
  CALL PDAF_set_iparam(6, type_trans, status_pdaf)       ! Type of ensemble transformation
  CALL PDAF_set_iparam(7, type_sqrt, status_pdaf)        ! Type of transform square-root (SEIK-sub4/ESTKF)
  CALL PDAF_set_iparam(8, observe_ens, status_pdaf)      ! Whether to apply observation operator to ensemble mean
  CALL PDAF_set_iparam(9, type_obs_init, status_pdaf)    ! Initialize observation before or after call to prepoststep


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

  CALL PDAF_init_forecast(next_observation_pdaf, distribute_state_pdaf, &
       prepoststep_pdaf, status_pdaf)


! ************************************************************************
! *** Set domain coordinate limits (for use with OMI's use_global_obs) ***
! ************************************************************************
  
    ! Get offset of local domain in global domain in x-direction
    off_nx = 0
    DO i = 1, mype_filter
       off_nx = off_nx + nx_p
    END DO

    lim_coords(1,1) = REAL(off_nx + 1)     ! West
    lim_coords(1,2) = REAL(off_nx + nx_p)  ! East
    lim_coords(2,1) = REAL(ny)             ! North
    lim_coords(2,2) = 1.0                  ! South

    CALL PDAFomi_set_domain_limits(lim_coords)

END SUBROUTINE init_pdaf
