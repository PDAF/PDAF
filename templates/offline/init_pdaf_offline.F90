!>  Interface routine to call initialization of PDAF
!!
!! This routine collects the initialization of variables for PDAF.
!! In addition, the initialization routine PDAF_init is called
!! to perform the internal initialization of PDAF.
!!
!! This variant is for the offline mode of PDAF.
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
       COMM_model, COMM_filter, COMM_couple, filterpe, abort_parallel
  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: dim_state_p, screen, filtertype, subtype, dim_ens, &
       type_forget, forget, &
       rank_ana_enkf, locweight, cradius, sradius, &
       type_trans, type_sqrt, &
       type_winf, limit_winf, pf_res_type, pf_noise_type, pf_noise_amp, &
       type_hyb, hyb_gamma, hyb_kappa, &
       observe_ens, type_obs_init, do_omi_obsstats
  USE obs_OBSTYPE_pdafomi, &      ! Variables for observation OBSTYPE
       ONLY: assim_OBSTYPE, rms_obs_OBSTYPE

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: filter_param_i(2) ! Integer parameter array for filter
  REAL    :: filter_param_r(1) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag

! *** External subroutines ***
  EXTERNAL :: init_ens_offline  ! Ensemble initialization
  

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

! *** Ensemble size ***
  dim_ens = 9        ! Size of ensemble for all ensemble filters

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

  !EnKF
  rank_ana_enkf = 0  ! EnKF: rank to be considered for inversion of HPH in analysis step

  ! NETF/LNETF/PF
  type_winf = 0      ! NETF/LNETF/PF: Type of weights inflation
  limit_winf = 0.0   ! NETF/LNETF/PF: Limit for weights inflation

  ! LKNETF
  type_hyb = 0       ! LKNETF: Type of hybrid weight
  hyb_gamma =  1.0   ! LKNETF: Hybrid filter weight for state (1.0: LETKF, 0.0: LNETF)
  hyb_kappa = 30.0   ! LKNETF: Hybrid norm for using skewness and kurtosis (type_hyb 3 or 4)

  ! PF
  pf_res_type = 1    ! PF: Resampling type for particle filter
  pf_noise_type = 0  ! PF: Type of pertubing noise
  pf_noise_amp = 0.0 ! PF: Noise amplitude for particle filter


! *********************************************************************
! ***   Settings for analysis steps  - used in call-back routines   ***
! *********************************************************************

! *** Which observation type to assimilate
  assim_OBSTYPE = .true.

! *** specifications for observations ***
  rms_obs_OBSTYPE = 0.5    ! Observation error standard deviation

! *** Localization settings
  locweight = 0     ! Type of localizating weighting
  cradius = 2.0     ! Cut-off radius for observation domain in local filters
  sradius = cradius ! Support radius for 5th-order polynomial
                    ! or radius for 1/e for exponential weighting


! ***********************************
! *** Some optional functionality ***
! ***********************************

! *** Parse command line options   ***
! *** This is optional, but useful ***

  call init_pdaf_parse()

! *** Possibly activate PDAF-OMI observation statistics ***

  IF (do_omi_obsstats) CALL PDAFomi_set_obs_diag(1)

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
       task_id, n_modeltasks, filterpe, init_ens_offline, &
       screen, status_pdaf)

  ! *** Additional parameter specifications ***
  ! *** -- These are all optional --        ***

  ! Generic settings
  CALL PDAF_set_iparam(5, type_forget, status_pdaf)      ! Type of forgetting factor
  CALL PDAF_set_iparam(6, type_trans, status_pdaf)       ! Type of ensemble transformation
  CALL PDAF_set_iparam(7, type_sqrt, status_pdaf)        ! Type of transform square-root (SEIK-sub4/ESTKF)
  CALL PDAF_set_iparam(8, observe_ens, status_pdaf)      ! Whether to apply observation operator to ensemble mean
  CALL PDAF_set_iparam(9, type_obs_init, status_pdaf)    ! Initialize observation before or after call to prepoststep

  ! Setting for EnKF and LEnKF
  IF (filtertype==PDAF_DA_ENKF .OR. filtertype==PDAF_DA_LENKF) THEN
     CALL PDAF_set_iparam(4, rank_ana_enkf, status_pdaf) ! Rank of EVD in EnKF (0 for direct inversion)
  END IF

  ! Settings for NETF and LNETF
  IF (filtertype==PDAF_DA_NETF .OR. filtertype==PDAF_DA_LNETF) THEN
     CALL PDAF_set_iparam(4, pf_noise_type, status_pdaf) ! Perturbation type
     CALL PDAF_set_iparam(7, type_winf, status_pdaf)     ! Type of weights inflation
     CALL PDAF_set_rparam(2, limit_winf, status_pdaf)    ! Limit for weights inflation
     CALL PDAF_set_rparam(3, pf_noise_amp, status_pdaf)  ! Noise amplitude
  END IF

  ! Settings for particle filter PF
  IF (filtertype==PDAF_DA_PF) THEN
     CALL PDAF_set_iparam(4, pf_noise_type, status_pdaf) ! Perturbation type
     CALL PDAF_set_iparam(6, pf_res_type, status_pdaf)   ! Resampling type
     CALL PDAF_set_iparam(7, type_winf, status_pdaf)     ! Type of weights inflation
     CALL PDAF_set_rparam(2, limit_winf, status_pdaf)    ! Limit for weights inflation
     CALL PDAF_set_rparam(3, pf_noise_amp, status_pdaf)  ! Noise amplitude
  END IF

  ! Settings for hybrid filter LKNETF
  IF (filtertype==PDAF_DA_LKNETF) THEN
     CALL PDAF_set_iparam(4, type_hyb, status_pdaf)      ! Choice of hybrid rule
     CALL PDAF_set_rparam(2, hyb_gamma, status_pdaf)     ! Hybrid filter weight for state
     CALL PDAF_set_rparam(3, hyb_kappa, status_pdaf)     ! Normalization factor for hybrid weight 
  END IF


! *** Check whether initialization of PDAF was successful ***
  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in initialization of PDAF - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF

END SUBROUTINE init_pdaf
