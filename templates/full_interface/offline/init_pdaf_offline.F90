!$Id$
!BOP
!
! !ROUTINE: init_pdaf - Interface routine to call initialization of PDAF
!
! !INTERFACE:
SUBROUTINE init_pdaf()

! !DESCRIPTION:
! This routine collects the initialization of variables for PDAF.
! In addition, the initialization routine PDAF_init is called
! such that the internal initialization of PDAF is performed.
! This variant is for the offline mode of PDAF.
!
! This routine is generic. However, it assumes a constant observation
! error (rms_obs). Further, with parallelization the local state
! dimension dim_state_p is used.
!
! !REVISION HISTORY:
! 2008-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF, &                     ! PDAF
       ONLY: PDAF_init
  USE mod_parallel, &             ! Parallelization variables
       ONLY: mype_world, n_modeltasks, task_id, &
       COMM_model, COMM_filter, COMM_couple, filterpe, abort_parallel
  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: dim_state_p, screen, filtertype, subtype, dim_ens, &
       rms_obs, incremental, type_forget, forget, &
       rank_ana_enkf, locweight, cradius, sradius, &
       filename, type_trans, type_sqrt, &
       type_winf, limit_winf, pf_res_type, pf_noise_type, pf_noise_amp, &
       type_hyb, hyb_gamma, hyb_kappa 

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: main
! Calls: init_pdaf_parse
! Calls: init_pdaf_info
! Calls: PDAF_init
!EOP

! Local variables
  INTEGER :: filter_param_i(7) ! Integer parameter array for filter
  REAL    :: filter_param_r(3) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag

  ! External subroutines
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
  screen      = 2    ! Write screen output (1) for output, (2) add timings

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
  incremental = 0    ! SEIK/LSEIK: (1) to perform incremental updating

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

! *** specifications for observations ***
  rms_obs = 0.5    ! Observation error standard deviation

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

! *** File names
  filename = 'output.dat'


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
! *** reduce this to selected filters.              ***
! ***                                               ***
! *** For all filters, first the arrays of integer  ***
! *** and real number parameters are initialized.   ***
! *** Subsequently, PDAF_init is called.            ***
! *****************************************************

  whichinit: IF (filtertype == 2) THEN
     ! *** EnKF with Monte Carlo init ***
     filter_param_i(1) = dim_state_p ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = rank_ana_enkf ! Rank of pseudo-inverse in analysis
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_i(5) = 0           ! Smoother lag (not implemented here)
     filter_param_r(1) = forget      ! Forgetting factor
     
     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 6,&
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_offline, &
          screen, status_pdaf)
  ELSEIF (filtertype == 9) THEN
     ! *** NETF ***
     filter_param_i(1) = dim_state_p ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = 0           ! Size of lag in smoother
     filter_param_i(4) = 0           ! Not used for NETF (Whether to perform incremental analysis)
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_i(7) = type_winf   ! Type of weights inflation
     filter_param_r(1) = forget      ! Forgetting factor
     filter_param_r(2) = limit_winf  ! Limit for weights inflation
     
     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 7, &
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_offline, &
          screen, status_pdaf)
  ELSEIF (filtertype == 10) THEN
     ! *** LNETF ***
     filter_param_i(1) = dim_state_p ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = 0           ! Size of lag in smoother
     filter_param_i(4) = 0           ! Not used for NETF (Whether to perform incremental analysis)
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_i(7) = type_winf   ! Type of weights inflation
     filter_param_r(1) = forget      ! Forgetting factor
     filter_param_r(2) = limit_winf  ! Limit for weights inflation
     
     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 7, &
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_offline, &
          screen, status_pdaf)
  ELSEIF (filtertype == 11) THEN
     ! *** Hybrid filter LKNETF                    ***
     filter_param_i(1) = dim_state_p ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = 0           ! Smoother lag (not implemented here)
     filter_param_i(4) = 0           ! Whether to perform incremental analysis (not implemented for LKNETF)
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_i(7) = type_hyb    ! Type of hybrid weight
     filter_param_r(1) = forget      ! Forgetting factor
     filter_param_r(2) = hyb_gamma   ! Hybrid filter weight for state
     filter_param_r(3) = hyb_kappa   ! Normalization factor for hybrid weight 
     
     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 7,&
          filter_param_r, 3, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_offline, &
          screen, status_pdaf)
  ELSEIF (filtertype == 12) THEN
     ! *** Particle Filter ***
     filter_param_i(1) = dim_state_p     ! State dimension
     filter_param_i(2) = dim_ens       ! Size of ensemble
     filter_param_r(1) = pf_noise_amp  ! Noise amplitude
     ! Optional parameters
     filter_param_i(3) = pf_res_type   ! Resampling type
     filter_param_i(4) = pf_noise_type ! Perturbation type
     filter_param_i(5) = type_forget   ! Type of forgetting factor
     filter_param_i(6) = type_winf     ! Type of weights inflation
     filter_param_r(2) = forget        ! Forgetting factor
     filter_param_r(3) = limit_winf    ! Limit for weights inflation

     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 6, &
          filter_param_r, 3, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_offline, &
          screen, status_pdaf)
  ELSE
     ! *** All other filters                       ***
     ! *** SEIK, LSEIK, ETKF, LETKF, ESTKF, LESTKF ***
     filter_param_i(1) = dim_state_p ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = 0           ! Smoother lag (not implemented here)
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
     filter_param_r(1) = forget      ! Forgetting factor
     
     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 7,&
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_offline, &
          screen, status_pdaf)
  END IF whichinit


! *** Check whether initialization of PDAF was successful ***
  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in initialization of PDAF - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF

END SUBROUTINE init_pdaf
