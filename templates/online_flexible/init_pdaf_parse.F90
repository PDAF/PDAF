!>  Parse command line options for PDAF
!!
!! This routine calls the command line parser to initialize
!! variables for the data assimilation with PDAF.
!!
!! Using the parser is optional and shows one possibility
!! to modify the variables of the compiled program. An 
!! alternative to this might be Fortran namelist files.
!!
!! __Revision history:__
!! * 2011-15 - Lars Nerger - Initial code extracted from init_pdaf
!! * Later revisions - see repository log
!!
SUBROUTINE init_pdaf_parse()

  USE PDAF, &             ! PDAF
       ONLY: PDAF_parse
  USE mod_assimilation, & ! Variables for assimilation
       ONLY: screen, filtertype, subtype, dim_ens, delt_obs, &
       model_error, model_err_amp, type_iau, steps_iau, type_forget, &
       forget, rank_ana_enkf, locweight, cradius, &
       sradius, type_trans, type_sqrt, dim_lag, type_hyb, &
       hyb_gamma, hyb_kappa, type_winf, limit_winf, &
       pf_res_type, pf_noise_type, pf_noise_amp
!  USE obs_OBSTYPE_pdafomi,   ! Variables for observation type OBSTYPE
!       ONLY: ...

  IMPLICIT NONE

! *** Local variables ***
  CHARACTER(len=32) :: handle  ! handle for command line parser


! **********************************
! *** Parse command line options ***
! **********************************

  ! Observation settings - particular for the implemented observation modules
!   handle = 'assim_OBSTYPE'           ! Whether to assimilation observation type OBSTYPE
!   CALL PDAF_parse(handle, assim_OBSTYPE)
!   handle = 'rms_obs_OBSTYPE'         ! Assumed uniform RMS error of observations OBSTYPE
!   CALL PDAF_parse(handle, rms_obs_OBSTYPE)

! The remaining PDAF_parse commands should be generic; usually no change necessary

  ! Observation settings
  handle = 'delt_obs'                ! Time step interval between filter analyses
  CALL PDAF_parse(handle, delt_obs)

  ! Settings for model and time stepping
  handle = 'model_error'             ! Control application of model error
  CALL PDAF_parse(handle, model_error)
  handle = 'model_err_amp'           ! Amplitude of model error
  CALL PDAF_parse(handle, model_err_amp)

  ! General settings for PDAF
  handle = 'screen'                  ! set verbosity of PDAF
  CALL PDAF_parse(handle, screen)
  handle = 'dim_ens'                 ! set ensemble size/rank of covar matrix
  CALL PDAF_parse(handle, dim_ens)
  handle = 'filtertype'              ! Choose filter algorithm
  CALL PDAF_parse(handle, filtertype)
  handle = 'subtype'                 ! Set subtype of filter
  CALL PDAF_parse(handle, subtype)
  handle = 'type_iau'                ! Set whether to use incremental updating
  CALL PDAF_parse(handle, type_iau)
  handle = 'steps_iau'               ! Number of time steps over which IAU is applied
  CALL PDAF_parse(handle, steps_iau)

  ! Settings for smoother
  handle = 'dim_lag'                 ! Size of lag in smoother
  CALL PDAF_parse(handle, dim_lag)

  ! Filter-specific settings
  handle = 'forget'                  ! Set forgetting factor
  CALL PDAF_parse(handle,forget)
  handle = 'type_forget'             ! Set type of forgetting factor
  CALL PDAF_parse(handle, type_forget)
  handle = 'type_trans'              ! Type of ensemble transformation in SEIK/ETKF/ESTKF/LSEIK/LETKF/LESTKF
  CALL PDAF_parse(handle, type_trans)
  handle = 'type_sqrt'               ! Set type of transformation square-root (SEIK-sub4, ESTKF)
  CALL PDAF_parse(handle, type_sqrt)
  handle = 'rank_ana_enkf'           ! Set rank for pseudo inverse in EnKF
  CALL PDAF_parse(handle, rank_ana_enkf)

  ! Settings for localization in LSEIK/LETKF
  handle = 'cradius'                 ! Set cut-off radius in grid points for observation domain
  CALL PDAF_parse(handle, cradius)
  handle = 'locweight'               ! Set type of localizating weighting
  CALL PDAF_parse(handle, locweight)
  sradius = cradius                  ! By default use cradius as support radius
  handle = 'sradius'                 ! Set support radius in grid points
             ! for 5th-order polynomial or radius for 1/e in exponential weighting
  CALL PDAF_parse(handle, sradius)

  ! Settings for nonlinear filters
  handle = 'pf_res_type'             ! Resampling type for particle filter
  CALL PDAF_parse(handle, pf_res_type)        
  handle = 'pf_noise_type'           ! Type of perturbing noise in PF
  CALL PDAF_parse(handle, pf_noise_type)        
  handle = 'pf_noise_amp'            ! Amplitude of perturbing noise in PF
  CALL PDAF_parse(handle, pf_noise_amp)        
  handle = 'type_winf'               ! Set type of weights inflation in NETF/LNETF
  CALL PDAF_parse(handle, type_winf)
  handle = 'limit_winf'              ! Set limit for weights inflation
  CALL PDAF_parse(handle, limit_winf)

  ! Hybrid weights for LKNETF
  handle = 'type_hyb'                ! Set type of hybrid weight
  CALL PDAF_parse(handle, type_hyb)
  handle = 'hyb_gamma'               ! Set hybrid filter weight for state (1.0 LETKF, 0.0 LNETF)
  CALL PDAF_parse(handle, hyb_gamma)
  handle = 'hyb_kappa'               ! Set hybrid norm (>1.0)
  CALL PDAF_parse(handle, hyb_kappa)

END SUBROUTINE init_pdaf_parse
