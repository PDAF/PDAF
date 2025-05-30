!$Id: init_pdaf_parse.F90 990 2022-03-24 16:21:30Z lnerger $
!BOP
!
! !ROUTINE: init_pdaf_parse - Parse command line options for PDAF
!
! !INTERFACE:
SUBROUTINE init_pdaf_parse()

! !DESCRIPTION:
! This routine calls the command line parser to initialize
! variables for the data assimilation with PDAF.
!
! Using the parser is optional and shows one possibility
! to modify the variables of the compiled program. An 
! alternative to this might be Fortran namelist files.
!
! !REVISION HISTORY:
! 2011-05 - Lars Nerger - Initial code extracted from init_pdaf
! Later revisions - see svn log
!
! !USES:
  USE parser, &           ! Parser function
       ONLY: parse
  USE mod_assimilation, &
       ONLY: screen, filtertype, subtype, dim_ens, delt_obs, &
       model_error, model_err_amp, type_forget, forget, &
       rank_ana_enkf, rms_obs, locweight, cradius, &
       sradius, type_trans, type_sqrt, dim_lag, type_hyb, &
       hyb_gamma, hyb_kappa, type_winf, limit_winf, &
       pf_res_type, pf_noise_type, pf_noise_amp

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: init_pdaf
! Calls: parse
!EOP

! Local variables
  CHARACTER(len=32) :: handle  ! handle for command line parser


! **********************************
! *** Parse command line options ***
! **********************************

  ! Settings for model and time stepping
  handle = 'model_error'             ! Control application of model error
  CALL parse(handle, model_error)
  handle = 'model_err_amp'           ! Amplitude of model error
  CALL parse(handle, model_err_amp)

  ! Observation settings
  handle = 'delt_obs'                ! Time step interval between filter analyses
  CALL parse(handle, delt_obs)
  handle = 'rms_obs'                 ! Assumed uniform RMS error of the observations
  CALL parse(handle, rms_obs)

  ! General settings for PDAF
  handle = 'screen'                  ! set verbosity of PDAF
  CALL parse(handle, screen)
  handle = 'dim_ens'                 ! set ensemble size/rank of covar matrix
  CALL parse(handle, dim_ens)
  handle = 'filtertype'              ! Choose filter algorithm
  CALL parse(handle, filtertype)
  handle = 'subtype'                 ! Set subtype of filter
  CALL parse(handle, subtype)

  ! Settings for smoother
  handle = 'dim_lag'                 ! Size of lag in smoother
  CALL parse(handle, dim_lag)

  ! Filter-specific settings
  handle = 'type_trans'              ! Type of ensemble transformation in SEIK/ETKF/LSEIK/LETKF
  CALL parse(handle, type_trans)
  handle = 'rank_ana_enkf'      ! Set rank for pseudo inverse in EnKF
  CALL parse(handle, rank_ana_enkf)
  handle = 'type_forget'             ! Set type of forgetting factor
  CALL parse(handle, type_forget)
  handle = 'forget'                  ! Set forgetting factor
  CALL parse(handle,forget)
  handle = 'type_sqrt'               ! Set type of transform square-root
  CALL parse(handle,type_sqrt)
  handle = 'pf_res_type'             ! Resampling type for particle filter
  CALL parse(handle, pf_res_type)        
  handle = 'pf_noise_type'           ! Type of perturbing noise in PF
  CALL parse(handle, pf_noise_type)        
  handle = 'pf_noise_amp'            ! Amplitude of perturbing noise in PF
  CALL parse(handle, pf_noise_amp)        
  handle = 'type_winf'               ! Set type of weights inflation in NETF/LNETF
  CALL parse(handle, type_winf)
  handle = 'limit_winf'              ! Set limit for weights inflation
  CALL parse(handle, limit_winf)

  ! Settings for localization in LSEIK/LETKF
  handle = 'cradius'                 ! Set cut-off radius in grid points for observation domain
  CALL parse(handle, cradius)
  handle = 'locweight'               ! Set type of localizating weighting
  CALL parse(handle, locweight)
  sradius = cradius                  ! By default use cradius as support radius
  handle = 'sradius'                 ! Set support radius in grid points
             ! for 5th-order polynomial or radius for 1/e in exponential weighting
  CALL parse(handle, sradius)

  ! Hybrid weights for LKNETF
  handle = 'type_hyb'                ! Set type of hybrid weight
  CALL parse(handle, type_hyb)
  handle = 'hyb_gamma'               ! Set hybrid filter weight for state (1.0 LETKF, 0.0 LNETF)
  CALL parse(handle, hyb_gamma)
  handle = 'hyb_kappa'               ! Set hybrid norm (>1.0)
  CALL parse(handle, hyb_kappa)

  ! Setting for file output

END SUBROUTINE init_pdaf_parse
