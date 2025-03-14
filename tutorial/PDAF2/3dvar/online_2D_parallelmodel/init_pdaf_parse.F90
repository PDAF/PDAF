!>  Parse command line options for PDAF
!!
!! This routine calls the command line parser to initialize
!! variables for the data assimilation with PDAF.
!!
!! Using the parser is optional and shows one possibility
!! to modify the variables of the compiled program. An 
!! alternative to this might be Fortran namelist files.
!!
!! This variant is adapted to the 3DVar schemes.
!!
!! __Revision history:__
!! * 2011-15 - Lars Nerger - Initial code extracted from init_pdaf
!! * Later revisions - see repository log
!!
SUBROUTINE init_pdaf_parse()

  USE parser, &           ! Parser function
       ONLY: parse
  USE mod_assimilation, & ! Variables for assimilation
       ONLY: screen, filtertype, subtype, dim_ens, delt_obs, &
       model_error, model_err_amp, incremental, type_forget, &
       forget, rank_ana_enkf, locweight, cradius, &
       sradius, type_trans, type_sqrt, dim_lag, &
       type_opt, mcols_cvec_ens, dim_cvec, beta_3dvar
  USE obs_A_pdafomi, &    ! Variables for observation type A
       ONLY: assim_A, rms_obs_A
  USE obs_B_pdafomi, &    ! Variables for observation type B
       ONLY: assim_B, rms_obs_B

  IMPLICIT NONE

! *** Local variables ***
  CHARACTER(len=32) :: handle  ! handle for command line parser


! **********************************
! *** Parse command line options ***
! **********************************

  ! Observation settings - particular for the implemented observation modules
  handle = 'assim_A'                 ! Whether to assimilation observation type A
  CALL parse(handle, assim_A)
  handle = 'assim_B'                 ! Whether to assimilation observation type B
  CALL parse(handle, assim_B)
  handle = 'rms_obs_A'               ! Assumed uniform RMS error of the observations type A
  CALL parse(handle, rms_obs_A)
  handle = 'rms_obs_B'               ! Assumed uniform RMS error of the observations type B
  CALL parse(handle, rms_obs_B)

! The remaining parse commands should be generic; usually no change necessary

  ! Observation settings
  handle = 'delt_obs'                ! Time step interval between filter analyses
  CALL parse(handle, delt_obs)

  ! Settings for model and time stepping
  handle = 'model_error'             ! Control application of model error
  CALL parse(handle, model_error)
  handle = 'model_err_amp'           ! Amplitude of model error
  CALL parse(handle, model_err_amp)

  ! General settings for PDAF
  handle = 'screen'                  ! set verbosity of PDAF
  CALL parse(handle, screen)
  handle = 'dim_ens'                 ! set ensemble size/rank of covar matrix
  CALL parse(handle, dim_ens)
  handle = 'filtertype'              ! Choose filter algorithm
  CALL parse(handle, filtertype)
  handle = 'subtype'                 ! Set subtype of filter
  CALL parse(handle, subtype)
  handle = 'incremental'             ! Set whether to use incremental updating
  CALL parse(handle, incremental)

  ! Settings for smoother
  handle = 'dim_lag'                 ! Size of lag in smoother
  CALL parse(handle, dim_lag)

  ! Filter-specific settings
  handle = 'forget'                  ! Set forgetting factor
  CALL parse(handle,forget)
  handle = 'type_forget'             ! Set type of forgetting factor
  CALL parse(handle, type_forget)
  handle = 'type_trans'              ! Type of ensemble transformation in SEIK/ETKF/ESTKF/LSEIK/LETKF/LESTKF
  CALL parse(handle, type_trans)
  handle = 'type_sqrt'               ! Set type of transformation square-root (SEIK-sub4, ESTKF)
  CALL parse(handle, type_sqrt)
  handle = 'rank_ana_enkf'           ! Set rank for pseudo inverse in EnKF
  CALL parse(handle, rank_ana_enkf)

  ! Settings for 3D-Var methods
  handle = 'type_opt'                ! Set solver type for 3D-Var
  CALL parse(handle, type_opt)
  dim_cvec = dim_ens
  handle = 'dim_cvec'                ! Size of control vector for 3D-Var
  CALL parse(handle, dim_cvec)
  handle = 'mcols_cvec_ens'          ! multiplication factor for dimension of ensemble control vector
  CALL parse(handle, mcols_cvec_ens)
  handle = 'beta_3dvar'              ! Hybrid weight for hybrid 3D-Var
  CALL parse(handle, beta_3dvar)

  ! Settings for localization in LSEIK/LETKF
  handle = 'cradius'                 ! Set cut-off radius in grid points for observation domain
  CALL parse(handle, cradius)
  handle = 'locweight'               ! Set type of localizating weighting
  CALL parse(handle, locweight)
  sradius = cradius                  ! By default use cradius as support radius
  handle = 'sradius'                 ! Set support radius in grid points
             ! for 5th-order polynomial or distance for 1/e in exponential weighting
  CALL parse(handle, sradius)

END SUBROUTINE init_pdaf_parse
