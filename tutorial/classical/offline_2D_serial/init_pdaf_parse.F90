!$Id: init_pdaf_parse.F90 1589 2015-06-12 11:57:58Z lnerger $
!BOP
!
! !ROUTINE: init_pdaf_parse - Parse command line options for PDAF
!
! !INTERFACE:
SUBROUTINE init_pdaf_parse()

! !DESCRIPTION:
! This routine calls the command line parser to initialize
! variables for the data assimilation with PDAF.
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
  USE mod_assimilation, & ! Variables for assimilation
       ONLY: screen, filtertype, subtype, dim_ens, delt_obs, &
       rms_obs, model_error, model_err_amp, incremental, type_forget, &
       forget, rank_ana_enkf, locweight, cradius, &
       sradius, filename, type_trans, dim_obs, &
       type_sqrt

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
  handle = 'dim_obs'                 ! Number of observations
  CALL parse(handle, dim_obs)

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

  ! Filter-specific settings
  handle = 'type_trans'              ! Type of ensemble transformation in SEIK/ETKF/LSEIK/LETKF
  CALL parse(handle, type_trans)
  handle = 'rank_ana_enkf'           ! Set rank for pseudo inverse in EnKF
  CALL parse(handle, rank_ana_enkf)
  handle = 'type_forget'             ! Set type of forgetting factor
  CALL parse(handle, type_forget)
  handle = 'forget'                  ! Set forgetting factor
  CALL parse(handle,forget)
  handle = 'type_sqrt'               ! Set type of transformation square-root (SEIK-sub4, ESTKF)
  CALL parse(handle, type_sqrt)

  ! Settings for localization in LSEIK/LETKF
  handle = 'cradius'                 ! Set cut-off radius in grid points for observation domain
  CALL parse(handle, cradius)
  handle = 'locweight'               ! Set type of localizating weighting
  CALL parse(handle, locweight)
  sradius = cradius                  ! By default use cradius as support radius
  handle = 'sradius'                 ! Set support radius in grid points
             ! for 5th-order polynomial or radius for 1/e in exponential weighting
  CALL parse(handle, sradius)

  ! Setting for file output
  handle = 'filename'                ! Set name of output file
  CALL parse(handle, filename)


END SUBROUTINE init_pdaf_parse
