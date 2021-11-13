!$Id$
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
! 2019-07 - Lars Nerger - Initial code copied from Lorenz-96 model
! Later revisions - see svn log
!
! !USES:
  USE parser, &           ! Parser function
       ONLY: parse
  USE mod_assimilation, &
       ONLY: screen, filtertype, subtype, dim_ens, delt_obs, &
       model_error, model_err_amp, incremental, &
       type_forget, forget, rank_analysis_enkf, &
       file_ini, type_ensinit, seedset, type_trans, &
       type_sqrt, stepnull_means, dim_lag, &
       twin_experiment, pf_res_type, init_dt, &
       pf_noise_type, pf_noise_amp, type_winf, limit_winf
  USE output_netcdf_asml, &
       ONLY: init_netcdf_asml, file_asml, delt_write_asml, write_states, &
       write_stats, write_ens
  USE obs_gp_pdafomi, &
       ONLY: rms_obs, file_obs, use_obs_mask, file_obs_mask, &
       use_maskfile, numobs, dx_obs, obs_err_type, file_syntobs

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
  handle = 'use_obs_mask'            ! Use a mask for observations with gaps
  CALL parse(handle, use_obs_mask)
  handle = 'use_maskfile'            ! Whether to read observation mask from file
  CALL parse(handle, use_maskfile)
  handle = 'numobs'                  ! Set number of observations to be used
  CALL parse(handle, numobs)         ! Used or observations 1 to numobs
  handle = 'dx_obs'                  ! Set number of observations to be used
  CALL parse(handle, dx_obs)         ! Grid point distance of obs
  handle = 'obs_err_type'            ! Set observation error type
  CALL parse(handle, obs_err_type)
  handle = 'twin_experiment'         ! Activate twin experiment
  CALL parse(handle, twin_experiment)

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

  ! General settings - external of PDAF
  handle = 'type_ensinit'            ! Define type of ensemble initialization
  CALL parse(handle, type_ensinit)   ! possible are 'eof' or 'rnd'
  handle = 'seedset'                 ! Choose set of seeds for init (only for 'rnd')
  CALL parse(handle, seedset)        ! valid values are 1 to 5
  handle = 'init_dt'                 ! Time step interval for random ensemble initialization
  CALL parse(handle, init_dt)        ! valid values are >=1

  ! Filter-specific settings
  handle = 'type_trans'              ! Type of ensemble transformation in (L)ESTKF/(L)SEIK/(L)ETKF
  CALL parse(handle, type_trans)
  handle = 'rank_analysis_enkf'      ! Set rank for pseudo inverse in EnKF
  CALL parse(handle, rank_analysis_enkf)
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

  ! Setting for file output
  handle = 'delt_write_asml'         ! Set write interval for output in assimilation cycles
  CALL parse(handle, delt_write_asml)
  handle = 'write_states'            ! Define, whether state information is written to file
  CALL parse(handle, write_states)
  handle = 'write_stats'             ! Define, whether to write higher-order ensemble statistics
  CALL parse(handle, write_stats)
  handle = 'write_ens'               ! Define, whether to write full ensemble
  CALL parse(handle, write_ens)
  handle = 'stepnull_means '         ! Step at which computation of time mean error is started
  CALL parse(handle, stepnull_means)
  handle = 'file_asml'               ! Set name of output file
  CALL parse(handle, file_asml)
  handle = 'file_ini'                ! Set name of initialization file
  CALL parse(handle, file_ini)
  handle = 'file_obs'                ! Set name of file holding observations
  CALL parse(handle, file_obs)
  handle = 'file_obs_mask'           ! Set name of file for observation mask
  CALL parse(handle, file_obs_mask)
  handle = 'file_syntobs'            ! Set name of file for synthetic observations generated by PDAF
  CALL parse(handle, file_syntobs)

END SUBROUTINE init_pdaf_parse
