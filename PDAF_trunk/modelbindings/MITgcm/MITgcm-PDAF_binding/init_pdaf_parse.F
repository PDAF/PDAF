C     $Id$
C!BOP
C
C     !ROUTINE: init_pdaf_parse - Parse command line options for PDAF
C
C     !INTERFACE:
      SUBROUTINE init_pdaf_parse()

C     !DESCRIPTION:
C     This routine calls the command line parser to initialize
C     variables for the data assimilation with PDAF.
C     Using the parser is optional and shows one possibility
C     to modify the variables of the compiled program. An 
C     alternative to this might be Fortran namelist files.
C
C     !REVISION HISTORY:
C     2011-05 - Lars Nerger - Initial code extracted from init_pdaf
C     Later revisions - see svn log
C
C     !USES:
      USE parser,               ! Parser function
     &     ONLY: parse
      USE mod_parallel_model,   ! Parallelization variables
     &     ONLY: mype_world
      USE mod_assimilation,     ! Variables for assimilation       
     &     ONLY: screen, filtertype, subtype, dim_ens,
     &     type_trans, type_sqrt, nsteps_between_obs, rms_obs,
     &     type_forget, forget, 
     &     locweight, local_range, srange, write_ens

      IMPLICIT NONE

C     !CALLING SEQUENCE:
C     Called by: init_pdaf
C     Calls: parse
C!EOP

C     Local variables
      CHARACTER(len=32) :: handle ! handle for command line parser


C     **********************************
C     *** Parse command line options ***
C     **********************************

C     Observation settings
      handle = 'nsteps_between_obs'      ! Time step interval between filter analyses
      CALL parse(handle, nsteps_between_obs)
      handle = 'rms_obs'                 ! Assumed uniform RMS error of the observations
      CALL parse(handle, rms_obs)

C     General settings for PDAF
      handle = 'screen'                  ! set verbosity of PDAF
      CALL parse(handle, screen)
      handle = 'dim_ens'                 ! set ensemble size/rank of covar matrix
      CALL parse(handle, dim_ens)
      handle = 'filtertype'              ! Choose filter algorithm
      CALL parse(handle, filtertype)
      handle = 'subtype'                 ! Set subtype of filter
      CALL parse(handle, subtype)

C     Filter-specific settings
      handle = 'type_trans'              ! Type of ensemble transformation in SEIK/ETKF/LSEIK/LETKF
      CALL parse(handle, type_trans)
      handle = 'type_forget'             ! Set type of forgetting factor
      CALL parse(handle, type_forget)
      handle = 'forget'                  ! Set forgetting factor
      CALL parse(handle,forget)
      handle = 'type_sqrt'               ! Set type of transformation square-root (SEIK-sub4, ESTKF)
      CALL parse(handle, type_sqrt)

C     Settings for localization in LSEIK/LETKF
      handle = 'local_range'             ! Set range in grid points for observation domain
      CALL parse(handle, local_range)
      handle = 'locweight'               ! Set type of localizating weighting
      CALL parse(handle, locweight)
      srange = local_range               ! By default use local_range as support range
      handle = 'srange'                  ! Set support range in grid points
             ! for 5th-order polynomial or range for 1/e in exponential weighting
      CALL parse(handle, srange)

      handle = 'write_ens'               ! Whether to write ensmeble files
      CALL parse(handle, write_ens)   

      
      END SUBROUTINE init_pdaf_parse
