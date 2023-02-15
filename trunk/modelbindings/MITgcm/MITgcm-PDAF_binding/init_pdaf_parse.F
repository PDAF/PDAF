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
      USE mod_assimilation,     ! Variables for assimilation
     &     ONLY: screen, filtertype, subtype, dim_ens,
     &     type_trans, type_sqrt, nsteps_between_obs,
     &     type_forget, forget,
     &     locweight, cradius, sradius, write_ens
      USE obs_etan_pdafomi,     ! Variables for SSH observations
     &     ONLY: assim_etan, rms_obs_etan

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
      handle = 'assim_etan'              ! Whether to assimilation observation type A
      CALL parse(handle, assim_etan)
      handle = 'rms_obs_etan'            ! Assumed uniform RMS error of the observations type A
      CALL parse(handle, rms_obs_etan)

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
      handle = 'cradius'                 ! Set cut-off radius in grid points for observation domain
      CALL parse(handle, cradius)
      handle = 'locweight'               ! Set type of localizating weighting
      CALL parse(handle, locweight)
      sradius = cradius                  ! By default use cradius as support radius
      handle = 'sradius'                 ! Set support radius in grid points
C             ! for 5th-order polynomial or radius for 1/e in exponential weighting
      CALL parse(handle, sradius)

      handle = 'write_ens'               ! Whether to write ensmeble files
      CALL parse(handle, write_ens)

      END SUBROUTINE init_pdaf_parse
