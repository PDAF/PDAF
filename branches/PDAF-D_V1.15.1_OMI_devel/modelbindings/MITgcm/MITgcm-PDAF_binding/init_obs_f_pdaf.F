C     $Id$
#include "CPP_OPTIONS.h"
C
C     !ROUTINE: init_obs_f_pdaf --- Initialize observation vector
C
C     !INTERFACE:
      SUBROUTINE init_obs_f_pdaf(step, dim_obs_f, observation_f)

C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C     Used in the filters: LSEIK/LETKF/LESTKF
C
C     The routine is called in PDAF\_lseik\_update
C     before the loop over all local analysis domains
C     is entered. It has to provide the full observation 
C     vector according to current time step (where 'full' 
C     means: all observations required for the localized 
C     analysis on the PE-local domain).  This routine 
C     is only used for LSEIK if a globally adaptive 
C     forgetting factor is requested, rather than an 
C     individual forgetting factor for each analysis 
C     domain. This routine has to be implemented 
C     consistently with the routines for the full 
C     observation dimension and the full observation 
C     operator. The forgetting factor will only be 
C     globally adaptive, if the full observation vector 
C     is the global observation vector.
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code 
C     Later revisions - see svn log
C
C     !USES:
      USE mod_assimilation,
     &     ONLY: obs

      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: step        ! Current time step
      INTEGER, INTENT(in) :: dim_obs_f   ! Dimension of full observation vector
      _RL, INTENT(out)   :: observation_f(dim_obs_f) ! Full observation vector

C     !CALLING SEQUENCE:
C     Called by: PDAF_lestkf_update  (as U_init_obs)
C     Called by: PDAF_letkf_update   (as U_init_obs)


C     ******************************************
C     *** Initialize full observation vector ***
C     ******************************************

      ! This is generic because it just using
      ! the vector obs filled in init_dim_obs_f
      
      observation_f(:) = obs(:)

      END SUBROUTINE init_obs_f_pdaf

