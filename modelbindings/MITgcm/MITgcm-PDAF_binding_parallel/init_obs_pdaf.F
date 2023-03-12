C     $Id$
#include "CPP_OPTIONS.h"
C
C     !ROUTINE: init_obs_pdaf --- Initialize observation vector
C
C     !INTERFACE:
      SUBROUTINE init_obs_pdaf(step, dim_obs_p, observation_p)

C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C     Used in the filters: SEEK/SEIK/EnKF/ETKF/ESTKF
C
C     The routine is called during the analysis step. 
C     It has to provide the PE-local observation vector 
C     for the current time step.
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !USES:
      USE mod_assimilation,
     &     ONLY: obs_p

      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: step                     ! Current time step
      INTEGER, INTENT(in) :: dim_obs_p                ! PE-local dimension of obs. vector
      _RL,  INTENT(out)   :: observation_p(dim_obs_p) ! PE-local observation vector

C     !CALLING SEQUENCE:
C     Called by: PDAF_etkf_analysis, PDAF_etkf_analysis_T
C     Called by: PDAF_estkf_analysis


C     ***************************************************************
C     *** Initialize observation vector for PE-local model domain ***
C     ***************************************************************
  
      ! This is generic because it is just using
      ! the vector obs filled in init_dim_obs_f
      
      observation_p(:) = obs_p(:)

      END SUBROUTINE init_obs_pdaf

