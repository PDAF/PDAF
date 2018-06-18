C     $Id$
#include "CPP_OPTIONS.h"
C
C     !ROUTINE: obs_op_pdaf --- Implementation of observation operator
C
C     !INTERFACE:
      SUBROUTINE obs_op_pdaf(step, dim_p, dim_obs_p, state_p, m_state_p)

C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C
C     The routine is called during the analysis step.
C     It has to perform the operation of the
C     observation operator acting on a state vector.
C     For domain decomposition, the action is on the
C     PE-local sub-domain of the state and has to 
C     provide the observed sub-state for the PE-local 
C     domain.
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !USES:
      USE mod_assimilation,
     &     ONLY: obs_index

      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: step                     ! Currrent time step
      INTEGER, INTENT(in) :: dim_p                    ! PE-local dimension of state
      INTEGER, INTENT(in) :: dim_obs_p                ! Dimension of observed state
      _RL,     INTENT(in)    :: state_p(dim_p)        ! PE-local model state
      _RL,     INTENT(out)   :: m_state_p(dim_obs_p)  ! PE-local observed state
      
C     !CALLING SEQUENCE:
C     Called by: PDAF_estkf_update  (as U_obs_op)
C     Called by: PDAF_etkf_update   (as U_obs_op)

C     *** local variables ***
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"

      INTEGER :: i 

      
C     *********************************************
C     *** Perform application of measurement    ***
C     *** operator H on vector or matrix column ***
C     *********************************************

      ! This is generic because it is just using
      ! obs_index filled in init_dim_obs_f.

      DO i = 1, dim_obs_p
         m_state_p(i) = state_p(obs_index(i))
      ENDDO

      END
