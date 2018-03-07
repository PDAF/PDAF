C     $Id$
#include "CPP_OPTIONS.h"
C
C     !ROUTINE: obs_op_f_pdaf --- Implementation of observation operator
C
C     !INTERFACE:
      SUBROUTINE obs_op_f_pdaf(step, dim_p, dim_obs_f,
     &      state_p, m_state_f)

C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C     Used in the filters: LSEIK/LETKF/LESTKF
C
C     The routine is called in PDAF\_X\_update
C     before the loop over all local analysis domains
C     is entered.  The routine has to perform the 
C     operation of the observation operator acting on 
C     a state vector.  The full vector of all 
C     observations required for the localized analysis
C     on the PE-local domain has to be initialized.
C     This is usually data on the PE-local domain plus 
C     some region surrounding the PE-local domain. 
C     This data is gathered by MPI operations. The 
C     gathering has to be done here, since in the loop 
C     through all local analysis domains, no global
C     MPI operations can be performed, because the 
C     number of local analysis domains can vary from 
C     PE to PE.
C
C     !REVISION HISTORY:
C     2013-09 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !USES:
      USE mod_assimilation,
     &     ONLY: obs_index_p, dim_obs_p
      
      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: step                   ! Current time step
      INTEGER, INTENT(in) :: dim_p                  ! PE-local dimension of state
      INTEGER, INTENT(in) :: dim_obs_f              ! Dimension of observed state
      REAL*8, INTENT(in)    :: state_p(dim_p)       ! PE-local model state
      REAL*8, INTENT(inout) :: m_state_f(dim_obs_f) ! PE-local observed state

C     !CALLING SEQUENCE:
C     Called by: PDAF_lestkf_update  (as U_obs_op)
C     Called by: PDAF_letkf_update   (as U_obs_op)

C     *** local variables ***
#include "SIZE.h"
      INTEGER :: i 
      INTEGER :: status        ! Status flag
      _RL, ALLOCATABLE :: m_state_p(:)   ! PE-local observed state


C     *********************************************
C     *** Perform application of measurement    ***
C     *** operator H on vector or matrix column ***
C     *********************************************

      ! This is generic because it is just using
      ! obs_index filled in init_dim_obs_f.
      
      ! Get process-local observed state
      ALLOCATE(m_state_p(dim_obs_p))

      DO i = 1, dim_obs_p
         m_state_p(i) = state_p(obs_index_p(i))
      ENDDO

      ! Get full observed state vector
      CALL PDAF_gather_obs_f(m_state_p, m_state_f, status)

      
C     *****************
C     *** Clean up  ***
C     *****************

      DEALLOCATE(m_state_p)

      END SUBROUTINE obs_op_f_pdaf
