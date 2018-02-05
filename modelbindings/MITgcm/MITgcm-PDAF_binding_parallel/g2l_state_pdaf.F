C     $Id$
C
C     !ROUTINE: g2l_state_pdaf --- Restrict a model state to a local analysis domain
C
C     !INTERFACE:
      SUBROUTINE g2l_state_pdaf(step, domain_p,
     &      dim_p, state_p, dim_l, state_l)

C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C     Used in the filters: LSEIK/LETKF/LESTKF
C
C     The routine is called during the loop over all
C     local analysis domains in PDAF\_lseik\_update
C     before the analysis on a single local analysis 
C     domain.  It has to project the full PE-local 
C     model state onto the current local analysis 
C     domain.
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !USES:
      USE mod_assimilation,
     &     ONLY: off_fields_p, nfields
      
      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: step           ! Current time step
      INTEGER, INTENT(in) :: domain_p       ! Current local analysis domain
      INTEGER, INTENT(in) :: dim_p          ! PE-local full state dimension
      INTEGER, INTENT(in) :: dim_l          ! Local state dimension
      REAL*8, INTENT(in)    :: state_p(dim_p) ! PE-local full state vector 
      REAL*8, INTENT(out)   :: state_l(dim_l) ! State vector on local analysis domain

C     !CALLING SEQUENCE:
C     Called by: PDAF_lseik_update    (as U_g2l_state)
C     Called by: PDAF_letkf_update    (as U_g2l_state)
C     Called by: PDAF_lestkf_update   (as U_g2l_state)

C     *** local variables ***
#include "SIZE.h"
      INTEGER :: i, lay, cnt  ! Counters
      
      
C     *************************************
C     *** Initialize local state vector ***
C     *************************************

C     Sea surface height
      state_l(1) = state_p(domain_p)

C     3D fields
      cnt = 2
      DO i = 2, nfields
         DO lay = 1, Nr
            state_l(cnt) = state_p(domain_p
     &           + off_fields_p(i) + (lay-1)*(sNx*sNy))
            cnt = cnt+1
         ENDDO
      ENDDO

      END SUBROUTINE g2l_state_pdaf
