C     $Id$
C
C     !ROUTINE: init_n_domains_pdaf --- Set number of local analysis domains
C
C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C     Used in the filters: LSEIK/LETKF/LESTKF
C
C     The routine is called in PDAF\_X\_update
C     at the beginning of the analysis step before
C     the loop through all local analysis domains.
C     It has to set the number of local analysis
C     domains for the PE-local domain.
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !INTERFACE:
      SUBROUTINE init_n_domains_pdaf(step, n_domains_p)
C     !USES:
      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in)  :: step        ! Current time step
      INTEGER, INTENT(out) :: n_domains_p ! PE-local number of analysis domains

C     !CALLING SEQUENCE:
C     Called by: PDAF_lseik_update   (as U_init_n_domains)
C     Called by: PDAF_lestkf_update  (as U_init_n_domains)
C     Called by: PDAF_letkf_update   (as U_init_n_domains)

C     *** local variables ***
#include "SIZE.h"


C     ************************************
C     *** Initialize number of domains ***
C     ************************************

      ! local analysis domains are vertical columns
      n_domains_p = sNx*SNy

      END SUBROUTINE init_n_domains_pdaf
