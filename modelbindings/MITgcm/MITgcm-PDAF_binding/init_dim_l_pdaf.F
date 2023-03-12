C     $Id$
C
C     !ROUTINE: init_dim_l_pdaf --- Set dimension of local model state
C
C     !INTERFACE:
      SUBROUTINE init_dim_l_pdaf(step, domain_p, dim_l)

C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C     Used in the filters: LSEIK/LETKF/LESTKF
C
C     The routine is called during analysis step
C     in the loop over all local analysis domain.
C     It has to set the dimension of local model 
C     state on the current analysis domain.
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !USES:
      USE mod_assimilation,
     &     ONLY: domain_index
      
      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in)  :: step     ! Current time step
      INTEGER, INTENT(in)  :: domain_p ! Current local analysis domain
      INTEGER, INTENT(out) :: dim_l    ! Local state dimension

C     !CALLING SEQUENCE:
C     Called by: PDAF_lestkf_update  (as U_init_dim_l)
C     Called by: PDAF_letkf_update   (as U_init_dim_l)

C     *** local variables ***
#include "SIZE.h"


C     ****************************************
C     *** Initialize local state dimension ***
C     ****************************************

      ! Local analysis domain is a vertical column
      ! State vector is one 2D and four 3D fields
      dim_l = 1 + 4*Nr

      ! Determine i,j-indices of local analysis domain
      domain_index(2) = ceiling(real(domain_p)/real(sNx))
      domain_index(1) = domain_p - (domain_index(2)-1)*sNx

      END SUBROUTINE init_dim_l_pdaf
