C     $Id$
C
C     !ROUTINE: init_dim_l_pdaf --- Set dimension of local model state
C
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
C     !INTERFACE:
      SUBROUTINE init_dim_l_pdaf(step, domain_p, dim_l)

C     !USES:
      USE mod_assimilation,
     &     ONLY: coords_l, id_lstate_in_pstate, n_fields, off_fields, id

      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in)  :: step     ! Current time step
      INTEGER, INTENT(in)  :: domain_p ! Current local analysis domain
      INTEGER, INTENT(out) :: dim_l    ! Local state dimension

C     !CALLING SEQUENCE:
C     Called by: PDAF_lestkf_update  (as U_init_dim_l)
C     Called by: PDAF_letkf_update   (as U_init_dim_l)

#include "CPP_OPTIONS.h"
#include "SIZE.h"
#include "GRID.h"

C     *** local variables ***
      INTEGER :: i, k, cnt       ! Counters
      INTEGER :: domain_index(2) ! i,j indices of local analysis domain


C     ****************************************
C     *** Initialize local state dimension ***
C     ****************************************

      ! Local analysis domain is a vertical column
      ! State vector is one 2D and four 3D fields
      dim_l = 1 + 4*Nr

C     **********************************************
C     *** Initialize coordinates of local domain ***
C     **********************************************

      ! Determine i,j-indices of local analysis domain
      domain_index(2) = ceiling(real(domain_p)/real(sNx))
      domain_index(1) = domain_p - (domain_index(2)-1)*sNx

      ! Determine coordinates of local analysis domain
      coords_l(1) = xC(domain_index(1), domain_index(2), 1, 1)
      coords_l(2) = yC(domain_index(1), domain_index(2), 1, 1)


C     ******************************************************
C     *** Initialize array of indices of the local state ***
C     ***  vector elements in the global state vector.   ***
C     ******************************************************

C     Allocate array
      IF (ALLOCATED(id_lstate_in_pstate))
     &     DEALLOCATE(id_lstate_in_pstate)
      ALLOCATE(id_lstate_in_pstate(dim_l))

C     EtaN is 2D
      id_lstate_in_pstate(1) = domain_p + off_fields(id%etan)

C     3D fields
C     This runs through all fields and excludes etaN.
C     This approach ensures its independence of the id values of etaN
      cnt = 2
      DO i = 1, n_fields
         IF (i /= id%etan) THEN
            DO k = 1, Nr
               id_lstate_in_pstate(cnt)
     &              = domain_p + off_fields(i) + (k-1)*(sNx*sNy)
               cnt = cnt+1
            ENDDO
         ENDIF
      ENDDO

      END SUBROUTINE init_dim_l_pdaf
