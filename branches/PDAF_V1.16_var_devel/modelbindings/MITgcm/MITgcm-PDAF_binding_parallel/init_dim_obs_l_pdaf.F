C     $Id$
#include "CPP_OPTIONS.h"
C
C     !ROUTINE: init_dim_obs_l_pdaf --- Set dimension of local observation vector
C
C     !INTERFACE:
      SUBROUTINE init_dim_obs_l_pdaf(domain_p, step,
     &      dim_obs_f, dim_obs_l)

C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C     Used in the filters: LSEIK/LETKF/LESTKF
C
C     The routine is called during the loop over
C     all local analysis domains. It has to set 
C     the dimension of the local observation vector 
C     for the current local analysis domain.
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !USES:
      USE mod_assimilation,
     &     ONLY: domain_index, coords_l, obs_coords_f,
     &     obs_index_l, distance_l, local_range
      
      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in)  :: domain_p   ! Current local analysis domain
      INTEGER, INTENT(in)  :: step       ! Current time step
      INTEGER, INTENT(in)  :: dim_obs_f  ! Full dimension of observation vector
      INTEGER, INTENT(out) :: dim_obs_l  ! Local dimension of observation vector

C     !CALLING SEQUENCE:
C     Called by: PDAF_lestkf_update  (as U_init_dim_obs_l)
C     Called by: PDAF_letkf_update   (as U_init_dim_obs_l)
C!EOP


C     *** local variables ***
#include "SIZE.h"
#include "GRID.h"
      INTEGER :: i, cnt
      _RL :: distance
      
C     **********************************************
C     *** Initialize local observation dimension ***
C     **********************************************

C     Here, we do the following
C     1. Determine coordinates of local analysis domain
C     2. Compute distance of each observation from local analysis domain
C        and count observation within localization radius local_range
C     3. Initialize vector of indices of local observations
C        in the full observation vector (used in g2l_obs)
C     4. Initialize vector of distances of local observations
C        (used in prodRinvA_l)
      
C     1. Coordinates of local analysis domain
      coords_l(1) = xC(domain_index(1), domain_index(2), 1, 1)
      coords_l(2) = yC(domain_index(1), domain_index(2), 1, 1)

C     2. Count observations within local_range
      dim_obs_l = 0
      DO i = 1, dim_obs_f
         distance = SQRT((coords_l(1) - obs_coords_f(1,i))**2 +
     &        (coords_l(2) - obs_coords_f(2,i))**2)
         IF (distance <= local_range) dim_obs_l = dim_obs_l + 1
      END DO

C     3. Initialize index array for local observations in full observed vector
C         and array of distances of local observations
      IF (ALLOCATED(obs_index_l)) DEALLOCATE(obs_index_l)
      IF (ALLOCATED(distance_l)) DEALLOCATE(distance_l)
      ALLOCATE(obs_index_l(dim_obs_l))
      ALLOCATE(distance_l(dim_obs_l))

      cnt = 0
      DO i = 1, dim_obs_f
         distance = SQRT((coords_l(1) - obs_coords_f(1,i))**2 +
     &        (coords_l(2) - obs_coords_f(2,i))**2)
         IF (distance <= local_range) THEN
            cnt = cnt + 1
            distance_l(cnt) = distance
            obs_index_l(cnt) = i
         END IF
      END DO

      END SUBROUTINE init_dim_obs_l_pdaf

