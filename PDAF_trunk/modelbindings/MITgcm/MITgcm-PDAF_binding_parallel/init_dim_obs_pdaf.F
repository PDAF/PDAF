C     $Id$
#include "CPP_OPTIONS.h"
C
C     !ROUTINE: init_dim_obs_pdaf --- Compute number of observations
C
C     !INTERFACE:
      SUBROUTINE init_dim_obs_pdaf(step, dim_obs_p)

C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C     Used in the filters: SEEK/SEIK/EnKF/ETKF/ESTKF
C
C     The routine is called at the beginning of each
C     analysis step.  It has to initialize the size of 
C     the observation vector according to the current 
C     time step for the PE-local domain.
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !USES:
      USE mod_assimilation,
     &     ONLY: obs_p, obs_index_p

      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in)  :: step        ! Current time step
      INTEGER, INTENT(out) :: dim_obs_p   ! Dimension of observation vector

C     !CALLING SEQUENCE:
C     Called by: PDAF_etkf_analysis, PDAF_etkf_analysis_T
C     Called by: PDAF_estkf_analysis, PDAF_estkf_analysis_fixed


C     *** Local variables
#include "SIZE.h"
#include "EEPARAMS.h"
      INTEGER :: i, j, cnt_obs, cnt_idx  ! Counters
      _RL, ALLOCATABLE :: ssh_obs(:,:)   ! 2D array of observations

      
C     ****************************************
C     *** Initialize observation dimension ***
C     ****************************************

C     We show here the usual procedure for the case that
C     there are gaps in the observations.
C     1. Read global observations from a file
C     2. If the observations could have gaps, count valid observations
C        on the process local sub-domain      
C     3. Store values and indices of valid observations
C        into arrays from module
      
C     1. Read global observations from a file      

      ! Here, we just set a value
      allocate(ssh_obs(Nx, Ny))
      ssh_obs = 0.5d0

      
C     2. Count valid observations

      ! Here we assume that observations with values
      ! between -2.0 and 2.0 are valid (actually all in our case)
      cnt_obs = 0
      DO i = myXGlobalLo, myXGlobalLo + sNx - 1
         DO j = myYGlobalLo, myYGlobalLo + sNy - 1
            IF (ssh_obs(i,j)>-2.0d0 .AND. ssh_obs(i,j)<2.0d0) THEN
               cnt_obs = cnt_obs + 1
            ENDIF
         ENDDO
      ENDDO

      ! Set number of observations
      dim_obs_p = cnt_obs

      
C     3. Store process-local observation values and coordinates

      ! Knowing the number of valid observation we can allocate
      ! arrays for the observations and indices of observations
      ! in the state vector.
      ! These arrays are declared in mod_assimilation and
      ! used in routines that are executed later during the
      ! analysis step.

      IF (ALLOCATED(obs_p)) DEALLOCATE(obs_p)
      IF (ALLOCATED(obs_index_p)) DEALLOCATE(obs_index_p)
      ALLOCATE(obs_p(dim_obs_p))          ! Vector of observations
      ALLOCATE(obs_index_p(dim_obs_p))    ! Index of observation in state vector

      cnt_obs = 0
      cnt_idx = 0
      DO i = myXGlobalLo, myXGlobalLo + sNx - 1
         DO j = myYGlobalLo, myYGlobalLo + sNy - 1
            cnt_idx = cnt_idx + 1
            IF (ssh_obs(i,j)>-2.0d0 .AND. ssh_obs(i,j)<2.0d0) THEN
               cnt_obs = cnt_obs + 1
               obs_index_p(cnt_obs) = cnt_idx
               obs_p(cnt_obs) = ssh_obs(i, j)
            ENDIF
         ENDDO
      ENDDO

      
C     *****************
C     *** Clean up  ***
C     *****************
      
      DEALLOCATE(ssh_obs)

      END

