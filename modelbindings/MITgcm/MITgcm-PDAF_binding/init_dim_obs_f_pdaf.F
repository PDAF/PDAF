C     $Id$
#include "CPP_OPTIONS.h"
C
C     !ROUTINE: init_dim_obs_f_pdaf --- Set full dimension of observations
C
C     !INTERFACE:
      SUBROUTINE init_dim_obs_f_pdaf(step, dim_obs_f)

C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C     Used in the filters: LSEIK/LETKF/LESTKF
C
C     The routine is called in PDAF\_lseik\_update 
C     at the beginning of the analysis step before 
C     the loop through all local analysis domains. 
C     It has to determine the dimension of the 
C     observation vector according to the current 
C     time step for all observations required for 
C     the analyses in the loop over all local 
C     analysis domains on the PE-local state domain.
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code
C     Later revisions - see svn log
c
C     !USES:
      USE mod_assimilation,
     &     ONLY: obs, obs_index, obs_coords
      
      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in)  :: step      ! Current time step
      INTEGER, INTENT(out) :: dim_obs_f ! Dimension of full observation vector

C     !CALLING SEQUENCE:
C     Called by: PDAF_lseik_update   (as U_init_dim_obs)
C     Called by: PDAF_lestkf_update  (as U_init_dim_obs)
C     Called by: PDAF_letkf_update   (as U_init_dim_obs)

C     *** local variables ***
#include "SIZE.h"
#include "GRID.h"
      INTEGER :: i, j, cnt_obs, cnt_idx  ! Counters
      _RL, ALLOCATABLE :: ssh_obs(:,:)   ! 2D array of observations

      
C     *********************************************
C     *** Initialize full observation dimension ***
C     *********************************************

C     We show here the usual procedure for the case that
C     there are gaps in the observations.
C     1. Read observations from a file
C     2. If the observations could have gaps, count valid observations
C     3. Store values and coordinates of valid observations
C        into arrays from module
      
C     1. Read observations from a file      

      ! Here, we just set a value
      allocate(ssh_obs(sNx, sNy))
      ssh_obs = 0.5d0

      
C     2. Count valid observations

      ! Here we assume that observations with values
      ! between -2.0 and 2.0 are valid (actually all in our case)

      cnt_obs = 0
      DO j = 1, sNy
         DO i = 1, sNx
            IF (ssh_obs(i,j)>-2.0d0 .AND. ssh_obs(i,j)<2.0d0) THEN
               cnt_obs = cnt_obs + 1
            ENDIF
         ENDDO
      ENDDO

      ! Set number of observations
      dim_obs_f = cnt_obs

      
C     3. Store observation values and coordinates

      ! Knowing the number of valid observation we can allocate
      ! arrays for the observations, indices of observations
      ! in the state vector and observation coordinates.
      ! These arrays are declared in mod_assimilation and
      ! used in routines that are executed later during the
      ! analysis step

      IF (ALLOCATED(obs)) DEALLOCATE(obs)
      IF (ALLOCATED(obs_index)) DEALLOCATE(obs_index)
      IF (ALLOCATED(obs_coords)) DEALLOCATE(obs_coords)
      ALLOCATE(obs(dim_obs_f))
      ALLOCATE(obs_index(dim_obs_f))
      ALLOCATE(obs_coords(2, dim_obs_f))

      cnt_obs = 0
      cnt_idx = 0
      DO j = 1, sNy
         DO i = 1, sNx
            cnt_idx = cnt_idx + 1
            IF (ssh_obs(i,j)>-2.0d0 .AND. ssh_obs(i,j)<2.0d0) THEN
               cnt_obs = cnt_obs + 1
               obs_index(cnt_obs) = cnt_idx
               obs(cnt_obs) = ssh_obs(i, j)
               obs_coords(1, cnt_obs) = xC(i, j, 1, 1)
               obs_coords(2, cnt_obs) = yC(i, j, 1, 1)
            ENDIF
         ENDDO
      ENDDO

      
C     *****************
C     *** Clean up  ***
C     *****************
      
      DEALLOCATE(ssh_obs)
      
      END SUBROUTINE init_dim_obs_f_pdaf

