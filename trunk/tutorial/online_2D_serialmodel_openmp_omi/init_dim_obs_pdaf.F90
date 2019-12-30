!$Id$
!>  Compute number of observations
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: SEEK/SEIK/EnKF/ETKF/ESTKF/NETF
!!
!! The routine is called at the beginning of each
!! analysis step.  It has to initialize the size of 
!! the observation vector according to the current 
!! time step for the PE-local domain.
!!
!! Implementation for the 2D online example
!! without parallelization.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_dim_obs_pdaf(step, dim_obs_p)

  USE mod_assimilation, &       ! Assimilation variables
       ONLY: obs_p, obs_index_p
  USE mod_model, &              ! Model variables
       ONLY: nx, ny

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step       !< Current time step
  INTEGER, INTENT(out) :: dim_obs_p  !< Dimension of observation vector

! *** Local variables ***
  INTEGER :: i, j                     ! Counters
  INTEGER :: cnt, cnt0                ! Counters
  REAL, ALLOCATABLE :: obs_field(:,:) ! Array for observation field read from file
  CHARACTER(len=2) :: stepstr         ! String for time step


! ****************************************
! *** Initialize observation dimension ***
! ****************************************

  ! Read observation field form file
  ALLOCATE(obs_field(ny, nx))

  IF (step<10) THEN
     WRITE (stepstr, '(i1)') step
  ELSE
     WRITE (stepstr, '(i2)') step
  END IF

  OPEN (12, file='../inputs_online/obs_step'//TRIM(stepstr)//'.txt', status='old')
  DO i = 1, ny
     READ (12, *) obs_field(i, :)
  END DO
  CLOSE (12)

  ! Count observations
  cnt = 0
  DO j = 1, nx
     DO i= 1, ny
        IF (obs_field(i,j) > -999.0) cnt = cnt + 1
     END DO
  END DO

  ! Set number of observations
  dim_obs_p = cnt

  ! Initialize vector of observations and index array
  IF (ALLOCATED(obs_index_p)) DEALLOCATE(obs_index_p)
  IF (ALLOCATED(obs_p)) DEALLOCATE(obs_p)
  ALLOCATE(obs_index_p(dim_obs_p))
  ALLOCATE(obs_p(dim_obs_p))

  cnt = 0
  cnt0 = 0
  DO j = 1, nx
     DO i= 1, ny
        cnt0 = cnt0 + 1
        IF (obs_field(i,j) > -999.0) THEN
           cnt = cnt + 1
           obs_index_p(cnt) = cnt0      ! Index of observation in state vector
           obs_p(cnt) = obs_field(i, j) ! Vector of observations
        END IF
     END DO
  END DO


! *** Clean up ***

  DEALLOCATE(obs_field)

END SUBROUTINE init_dim_obs_pdaf

