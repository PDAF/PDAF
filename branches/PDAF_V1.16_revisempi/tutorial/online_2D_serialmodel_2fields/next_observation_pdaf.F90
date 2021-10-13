!$Id$
!>  Initialize information on next observation
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all filters
!!
!! The subroutine is called before each forecast phase
!! by PDAF_get_state. It has to initialize the number 
!! of time steps until the next available observation 
!! (nsteps) and the current model time (time). In 
!! addition the exit flag (exit) has to be initialized.
!! It indicates if the data assimilation process is 
!! completed such that the ensemble loop in the model 
!! routine can be exited.
!!
!! The routine is called from PDAF_get_state by all processes
!!         
!! Version for the 2D tutorial model.
!!
!! __Revision history:__
!! * 2004-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE next_observation_pdaf(stepnow, nsteps, doexit, time)

  USE mod_assimilation, &     ! Assimilation variables
       ONLY: delt_obs
  USE mod_parallel_pdaf, &    ! Parallelization variables
       ONLY: mype_world
  USE mod_model, &            ! Model variables
       ONLY: total_steps

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: stepnow  !< Number of the current time step
  INTEGER, INTENT(out) :: nsteps   !< Number of time steps until next obs
  INTEGER, INTENT(out) :: doexit   !< Whether to exit forecasting (1 for exit)
  REAL, INTENT(out)    :: time     !< Current model (physical) time


! *******************************************************
! *** Set number of time steps until next observation ***
! *******************************************************

  time = 0.0          ! Not used in this implementation

  IF (stepnow + nsteps <= total_steps) THEN
     ! *** During the assimilation process ***
     nsteps = delt_obs   ! This assumes a constant time step interval
     doexit = 0          ! Do not exit assimilation

     IF (mype_world == 0) WRITE (*, '(i7, 3x, a, i7)') &
          stepnow, 'Next observation at time step', stepnow + nsteps
  ELSE
     ! *** End of assimilation process ***
     nsteps = 0          ! No more steps
     doexit = 1          ! Exit assimilation

     IF (mype_world == 0) WRITE (*, '(i7, 3x, a)') &
          stepnow, 'No more observations - end assimilation'
  END IF

END SUBROUTINE next_observation_pdaf
