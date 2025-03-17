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
!! __Revision history:__
!! * 2004-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE next_observation_pdaf(stepnow, nsteps, doexit, time)

  USE mod_assimilation, &        ! Assimilation variables
       ONLY: delt_obs, mod_time => time
  USE mod_model, &               ! Module provided by model code
       ONLY: dt, step_final

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: stepnow  !< Number of the current time step
  INTEGER, INTENT(out) :: nsteps   !< Number of time steps until next obs
  INTEGER, INTENT(out) :: doexit   !< Whether to exit forecasting (1 for exit)
  REAL, INTENT(out)    :: time     !< Current model (physical) time


! *************************************************************
! *** Determine number of time steps until next observation ***
! *************************************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE next_observation_pdaf.F90: Set number of time steps in forecast!'

  IF (stepnow < step_final) THEN
     nsteps = delt_obs   ! This assumes a constant time step interval
  ELSE
     nsteps = 0
  END IF


! *********************************
! *** Set current physical time ***
! *********************************

  time = mod_time


! *********************
! *** Set exit flag ***
! *********************

  ! Below is the usual logic distinguishing the different cases

  setexit: IF (stepnow == step_final) THEN
    ! Already at final time step
     WRITE (*, '(3x, a,i7, 3x, a)') &
          'PDAFuser', stepnow,'No more observations, exit filtering'
     doexit = 1

  ELSE IF (stepnow + nsteps < step_final) THEN setexit
     ! Next observation ahead
     WRITE (*, '(3x, a,i7, 3x, a, i7)') &
         'PDAFuser', stepnow, 'Next observation at time step', stepnow + nsteps
     doexit = 0

  ELSE IF (stepnow + nsteps == step_final) THEN setexit
     ! Final observation ahead
     WRITE (*, '(3x, a, i7, 3x, a, i7)') &
         'PDAFuser',stepnow, 'Final observation at time step', stepnow + nsteps
     doexit = 0

  ELSE IF (stepnow < step_final) THEN setexit
     ! Only forecasting requested
     ! reset time steps and MOD_TIME
     nsteps = step_final - stepnow
     mod_time = mod_time - REAL(nsteps) * dt + REAL(step_final - stepnow) * dt
     doexit = 0
     WRITE (*, '(3x, a, i7, 3x, a, i7)') &
         'PDAFuser',stepnow, 'No more observations, evolve up to time step', stepnow + nsteps

  END IF setexit

END SUBROUTINE next_observation_pdaf
