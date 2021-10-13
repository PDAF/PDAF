!$Id$
!>  Routine to initialize information on next observation
!!
!! User-supplied call-back routine for PDAF.
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
!! The routine is called by all filter processes. 
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
SUBROUTINE next_observation_pdaf(stepnow, nsteps, doexit, time)

  USE mod_parallel_pdaf, ONLY: mype_filter_echam, task_id
  USE mod_assim_pdaf, ONLY: step_null
  USE mod_assim_atm_pdaf, ONLY: delt_obs_atm, delt_obs_atm_offset, dp

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: stepnow  !< Number of the current time step
  INTEGER, INTENT(out) :: nsteps   !< Number of time steps until next obs
  INTEGER, INTENT(out) :: doexit   !< Whether to exit forecasting (1 for exit)
  REAL(dp), INTENT(out):: time     !< Current model (physical) time


! *************************************************************
! *** Determine number of time steps until next observation ***
! *************************************************************

  IF (stepnow==0) THEN
     nsteps=delt_obs_atm_offset + delt_obs_atm
  ELSE
     nsteps=delt_obs_atm
  END IF

  IF (mype_filter_echam==0 .AND. task_id==1) THEN
     WRITE (*,'(a,i8,a)') 'ECHAM-PDAF: Next observation after ', nsteps ,' time steps'
  ENDIF

  
! *********************************
! *** Set current physical time ***
! *********************************

  time = 1.0

! *********************
! *** Set exit flag ***
! *********************

  doexit = 0

END SUBROUTINE next_observation_pdaf
