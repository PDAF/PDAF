!$Id: next_observation_pdaf.f90 2135 2019-11-22 18:56:29Z lnerger $
!BOP
!
! !ROUTINE: next_observation_pdaf --- Initialize information on next observation
!
! !INTERFACE:
SUBROUTINE next_observation_pdaf(stepnow, nsteps, doexit, time)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! The subroutine is called before each forecast phase
! by PDAF\_get\_state. It has to initialize the number 
! of time steps until the next available observation 
! (nsteps) and the current model time (time). In 
! addition the exit flag (exit) has to be initialized.
! It indicates if the data assimilation process is 
! completed such that the ensemble loop in the model 
! routine can be exited.
!
! The routine is called by all filter processes. 
!
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, ONLY: mype_world
  USE mod_assim_pdaf, ONLY: delt_obs_atm
  USE mo_kind_pdaf, ONLY: dp

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: stepnow  ! Number of the current time step
  INTEGER, INTENT(out) :: nsteps   ! Number of time steps until next obs
  INTEGER, INTENT(out) :: doexit   ! Whether to exit forecasting (1 for exit)
  REAL(dp), INTENT(out):: time     ! Current model (physical) time

! !CALLING SEQUENCE:
! Called by: PDAF_get_state   (as U_next_obs)
!EOP


! *************************************************************
! *** Determine number of time steps until next observation ***
! *************************************************************

  nsteps=delt_obs_atm
  
! *********************************
! *** Set current physical time ***
! *********************************

  time = 1.0

! *********************
! *** Set exit flag ***
! *********************

  doexit = 0

END SUBROUTINE next_observation_pdaf
