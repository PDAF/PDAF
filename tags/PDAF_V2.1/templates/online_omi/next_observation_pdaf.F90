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
!! __Revision history:__
!! * 2004-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE next_observation_pdaf(stepnow, nsteps, doexit, time)

  USE mod_assimilation, &        ! Assimilation variables
       ONLY: delt_obs

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

  nsteps = delt_obs

! *********************************
! *** Set current physical time ***
! *********************************

!   time = ??

! *********************
! *** Set exit flag ***
! *********************

!   doexit = ??

END SUBROUTINE next_observation_pdaf
