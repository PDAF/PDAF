!$Id: next_observation.F90 1590 2015-06-12 11:59:35Z lnerger $
!BOP
!
! !ROUTINE: next_observation --- Initialize information on next observation
!
! !INTERFACE:
SUBROUTINE next_observation(stepnow, nsteps, doexit, time)

! !DESCRIPTION:
! User-supplied routine for PDAF (all filters):
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
! This version is for the Lorenz96 model
! without parallelization.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: delt_obs, mod_time => time, have_obs
  USE mod_model, &
       ONLY: dt, step_final
  USE output_netcdf_asml, &
       ONLY: close_netcdf_asml

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: stepnow  ! Number of the current time step
  INTEGER, INTENT(out) :: nsteps   ! Number of time steps until next obs
  INTEGER, INTENT(out) :: doexit   ! Whether to exit forecasting (1 for exit)
  REAL, INTENT(out)    :: time     ! Current model (physical) time

! !CALLING SEQUENCE:
! Called by: PDAF_get_state   (as U_next_obs)
!EOP


! *************************************************************
! *** Determine number of time steps until next observation ***
! *************************************************************

  IF (stepnow < step_final) THEN
     nsteps = delt_obs   ! This assumes a constant time step interval
  ELSE
     nsteps = 0
  END IF

  ! *** Set current physical time ***
  time = mod_time  

  ! *** Set time index of next observation ***
  ! *** With MOD_TIME we pre-set the time  ***
  ! *** index for use for IO via module    ***
  ! *** MOD_ASSIMILATION such that the     ***
  ! *** forecast time is correct in IO     ***
  mod_time = mod_time + REAL(nsteps) * dt


! *************************************
! *** Dimension of next observation ***
! *************************************

  setexit: IF (stepnow == step_final) THEN
    ! Already at final time step
     WRITE (*, '(i7, 3x, a)') &
          stepnow,'No more observations, exit filtering'
     doexit = 1
     have_obs = .FALSE.

     ! Close NetCDF output file
     CALL close_netcdf_asml()

     ! Close NetCDF file holding true states
     CALL close_netcdf_state()

  ELSE IF (stepnow + nsteps < step_final) THEN setexit
     ! Next observation ahead
     WRITE (*, '(i7, 3x, a, i7)') &
         stepnow, 'Next observation at time step', stepnow + nsteps
     doexit = 0
     have_obs = .TRUE.

  ELSE IF (stepnow + nsteps == step_final) THEN setexit
     ! Final observation ahead
     WRITE (*, '(i7, 3x, a, i7)') &
         stepnow, 'Final observation at time step', stepnow + nsteps
     doexit = 0
     have_obs = .TRUE.

  ELSE IF (stepnow < step_final) THEN setexit
     ! Only forecasting requested
     ! reset time steps and MOD_TIME
     nsteps = step_final - stepnow
     mod_time = mod_time - REAL(nsteps) * dt + REAL(step_final - stepnow) * dt
     doexit = 0
     have_obs = .FALSE.
     WRITE (*, '(i7, 3x, a, i7)') &
         stepnow, 'No more observations, evolve up to time step', stepnow + nsteps
  END IF setexit

END SUBROUTINE next_observation
