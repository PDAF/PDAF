!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz, Guowei He and Mukund Pondkule (Forschungszentrum Juelich GmbH)
!
!This file is part of TSMP-PDAF
!
!TSMP-PDAF is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TSMP-PDAF is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU LesserGeneral Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public License
!along with TSMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------
!next_observation_pdaf.F90: TSMP-PDAF implementation of routine
!                           'next_observation_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: next_observation_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: next_observation_pdaf --- Initialize information on next observation
!
! !INTERFACE:
SUBROUTINE next_observation_pdaf(stepnow, nsteps, doexit, time)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
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
! The routine is called by all processes. 
!
! !REVISION HISTORY:
! 2013-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: delt_obs, toffset, screen
  USE mod_parallel_pdaf, &
       ONLY: mype_world
  USE mod_tsmp, &
       ONLY: total_steps
  USE mod_tsmp, ONLY: da_interval
  USE mod_tsmp, ONLY: da_interval_final
  USE mod_tsmp, ONLY: flexible_da_interval
  USE mod_assimilation, &
       ONLY: obs_filename
  use mod_read_obs, &
       only: check_n_observationfile
  use mod_read_obs, ONLY: check_n_observationfile_da_interval
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: stepnow  ! Number of the current time step
  INTEGER, INTENT(out) :: nsteps   ! Number of time steps until next obs
  INTEGER, INTENT(out) :: doexit   ! Whether to exit forecasting (1 for exit)
  REAL, INTENT(out)    :: time     ! Current model (physical) time

! !CALLING SEQUENCE:
! Called by: PDAF_get_state   (as U_next_obs)
!EOP

  !kuw: local variables
  integer :: counter
  integer :: no_obs=0
  character (len = 110) :: fn
  !kuw end

  REAL :: da_interval_new
  
  time = 0.0    ! Not used in fully-parallel implementation variant
  doexit = 0

  !kuw: implementation for at least 1 existing observation per observation file
  !!print *, "stepnow", stepnow
  !write(*,*)'stepnow (in next_observation_pdaf):',stepnow
  !nsteps = delt_obs
  !kuw end

  !kuw: check, for observation file with at least 1 observation
!  counter = stepnow 
  counter = stepnow 
  !nsteps  = 0

  if (mype_world==0 .and. screen > 2) then
      write(*,*) 'TSMP-PDAF (in next_observation_pdaf.F90) total_steps: ',total_steps
  end if

  do
    !nsteps  = nsteps  + delt_obs 
    counter = counter + delt_obs

    ! Exit if past last observation file
    !if(counter>total_steps) exit
    if(counter>(total_steps+toffset)) then
      exit
    end if

    ! Check observation file #counter for observations
    write(fn, '(a, i5.5)') trim(obs_filename)//'.', counter
    call check_n_observationfile(fn,no_obs)

    ! Exit loop if observation file contains observations
    if(no_obs>0) then
      exit
    end if

  end do

  ! Set number of steps for PDAF
  nsteps = counter - stepnow

  ! flexible_da_interval should be input (0/1)
  if(flexible_da_interval.eq.1) then

#ifdef PDAF_DEBUG
    ! Error Check: delt_obs must be one for flexible time stepping
    if (delt_obs /= 1) then
      write(*,'(a,i10)') "delt_obs = ", delt_obs
      write(*,'(a)') "delt_obs must be one for flexible time stepping"
      stop "Stopped from incorrect delt_obs"
    end if

    ! Warning: nsteps should be one
    if(nsteps > 1) then
      write(*,'(a,i10)') "WARNING: nsteps = ", nsteps
      write(*,'(a)') "WARNING: nsteps should be one for flexible time stepping"
      write(*,'(a)') "WARNING: Any time differences can be encoded in observation files"
      write(*,'(a)') "WARNING: using the variable da_interval."
    end if
#endif

    ! Initialize da_interval_variable as zero
    da_interval_new = 0.0

    if(counter>(total_steps+toffset)) then
      ! Set da_interval to da_interval_final from EnKF input file
      da_interval_new = da_interval_final
    else
      ! Set da_interval from observation files
      call check_n_observationfile_da_interval(fn,da_interval_new)
    end if

#ifdef PDAF_DEBUG
    ! Error Check: da_interval_new should be set to at least one
    if(da_interval_new < 1.0) then
      write(*,'(a,es22.15)') "da_interval_new = ", da_interval_new
      write(*,'(a)') "da_interval_new is too small, should be minimum of one"
      stop "Stopped from incorrect da_interval_new"
    end if
#endif

    ! Update da_interval
    da_interval = da_interval_new

    if (mype_world==0 .and. screen > 2) then
      write(*,'(a,es22.15)')'TSMP-PDAF (next_observation_pdaf.F90) da_interval: ', da_interval
    end if

  end if

  if (mype_world==0 .and. screen > 2) then
      write(*,*)'TSMP-PDAF (next_observation_pdaf.F90) stepnow: ',stepnow
      write(*,*)'TSMP-PDAF (next_observation_pdaf.F90) no_obs, nsteps, counter: ',no_obs,nsteps,counter
  end if
  !kuw end




!  IF (stepnow + nsteps <= total_steps) THEN
!   if (2<1) then
!    ! *** During the assimilation process ***
!    nsteps = delt_obs   ! This assumes a constant time step interval
!    doexit = 0          ! Do not exit assimilation
!    IF (mype_world == 0) WRITE (*, '(i7, 3x, a, i7)') &
!         stepnow, 'Next observation at time step', stepnow + nsteps
! ELSE
!    ! *** End of assimilation process ***
!    nsteps = 0          ! No more steps
!    doexit = 1          ! Exit assimilation
!    IF (mype_world == 0) WRITE (*, '(i7, 3x, a)') &
!         stepnow, 'No more observations - end assimilation'
! END IF
! *******************************************************
! *** Set number of time steps until next observation ***
! *******************************************************

!   nsteps = ???

! *********************
! *** Set exit flag ***
! *********************

!   doexit = ??
  !print *, "next_observation_pdaf finished"

END SUBROUTINE next_observation_pdaf


