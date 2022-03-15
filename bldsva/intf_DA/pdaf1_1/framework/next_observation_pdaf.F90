!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz, Guowei He and Mukund Pondkule (Forschungszentrum Juelich GmbH)
!
!This file is part of TerrSysMP-PDAF
!
!TerrSysMP-PDAF is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TerrSysMP-PDAF is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU LesserGeneral Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public License
!along with TerrSysMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------
!next_observation_pdaf.F90: TerrSysMP-PDAF implementation of routine
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
! The routine is called by all processes. 
!
! !REVISION HISTORY:
! 2013-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: delt_obs, toffset, screen
  USE mod_parallel_model, &
       ONLY: mype_world, total_steps
  USE mod_assimilation, &
       ONLY: obs_filename
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
  integer :: counter,no_obs=0
  character (len = 110) :: fn
  !kuw end
  
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
      write(*,*) 'PDAF (in next_observation_pdaf.F90) total_steps: ',total_steps
  end if
  do
    !nsteps  = nsteps  + delt_obs 
    counter = counter + delt_obs
    !if(counter>total_steps) exit
    if(counter>(total_steps+toffset)) exit
    write(fn, '(a, i5.5)') trim(obs_filename)//'.', counter
    call check_n_observationfile(fn,no_obs)
    if(no_obs>0) exit
  end do
  nsteps = counter - stepnow
  if (mype_world==0 .and. screen > 2) then
      write(*,*)'PDAF (next_observation_pdaf.F90) stepnow: ',stepnow
      write(*,*)'PDAF (next_observation_pdaf.F90) no_obs, nsteps, counter: ',no_obs,nsteps,counter
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

subroutine check_n_observationfile(fn,nn)
  use netcdf
  implicit none
  character(len=*),intent(in) :: fn
  integer, intent(out)        :: nn

  integer :: ncid, varid, status !,dimid
  !character (len = *), parameter :: dim_name = "dim_obs"
  character (len = *), parameter :: varname = "no_obs"
  character(len = nf90_max_name) :: recorddimname

  call check(nf90_open(fn, nf90_nowrite, ncid))
  !call check(nf90_inq_dimid(ncid, dim_name, dimid))
  !call check(nf90_inquire_dimension(ncid, dimid, recorddimname, nn))
  call check( nf90_inq_varid(ncid, varname, varid) )
  call check( nf90_get_var(ncid, varid, nn) )
  call check(nf90_close(ncid))

end subroutine
