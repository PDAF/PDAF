! Copyright (c) 2004-2023 Lars Nerger
!
! This file is part of PDAF.
!
! PDAF is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version
! 3 of the License, or (at your option) any later version.
!
! PDAF is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with PDAF.  If not, see <http://www.gnu.org/licenses/>.
!
!$Id: PDAF_set_comm_pdaf.F90 918 2021-12-03 07:42:19Z lnerger $


!> Set debugging flag
!!
!! This routine set the debug flag for PDAF. 
!! One can set the flag dependent on the local analysis
!! domain, the MPI rank, or the OpenMP thread ID, or
!! and combination of them.
!!
!! For debugval>0 additional information is written by
!! the PDAF routine to stdout. One should activate the 
!! debugging before calling some selected routine(s) and
!! deactivate it with debugval=0 afterwards. This allows 
!! for a targeted checking of the functionality.
!!
!! Note: The debugging of PDAF is independent of that 
!! for PDAF-OMI.
!!
!!  This is a core routine of PDAF and
!!  should not be changed by the user   !
!!
!! __Revision history:__
!! * 2022-07 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE PDAF_set_debug_flag(debugval)

  USE PDAF_mod_filter, &
       ONLY: debug
  USE PDAF_mod_filtermpi, &
       ONLY: mype, filterpe, mype_world

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: debugval          !< Value for debugging flag

! *** Local variables ***
  INTEGER, SAVE :: debug_save=0            ! Store previous debug value to indicate deactivation


! *** Set debugging flag ***

  debug = debugval

  ! Print debug information
  IF (debug>0) THEN
     IF (filterpe) THEN
        WRITE (*,*) '++ PDAF-debug set_debug_flag: mype_filter', mype, 'activate', debug
     ELSE
        WRITE (*,*) '++ PDAF-debug set_debug_flag: mype_world', mype_world, 'activate', debug
     END IF
  ELSE 
     IF (debug_save>0 .AND. debug==0) THEN
        IF (filterpe) THEN
           WRITE (*,*) '++ PDAF-debug set_debug_flag: mype_filter', mype, 'deactivate'
        ELSE
           WRITE (*,*) '++ PDAF-debug set_debug_flag: mype_world', mype_world, 'deactivate'
        END IF
     END IF
  END IF

  ! Save current value of debug
  debug_save = debug

END SUBROUTINE PDAF_set_debug_flag
