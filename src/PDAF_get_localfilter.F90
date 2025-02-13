! Copyright (c) 2004-2025 Lars Nerger
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
!
!> Query whther chosen filter is domain-localized
!!
!! Routine to return the information whether the current filter
!! is domain-localized. The valu eof localfilter is set in
!! the initialization routine of a filter.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2020-03 - Lars Nerger - Initial code
!! * Later revisions - see svn log
!!
SUBROUTINE PDAF_get_localfilter(localfilter_out)

  USE PDAF_mod_filter, &
       ONLY: localfilter

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(out) :: localfilter_out  !< Whether the filter is domain-localized

  
! ***********************
! *** Set localfilter ***
! ***********************

  localfilter_out = localfilter
  
END SUBROUTINE PDAF_get_localfilter
