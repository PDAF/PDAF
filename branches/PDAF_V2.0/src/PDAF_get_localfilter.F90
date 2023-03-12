! Copyright (c) 2004-2021 Lars Nerger
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
!$Id$
!BOP
!
! !ROUTINE: PDAF_get_localfilter --- Query whther chosen filter is domain-localized
SUBROUTINE PDAF_get_localfilter(lfilter)

! !DESCRIPTION:
! Routine to return the information whther the current filter is domain-localized.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2020-03 - Lars Nerger - Initial code
! Later revisions - see svn log

  USE PDAF_mod_filter, &
       ONLY: localfilter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(out) :: lfilter  ! Whether the filter is domain-localized

! !CALLING SEQUENCE:
! Called by: model code
!EOP

  
! ***********************
! *** Set localfilter ***
! ***********************

  lfilter = localfilter
  
END SUBROUTINE PDAF_get_localfilter
