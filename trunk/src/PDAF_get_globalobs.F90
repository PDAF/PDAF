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
! !ROUTINE: PDAF_get_globalobs --- Query whether chosen filter needs global observations
SUBROUTINE PDAF_get_globalobs(gobs)

! !DESCRIPTION:
! Routine to return the information whether the current filter uses global observations.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2022-09 - Lars Nerger - Initial code
! Later revisions - see svn log

  USE PDAF_mod_filter, &
       ONLY: globalobs

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(out) :: gobs ! Whether the filter uses global obs.

! !CALLING SEQUENCE:
! Called by: PDAFomi
! Called by: model code
!EOP

  
! ***********************
! *** Set localfilter ***
! ***********************

  gobs = globalobs
  
END SUBROUTINE PDAF_get_globalobs
