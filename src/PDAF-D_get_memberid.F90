! Copyright (c) 2004-2016 Lars Nerger
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
!$Id: PDAF-D_get_memberid.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_get_memberid --- Query ensemble index of the current member
!
! !INTERFACE:
SUBROUTINE PDAF_get_memberid(memberid)

! !DESCRIPTION:
! Helper routine for PDAF.
! The routine allows to query the member index of the ensemble
! state that is currently integrated.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2012-03 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filter, &
       ONLY: member

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER,INTENT(inout) :: memberid    ! Index in the local ensemble
!EOP

! *** Set ensemble member ***

  memberid = member

END SUBROUTINE PDAF_get_memberid
