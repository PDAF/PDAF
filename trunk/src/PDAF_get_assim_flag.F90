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
! !ROUTINE: PDAF_get_assim_flag --- Query whether assimilation was performed at current time step
!
! !INTERFACE:
SUBROUTINE PDAF_get_assim_flag(did_assim)

! !DESCRIPTION:
! Helper routine for PDAF.
! The routine allows to query whether observations were assimilated 
! at the most recent to to PDAF_assimilate_X.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2018-03 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filter, &
       ONLY: assim_flag

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER,INTENT(out) :: did_assim    ! Flag: (1) for assimilation; (0) else
!EOP

! *** Set ensemble member ***

  did_assim = assim_flag

END SUBROUTINE PDAF_get_assim_flag
