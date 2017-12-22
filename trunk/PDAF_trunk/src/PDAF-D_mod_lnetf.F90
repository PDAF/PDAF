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
!$Id: PDAF-D_mod_lnetf.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !MODULE:
MODULE PDAF_mod_lnetf
  
! !DESCRIPTION:
! This module provides variables shared between the
! subroutines of of LNETF in PDAF.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2016-06 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:

  USE PDAF_mod_filter   ! Include general filter variables

  IMPLICIT NONE
  SAVE

! !PUBLIC DATA MEMBERS:

  ! *** Filter fields ***
  REAL, ALLOCATABLE :: weights(:)          ! Vector of weights
!EOP

END MODULE PDAF_mod_lnetf
