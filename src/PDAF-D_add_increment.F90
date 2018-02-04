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
!$Id: PDAF-D_add_increment.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_add_increment --- Add analysis increment to state vector
!
! !INTERFACE:
SUBROUTINE PDAF_add_increment(dim_p, state_p)

! !DESCRIPTION:
! Interface routine for PDAF with incremental analysis updating. 
! The routine adds the full analysis increment to the given state vector .
! This is required, e.g. for the PREPOSTSTEP routines to test the
! analysis state before the next forecast.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2006-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filter, &
       ONLY: state_inc

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER,INTENT(in) :: dim_p          ! State dimension
  REAL,INTENT(inout) :: state_p(dim_p) ! State vector
!EOP

! *** Add increment ***

  state_p = state_p + state_inc

END SUBROUTINE PDAF_add_increment
