! Copyright (c) 2004-2020 Lars Nerger
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
! !ROUTINE: PDAF_print_info --- Print information for PDAF (timing and memory) to screen
!
! !INTERFACE:
SUBROUTINE PDAF_print_info(printtype)

! !DESCRIPTION:
! This routine displays the information from PDAF.
! Possible are to display the timing information and
! allocated memory.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2008-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filter, &
       ONLY: filterstr

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(inout) :: printtype    ! Type of screen output:  
                                         ! (1) timings, (2) memory
!EOP


! ********************************
! *** Print screen information ***
! ********************************

  IF (TRIM(filterstr) == 'ESTKF') THEN
     CALL PDAF_estkf_memtime(printtype)
  ELSE IF (TRIM(filterstr) == 'LESTKF') THEN
     CALL PDAF_lestkf_memtime(printtype)
  END IF

END SUBROUTINE PDAF_print_info
