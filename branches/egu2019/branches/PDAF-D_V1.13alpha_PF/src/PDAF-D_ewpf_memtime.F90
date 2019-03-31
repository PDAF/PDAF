! Copyright (c) 2014-2018 Paul Kirchgessner
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
! !ROUTINE: PDAF_ewpf_memtime --- Display timing and memory information for EWPF
!
! !INTERFACE:
SUBROUTINE PDAF_ewpf_memtime(printtype)

! !DESCRIPTION:
! This routine displays the PDAF-internal timing and
! memory information for the ETKF.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2014-05 - Paul Kirchgessner - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_timer, &
       ONLY: PDAF_time_tot
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount_get
  USE PDAF_mod_ewpf, &
       ONLY: subtype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: printtype    ! Type of screen output:  
                                      ! (1) timings, (2) memory
!EOP

! *** Local variables ***
  INTEGER :: i   ! Counter


! ********************************
! *** Print screen information ***
! ********************************
  ptype: IF (printtype == 1) THEN

! **************************************
! *** Print basic timing information ***
! **************************************

     ! Generic part
     WRITE (*, '(//21x, a)') 'PDAF Timing information'
     WRITE (*, '(10x, 45a)') ('-', i=1, 45)
     WRITE (*, '(12x, a, F11.3, 1x, a)') &
          'Generate state ensemble:', pdaf_time_tot(1), 's'
     IF (subtype_filter /= 5) THEN
        WRITE (*, '(18x, a, F11.3, 1x, a)') 'Time of forecasts:', pdaf_time_tot(2), 's'
     END IF

     ! Filter-specific part
     WRITE (*, '(14x, a, F11.3, 1x, a)') 'Time for assimilation:', pdaf_time_tot(3), 's'

     ! Generic part A
     WRITE (*, '(16x, a, F11.3, 1x, a)') 'Time of prepoststep:', pdaf_time_tot(5), 's'

     ! Generic part B
     WRITE (*, '(16x, a, F11.3, 1x, a)') 'Time of proposal step:', pdaf_time_tot(6), 's'

     ! Generic part C
     WRITE (*, '(16x, a, F11.3, 1x, a)') 'Time of put state EWPF:', pdaf_time_tot(7), 's'

     ! Generic part C
     WRITE (*, '(16x, a, F11.3, 1x, a)') 'Time of equal weights:', pdaf_time_tot(8), 's'

  ELSE IF (printtype == 2) THEN ptype

! *******************************
! *** Print allocated memory  ***
! *******************************

     WRITE (*, '(/23x, a)') 'PDAF Memory overview'
     WRITE (*, '(10x, 45a)') ('-', i=1, 45)
     WRITE (*, '(21x, a, f10.3, a)') 'Allocated memory  (MB)'
     WRITE (*, '(14x, a, f10.5, a)') &
          'state and A:', pdaf_memcount_get(1, 'M'), ' MB (persistent)'
     WRITE (*, '(11x, a, f10.5, a)') &
          'ensemble array:', pdaf_memcount_get(2, 'M'), ' MB (persistent)'
     WRITE (*, '(12x, a, f10.5, a)') &
          'analysis step:', pdaf_memcount_get(3, 'M'), ' MB (temporary)'

  ENDIF ptype

END SUBROUTINE PDAF_ewpf_memtime
