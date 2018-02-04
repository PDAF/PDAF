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
!$Id: PDAF-D_seek_memtime.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_seek_memtime --- Display timing and memory information for SEEK
!
! !INTERFACE:
SUBROUTINE PDAF_seek_memtime(printtype)

! !DESCRIPTION:
! This routine displays the PDAF-internal timing and
! memory information for the SEEK filter.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2008-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_timer, &
       ONLY: PDAF_time_tot
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount_get
  USE PDAF_mod_filter, &
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
     WRITE (*, '(//a, 21x, a)') 'PDAF', 'PDAF Timing information'
     WRITE (*, '(a, 10x, 45a)') 'PDAF', ('-', i=1, 45)
     WRITE (*, '(a, 17x, a, F11.3, 1x, a)') &
             'PDAF', 'EOF initialization:', pdaf_time_tot(1), 's'
     IF (subtype_filter /= 5) THEN
        WRITE (*, '(a, 18x, a, F11.3, 1x, a)') 'PDAF', 'Time of forecasts:', pdaf_time_tot(2), 's'
     END IF

     ! Filter-specific part
     WRITE (*, '(a, 14x, a, F11.3, 1x, a)') 'PDAF', 'Time of assimilations:', pdaf_time_tot(3), 's'
     WRITE (*, '(a, 15x, a, F11.3, 1x, a)') 'PDAF', 're-diagonalize covar:', pdaf_time_tot(4), 's'

     ! Generic part B
     WRITE (*, '(a, 16x, a, F11.3, 1x, a)') 'PDAF', 'Time of prepoststep:', pdaf_time_tot(5), 's'

  ELSE IF (printtype == 2) THEN ptype

! *******************************
! *** Print allocated memory  ***
! *******************************

     WRITE (*, '(//a, 23x, a)') 'PDAF', 'PDAF Memory overview'
     WRITE (*, '(a, 10x, 45a)') 'PDAF', ('-', i=1, 45)
     WRITE (*, '(a, 21x, a, f10.3, a)') 'PDAF', 'Allocated memory  (MB)'
     WRITE (*, '(a, 14x, a, f10.5, a)') &
          'PDAF', 'state and U:', pdaf_memcount_get(1, 'M'), ' MB (persistent)'
     WRITE (*, '(a, 9x, a, f10.5, a)') &
          'PDAF', 'covariance modes:', pdaf_memcount_get(2, 'M'), ' MB (persistent)'
     WRITE (*, '(a, 12x, a, f10.5, a)') &
          'PDAF', 'analysis step:', pdaf_memcount_get(3, 'M'), ' MB (temporary)'
     WRITE (*, '(a, 9x, a, f10.5, a)') &
          'PDAF', 'reinitialization:', pdaf_memcount_get(4, 'M'), ' MB (temporary)'

  ELSE IF (printtype == 3) THEN ptype

! *********************************************
! *** Print second-level timing information ***
! *********************************************

     ! Generic part
     WRITE (*, '(//a, 21x, a)') 'PDAF', 'PDAF Timing information'
     WRITE (*, '(a, 10x, 45a)') 'PDAF', ('-', i=1, 45)
     WRITE (*, '(a, 18x, a, F11.3, 1x, a)') &
             'PDAF', 'EOF initialization (1):', pdaf_time_tot(1), 's'
     IF (subtype_filter /= 5) THEN
        WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'Time of forecasts (2):', pdaf_time_tot(2), 's'
     END IF

     ! Filter-specific part
     WRITE (*, '(a, 15x, a, F11.3, 1x, a)') 'PDAF', 'Time of assimilations (3):', pdaf_time_tot(3), 's'
     WRITE (*, '(a, 25x, a, F11.3, 1x, a)') 'PDAF', 'get residual (12):', pdaf_time_tot(12), 's'
     WRITE (*, '(a, 24x, a, F11.3, 1x, a)') 'PDAF', 'compute new U (10):', pdaf_time_tot(10), 's'
     WRITE (*, '(a, 23x, a, F11.3, 1x, a)') 'PDAF', 'solve for gain (13):', pdaf_time_tot(13), 's'
     WRITE (*, '(a, 25x, a, F11.3, 1x, a)') 'PDAF', 'update state (14):', pdaf_time_tot(14), 's'
     WRITE (*, '(a, 16x, a, F11.3, 1x, a)') 'PDAF', 're-diagonalize covar (4):', pdaf_time_tot(4), 's'
     WRITE (*, '(a, 17x, a, F11.3, 1x, a)') 'PDAF', 'prepare mode weights (20):', pdaf_time_tot(20), 's'
     WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'gather mode matrix (21):', pdaf_time_tot(21), 's'
     WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'update mode matrix (22):', pdaf_time_tot(22), 's'

     ! Generic part B
     WRITE (*, '(a, 17x, a, F11.3, 1x, a)') 'PDAF', 'Time of prepoststep (5):', pdaf_time_tot(5), 's'

  ELSE IF (printtype == 4) THEN ptype

! *****************************************
! *** Print detailed timing information ***
! *****************************************

     ! Generic part
     WRITE (*, '(//a, 21x, a)') 'PDAF', 'PDAF Timing information'
     WRITE (*, '(a, 10x, 45a)') 'PDAF', ('-', i=1, 45)
     WRITE (*, '(a, 18x, a, F11.3, 1x, a)') &
             'PDAF', 'EOF initialization (1):', pdaf_time_tot(1), 's'
     IF (subtype_filter /= 5) THEN
        WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'Time of forecasts (2):', pdaf_time_tot(2), 's'
     END IF

     ! Filter-specific part
     WRITE (*, '(a, 15x, a, F11.3, 1x, a)') 'PDAF', 'Time of assimilations (3):', pdaf_time_tot(3), 's'
     WRITE (*, '(a, 25x, a, F11.3, 1x, a)') 'PDAF', 'get residual (12):', pdaf_time_tot(12), 's'
     WRITE (*, '(a, 24x, a, F11.3, 1x, a)') 'PDAF', 'compute new U (10):', pdaf_time_tot(10), 's'
     WRITE (*, '(a, 34x, a, F11.3, 1x, a)') 'PDAF', 'H V_p (30):', pdaf_time_tot(30), 's'
     WRITE (*, '(a, 26x, a, F11.3, 1x, a)') 'PDAF', 'complete Uinv (31):', pdaf_time_tot(31), 's'
     WRITE (*, '(a, 23x, a, F11.3, 1x, a)') 'PDAF', 'solve for gain (13):', pdaf_time_tot(13), 's'
     WRITE (*, '(a, 25x, a, F11.3, 1x, a)') 'PDAF', 'update state (14):', pdaf_time_tot(14), 's'
     WRITE (*, '(a, 16x, a, F11.3, 1x, a)') 'PDAF', 're-diagonalize covar (4):', pdaf_time_tot(4), 's'
     WRITE (*, '(a, 17x, a, F11.3, 1x, a)') 'PDAF', 'prepare mode weights (20):', pdaf_time_tot(20), 's'
     WRITE (*, '(a, 28x, a, F11.3, 1x, a)') 'PDAF', 'invert Uinv (32):', pdaf_time_tot(32), 's'
     WRITE (*, '(a, 32x, a, F11.3, 1x, a)') 'PDAF', 'SQRT(U) (33):', pdaf_time_tot(33), 's'
     WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'gather mode matrix (21):', pdaf_time_tot(21), 's'
     WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'update mode matrix (22):', pdaf_time_tot(22), 's'

     ! Generic part B
     WRITE (*, '(a, 17x, a, F11.3, 1x, a)') 'PDAF', 'Time of prepoststep (5):', pdaf_time_tot(5), 's'

  END IF ptype


END SUBROUTINE PDAF_seek_memtime
