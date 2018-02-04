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
!$Id: PDAF-D_etkf_memtime.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_etkf_memtime --- Display timing and memory information for ETKF
!
! !INTERFACE:
SUBROUTINE PDAF_etkf_memtime(printtype)

! !DESCRIPTION:
! This routine displays the PDAF-internal timing and
! memory information for the ETKF.
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
       ONLY: subtype_filter, dim_lag

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
     WRITE (*, '(a, 12x, a, F11.3, 1x, a)') &
          'PDAF', 'Generate state ensemble:', pdaf_time_tot(1), 's'
     IF (subtype_filter /= 5) THEN
        WRITE (*, '(a, 18x, a, F11.3, 1x, a)') 'PDAF', 'Time of forecasts:', pdaf_time_tot(2), 's'
     END IF

     ! Filter-specific part
     WRITE (*, '(a, 14x, a, F11.3, 1x, a)') 'PDAF', 'Time for assimilation:', pdaf_time_tot(3), 's'

     ! Generic part B
     WRITE (*, '(a, 16x, a, F11.3, 1x, a)') 'PDAF', 'Time of prepoststep:', pdaf_time_tot(5), 's'

  ELSE IF (printtype == 2) THEN ptype

! *******************************
! *** Print allocated memory  ***
! *******************************

     WRITE (*, '(/a, 23x, a)') 'PDAF', 'PDAF Memory overview'
     WRITE (*, '(a, 10x, 45a)') 'PDAF', ('-', i=1, 45)
     WRITE (*, '(a, 21x, a, f10.3, a)') 'PDAF', 'Allocated memory  (MB)'
     WRITE (*, '(a, 14x, a, f10.5, a)') 'PDAF', &
          'state and A:', pdaf_memcount_get(1, 'M'), ' MB (persistent)'
     WRITE (*, '(a, 11x, a, f10.5, a)') 'PDAF', &
          'ensemble array:', pdaf_memcount_get(2, 'M'), ' MB (persistent)'
     WRITE (*, '(a, 12x, a, f10.5, a)') 'PDAF', &
          'analysis step:', pdaf_memcount_get(3, 'M'), ' MB (temporary)'

  ELSE IF (printtype == 3) THEN ptype

! *********************************************
! *** Print second-level timing information ***
! *********************************************

     ! Generic part
     WRITE (*, '(//a, 21x, a)') 'PDAF', 'PDAF Timing information'
     WRITE (*, '(a, 10x, 45a)') 'PDAF', ('-', i=1, 45)
     WRITE (*, '(a, 13x, a, F11.3, 1x, a)') &
          'PDAF', 'Generate state ensemble (1):', pdaf_time_tot(1), 's'
     IF (subtype_filter /= 5) THEN
        WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'Time of forecasts (2):', pdaf_time_tot(2), 's'
     END IF

     ! Filter-specific part
     WRITE (*, '(a, 15x, a, F11.3, 1x, a)') 'PDAF', 'Time for assimilation (3):', pdaf_time_tot(3), 's'
     WRITE (*, '(a, 23x, a, F11.3, 1x, a)') 'PDAF', 'get mean state (11):', pdaf_time_tot(11), 's'
     WRITE (*, '(a, 11x, a, F11.3, 1x, a)') 'PDAF', 'init observation dimension (15):', pdaf_time_tot(15), 's'
     WRITE (*, '(a, 24x, a, F11.3, 1x, a)') 'PDAF', 'init residual (12):', pdaf_time_tot(12), 's'
     WRITE (*, '(a, 16x, a, F11.3, 1x, a)') 'PDAF', 'compute new inverse U (10):', pdaf_time_tot(10), 's'
     WRITE (*, '(a, 14x, a, F11.3, 1x, a)') 'PDAF', 'get state weight vector (13):', pdaf_time_tot(13), 's'
     IF (subtype_filter /=3) THEN
        WRITE (*, '(a, 13x, a, F11.3, 1x, a)') 'PDAF', 'prepare ensemble weights (20):', pdaf_time_tot(20), 's'
        WRITE (*, '(a, 16x, a, F11.3, 1x, a)') 'PDAF', 'store ensemble matrix (21):', pdaf_time_tot(21), 's'
        WRITE (*, '(a, 22x, a, F11.3, 1x, a)') 'PDAF', 'update ensemble (22):', pdaf_time_tot(22), 's'
        IF (dim_lag >0) &
             WRITE (*, '(a, 20x, a, F11.3, 1x, a)') 'PDAF', 'perform smoothing (17):', pdaf_time_tot(17), 's'
     ELSE
        WRITE (*, '(a, 12x, a, F11.3, 1x, a)') 'PDAF', 'update state and ensemble (14):', pdaf_time_tot(14), 's'
     END IF
     ! Generic part B
     WRITE (*, '(a, 17x, a, F11.3, 1x, a)') 'PDAF', 'Time of prepoststep (5):', pdaf_time_tot(5), 's'

  ELSE IF (printtype == 4) THEN ptype

! *****************************************
! *** Print detailed timing information ***
! *****************************************

     ! Generic part
     WRITE (*, '(//a, 21x, a)') 'PDAF', 'PDAF Timing information'
     WRITE (*, '(a, 10x, 45a)') 'PDAF', ('-', i=1, 45)
     WRITE (*, '(a, 13x, a, F11.3, 1x, a)') &
          'PDAF', 'Generate state ensemble (1):', pdaf_time_tot(1), 's'
     IF (subtype_filter /= 5) THEN
        WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'Time of forecasts (2):', pdaf_time_tot(2), 's'
     END IF

     ! Filter-specific part
     WRITE (*, '(a, 15x, a, F11.3, 1x, a)') 'PDAF', 'Time for assimilation (3):', pdaf_time_tot(3), 's'
     WRITE (*, '(a, 23x, a, F11.3, 1x, a)') 'PDAF', 'get mean state (11):', pdaf_time_tot(11), 's'
     WRITE (*, '(a, 11x, a, F11.3, 1x, a)') 'PDAF', 'init observation dimension (15):', pdaf_time_tot(15), 's'
     WRITE (*, '(a, 24x, a, F11.3, 1x, a)') 'PDAF', 'init residual (12):', pdaf_time_tot(12), 's'
     WRITE (*, '(a, 16x, a, F11.3, 1x, a)') 'PDAF', 'compute new inverse U (10):', pdaf_time_tot(10), 's'
     WRITE (*, '(a, 35x, a, F11.3, 1x, a)') 'PDAF', 'HZ_p (30):', pdaf_time_tot(30), 's'
     WRITE (*, '(a, 26x, a, F11.3, 1x, a)') 'PDAF', 'complete Uinv (31):', pdaf_time_tot(31), 's'
     WRITE (*, '(a, 14x, a, F11.3, 1x, a)') 'PDAF', 'get state weight vector (13):', pdaf_time_tot(13), 's'
     IF (subtype_filter /= 3) THEN
        WRITE (*, '(a, 13x, a, F11.3, 1x, a)') 'PDAF', 'prepare ensemble weights (20):', pdaf_time_tot(20), 's'
        WRITE (*, '(a, 16x, a, F11.3, 1x, a)') 'PDAF', 'store ensemble matrix (21):', pdaf_time_tot(21), 's'
        WRITE (*, '(a, 22x, a, F11.3, 1x, a)') 'PDAF', 'update ensemble (22):', pdaf_time_tot(22), 's'
     IF (dim_lag >0) &
          WRITE (*, '(a, 20x, a, F11.3, 1x, a)') 'PDAF', 'perform smoothing (17):', pdaf_time_tot(17), 's'
     ELSE
        WRITE (*, '(a, 12x, a, F11.3, 1x, a)') 'PDAF', 'update state and ensemble (14):', pdaf_time_tot(14), 's'
     END IF

     ! Generic part
     WRITE (*, '(a, 17x, a, F11.3, 1x, a)') 'PDAF', 'Time of prepoststep (5):', pdaf_time_tot(5), 's'
  END IF ptype


END SUBROUTINE PDAF_etkf_memtime
