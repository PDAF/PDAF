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
!$Id: PDAF-D_enkf_memtime.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_enkf_memtime --- Display timing and memory information for EnKF
!
! !INTERFACE:
SUBROUTINE PDAF_enkf_memtime(printtype)

! !DESCRIPTION:
! This routine displays the PDAF-internal timing and
! memory information for the EnKF.
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
       ONLY: subtype_filter, rank_ana_enkf

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
     WRITE (*, '(a, 14x, a, F11.3, 1x, a)') 'PDAF', 'Time of assimilations:', pdaf_time_tot(3), 's'

     ! Generic part
     WRITE (*, '(a, 16x, a, F11.3, 1x, a)') 'PDAF', 'Time of prepoststep:', pdaf_time_tot(5), 's'

  ELSE IF (printtype == 2) THEN ptype

! *******************************
! *** Print allocated memory  ***
! *******************************

     WRITE (*, '(/a, 23x, a)') 'PDAF', 'PDAF Memory overview'
     WRITE (*, '(a, 10x, 45a)') 'PDAF', ('-', i=1, 45)
     WRITE (*, '(a, 21x, a, f10.3, a)') 'PDAF', 'Allocated memory  (MB)'
     WRITE (*, '(a, 20x, a, f10.5, a)') &
          'PDAF', 'state:', pdaf_memcount_get(1, 'M'), ' MB (persistent)'
     WRITE (*, '(a, 11x, a, f10.5, a)') &
          'PDAF', 'ensemble array:', pdaf_memcount_get(2, 'M'), ' MB (persistent)'
     WRITE (*, '(a, 12x, a, f10.5, a)') &
          'PDAF', 'analysis step:', pdaf_memcount_get(3, 'M'), ' MB (temporary)'

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
     WRITE (*, '(a, 15x, a, F11.3, 1x, a)') 'PDAF', 'Time of assimilations (3):', pdaf_time_tot(3), 's'
     WRITE (*, '(a, 23x, a, F11.3, 1x, a)') 'PDAF', 'get mean state (11):', pdaf_time_tot(11), 's'
     IF (subtype_filter == 1) THEN
        WRITE (*, '(a, 17x, a, F11.3, 1x, a)') 'PDAF', 'compute HPH+R and HP (10):', pdaf_time_tot(10), 's'
     ELSE
        WRITE (*, '(a, 24x, a, F11.3, 1x, a)') 'PDAF', 'compute HPH+R (10):', pdaf_time_tot(10), 's'
     END IF
     WRITE (*, '(a, 15x, a, F11.3, 1x, a)') 'PDAF', 'generate obs. ensemble (15):', pdaf_time_tot(15), 's'
     WRITE (*, '(a, 24x, a, F11.3, 1x, a)') 'PDAF', 'init residual (12):', pdaf_time_tot(12), 's'
     WRITE (*, '(a, 14x, a, F11.3, 1x, a)') 'PDAF', 'ensemble transformation (14):', pdaf_time_tot(14), 's'

     ! Generic part
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
     WRITE (*, '(a, 15x, a, F11.3, 1x, a)') 'PDAF', 'Time of assimilations (3):', pdaf_time_tot(3), 's'
     WRITE (*, '(a, 23x, a, F11.3, 1x, a)') 'PDAF', 'get mean state (11):', pdaf_time_tot(11), 's'
     IF (subtype_filter == 1) THEN
        WRITE (*, '(a, 17x, a, F11.3, 1x, a)') 'PDAF', 'compute HPH+R and HP (10):', pdaf_time_tot(10), 's'
     ELSE
        WRITE (*, '(a, 24x, a, F11.3, 1x, a)') 'PDAF', 'compute HPH+R (10):', pdaf_time_tot(10), 's'
     END IF
     WRITE (*, '(a, 31x, a, F11.3, 1x, a)') 'PDAF', 'HXmean_p (33):', pdaf_time_tot(33), 's'
     WRITE (*, '(a, 33x, a, F11.3, 1x, a)') 'PDAF', 'HXpert (30):', pdaf_time_tot(30), 's'
     IF (subtype_filter == 1) THEN
        WRITE (*, '(a, 35x, a, F11.3, 1x, a)') 'PDAF', 'HP_p (31):', pdaf_time_tot(31), 's'
     END IF
     WRITE (*, '(a, 36x, a, F11.3, 1x, a)') 'PDAF', 'HPH (32):', pdaf_time_tot(32), 's'
     WRITE (*, '(a, 15x, a, F11.3, 1x, a)') 'PDAF', 'generate obs. ensemble (15):', pdaf_time_tot(15), 's'
     WRITE (*, '(a, 24x, a, F11.3, 1x, a)') 'PDAF', 'init residual (12):', pdaf_time_tot(12), 's'
     WRITE (*, '(a, 14x, a, F11.3, 1x, a)') 'PDAF', 'ensemble transformation (14):', pdaf_time_tot(14), 's'
     IF (rank_ana_enkf > 0) THEN
        WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'compute representers (13):', pdaf_time_tot(13), 's'
        WRITE (*, '(a, 27x, a, F11.3, 1x, a)') 'PDAF', 'pseudo inverse (36):', pdaf_time_tot(36), 's'
        WRITE (*, '(a, 21x, a, F11.3, 1x, a)') 'PDAF', 'compute representers (37):', pdaf_time_tot(37), 's'
     ELSE
        WRITE (*, '(a, 17x, a, F11.3, 1x, a)') 'PDAF', 'solve for representers (13):', pdaf_time_tot(13), 's'
     END IF
     IF (subtype_filter == 0) THEN
        WRITE (*, '(a, 16x, a, F11.3, 1x, a)') 'PDAF', 'blocked ensemble update (16):', pdaf_time_tot(16), 's'
        WRITE (*, '(a, 27x, a, F11.3, 1x, a)') 'PDAF', 'HX*representer (31):', pdaf_time_tot(31), 's'
        WRITE (*, '(a, 29x, a, F11.3, 1x, a)') 'PDAF', 'blocked X-Xm (35):', pdaf_time_tot(35), 's'
        WRITE (*, '(a, 26x, a, F11.3, 1x, a)') 'PDAF', 'update ensemble (38):', pdaf_time_tot(38), 's'
     ELSE
        WRITE (*, '(a, 24x, a, F11.3, 1x, a)') 'PDAF', 'ensemble update (16):', pdaf_time_tot(16), 's'
     END IF

     ! Generic part
     WRITE (*, '(a, 17x, a, F11.3, 1x, a)') 'PDAF', 'Time of prepoststep (5):', pdaf_time_tot(5), 's'
  END IF ptype

  
END SUBROUTINE PDAF_enkf_memtime
