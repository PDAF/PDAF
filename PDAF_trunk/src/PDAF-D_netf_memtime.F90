! Copyright (c) 2004-2018 Lars Nerger
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
! !ROUTINE: PDAF_netf_memtime --- Display timing and memory information for NETF
!
! !INTERFACE:
SUBROUTINE PDAF_netf_memtime(printtype)

! !DESCRIPTION:
! This routine displays the PDAF-internal timing and
! memory information for the SEIK filter.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2011-09 - Lars Nerger - Initial code
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
     WRITE (*, '(a, 21x, a)') 'PDAF', 'Allocated memory  (MiB)'
     WRITE (*, '(a, 14x, a, 1x, f10.3, a)') &
          'PDAF', 'state and A:', pdaf_memcount_get(1, 'M'), ' MiB (persistent)'
     WRITE (*, '(a, 11x, a, 1x, f10.3, a)') &
          'PDAF', 'ensemble array:', pdaf_memcount_get(2, 'M'), ' MiB (persistent)'
     WRITE (*, '(a, 12x, a, 1x, f10.3, a)') &
          'PDAF', 'analysis step:', pdaf_memcount_get(3, 'M'), ' MiB (temporary)'

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
        WRITE (*, '(a, 7x, a, F11.3, 1x, a)') 'PDAF', 'Time to collect/distribute ens (19):', pdaf_time_tot(19), 's'
     END IF

     ! Filter-specific part
     WRITE (*, '(a, 15x, a, F11.3, 1x, a)') 'PDAF', 'Time for assimilation (3):', pdaf_time_tot(3), 's'
     WRITE (*, '(a, 11x, a, F11.3, 1x, a)') 'PDAF', 'init observation dimension (15):', pdaf_time_tot(15), 's'
     WRITE (*, '(a, 15x, a, F11.3, 1x, a)') 'PDAF', 'compute filter weights (12):', pdaf_time_tot(12), 's'
     WRITE (*, '(a, 21x, a, F11.3, 1x, a)') 'PDAF', 'compute matrix A (10):', pdaf_time_tot(10), 's'
     WRITE (*, '(a, 13x, a, F11.3, 1x, a)') 'PDAF', 'compute transform matrix (13):', pdaf_time_tot(13), 's'
     WRITE (*, '(a, 16x, a, F11.3, 1x, a)') 'PDAF', 'store ensemble matrix (21):', pdaf_time_tot(21), 's'
     WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'transform ensemble (22):', pdaf_time_tot(22), 's'
     IF (dim_lag >0) THEN
        WRITE (*, '(a, 11x, a, F11.3, 1x, a)') 'PDAF', 'compute smoother transform (16):', pdaf_time_tot(16), 's'
        WRITE (*, '(a, 20x, a, F11.3, 1x, a)') 'PDAF', 'perform smoothing (17):', pdaf_time_tot(17), 's'
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
        WRITE (*, '(a, 7x, a, F11.3, 1x, a)') 'PDAF', 'Time to collect/distribute ens (19):', pdaf_time_tot(19), 's'
     END IF

     ! Filter-specific part
     WRITE (*, '(a, 15x, a, F11.3, 1x, a)') 'PDAF', 'Time for assimilation (3):', pdaf_time_tot(3), 's'
     WRITE (*, '(a, 11x, a, F11.3, 1x, a)') 'PDAF', 'init observation dimension (15):', pdaf_time_tot(15), 's'
     WRITE (*, '(a, 23x, a, F11.3, 1x, a)') 'PDAF', 'inflate ensemble (34):', pdaf_time_tot(34), 's'
     WRITE (*, '(a, 15x, a, F11.3, 1x, a)') 'PDAF', 'compute filter weights (12):', pdaf_time_tot(12), 's'
     WRITE (*, '(a, 21x, a, F11.3, 1x, a)') 'PDAF', 'compute matrix A (10):', pdaf_time_tot(10), 's'
     WRITE (*, '(a, 13x, a, F11.3, 1x, a)') 'PDAF', 'compute transform matrix (13):', pdaf_time_tot(13), 's'
     WRITE (*, '(a, 16x, a, F11.3, 1x, a)') 'PDAF', 'store ensemble matrix (21):', pdaf_time_tot(21), 's'
     WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'transform ensemble (22):', pdaf_time_tot(22), 's'
     IF (dim_lag >0) THEN
        WRITE (*, '(a, 11x, a, F11.3, 1x, a)') 'PDAF', 'compute smoother transform (16):', pdaf_time_tot(16), 's'
        WRITE (*, '(a, 20x, a, F11.3, 1x, a)') 'PDAF', 'perform smoothing (17):', pdaf_time_tot(17), 's'
     END IF

     ! Generic part B
     WRITE (*, '(a, 17x, a, F11.3, 1x, a)') 'PDAF', 'Time of prepoststep (5):', pdaf_time_tot(5), 's'
  END IF ptype


END SUBROUTINE PDAF_netf_memtime
