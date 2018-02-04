! Copyright (c) 2014-2016 Paul Kirchgessner
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
!$Id: PDAF-D_lnetf_memtime.F90 1672 2016-12-02 09:09:38Z lnerger $
!BOP
!
! !ROUTINE: PDAF_lnetf_memtime --- Display timing and memory information for LNETF
!
! !INTERFACE:
SUBROUTINE PDAF_lnetf_memtime(printtype)

! !DESCRIPTION:
! This routine displays the PDAF-internal timing and
! memory information for the LNETF.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2014-05 - Paul Kirchgessner - Initial code based on LETKF code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_timer, &
       ONLY: PDAF_time_tot
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount_get
  USE PDAF_mod_filter, &
       ONLY: subtype_filter, dim_lag, type_forget
  USE PDAF_mod_filterMPI, &
       ONLY: filterpe
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

  ptype: IF (printtype == 1 .and. filterpe) THEN

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
     WRITE (*, '(8x, a, F11.3, 1x, a)') 'Time for analysis/transform:', pdaf_time_tot(3), 's'
     WRITE (*, '(18x, a, F11.3, 1x, a)') 'global preparations:', pdaf_time_tot(4), 's'
     WRITE (*, '(18x, a, F11.3, 1x, a)') 'local analysis loop:', pdaf_time_tot(6), 's'

     ! Generic part B
     WRITE (*, '(16x, a, F11.3, 1x, a)') 'Time of prepoststep:', pdaf_time_tot(5), 's'

  ELSE IF (printtype == 2 .and. filterpe) THEN ptype

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

  ELSE IF (printtype == 3 .and. filterpe) THEN ptype

! *********************************************
! *** Print second-level timing information ***
! *********************************************

     ! Generic part
     WRITE (*, '(//21x, a)') 'PDAF Timing information'
     WRITE (*, '(10x, 45a)') ('-', i=1, 45)
     WRITE (*, '(13x, a, F11.3, 1x, a)') &
          'Generate state ensemble (1):', pdaf_time_tot(1), 's'
     IF (subtype_filter /= 5) THEN
        WRITE (*, '(19x, a, F11.3, 1x, a)') 'Time of forecasts (2):', pdaf_time_tot(2), 's'
     END IF

     ! Filter-specific part
     WRITE (*, '(15x, a, F11.3, 1x, a)') 'Time for assimilation (3):', pdaf_time_tot(3), 's'
     WRITE (*, '(19x, a, F11.3, 1x, a)') 'global preparations (4):', pdaf_time_tot(4), 's'
     WRITE (*, '(25x, a, F11.3, 1x, a)') 'get mean state (11):', pdaf_time_tot(11), 's'
     WRITE (*, '(19x, a, F11.3, 1x, a)') 'local analysis loop (6):', pdaf_time_tot(6), 's'
     WRITE (*, '(16x, a, F11.3, 1x, a)') 'search local obs. domain (9):', pdaf_time_tot(9), 's'
     WRITE (*, '(24x, a, F11.3, 1x, a)') 'global to local (15):', pdaf_time_tot(15), 's'
     WRITE (*, '(26x, a, F11.3, 1x, a)') 'local analysis (7):', pdaf_time_tot(7), 's'
     WRITE (*, '(24x, a, F11.3, 1x, a)') 'local to global (16):', pdaf_time_tot(16), 's'
     IF (dim_lag >0) THEN
        WRITE (*, '(13x, a, F11.3, 1x, a)') 'compute smoother transform (17):', pdaf_time_tot(17), 's'
        WRITE (*, '(22x, a, F11.3, 1x, a)') 'perform smoothing (18):', pdaf_time_tot(18), 's'
     END IF

     ! Generic part B
     WRITE (*, '(17x, a, F11.3, 1x, a)') 'Time of prepoststep (5):', pdaf_time_tot(5), 's'

  ELSE IF (printtype == 4 .and. filterpe) THEN ptype

! *****************************************
! *** Print detailed timing information ***
! *****************************************

     ! Generic part
     WRITE (*, '(//21x, a)') 'PDAF Timing information'
     WRITE (*, '(10x, 45a)') ('-', i=1, 45)
     WRITE (*, '(13x, a, F11.3, 1x, a)') &
          'Generate state ensemble (1):', pdaf_time_tot(1), 's'
     IF (subtype_filter /= 5) THEN
        WRITE (*, '(19x, a, F11.3, 1x, a)') 'Time of forecasts (2):', pdaf_time_tot(2), 's'
     END IF

     ! Filter-specific part
     WRITE (*, '(15x, a, F11.3, 1x, a)') 'Time for assimilation (3):', pdaf_time_tot(3), 's'
     WRITE (*, '(19x, a, F11.3, 1x, a)') 'global preparations (4):', pdaf_time_tot(4), 's'
     IF (dim_lag>0) &
          WRITE (*, '(16x, a, F11.3, 1x, a)')  'compute HX for smoother (27):', pdaf_time_tot(27), 's'
     IF (type_forget==0 .OR. type_forget==1) &
          WRITE (*, '(16x, a, F11.3, 1x, a)') 'apply forgetting factor (14):', pdaf_time_tot(14), 's'
     WRITE (*, '(18x, a, F11.3, 1x, a)')  'compute HX for filter (12):', pdaf_time_tot(12), 's'
     WRITE (*, '(25x, a, F11.3, 1x, a)') 'get mean state (11):', pdaf_time_tot(11), 's'
     WRITE (*, '(17x, a, F11.3, 1x, a)') 'generate random matrix (13):', pdaf_time_tot(13), 's'
     WRITE (*, '(19x, a, F11.3, 1x, a)') 'local analysis loop (6):', pdaf_time_tot(6), 's'
     WRITE (*, '(16x, a, F11.3, 1x, a)') 'search local obs. domain (9):', pdaf_time_tot(9), 's'
     WRITE (*, '(24x, a, F11.3, 1x, a)') 'global to local (15):', pdaf_time_tot(15), 's'
     WRITE (*, '(26x, a, F11.3, 1x, a)') 'local analysis (7):', pdaf_time_tot(7), 's'
     WRITE (*, '(27x, a, F11.3, 1x, a)') 'init local obs (21):', pdaf_time_tot(21), 's'
     WRITE (*, '(17x, a, F11.3, 1x, a)') 'compute particle weights (22):', pdaf_time_tot(22), 's'
     WRITE (*, '(32x, a, F11.3, 1x, a)') 'compute A (23):', pdaf_time_tot(23), 's'
     WRITE (*, '(17x, a, F11.3, 1x, a)') 'compute transform matrix (24):', pdaf_time_tot(24), 's'
       WRITE (*, '(16x, a, F11.3, 1x, a)') 'compute SVD of W-ww^t and T (31):', pdaf_time_tot(31), 's'
       WRITE (*, '(17x, a, F11.3, 1x, a)') 'compute product of T and A (32):', pdaf_time_tot(32), 's'
       WRITE (*, '(24x, a, F11.3, 1x, a)') 'multiply with svals (33):', pdaf_time_tot(33), 's'
       IF (type_forget==2 .OR. type_forget==3) &
            WRITE (*, '(20x, a, F11.3, 1x, a)') 'apply forgetting factor (34):', pdaf_time_tot(34), 's'
       WRITE (*, '(24x, a, F11.3, 1x, a)') 'apply random matrix (35):', pdaf_time_tot(35), 's'
     WRITE (*, '(23x, a, F11.3, 1x, a)') 'transform ensemble (25):', pdaf_time_tot(25), 's'
       WRITE (*, '(22x, a, F11.3, 1x, a)') 'store ensemble matrix (36):', pdaf_time_tot(36), 's'
       WRITE (*, '(28x, a, F11.3, 1x, a)') 'update ensemble (37):', pdaf_time_tot(37), 's'

     WRITE (*, '(24x, a, F11.3, 1x, a)') 'local to global (16):', pdaf_time_tot(16), 's'
     IF (dim_lag >0) THEN
        WRITE (*, '(13x, a, F11.3, 1x, a)') 'compute smoother transform (17):', pdaf_time_tot(17), 's'
        WRITE (*, '(22x, a, F11.3, 1x, a)') 'perform smoothing (18):', pdaf_time_tot(18), 's'
     END IF

     ! Generic part B
     WRITE (*, '(17x, a, F11.3, 1x, a)') 'Time of prepoststep (5):', pdaf_time_tot(5), 's'

  END IF ptype


END SUBROUTINE PDAF_lnetf_memtime
