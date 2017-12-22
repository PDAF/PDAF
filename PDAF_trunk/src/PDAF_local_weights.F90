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
!$Id: PDAF_local_weights.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_local_weights --- Compute weight functions for localization
!
! !INTERFACE:
SUBROUTINE PDAF_local_weights(wtype, cradius, sradius, dim, distance, &
     weight, verbose)

! !DESCRIPTION:
! This routine initializates a vector hold weight coefficients
! for localization. The weights can be applied to localization
! on the state covariance matrix of the observation error 
! covariance matrix.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2010-06 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: wtype          ! Type of weight function
          ! (0): unit weight (=1 up to distance=cradius)
          ! (1): exponential decrease (1/e at distance=sradius; 0 for distance>cradius)
          ! (2): 5th order polynomial (Gaspari&Cohn 1999; 0 for distance>cradius)
  REAL, INTENT(in)    :: cradius        ! Parameter for cut-off
  REAL, INTENT(in)    :: sradius        ! Support radius 
  INTEGER, INTENT(in) :: dim            ! Size of distance and weight arrays
  REAL, INTENT(in)    :: distance(dim)  ! Array holding distances
  REAL, INTENT(out)   :: weight(dim)    ! Array for weights
  INTEGER, INTENT(in) :: verbose        ! Verbosity flag
  

! *** Local variables ***
  INTEGER :: i   ! Counter
  INTEGER :: verbose_s  ! verbosity flag
  REAL    :: distance_s ! Distance for single observation
  REAL    :: weight_s   ! Weight for single observation


! *******************************
! *** Initialize weight array ***
! *******************************

  DO i = 1, dim

     IF (verbose >= 1 .AND. i == 1) THEN
        verbose_s = 1
     ELSE
        verbose_s = 0
     END IF

     ! Set distance
     distance_s = distance(i)

     ! Compute weight
     CALL PDAF_local_weight(wtype, 0, cradius, sradius, distance_s, &
     1, 1, 1.0, 0.0, weight_s, verbose_s)

     ! Initialize element of weight array
     weight(i) = weight_s

  END DO

END SUBROUTINE PDAF_local_weights
