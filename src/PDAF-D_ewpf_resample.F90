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
! !ROUTINE: PDAF_ewpf_resample --- resampling indices
!
! !INTERFACE:
SUBROUTINE PDAF_ewpf_resample(dim_ens, weights, index_out)

! !DESCRIPTION:
! Compute resampling indices
!
! Variant for ETKF with domain decomposition.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2014-05 - Paul Kirchgessner - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

  INTEGER, INTENT(in) :: dim_ens		!Number of particles
  REAL, DIMENSION(dim_ens), INTENT(inout) :: weights !weights of particles
  INTEGER, INTENT(out), DIMENSION(dim_ens) :: index_out

  ! LOCAL VARIALBES
  INTEGER :: old, particle, k !index 
  INTEGER, SAVE :: iseed(4) !random seed
  INTEGER :: dupepos, dupe
  INTEGER :: i,j
  REAL :: draw, point
  REAL, ALLOCATABLE :: cumweights(:) !cumulative weight
  INTEGER, ALLOCATABLE :: NEW(:)
  INTEGER, ALLOCATABLE :: tobereplaced(:)
  LOGICAL :: in

  ALLOCATE(cumweights(dim_ens))
  ALLOCATE(NEW(dim_ens))
  ALLOCATE(tobereplaced(dim_ens))

  ! Normalize weights
  weights = EXP(-weights + MINVAL(weights))
  weights = weights/SUM(weights)

  ! Calculate cumulative weithgs
  cumweights = 0.0

  cumweights(1) = weights(1)
  DO i = 2,dim_ens
     cumweights(i) = cumweights(i-1) + weights(i)
  END DO

!LN  CALL PDAF_randomseed(iseed)
  iseed(1) = 10
  iseed(2) = 102
  iseed(3) = 230
  iseed(4) = 17
  

  CALL dlarnv(1,iseed,1,draw)
  draw = draw/dim_ens

  NEW = 0
  old = 1

  DO particle = 1,dim_ens
     point = REAL(particle-1)/dim_ens + draw

     DO
        IF (point <= cumweights(old)) THEN
           NEW(particle) = old
           EXIT
        ELSE
           old = old + 1
        END IF
     END DO
  END DO

  tobereplaced = -1


  DO i = 1,dim_ens
     in = .FALSE.
     DO j= 1,dim_ens
        IF(i == NEW(j)) THEN
           in = .TRUE.
           EXIT
        END IF
     END DO
     IF(.NOT. in) tobereplaced(i) = i
  END DO

  dupepos = 1
  index_out = (/ (i,i=1,dim_ens) /)
  DO i = 1, dim_ens
     IF(tobereplaced(i) /= -1) THEN
        DO k = dupepos, dim_ens
           IF(NEW(k) == NEW(k+1)) THEN
              dupe = NEW(k+1)
              dupepos = k+1
              EXIT
           END IF
        END DO
        index_out(tobereplaced(i)) = dupe
     END IF
  END DO

  weights = REAL(1.0)/REAL(dim_ens)

  weights = -LOG(weights)

  DEALLOCATE(cumweights, NEW, tobereplaced)

END SUBROUTINE PDAF_ewpf_resample

!-------------------------------------------

SUBROUTINE PDAF_randomseed(seed)

  IMPLICIT NONE

  INTEGER, INTENT(inout) :: seed(4)
  INTEGER :: clock

  ! I do not really recommend to use this routine,
  ! because the seed is not reproducible

  !generate random seed for the use of blas 
  CALL SYSTEM_CLOCK(count=clock)

  seed(1) = MOD(clock, 4095)+2
  seed(2) = MOD(clock-1431,4095)+17
  seed(3) = MOD(clock+10,3087)+13

  CALL SYSTEM_CLOCK(count= clock)
  seed(4) = 2*MOD(clock,2045)+1 !needs to be odd 

END SUBROUTINE PDAF_randomseed
