! Copyright (c) 2019-2025 Lars Nerger
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
!
!> Inflation for particle filters
!!
!! This routine compute an adaptive inflation using the effective sample
!! size N_eff according to
!!      N_eff / N >= alpha
!! whether N_eff is itertively computed with increasing inflation of the
!! observation error variance.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2019-08 - Lars Nerger
!! * Later revisions - see svn log
!!
SUBROUTINE PDAF_inflate_weights(screen, dim_ens, alpha, weights)

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: screen            !< verbosity flag
  INTEGER, INTENT(in) :: dim_ens           !< Ensemble size
  REAL, INTENT(inout) :: weights(dim_ens)  !< weights (before and after inflation)
  REAL, INTENT(in) :: alpha                !< Minimum limit of n_eff / N
       
! Local variables
  INTEGER :: i
  REAL :: alpha_iter
  REAL, ALLOCATABLE :: logw(:)      ! Logarith of particle weights
  REAL, ALLOCATABLE :: aweights(:)  ! temporary weights
  REAL :: a_step, tot_weight
  REAL :: alpha_lim, n_eff


! ******************
! *** Initialize ***
! ******************

  ALLOCATE(logw(dim_ens))
  ALLOCATE(aweights(dim_ens))

  ! Get logarithm of weights
  DO i=1, dim_ens
     logw(i) = LOG(weights(i))
  END DO

  ! Store initial weigts
  aweights = weights

  IF (screen>0) THEN
     WRITE (*,'(a, 5x, a, F10.3)') &
          'PDAF','--- Inflate weights according to N_eff/N > ', alpha
  END IF
  

! **********************************************
! *** Determine inflation according to alpha ***
! **********************************************

  alpha_iter = 0.0
  a_step = 0.05

  ! Set limit
  alpha_lim = alpha * REAL(dim_ens)

  aloop: DO

     IF (alpha_iter >= 1.0) THEN
        alpha_iter = 1.0
        EXIT aloop
     END IF

     ! scale 
     DO i = 1, dim_ens
        aweights(i) = EXP(logw(i) * (1.0-alpha_iter))
     END DO
  
     ! Normalize weights
     tot_weight = 0.0
     DO i = 1, dim_ens
        tot_weight = tot_weight + aweights(i)
     END DO
     IF (tot_weight /= 0.0) THEN
        aweights = aweights / tot_weight

        ! Compute effective ensemble size
        CALL PDAF_diag_effsample(dim_ens, aweights, n_eff)

        ! If limiting condition is fullfileed, exit loop
        IF (REAL(n_eff) > alpha_lim) EXIT aloop
     END IF

     alpha_iter = alpha_iter + a_step

  END DO aloop

  ! Store final inflated weigts
  weights = aweights


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(aweights, logw)

END SUBROUTINE PDAF_inflate_weights
