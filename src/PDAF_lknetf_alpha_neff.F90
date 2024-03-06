! Copyright (c) 2019-2024 Lars Nerger
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
!BOP
!
! !ROUTINE: PDAF_lknetf_alpha_neff --- get tempering weight according to N_eff
!
! !INTERFACE:
SUBROUTINE PDAF_lknetf_alpha_neff(dim_ens, weights, hlimit, alpha)

! !DESCRIPTION:
! Routine to compute an adaptive tempering factor alpha
! according to the effective sample size.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2019-08 - Lars Nerger
! Later revisions - see svn log
!
! !USES:

  IMPLICIT NONE

! !ARGUMENTS:
! ! Variable naming scheme:
! !   suffix _p: Denotes a full variable on the PE-local domain
! !   suffix _l: Denotes a local variable on the current analysis domain
  INTEGER, INTENT(in) :: dim_ens        ! Size of ensemble 
  REAL, INTENT(in) :: weights(dim_ens)  ! Weights
  REAL, INTENT(in) :: hlimit            ! Minimum of n_eff / N
  REAL, INTENT(inout) :: alpha          ! hybrid weight
  
! !CALLING SEQUENCE:
! Called by: PDAF_lknetf_analysis_T
! Calls: PDAF_memcount
!EOP
       
! Local variables
  INTEGER :: i
  REAL, ALLOCATABLE :: locw(:)
  REAL, ALLOCATABLE :: hweights(:)
  REAL :: a_step, tot_weight
  REAL :: nhlim, n_eff


! *****************************************
! *** Iteratively compute alpha so that ***
! ***      N_eff/dim_ens > hlimit       ***
! *****************************************


  ALLOCATE(locw(dim_ens))
  ALLOCATE(hweights(dim_ens))

  ! Get logarithm of weights
  DO i=1, dim_ens
     locw(i) = LOG(weights(i))
  END DO

  ! Initialize iterations
  alpha = 0.0
  a_step = 0.05
  nhlim = ABS(hlimit * REAL(dim_ens))

  aloop: DO

     ! scale 
     DO i = 1, dim_ens
        hweights(i) = EXP(locw(i) * (1-alpha))
     END DO
  
     ! Normalize weights
     tot_weight = 0.0
     DO i = 1, dim_ens
        tot_weight = tot_weight + hweights(i)
     END DO
     IF (tot_weight /= 0.0) THEN
        hweights = hweights / tot_weight

        ! Compute effective ensemble size
        CALL PDAF_diag_effsample(dim_ens, hweights, n_eff)

        IF (REAL(n_eff) >= nhlim) EXIT aloop
     END IF

     IF (alpha>=1.0) THEN
        alpha = 1.0
        EXIT aloop
     END IF

     alpha = alpha + a_step

  END DO aloop


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(hweights, locw)

END SUBROUTINE PDAF_lknetf_alpha_neff
