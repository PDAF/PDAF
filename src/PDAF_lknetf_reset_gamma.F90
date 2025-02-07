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
!> reset gamma value of LKNETF
!!
!! This routine resets the hybrid weight value
!! gamma in the LKNETF.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Resision history
!! * 2019-09 - Lars Nerger - initial code
!! * Later revisions - see repository log
!!
SUBROUTINE PDAF_lknetf_reset_gamma(gamma_in)

  USE PDAF_lknetf, &
       ONLY: hyb_g

  IMPLICIT NONE

! *** Argument ***
  REAL, INTENT(in) :: gamma_in     !< Prescribed hybrid weight


! *** Set hybrid weights ***
  hyb_g = gamma_in

END SUBROUTINE PDAF_lknetf_reset_gamma
