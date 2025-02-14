! Copyright (c) 2004-2025 Lars Nerger
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
!> PDAF_set_iparam --- Set integer parameter for PDAF
!!
!! This routine lets the user set the value of a
!! method-specific integer parameter. The routine
!! simply calls PDAF_set_iparam_filters, which 
!! includes the method-specific calls.
!!
!!    ! This is a core routine of PDAF and !
!!    ! should not be changed by the user  !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Other revisions - see repository log
!!
SUBROUTINE PDAF_set_iparam(id, value, flag)

  USE PDAF_utils_filters, &
       ONLY: PDAF_set_iparam_filters

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: id       !< Index of parameter
  INTEGER, INTENT(in) :: value    !< Parameter value
  INTEGER, INTENT(out) :: flag    !< Status flag: 0 for no error


! ********************************
! *** Print screen information ***
! ********************************

  CALL PDAF_set_iparam(id, value, flag)

END SUBROUTINE PDAF_set_iparam
