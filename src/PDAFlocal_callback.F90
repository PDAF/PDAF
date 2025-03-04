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
!
!> Project global to local vector according to index array
!!
!! Project a global to a local state vector for the localized filters.
!! This is the full callback function to be used internally. The mapping 
!! is done using the index vector id_lstate_in_pstate that is initialize
!! in PDAF_local_set_index.
!!
!! __Revision history:__
!! 2024-08 - Lars Nerger - Initial code
!! Other revisions - see repository log
!!
SUBROUTINE PDAFlocal_g2l_cb(step, domain_p, dim_p, state_p, dim_l, state_l)

  USE PDAFlocal, &
       ONLY: id_lstate_in_pstate, PDAFlocal_was_used

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step           !< Current time step
  INTEGER, INTENT(in) :: domain_p       !< Current local analysis domain
  INTEGER, INTENT(in) :: dim_p          !< PE-local full state dimension
  INTEGER, INTENT(in) :: dim_l          !< Local state dimension
  REAL, INTENT(in)    :: state_p(dim_p) !< PE-local full state vector 
  REAL, INTENT(out)   :: state_l(dim_l) !< State vector on local analysis domain

! *** local variables ***
  INTEGER :: i                  ! Counter
  INTEGER :: dummy              ! Dummy integer variable


! *************************************
! *** Initialize local state vector ***
! *************************************

  ! Dummy initializations to prevent compiler warnings
  dummy = domain_p
  dummy = step

  ! Set flag that PDAFlocal was used
  PDAFlocal_was_used = .TRUE.

  DO i = 1, dim_l
     state_l(i) = state_p(id_lstate_in_pstate(i))
  END DO
   
END SUBROUTINE PDAFlocal_g2l_cb

!--------------------------------------------------------------------------
!>Initialize global vector elements from local state vector
!!
!! Initialize elements of a global state vector from a local state vector.
!! This is the full callback function to be used internally. The mapping 
!! is done using the index vector id_lstate_in_pstate that is initialize
!! in PDAF_local_set_index.
!!
!! To exclude any element of the local state vector from the initialization
!! one can set the corresponding index value to 0.
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAFlocal_l2g_cb(step, domain_p, dim_l, state_l, dim_p, state_p)

  USE PDAFlocal, &
       ONLY: id_lstate_in_pstate, l2g_weights

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step           !< Current time step
  INTEGER, INTENT(in) :: domain_p       !< Current local analysis domain
  INTEGER, INTENT(in) :: dim_l          !< Local state dimension
  INTEGER, INTENT(in) :: dim_p          !< PE-local full state dimension
  REAL, INTENT(in)    :: state_l(dim_l) !< State vector on local analysis domain
  REAL, INTENT(inout) :: state_p(dim_p) !< PE-local full state vector 
  
! *** local variables ***
  INTEGER :: i                  ! Counter
  INTEGER :: dummy              ! Dummy integer variable


! **************************************************
! *** Initialize elements of global state vector ***
! **************************************************

  ! Dummy initializations to prevent compiler warnings
  dummy = domain_p
  dummy = step

  IF (.NOT.ALLOCATED(l2g_weights)) THEN
     ! Initialize global state vector with full updated local state
     DO i = 1, dim_l
        state_p(id_lstate_in_pstate(i)) = state_l(i)
     END DO
  ELSE
     ! Apply increment weight when initializaing global state vector from local state vector
     DO i = 1, dim_l
        state_p(id_lstate_in_pstate(i)) = state_p(id_lstate_in_pstate(i)) &
             + l2g_weights(i) * (state_l(i) - state_p(id_lstate_in_pstate(i)))
     END DO
  END IF

END SUBROUTINE PDAFlocal_l2g_cb
