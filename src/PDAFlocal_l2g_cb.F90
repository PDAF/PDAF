! Copyright (c) 2004-2024 Lars Nerger
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
! !ROUTINE: PDAFlocal_l2g_cb - Initialize global vector elements from local state vector
!
! !INTERFACE:
SUBROUTINE PDAFlocal_l2g_cb(step, domain_p, dim_l, state_l, dim_p, state_p)

! !DESCRIPTION:
! Initialize elements of a global state vector from a local state vector.
! This is the full callback function to be used internally. The mapping 
! is done using the index vector id_lstate_in_pstate that is initialize
! in PDAF_local_set_index.
!
! To exclude any element of the local state vector from the initialization
! one can set the corresponding index value to 0.
!
! !REVISION HISTORY:
! 2024-08 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAFlocal, &
       ONLY: id_lstate_in_pstate, l2g_weights

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step           !< Current time step
  INTEGER, INTENT(in) :: domain_p       !< Current local analysis domain
  INTEGER, INTENT(in) :: dim_l          !< Local state dimension
  INTEGER, INTENT(in) :: dim_p          !< PE-local full state dimension
  REAL, INTENT(in)    :: state_l(dim_l) !< State vector on local analysis domain
  REAL, INTENT(inout) :: state_p(dim_p) !< PE-local full state vector 

! !CALLING SEQUENCE:
! Called by filter routine
!EOP
  
! *** local variables ***
  INTEGER :: i                  ! Counter


! **************************************************
! *** Initialize elements of global state vector ***
! **************************************************

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
