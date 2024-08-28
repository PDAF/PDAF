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
! !ROUTINE: PDAFlocal_g2l_cb - Project global to local vector according to index array
!
! !INTERFACE:
SUBROUTINE PDAFlocal_g2l_cb(step, domain_p, dim_p, state_p, dim_l, state_l)

! !DESCRIPTION:
! Project a global to a local state vector for the localized filters.
! This is the full callback function to be used internally. The mapping 
! is done using the index vector id_lstate_in_pstate that is initialize
! in PDAF_local_set_index.
!
! !REVISION HISTORY:
! 2024-08 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAFlocal, &
       ONLY: id_lstate_in_pstate, PDAFlocal_was_used

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step           !< Current time step
  INTEGER, INTENT(in) :: domain_p       !< Current local analysis domain
  INTEGER, INTENT(in) :: dim_p          !< PE-local full state dimension
  INTEGER, INTENT(in) :: dim_l          !< Local state dimension
  REAL, INTENT(in)    :: state_p(dim_p) !< PE-local full state vector 
  REAL, INTENT(out)   :: state_l(dim_l) !< State vector on local analysis domain

! !CALLING SEQUENCE:
! Called by filter routine
!EOP
  
! *** local variables ***
  INTEGER :: i                  ! Counter


! *************************************
! *** Initialize local state vector ***
! *************************************

  ! Set flag that PDAFlocal was used
  PDAFlocal_was_used = .TRUE.

  DO i = 1, dim_l
     state_l(i) = state_p(id_lstate_in_pstate(i))
  END DO
   
END SUBROUTINE PDAFlocal_g2l_cb
