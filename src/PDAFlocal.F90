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

!> PDAF-LOCAL routines handling localization
!!
!! This module contains subroutines that handle the localization of
!! the state vector
!!
!! * PDAFlocal_set_indices \n
!!        Set indices of elements of lcoal state vector in global state vector
!! * PDAFlocal_set_increment_weights \n
!!        Set optional increment weights applied when upating the lobal state vector
!!        from the local analysis state vector
!! * PDAFlocal_clear_increment_weights \n
!!        Deallocate vector of increment weights. Afterwards l2g_state is applied without weights
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
MODULE PDAFlocal

  USE PDAFlocal_interfaces    ! Interface defintions for put_state and assimilate routines

  IMPLICIT NONE
  SAVE

  INTEGER, ALLOCATABLE :: id_lstate_in_pstate(:) !< Indices of local state vector in PE-local global state vector
  REAL, ALLOCATABLE :: l2g_weights(:)            !< Increment weights applied in l2g_state
  LOGICAL :: PDAFlocal_was_used = .FALSE.        !< Flag whether PDAFlocal was used (set in PDAFlocal_g2l_cb)

!$OMP THREADPRIVATE(id_lstate_in_pstate, l2g_weights)

!-------------------------------------------------------------------------------
  
  INTERFACE
     SUBROUTINE PDAFlocal_g2l_cb(step, domain_p, dim_p, state_p, dim_l, state_l)
       INTEGER, INTENT(in) :: step           !< Current time step
       INTEGER, INTENT(in) :: domain_p       !< Current local analysis domain
       INTEGER, INTENT(in) :: dim_p          !< PE-local full state dimension
       INTEGER, INTENT(in) :: dim_l          !< Local state dimension
       REAL, INTENT(in)    :: state_p(dim_p) !< PE-local full state vector 
       REAL, INTENT(out)   :: state_l(dim_l) !< State vector on local analysis domain
     END SUBROUTINE PDAFlocal_g2l_cb
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_l2g_cb(step, domain_p, dim_l, state_l, dim_p, state_p)
       INTEGER, INTENT(in) :: step           !< Current time step
       INTEGER, INTENT(in) :: domain_p       !< Current local analysis domain
       INTEGER, INTENT(in) :: dim_l          !< Local state dimension
       INTEGER, INTENT(in) :: dim_p          !< PE-local full state dimension
       REAL, INTENT(in)    :: state_l(dim_l) !< State vector on local analysis domain
       REAL, INTENT(inout) :: state_p(dim_p) !< PE-local full state vector 
     END SUBROUTINE PDAFlocal_l2g_cb
  END INTERFACE

END MODULE PDAFlocal
