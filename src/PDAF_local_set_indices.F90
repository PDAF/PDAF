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
!$Id: PDAF_set_comm_pdaf.F90 918 2021-12-03 07:42:19Z lnerger $


!> Get value of a correlation function
!!
!! This routine initializes a PDAF_internal local index array
!! for the mapping between the global and local state vectors
!!
!!  This is a core routine of PDAF and
!!  should not be changed by the user   !
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE PDAF_local_set_indices(dim_l, map)

  USE PDAF_mod_filter, &
       ONLY: id_lstate_in_pstate

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_l          !< Dimension of local state vector
  INTEGER, INTENT(in) :: map(dim_l)     !< Index array for mapping

! *** Local variables ***
  REAL :: shalf                         ! Half of support distance for Gaspari-Cohn


! ********************************************
! *** Initialize PDAF_internal index array ***
! ********************************************

  IF (ALLOCATED(id_lstate_in_pstate)) DEALLOCATE(id_lstate_in_pstate)
  ALLOCATE(id_lstate_in_pstate(dim_l))

  id_lstate_in_pstate(:) = map(:)

END SUBROUTINE PDAF_local_set_indices
