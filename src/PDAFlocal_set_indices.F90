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

!> Set index vector to map between global and local state vectors
!!
!! This routine initializes a PDAF_internal local index array
!! for the mapping between the global and local state vectors
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
SUBROUTINE PDAFlocal_set_indices(dim_l, map)

  USE PDAF_mod_filter, &
       ONLY: debug
  USE PDAFlocal, &
       ONLY: id_lstate_in_pstate

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_l          !< Dimension of local state vector
  INTEGER, INTENT(in) :: map(dim_l)     !< Index array for mapping


! ********************************************
! *** Initialize PDAF_internal index array ***
! ********************************************

  IF (ALLOCATED(id_lstate_in_pstate)) DEALLOCATE(id_lstate_in_pstate)
  ALLOCATE(id_lstate_in_pstate(dim_l))

  id_lstate_in_pstate(:) = map(:)

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAFlocal_set_indices:', debug, 'indices', id_lstate_in_pstate(1:dim_l)
  END IF

END SUBROUTINE PDAFlocal_set_indices
