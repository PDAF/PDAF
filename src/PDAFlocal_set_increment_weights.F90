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

!> Set vector of local increment weights
!!
!! This routine initializes a PDAF_internal local array
!! of increment weights. The weights are applied in 
!! in PDAF_local_l2g_cb, when the global state vector
!! is initialized from the local state vector. These can
!! e.g. be used to apply a vertical localization.
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
SUBROUTINE PDAFlocal_set_increment_weights(dim_l, weights)

  USE PDAF_mod_filter, &
       ONLY: debug
  USE PDAFlocal, &
       ONLY: l2g_weights

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_l          !< Dimension of local state vector
  REAL, INTENT(in) :: weights(dim_l)    !< Weights array


! ********************************************
! *** Initialize PDAF_internal index array ***
! ********************************************

  IF (ALLOCATED(l2g_weights)) DEALLOCATE(l2g_weights)
  ALLOCATE(l2g_weights(dim_l))

  l2g_weights(:) = weights(:)

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAFlocal_set_increment_weights -- Set local increment weights'
     WRITE (*,*) '++ PDAF-debug PDAFlocal_set_increment_weights:', debug, 'weights', l2g_weights(1:dim_l)
  END IF

END SUBROUTINE PDAFlocal_set_increment_weights
