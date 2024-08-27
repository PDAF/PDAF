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

!> Deallocate vector of local increment weights
!!
!! This routine simply deallocates the local increment
!! weight vector if it is allocated.
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
SUBROUTINE PDAFlocal_clear_increment_weights()

  USE PDAF_mod_filter, &
       ONLY: debug
  USE PDAFlocal, &
       ONLY: l2g_weights

  IMPLICIT NONE


! *****************************************
! *** Deallocate increment weight array ***
! *****************************************

  IF (ALLOCATED(l2g_weights)) DEALLOCATE(l2g_weights)

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAFlocal_free_increment_weights -- Unset local increment weights'
  END IF

END SUBROUTINE PDAFlocal_clear_increment_weights
