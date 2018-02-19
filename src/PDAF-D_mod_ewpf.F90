! Copyright (c) 2004-2018 Lars Nerger
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
! !MODULE:
MODULE PDAF_mod_ewpf
  
! !DESCRIPTION:
! This module provides variables shared between the
! subroutines of PDAF.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2018-02 - Paul Kirchgessner - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filter

  IMPLICIT NONE
  SAVE

  LOGICAL :: use_model_error = .false.
  LOGICAL :: assimilation=.true.     ! Flag is assimilation has been performed
  LOGICAL :: isnudging = .false.

  INTEGER :: type_nudging = 1
  REAL :: bt = 0.3
  REAL :: start_nudging = 0.5
  REAL :: keep = 0.9 ! Fraction of particles keep at assimilation time


  REAL :: weight  ! use in parallel formulation of ewpf
  REAL, ALLOCATABLE :: weights(:)   ! PF weights
  REAL, ALLOCATABLE :: observation(:) !Store current observation for multiple
                                      ! use of them in the proposal filter.
  REAL, ALLOCATABLE :: x_last(:,:)    ! Ensemble matrix in case of a proposal filter
                                      !    or matrix of eigenvectors from EOF computation
  REAL, ALLOCATABLE :: HX_last(:)  ! store HX-obs from last timestep
  REAL,ALLOCATABLE :: state_last(:)
 
!EOP

END MODULE PDAF_mod_ewpf
