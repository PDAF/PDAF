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
!> Set seedset for random number generation
!!
!! Helper routine for PDAF.
!! The routine allows to set the seedset index that
!! is used in PDAF_generate_rndmat. Values between
!! 1 and 20 are allowed.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! 2021-02 - Lars Nerger - Initial code
!! Later revisions - see svn log
!!
SUBROUTINE PDAF_set_seedset(seedset_in)

  USE PDAF_mod_filter, &
       ONLY: seedset, new_seedset

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER,INTENT(in) :: seedset_in    !< Seedset index (1-20)


! *** Set ensemble member ***

  IF (seedset_in>0 .AND. seedset_in<21) THEN
     seedset = seedset_in
     new_seedset = .TRUE.
  ELSE
     write (*,*) 'PDAF-ERROR: PDAF_set_seedset - Invalid value for seedset'
  END IF

END SUBROUTINE PDAF_set_seedset
