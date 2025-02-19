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
!> Module providing shared variables for DA-methods
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAF_DA
  
  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: PDAF_DA_SEIK = 1
  INTEGER, PARAMETER :: PDAF_DA_LSEIK = 3
  INTEGER, PARAMETER :: PDAF_DA_ENKF = 2
  INTEGER, PARAMETER :: PDAF_DA_LENKF = 8
  INTEGER, PARAMETER :: PDAF_DA_ETKF = 4
  INTEGER, PARAMETER :: PDAF_DA_LETKF = 5
  INTEGER, PARAMETER :: PDAF_DA_ESTKF = 6
  INTEGER, PARAMETER :: PDAF_DA_LESTKF = 7
  INTEGER, PARAMETER :: PDAF_DA_NETF = 9
  INTEGER, PARAMETER :: PDAF_DA_LNETF = 10
  INTEGER, PARAMETER :: PDAF_DA_LKNETF = 11
  INTEGER, PARAMETER :: PDAF_DA_PF = 12
  INTEGER, PARAMETER :: PDAF_DA_GENOBS = 100
  INTEGER, PARAMETER :: PDAF_DA_3DVAR = 200

END MODULE PDAF_DA
