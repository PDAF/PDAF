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
!> Interface definitions for PDAF
!!
!! Module providing interface definition of the PDAF routines that
!! are called from the model code.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code by copying PDAF_interfaces_module
!! * Other revisions - see repository log
MODULE PDAF

  USE PDAF_DA
  USE PDAF_IAU, &
       ONLY: PDAF_iau_init, PDAF_iau_init_inc, PDAF_iau_reset
  USE PDAF_analysis_utils, ONLY: PDAF_seik_Omega, PDAF_ens_Omega, &
       PDAF_subtract_rowmean, PDAF_subtract_colmean, PDAF_set_forget, &
       PDAF_inflate_ens, PDAF_add_particle_noise
  USE PDAFlocal_interfaces
  USE PDAF_assim_interfaces
  USE PDAF_utils_interfaces
  USE PDAFomi_interfaces

END MODULE PDAF
