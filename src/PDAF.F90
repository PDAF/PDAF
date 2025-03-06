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
  USE PDAF_set
  USE PDAF_get
  USe PDAF_utils
  USE PDAF_diag
  USE PDAF_sample
  USE PDAF_comm_obs
  USE PDAF_IAU, &
       ONLY: PDAF_iau_init, PDAF_iau_init_inc, PDAF_iau_reset, &
       PDAF_iau_set_weights, PDAF_iau_set_pointer, PDAF_iau_add_inc
  USE PDAF_analysis_utils, ONLY: PDAF_seik_Omega, PDAF_ens_Omega, &
       PDAF_seik_TtimesA, PDAF_seik_matrixT, &
       PDAF_subtract_rowmean, PDAF_subtract_colmean, &
       PDAF_inflate_ens, PDAF_add_particle_noise, PDAF_inflate_weights, &
       PDAF_generate_rndmat, PDAF_local_weight, PDAF_local_weights
  USE PDAF_info, &
       ONLY: PDAF_print_info
  USE PDAFlocal

  USE PDAF3_init
  USE PDAF_assimilate_ens
  USE PDAF_put_state_ens

  USE PDAF_assimilate_ens_nondiagR
  USE PDAF_put_state_ens_nondiagR

  USE PDAF_assimilate_3dvars
  USE PDAF_put_state_3dvars
  USE PDAF_assimilate_3dvars_nondiagR
  USE PDAF_put_state_3dvars_nondiagR

  USE PDAFomi_assimilate_ens
  USE PDAFomi_put_state_ens

  USE PDAFomi_assimilate_ens_nondiagR
  USE PDAFomi_put_state_ens_nondiagR

  USE PDAFomi_assimilate_3dvars
  USE PDAFomi_put_state_3dvars

  USE PDAFlocalomi_assimilate_ens
  USE PDAFlocalomi_put_state_ens
  USE PDAFlocalomi_assimilate_3dvars
  USE PDAFlocalomi_put_state_3dvars

  USE PDAFlocal_assimilate_ens
  USE PDAFlocal_put_state_ens
  USE PDAFlocal_assimilate_3dvars
  USE PDAFlocal_put_state_3dvars

END MODULE PDAF
