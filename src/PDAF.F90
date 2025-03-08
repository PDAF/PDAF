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

  USE PDAF_lknetf, &
       ONLY: PDAF_lknetf_reset_gamma

  ! Interfaces for init and get_state
  USE PDAFinit
  USE PDAF_init_mod
  USE PDAFget_state

  ! Full put_state interfaces
  USE PDAFput_state_enkf
  USE PDAFput_state_ensrf
  USE PDAFput_state_estkf
  USE PDAFput_state_etkf
  USE PDAFput_state_lenkf
  USE PDAFput_state_lestkf
  USE PDAFput_state_letkf
  USE PDAFput_state_lknetf
  USE PDAFput_state_lnetf
  USE PDAFput_state_lseik
  USE PDAFput_state_netf
  USE PDAFput_state_pf
  USE PDAFput_state_prepost
  USE PDAFput_state_seik

  USE PDAFput_state_generate_obs

  USE PDAFput_state_3dvar
  USE PDAFput_state_en3dvar_estkf
  USE PDAFput_state_en3dvar_lestkf
  USE PDAFput_state_hyb3dvar_estkf
  USE PDAFput_state_hyb3dvar_lestkf

  ! Full assimilate interfaces
  USE PDAFassimilate_enkf
  USE PDAFassimilate_ensrf
  USE PDAFassimilate_estkf
  USE PDAFassimilate_etkf
  USE PDAFassimilate_lenkf
  USE PDAFassimilate_lestkf
  USE PDAFassimilate_letkf
  USE PDAFassimilate_lknetf
  USE PDAFassimilate_lnetf
  USE PDAFassimilate_lseik
  USE PDAFassimilate_netf
  USE PDAFassimilate_pf
  USE PDAFassimilate_seik

  USE PDAFassimilate_prepost
  USE PDAFprepost

  USE PDAFgenerate_obs

  USE PDAFassimilate_3dvar
  USE PDAFassimilate_en3dvar_estkf
  USE PDAFassimilate_en3dvar_lestkf
  USE PDAFassimilate_hyb3dvar_estkf
  USE PDAFassimilate_hyb3dvar_lestkf

  ! PDAF-3 advanced interfaces
  USE PDAF_assimilate_ens
  USE PDAF_put_state_ens
  USE PDAF_assimilate_ens_nondiagR
  USE PDAF_put_state_ens_nondiagR

  USE PDAF_assimilate_3dvars
  USE PDAF_put_state_3dvars
  USE PDAF_assimilate_3dvars_nondiagR
  USE PDAF_put_state_3dvars_nondiagR


  ! PDAF-2 OMI interfaces
  USE PDAFomi_assimilate_ens
  USE PDAFomi_put_state_ens
  USE PDAFomi_assimilate_ens_nondiagR
  USE PDAFomi_put_state_ens_nondiagR

  USE PDAFomi_assimilate_3dvars
  USE PDAFomi_put_state_3dvars

  ! PDAF-2 LOCALOMI interfaces
  USE PDAFlocalomi_assimilate_ens
  USE PDAFlocalomi_put_state_ens
  USE PDAFlocalomi_assimilate_3dvars
  USE PDAFlocalomi_put_state_3dvars

  ! PDAF-2 LOCAL interfaces
  USE PDAFlocal_assimilate_ens
  USE PDAFlocal_put_state_ens
  USE PDAFlocal_assimilate_3dvars
  USE PDAFlocal_put_state_3dvars

  ! PDAF-OMI interfaces
  USE PDAFomi

END MODULE PDAF
