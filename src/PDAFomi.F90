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
!$Id: PDAFomi.F90 1147 2023-03-12 16:14:34Z lnerger $

!> PDAF-OMI main module
!!
!! This module combines the different PDAFomi modules
!! so that a user only needs to use-include from a
!! single module.
!!
!! __Revision history:__
!! * 2020-05 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAFomi

  USE PDAFomi_obs_f, &
       ONLY: obs_f, PDAFomi_gather_obs, PDAFomi_gather_obsstate, PDAFomi_init_obs_f, &
       PDAFomi_init_obsvars_f, PDAFomi_init_obsvar_f, PDAFomi_prodRinvA, PDAFomi_likelihood, &
       PDAFomi_add_obs_error, PDAFomi_init_obscovar, PDAFomi_init_obserr_f, PDAFomi_set_domain_limits, &
       PDAFomi_get_domain_limits_unstr, PDAFomi_get_local_ids_obs_f, PDAFomi_limit_obs_f, &
       PDAFomi_gather_dim_obs_f, PDAFomi_gather_obs_f_flex, PDAFomi_gather_obs_f2_flex, &
       PDAFomi_omit_by_inno, PDAFomi_obsstats, PDAFomi_gather_obsdims, PDAFomi_check_error, &
       PDAFomi_set_doassim, PDAFomi_set_disttype, PDAFomi_set_ncoord, PDAFomi_set_obs_err_type, &
       PDAFomi_set_use_global_obs, PDAFomi_set_inno_omit, PDAFomi_set_inno_omit_ivar, &
       PDAFomi_set_id_obs_p, PDAFomi_set_icoeff_p, PDAFomi_set_domainsize, PDAFomi_set_globalobs
  USE PDAFomi_obs_l
  USE PDAFomi_dim_obs_l
  USE PDAFomi_obs_op
  USE PDAFomi_obs_diag

  USE PDAFomi_assimilate_ens
  USE PDAFomi_assimilate_ens_nondiagR
  USE PDAFomi_assimilate_3dvars

  USE PDAFomi_put_state_ens
  USE PDAFomi_put_state_ens_nondiagR
  USE PDAFomi_put_state_3dvars

END MODULE PDAFomi

