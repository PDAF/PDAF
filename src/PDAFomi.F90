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
!$Id: PDAFomi.F90 1147 2023-03-12 16:14:34Z lnerger $

!> PDAF-OMI main module
!!
!! This module combines the different PDAFomi modules
!! so that a user only needs to use-include from a
!! single module.
!!
!! __Revision history:__
!! * 2020-05 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE PDAFomi

  USE PDAFomi_obs_f
  USE PDAFomi_obs_l
  USE PDAFomi_dim_obs_l
  USE PDAFomi_obs_op

END MODULE PDAFomi

