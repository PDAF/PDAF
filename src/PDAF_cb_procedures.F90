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
!> Abstract interfaces for call-back routines
!!
!! This module provides the abstract interfaces for the 
!! user-supplied call-back routines. Together with using
!! procedure declarations they allow to obtain an
!! explicit interface for the call-back routines
MODULE PDAF_cb_procedures

  IMPLICIT NONE

  ABSTRACT INTERFACE 
     SUBROUTINE init_ens_cb(filtertype, dim_p, dim_ens, &
          state_p, Uinv, ens_p, flag)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: filtertype      !< Type of filter to initialize
       INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
       INTEGER, INTENT(in) :: dim_ens         !< Size of ensemble
       REAL, INTENT(inout) :: state_p(dim_p)              !< PE-local model state
       REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1)   !< Array not referenced for ensemble filters
       REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)       !< PE-local state ensemble
       INTEGER, INTENT(inout) :: flag         !< PDAF status flag
     END SUBROUTINE init_ens_cb
  END INTERFACE


!!   U_init_ens
!   U_collect_state
!   U_distribute_state
!   U_next_observation
!   U_prodRinvA
! 
!   U_init_dim_obs
!   U_obs_op
!   U_init_obsvar
!   U_init_obsvars
!   U_init_obs
!   U_prepoststep
!   U_localize
!   U_localize_covar_serial
!   U_add_obs_err
!   U_init_obs_covar
!   U_init_obserr_f
!   U_get_obs_f
! 
!   U_init_n_domains_p
!   U_init_dim_l
!   U_g2l_state
!   U_l2g_state
! 
!   U_init_obs_l
!   U_prodRinvA_l
!   U_likelihood_l
!   U_init_dim_obs_l
!   U_g2l_obs
!   U_init_obsvar_l
! 
!   U_prodRinvA_hyb_l
!   U_likelihood_hyb_l

END MODULE PDAF_cb_procedures
