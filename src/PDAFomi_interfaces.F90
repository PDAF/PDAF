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
!! * 2025-02 - Lars Nerger - Initial code split from PDAF_interfaces_module
!! * Later revisions - see svn log
MODULE PDAFomi_interfaces

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
!#include "typedefs.h"

  INTERFACE
     SUBROUTINE PDAFomi_put_state_global(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_prepoststep, flag)
       INTEGER, INTENT(out) :: flag   ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &         ! Initialize dimension of observation vector
            U_obs_op, &               ! Observation operator
            U_prepoststep             ! User supplied pre/poststep routine
     END SUBROUTINE PDAFomi_put_state_global
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_global_nondiagR(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_prodRinvA, U_prepoststep, flag)
       INTEGER, INTENT(out) :: flag   ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &         ! Initialize dimension of observation vector
            U_obs_op, &               ! Observation operator
            U_prodRinvA, &            ! Provide product R^-1 A
            U_prepoststep             ! User supplied pre/poststep routine
     END SUBROUTINE PDAFomi_put_state_global_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_nonlin_nondiagR(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_likelihood, U_prepoststep, flag)
       INTEGER, INTENT(out) :: flag   ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &         ! Initialize dimension of observation vector
            U_obs_op, &               ! Observation operator
            U_likelihood, &           ! Compute likelihood
            U_prepoststep             ! User supplied pre/poststep routine
     END SUBROUTINE PDAFomi_put_state_nonlin_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_enkf_nondiagR(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_prepoststep, U_add_obs_err, U_init_obs_covar, flag)
       INTEGER, INTENT(out) :: flag   ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &         ! Initialize dimension of observation vector
            U_obs_op, &               ! Observation operator
            U_add_obs_err, &          ! Add obs error covariance R to HPH in EnKF
            U_init_obs_covar, &       ! Initialize obs. error cov. matrix R in EnKF
            U_prepoststep             ! User supplied pre/poststep routine
     END SUBROUTINE PDAFomi_put_state_enkf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_lenkf_nondiagR(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_prepoststep, U_localize, U_add_obs_err, U_init_obs_covar, flag)
       INTEGER, INTENT(out) :: flag   ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &         ! Initialize dimension of observation vector
            U_obs_op, &               ! Observation operator
            U_add_obs_err, &          ! Add obs error covariance R to HPH in EnKF
            U_init_obs_covar, &       ! Initialize obs. error cov. matrix R in EnKF
            U_localize, &             ! Apply localization to HP and HPH^T
            U_prepoststep             ! User supplied pre/poststep routine
     END SUBROUTINE PDAFomi_put_state_lenkf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_global(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_prepoststep, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &      ! Routine to distribute a state vector
            U_init_dim_obs, &          ! Initialize dimension of observation vector
            U_obs_op, &                ! Observation operator
            U_prepoststep, &           ! User supplied pre/poststep routine
            U_next_observation         ! Provide time step and time of next observation
     END SUBROUTINE PDAFomi_assimilate_global
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_local(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_prepoststep, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
          U_g2l_state, U_l2g_state, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_obs_op, &                ! Observation operator
            U_init_n_domains_p, &      ! Provide number of local analysis domains
            U_init_dim_l, &            ! Init state dimension for local ana. domain
            U_init_dim_obs, &          ! Initialize dimension of observation vector
            U_init_dim_obs_l, &        ! Initialize dim. of obs. vector for local ana. domain
            U_g2l_state, &             ! Get state on local ana. domain from full state
            U_l2g_state, &             ! Init full state from state on local analysis domain
            U_prepoststep              ! User supplied pre/poststep routine
     END SUBROUTINE PDAFomi_put_state_local
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_local_nondiagR(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_prepoststep, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_prodRinvA_l, &
          U_g2l_state, U_l2g_state, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_obs_op, &                ! Observation operator
            U_init_n_domains_p, &      ! Provide number of local analysis domains
            U_init_dim_l, &            ! Init state dimension for local ana. domain
            U_init_dim_obs, &          ! Initialize dimension of observation vector
            U_init_dim_obs_l, &        ! Initialize dim. of obs. vector for local ana. domain
            U_prodRinvA_l, &           ! Provide product R^-1 A on local analysis domain
            U_g2l_state, &             ! Get state on local ana. domain from full state
            U_l2g_state, &             ! Init full state from state on local analysis domain
            U_prepoststep              ! User supplied pre/poststep routine
     END SUBROUTINE PDAFomi_put_state_local_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_lnetf_nondiagR(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_prepoststep, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_likelihood_l, &
          U_g2l_state, U_l2g_state, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_obs_op, &                ! Observation operator
            U_init_n_domains_p, &      ! Provide number of local analysis domains
            U_init_dim_l, &            ! Init state dimension for local ana. domain
            U_init_dim_obs, &          ! Initialize dimension of observation vector
            U_init_dim_obs_l, &        ! Initialize dim. of obs. vector for local ana. domain
            U_likelihood_l, &          ! Compute likelihood and apply localization
            U_g2l_state, &             ! Get state on local ana. domain from full state
            U_l2g_state, &             ! Init full state from state on local analysis domain
            U_prepoststep              ! User supplied pre/poststep routine
     END SUBROUTINE PDAFomi_put_state_lnetf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_lknetf_nondiagR(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_prepoststep, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
          U_prodRinvA_l, U_prodRinvA_hyb_l, U_likelihood_l, U_likelihood_hyb_l, &
          U_g2l_state, U_l2g_state, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_obs_op, &                ! Observation operator
            U_init_n_domains_p, &      ! Provide number of local analysis domains
            U_init_dim_l, &            ! Init state dimension for local ana. domain
            U_init_dim_obs, &          ! Initialize dimension of observation vector
            U_init_dim_obs_l, &        ! Initialize dim. of obs. vector for local ana. domain
            U_prodRinvA_l, &           ! Provide product R^-1 A on local analysis domain
            U_prodRinvA_hyb_l, &       ! Product R^-1 A on local analysis domain with hybrid weight
            U_likelihood_l, &          ! Compute likelihood and apply localization
            U_likelihood_hyb_l, &      ! Compute likelihood and apply localization with tempering
            U_g2l_state, &             ! Get state on local ana. domain from full state
            U_l2g_state, &             ! Init full state from state on local analysis domain
            U_prepoststep              ! User supplied pre/poststep routine
     END SUBROUTINE PDAFomi_put_state_lknetf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_local(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_prepoststep, U_init_n_domains_p, U_init_dim_l, &
          U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &      ! Routine to distribute a state vector
            U_obs_op, &                ! Observation operator
            U_init_n_domains_p, &      ! Provide number of local analysis domains
            U_init_dim_l, &            ! Init state dimension for local ana. domain
            U_init_dim_obs, &          ! Initialize dimension of observation vector
            U_init_dim_obs_l, &        ! Initialize dim. of obs. vector for local ana. domain
            U_g2l_state, &             ! Get state on local ana. domain from full state
            U_l2g_state, &             ! Init full state from state on local analysis domain
            U_prepoststep, &           ! User supplied pre/poststep routine
            U_next_observation         ! Provide time step and time of next observation
     END SUBROUTINE PDAFomi_assimilate_local
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_local_nondiagR(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_prepoststep, U_init_n_domains_p, U_init_dim_l, &
          U_init_dim_obs_l, U_prodRinvA_l, U_g2l_state, U_l2g_state, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &      ! Routine to distribute a state vector
            U_obs_op, &                ! Observation operator
            U_init_n_domains_p, &      ! Provide number of local analysis domains
            U_init_dim_l, &            ! Init state dimension for local ana. domain
            U_init_dim_obs, &          ! Initialize dimension of observation vector
            U_init_dim_obs_l, &        ! Initialize dim. of obs. vector for local ana. domain
            U_prodRinvA_l, &           ! Provide product R^-1 A on local analysis domain
            U_g2l_state, &             ! Get state on local ana. domain from full state
            U_l2g_state, &             ! Init full state from state on local analysis domain
            U_prepoststep, &           ! User supplied pre/poststep routine
            U_next_observation         ! Provide time step and time of next observation
     END SUBROUTINE PDAFomi_assimilate_local_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_lnetf_nondiagR(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_prepoststep, U_init_n_domains_p, U_init_dim_l, &
          U_init_dim_obs_l, U_likelihood_l, U_g2l_state, U_l2g_state, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &      ! Routine to distribute a state vector
            U_obs_op, &                ! Observation operator
            U_init_n_domains_p, &      ! Provide number of local analysis domains
            U_init_dim_l, &            ! Init state dimension for local ana. domain
            U_init_dim_obs, &          ! Initialize dimension of observation vector
            U_init_dim_obs_l, &        ! Initialize dim. of obs. vector for local ana. domain
            U_likelihood_l, &          ! Compute likelihood and apply localization
            U_g2l_state, &             ! Get state on local ana. domain from full state
            U_l2g_state, &             ! Init full state from state on local analysis domain
            U_prepoststep, &           ! User supplied pre/poststep routine
            U_next_observation         ! Provide time step and time of next observation
     END SUBROUTINE PDAFomi_assimilate_lnetf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_lknetf_nondiagR(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_prepoststep, U_init_n_domains_p, U_init_dim_l, &
          U_init_dim_obs_l, U_prodRinvA_l, U_prodRinvA_hyb_l, U_likelihood_l, U_likelihood_hyb_l, &
          U_g2l_state, U_l2g_state, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &      ! Routine to distribute a state vector
            U_obs_op, &                ! Observation operator
            U_init_n_domains_p, &      ! Provide number of local analysis domains
            U_init_dim_l, &            ! Init state dimension for local ana. domain
            U_init_dim_obs, &          ! Initialize dimension of observation vector
            U_init_dim_obs_l, &        ! Initialize dim. of obs. vector for local ana. domain
            U_prodRinvA_l, &           ! Provide product R^-1 A on local analysis domain
            U_prodRinvA_hyb_l, &       ! Product R^-1 A on local analysis domain with hybrid weight
            U_likelihood_l, &          ! Compute likelihood and apply localization
            U_likelihood_hyb_l, &      ! Compute likelihood and apply localization with tempering
            U_g2l_state, &             ! Get state on local ana. domain from full state
            U_l2g_state, &             ! Init full state from state on local analysis domain
            U_prepoststep, &           ! User supplied pre/poststep routine
            U_next_observation         ! Provide time step and time of next observation
     END SUBROUTINE PDAFomi_assimilate_lknetf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_global_nondiagR(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_prodRinvA, U_prepoststep, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &      ! Routine to distribute a state vector
            U_init_dim_obs, &          ! Initialize dimension of observation vector
            U_obs_op, &                ! Observation operator
            U_prodRinvA, &             ! Provide product R^-1 A
            U_prepoststep, &           ! User supplied pre/poststep routine
            U_next_observation         ! Provide time step and time of next observation
     END SUBROUTINE PDAFomi_assimilate_global_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_nonlin_nondiagR(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_likelihood, U_prepoststep, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &      ! Routine to distribute a state vector
            U_init_dim_obs, &          ! Initialize dimension of observation vector
            U_obs_op, &                ! Observation operator
            U_likelihood, &            ! Compute likelihood
            U_prepoststep, &           ! User supplied pre/poststep routine
            U_next_observation         ! Provide time step and time of next observation
     END SUBROUTINE PDAFomi_assimilate_nonlin_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_enkf_nondiagR(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_add_obs_err, U_init_obs_covar, &
          U_prepoststep, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &      ! Routine to distribute a state vector
            U_init_dim_obs, &          ! Initialize dimension of observation vector
            U_obs_op, &                ! Observation operator
            U_add_obs_err, &           ! Add obs error covariance R to HPH in EnKF
            U_init_obs_covar, &        ! Initialize obs. error cov. matrix R in EnKF
            U_prepoststep, &           ! User supplied pre/poststep routine
            U_next_observation         ! Provide time step and time of next observation
     END SUBROUTINE PDAFomi_assimilate_enkf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_lenkf_nondiagR(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_prepoststep, U_localize, &
          U_add_obs_error, U_init_obs_covar, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_distribute_state, &     ! Routine to distribute a state vector
            U_init_dim_obs, &      ! Initialize dimension of observation vector
            U_obs_op, &            ! Observation operator
            U_init_obs_covar, &    ! Initialize obs. error cov. matrix R in EnKF
            U_prepoststep, &       ! User supplied pre/poststep routine
            U_localize, &          ! Apply localization to HP and HPH^T
            U_add_obs_error, &     ! Add obs error covariance R to HPH in EnKF
            U_next_observation     ! Provide time step and time of next observation
     END SUBROUTINE PDAFomi_assimilate_lenkf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_local_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFomi_put_state_local_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_local_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFomi_assimilate_local_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_local_nondiagR_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFomi_assimilate_local_nondiagR_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_global_nondiagR_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFomi_assimilate_global_nondiagR_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_enkf_nondiagR_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFomi_assimilate_enkf_nondiagR_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_nonlin_nondiagR_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFomi_assimilate_nonlin_nondiagR_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_lenkf_nondiagR_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFomi_assimilate_lenkf_nondiagR_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_lnetf_nondiagR_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFomi_assimilate_lnetf_nondiagR_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_lknetf_nondiagR_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFomi_assimilate_lknetf_nondiagR_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_global_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFomi_put_state_global_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_global_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFomi_assimilate_global_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_lenkf_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFomi_assimilate_lenkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_lenkf_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFomi_put_state_lenkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_local_nondiagR_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFomi_put_state_local_nondiagR_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_global_nondiagR_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFomi_put_state_global_nondiagR_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_enkf_nondiagR_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFomi_put_state_enkf_nondiagR_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_nonlin_nondiagR_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFomi_put_state_nonlin_nondiagR_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_lenkf_nondiagR_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFomi_put_state_lenkf_nondiagR_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_lknetf_nondiagR_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFomi_put_state_lknetf_nondiagR_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_lnetf_nondiagR_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFomi_put_state_lnetf_nondiagR_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_lenkf(U_collect_state, U_init_dim_obs, U_obs_op,  &
          U_prepoststep, U_localize, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_init_dim_obs, &          ! Initialize dimension of observation vector
            U_obs_op, &                ! Observation operator
            U_prepoststep, &           ! User supplied pre/poststep routine
            U_localize                 ! Apply localization to HP and HPH^T
     END SUBROUTINE PDAFomi_put_state_lenkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_lenkf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_prepoststep, U_localize, &
          U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &      ! Routine to distribute a state vector
            U_init_dim_obs, &          ! Initialize dimension of observation vector
            U_obs_op, &                ! Observation operator
            U_prepoststep, &           ! User supplied pre/poststep routine
            U_localize, &              ! Apply localization to HP and HPH^T
            U_next_observation         ! Provide time step and time of next observation
     END SUBROUTINE PDAFomi_assimilate_lenkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_generate_obs(U_collect_state, U_init_dim_obs_f, U_obs_op_f, &
          U_get_obs_f, U_prepoststep, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_init_dim_obs_f, &        ! Initialize dimension of observation vector
            U_obs_op_f, &              ! Observation operator
            U_get_obs_f, &             ! Provide observation vector to user
            U_prepoststep              ! User supplied pre/poststep routine
     END SUBROUTINE PDAFomi_put_state_generate_obs
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_generate_obs(U_collect_state, U_distribute_state, &
          U_init_dim_obs_f, U_obs_op_f, U_get_obs_f, U_prepoststep, &
          U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &      ! Routine to distribute a state vector
            U_init_dim_obs_f, &        ! Initialize dimension of observation vector
            U_obs_op_f, &              ! Observation operator
            U_get_obs_f, &             ! Provide observation vector to user
            U_prepoststep, &           ! User supplied pre/poststep routine
            U_next_observation         ! Provide time step and time of next observation
     END SUBROUTINE PDAFomi_generate_obs
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_3dvar(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, &
          cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          prepoststep_pdaf, next_observation_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
            distribute_state_pdaf, &        ! Routine to distribute a state vector
            next_observation_pdaf, &        ! Provide time step, time and dimension of next observation
            prepoststep_pdaf                ! User supplied pre/poststep routine
       EXTERNAL :: init_dim_obs_pdaf, &     ! Initialize dimension of observation vector
            obs_op_pdaf, &                  ! Observation operator
            cvt_pdaf, &                     ! Apply control vector transform matrix to control vector
            cvt_adj_pdaf, &                 ! Apply adjoint control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
       EXTERNAL :: PDAFomi_init_obs_f_cb, & ! Initialize observation vector
            PDAFomi_prodRinvA_cb            ! Provide product R^-1 A
     END SUBROUTINE PDAFomi_assimilate_3dvar
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_3dvar_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, prodRinvA_pdaf, &
          cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          prepoststep_pdaf, next_observation_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
            distribute_state_pdaf, &        ! Routine to distribute a state vector
            next_observation_pdaf, &        ! Provide time step, time and dimension of next observation
            prepoststep_pdaf                ! User supplied pre/poststep routine
       EXTERNAL :: init_dim_obs_pdaf, &     ! Initialize dimension of observation vector
            obs_op_pdaf, &                  ! Observation operator
            prodRinvA_pdaf, &               ! Provide product R^-1 A
            cvt_pdaf, &                     ! Apply control vector transform matrix to control vector
            cvt_adj_pdaf, &                 ! Apply adjoint control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
     END SUBROUTINE PDAFomi_assimilate_3dvar_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_en3dvar_estkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          prepoststep_pdaf, next_observation_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
            distribute_state_pdaf, &        ! Routine to distribute a state vector
            next_observation_pdaf, &        ! Provide time step, time and dimension of next observation
            prepoststep_pdaf                ! User supplied pre/poststep routine
       EXTERNAL :: init_dim_obs_pdaf, &     ! Initialize dimension of observation vector
            obs_op_pdaf, &                  ! Observation operator
            cvt_ens_pdaf, &                 ! Apply control vector transform matrix to control vector
            cvt_adj_ens_pdaf, &             ! Apply adjoint control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
     END SUBROUTINE PDAFomi_assimilate_en3dvar_estkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_en3dvar_estkf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, prodRinvA_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          prepoststep_pdaf, next_observation_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
            distribute_state_pdaf, &        ! Routine to distribute a state vector
            next_observation_pdaf, &        ! Provide time step, time and dimension of next observation
            prepoststep_pdaf                ! User supplied pre/poststep routine
       EXTERNAL :: init_dim_obs_pdaf, &     ! Initialize dimension of observation vector
            obs_op_pdaf, &                  ! Observation operator
            prodRinvA_pdaf, &               ! Provide product R^-1 A
            cvt_ens_pdaf, &                 ! Apply control vector transform matrix to control vector
            cvt_adj_ens_pdaf, &             ! Apply adjoint control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
     END SUBROUTINE PDAFomi_assimilate_en3dvar_estkf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_en3dvar_lestkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, next_observation_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
            distribute_state_pdaf, &        ! Routine to distribute a state vector
            next_observation_pdaf, &        ! Provide time step, time and dimension of next observation
            prepoststep_pdaf                ! User supplied pre/poststep routine
       EXTERNAL :: cvt_ens_pdaf, &          ! Apply control vector transform matrix to control vector
            cvt_adj_ens_pdaf, &             ! Apply adjoint control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
       EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
            init_dim_l_pdaf, &              ! Init state dimension for local ana. domain
            g2l_state_pdaf, &               ! Get state on local ana. domain from full state
            l2g_state_pdaf, &               ! Init full state from local state
            init_dim_obs_f_pdaf, &          ! Initialize dimension of full observation vector
            obs_op_f_pdaf, &                ! Full observation operator
            init_dim_obs_l_pdaf             ! Initialize local dimimension of obs. vector
     END SUBROUTINE PDAFomi_assimilate_en3dvar_lestkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_en3dvar_lestkf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, prodRinvA_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, next_observation_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
            distribute_state_pdaf, &        ! Routine to distribute a state vector
            next_observation_pdaf, &        ! Provide time step, time and dimension of next observation
            prepoststep_pdaf                ! User supplied pre/poststep routine
       EXTERNAL :: cvt_ens_pdaf, &          ! Apply control vector transform matrix to control vector
            cvt_adj_ens_pdaf, &             ! Apply adjoint control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
       EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
            prodRinvA_pdaf, &               ! Provide product R^-1 A
            init_dim_l_pdaf, &              ! Init state dimension for local ana. domain
            g2l_state_pdaf, &               ! Get state on local ana. domain from full state
            l2g_state_pdaf, &               ! Init full state from local state
            init_dim_obs_f_pdaf, &          ! Initialize dimension of full observation vector
            obs_op_f_pdaf, &                ! Full observation operator
            init_dim_obs_l_pdaf, &          ! Initialize local dimimension of obs. vector
            prodRinvA_l_pdaf                ! Provide product R^-1 A with localization
     END SUBROUTINE PDAFomi_assimilate_en3dvar_lestkf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_hyb3dvar_estkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
          obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf, next_observation_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
            distribute_state_pdaf, &        ! Routine to distribute a state vector
            next_observation_pdaf, &        ! Provide time step, time and dimension of next observation
            prepoststep_pdaf                ! User supplied pre/poststep routine
       EXTERNAL :: init_dim_obs_pdaf, &     ! Initialize dimension of observation vector
            obs_op_pdaf, &                  ! Observation operator
            cvt_pdaf, &                     ! Apply control vector transform matrix to control vector
            cvt_adj_pdaf, &                 ! Apply adjoint control vector transform matrix
            cvt_ens_pdaf, &                 ! Apply ensemble control vector transform matrix to control vector
            cvt_adj_ens_pdaf, &             ! Apply adjoint ensemble control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
     END SUBROUTINE PDAFomi_assimilate_hyb3dvar_estkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_hyb3dvar_estkf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, prodRinvA_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
          obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf, next_observation_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
            distribute_state_pdaf, &        ! Routine to distribute a state vector
            next_observation_pdaf, &        ! Provide time step, time and dimension of next observation
            prepoststep_pdaf                ! User supplied pre/poststep routine
       EXTERNAL :: init_dim_obs_pdaf, &     ! Initialize dimension of observation vector
            obs_op_pdaf, &                  ! Observation operator
            prodRinvA_pdaf, &               ! Provide product R^-1 A
            cvt_pdaf, &                     ! Apply control vector transform matrix to control vector
            cvt_adj_pdaf, &                 ! Apply adjoint control vector transform matrix
            cvt_ens_pdaf, &                 ! Apply control vector transform matrix to control vector
            cvt_adj_ens_pdaf, &             ! Apply adjoint control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
     END SUBROUTINE PDAFomi_assimilate_hyb3dvar_estkf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_hyb3dvar_lestkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, next_observation_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
            distribute_state_pdaf, &        ! Routine to distribute a state vector
            next_observation_pdaf, &        ! Provide time step, time and dimension of next observation
            prepoststep_pdaf                ! User supplied pre/poststep routine
       EXTERNAL :: cvt_ens_pdaf, &          ! Apply control vector transform matrix to control vector
            cvt_adj_ens_pdaf, &             ! Apply adjoint control vector transform matrix
            cvt_pdaf, &                     ! Apply control vector transform matrix to control vector
            cvt_adj_pdaf, &                 ! Apply adjoint control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
       EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
            init_dim_l_pdaf, &              ! Init state dimension for local ana. domain
            g2l_state_pdaf, &               ! Get state on local ana. domain from full state
            l2g_state_pdaf, &               ! Init full state from local state
            init_dim_obs_f_pdaf, &          ! Initialize dimension of full observation vector
            obs_op_f_pdaf, &                ! Full observation operator
            init_dim_obs_l_pdaf             ! Initialize local dimimension of obs. vector
     END SUBROUTINE PDAFomi_assimilate_hyb3dvar_lestkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_assimilate_hyb3dvar_lestkf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, prodRinvA_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, next_observation_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
            distribute_state_pdaf, &        ! Routine to distribute a state vector
            next_observation_pdaf, &        ! Provide time step, time and dimension of next observation
            prepoststep_pdaf                ! User supplied pre/poststep routine
       EXTERNAL :: cvt_ens_pdaf, &          ! Apply control vector transform matrix to control vector
            cvt_adj_ens_pdaf, &             ! Apply adjoint control vector transform matrix
            cvt_pdaf, &                     ! Apply control vector transform matrix to control vector
            cvt_adj_pdaf, &                 ! Apply adjoint control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
       EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
            prodRinvA_pdaf, &               ! Provide product R^-1 A
            init_dim_l_pdaf, &              ! Init state dimension for local ana. domain
            g2l_state_pdaf, &               ! Get state on local ana. domain from full state
            l2g_state_pdaf, &               ! Init full state from local state
            init_dim_obs_f_pdaf, &          ! Initialize dimension of full observation vector
            obs_op_f_pdaf, &                ! Full observation operator
            init_dim_obs_l_pdaf, &          ! Initialize local dimimension of obs. vector
            prodRinvA_l_pdaf                ! Provide product R^-1 A with localization
     END SUBROUTINE PDAFomi_assimilate_hyb3dvar_lestkf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_3dvar(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
          cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
            prepoststep_pdaf                ! User supplied pre/poststep routine
       EXTERNAL :: init_dim_obs_pdaf, &     ! Initialize dimension of observation vector
            obs_op_pdaf, &                  ! Observation operator
            cvt_pdaf, &                     ! Apply control vector transform matrix to control vector
            cvt_adj_pdaf, &                 ! Apply adjoint control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
     END SUBROUTINE PDAFomi_put_state_3dvar
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_3dvar_nondiagR(collect_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, prodRinvA_pdaf, &
          cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
            prepoststep_pdaf                ! User supplied pre/poststep routine
       EXTERNAL :: init_dim_obs_pdaf, &     ! Initialize dimension of observation vector
            obs_op_pdaf, &                  ! Observation operator
            prodRinvA_pdaf, &               ! Provide product R^-1 A
            cvt_pdaf, &                     ! Apply control vector transform matrix to control vector
            cvt_adj_pdaf, &                 ! Apply adjoint control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
     END SUBROUTINE PDAFomi_put_state_3dvar_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_en3dvar_estkf(collect_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          prepoststep_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
            prepoststep_pdaf                ! User supplied pre/poststep routine
       EXTERNAL :: init_dim_obs_pdaf, &     ! Initialize dimension of observation vector
            obs_op_pdaf, &                  ! Observation operator
            cvt_ens_pdaf, &                 ! Apply control vector transform matrix to control vector
            cvt_adj_ens_pdaf, &             ! Apply adjoint control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
     END SUBROUTINE PDAFomi_put_state_en3dvar_estkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_en3dvar_estkf_nondiagR(collect_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, prodRinvA_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          prepoststep_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
            prepoststep_pdaf                ! User supplied pre/poststep routine
       EXTERNAL :: init_dim_obs_pdaf, &     ! Initialize dimension of observation vector
            obs_op_pdaf, &                  ! Observation operator
            prodRinvA_pdaf, &               ! Provide product R^-1 A
            cvt_ens_pdaf, &                 ! Apply control vector transform matrix to control vector
            cvt_adj_ens_pdaf, &             ! Apply adjoint control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
     END SUBROUTINE PDAFomi_put_state_en3dvar_estkf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_en3dvar_lestkf(collect_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
       prepoststep_pdaf                ! User supplied pre/poststep routine
       EXTERNAL :: cvt_ens_pdaf, &          ! Apply control vector transform matrix to control vector
            cvt_adj_ens_pdaf, &             ! Apply adjoint control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
       EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
            init_dim_l_pdaf, &              ! Init state dimension for local ana. domain
            g2l_state_pdaf, &               ! Get state on local ana. domain from full state
            l2g_state_pdaf, &               ! Init full state from local state
            init_dim_obs_f_pdaf, &          ! Initialize dimension of full observation vector
            obs_op_f_pdaf, &                ! Full observation operator
            init_dim_obs_l_pdaf             ! Initialize local dimimension of obs. vector
     END SUBROUTINE PDAFomi_put_state_en3dvar_lestkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_en3dvar_lestkf_nondiagR(collect_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, prodRinvA_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag    ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
       prepoststep_pdaf                     ! User supplied pre/poststep routine
       EXTERNAL :: cvt_ens_pdaf, &          ! Apply control vector transform matrix to control vector
            cvt_adj_ens_pdaf, &             ! Apply adjoint control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
       EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
            init_dim_l_pdaf, &              ! Init state dimension for local ana. domain
            prodRinvA_pdaf, &               ! Provide product R^-1 A
            g2l_state_pdaf, &               ! Get state on local ana. domain from full state
            l2g_state_pdaf, &               ! Init full state from local state
            init_dim_obs_f_pdaf, &          ! Initialize dimension of full observation vector
            obs_op_f_pdaf, &                ! Full observation operator
            init_dim_obs_l_pdaf, &          ! Initialize local dimimension of obs. vector
            prodRinvA_l_pdaf                ! Provide product R^-1 A with localization
     END SUBROUTINE PDAFomi_put_state_en3dvar_lestkf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_hyb3dvar_estkf(collect_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
          obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
            prepoststep_pdaf                ! User supplied pre/poststep routine
       EXTERNAL :: init_dim_obs_pdaf, &     ! Initialize dimension of observation vector
            obs_op_pdaf, &                  ! Observation operator
            cvt_pdaf, &                     ! Apply control vector transform matrix to control vector
            cvt_adj_pdaf, &                 ! Apply adjoint control vector transform matrix
            cvt_ens_pdaf, &                 ! Apply ensemble control vector transform matrix to control vector
            cvt_adj_ens_pdaf, &             ! Apply adjoint ensemble control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
     END SUBROUTINE PDAFomi_put_state_hyb3dvar_estkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_hyb3dvar_estkf_nondiagR(collect_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, prodRinvA_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
          obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
            prepoststep_pdaf                ! User supplied pre/poststep routine
       EXTERNAL :: init_dim_obs_pdaf, &     ! Initialize dimension of observation vector
            obs_op_pdaf, &                  ! Observation operator
            prodRinvA_pdaf, &               ! Provide product R^-1 A
            cvt_pdaf, &                     ! Apply control vector transform matrix to control vector
            cvt_adj_pdaf, &                 ! Apply adjoint control vector transform matrix
            cvt_ens_pdaf, &                 ! Apply control vector transform matrix to control vector
            cvt_adj_ens_pdaf, &             ! Apply adjoint control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
     END SUBROUTINE PDAFomi_put_state_hyb3dvar_estkf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_hyb3dvar_lestkf(collect_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
            prepoststep_pdaf                ! User supplied pre/poststep routine
       EXTERNAL :: cvt_ens_pdaf, &          ! Apply control vector transform matrix to control vector
            cvt_adj_ens_pdaf, &             ! Apply adjoint control vector transform matrix
            cvt_pdaf, &                     ! Apply control vector transform matrix to control vector
            cvt_adj_pdaf, &                 ! Apply adjoint control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
       EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
            init_dim_l_pdaf, &              ! Init state dimension for local ana. domain
            g2l_state_pdaf, &               ! Get state on local ana. domain from full state
            l2g_state_pdaf, &               ! Init full state from local state
            init_dim_obs_f_pdaf, &          ! Initialize dimension of full observation vector
            obs_op_f_pdaf, &                ! Full observation operator
            init_dim_obs_l_pdaf             ! Initialize local dimimension of obs. vector
     END SUBROUTINE PDAFomi_put_state_hyb3dvar_lestkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFomi_put_state_hyb3dvar_lestkf_nondiagR(collect_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, prodRinvA_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
       prepoststep_pdaf                ! User supplied pre/poststep routine
       EXTERNAL :: cvt_ens_pdaf, &          ! Apply control vector transform matrix to control vector
            cvt_adj_ens_pdaf, &             ! Apply adjoint control vector transform matrix
            cvt_pdaf, &                     ! Apply control vector transform matrix to control vector
            cvt_adj_pdaf, &                 ! Apply adjoint control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
       EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
            init_dim_l_pdaf, &              ! Init state dimension for local ana. domain
            prodRinvA_pdaf, &               ! Provide product R^-1 A
            g2l_state_pdaf, &               ! Get state on local ana. domain from full state
            l2g_state_pdaf, &               ! Init full state from local state
            init_dim_obs_f_pdaf, &          ! Initialize dimension of full observation vector
            obs_op_f_pdaf, &                ! Full observation operator
            init_dim_obs_l_pdaf, &          ! Initialize local dimimension of obs. vector
            prodRinvA_l_pdaf                ! Provide product R^-1 A with localization
     END SUBROUTINE PDAFomi_put_state_hyb3dvar_lestkf_nondiagR
  END INTERFACE

END MODULE PDAFomi_interfaces
