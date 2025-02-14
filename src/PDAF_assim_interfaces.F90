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
!! * 2012-05 - Lars Nerger - Initial code
!! * 2025-02 - Lars Nerger - Split for better coding overview
!! * Other revisions - see repository log
MODULE PDAF_assim_interfaces

  INTERFACE 
     SUBROUTINE PDAF_init(filtertype, subtype, stepnull, param_int, dim_pint, &
          param_real, dim_preal, COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, in_filterpe, U_init_ens, in_screen, &
          flag)
       INTEGER, INTENT(in) :: filtertype     ! Type of filter
       INTEGER, INTENT(in) :: subtype        ! Sub-type of filter
       INTEGER, INTENT(in) :: stepnull       ! Initial time step of assimilation
       INTEGER, INTENT(in) :: dim_pint       ! Number of integer parameters
       INTEGER, INTENT(inout) :: param_int(dim_pint) ! Integer parameter array
       INTEGER, INTENT(in) :: dim_preal      ! Number of real parameter 
       REAL, INTENT(inout) :: param_real(dim_preal) ! Real parameter array
       INTEGER, INTENT(in) :: COMM_model     ! Model communicator
       INTEGER, INTENT(in) :: COMM_couple    ! Coupling communicator
       INTEGER, INTENT(in) :: COMM_filter    ! Filter communicator
       INTEGER, INTENT(in) :: task_id        ! Id of my ensemble task
       INTEGER, INTENT(in) :: n_modeltasks   ! Number of parallel model tasks
       LOGICAL, INTENT(in) :: in_filterpe    ! Is my PE a filter-PE?
       INTEGER, INTENT(in) :: in_screen      ! Control screen output:
       INTEGER, INTENT(out):: flag           ! Status flag, 0: no error, error codes:
       EXTERNAL :: U_init_ens  ! User-supplied routine for ensemble initialization
     END SUBROUTINE PDAF_init
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_get_state_si(nsteps, time, doexit, flag)
       INTEGER, INTENT(inout) :: nsteps  ! Flag and number of time steps
       REAL, INTENT(out)      :: time    ! current model time
       INTEGER, INTENT(inout) :: doexit  ! Whether to exit from forecasts
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_get_state_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_get_state(steps, time, doexit, U_next_observation, U_distribute_state, &
          U_prepoststep, flag)
       INTEGER, INTENT(inout) :: steps   ! Flag and number of time steps
       REAL, INTENT(out)      :: time    ! current model time
       INTEGER, INTENT(inout) :: doexit  ! Whether to exit from forecasts
       INTEGER, INTENT(inout) :: flag    ! Status flag
       EXTERNAL :: U_next_observation, & ! Provide time step and time of next observation
            U_distribute_state, &        ! Routine to distribute a state vector
            U_prepoststep                ! User supplied pre/poststep routine
     END SUBROUTINE PDAF_get_state
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_seek(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_prepoststep, U_prodRinvA, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_obs_op, &             ! Observation operator
            U_init_obs, &           ! Initialize observation vector
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_prodRinvA             ! Provide product R^-1 HV
     END SUBROUTINE PDAF_put_state_seek
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_seek_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_seek_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_seek(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_prodRinvA, &
          U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &   ! Routine to distribute a state vector
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_obs_op, &             ! Observation operator
            U_init_obs, &           ! Initialize observation vector
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_prodRinvA, &          ! Provide product R^-1 HV
            U_next_observation      ! Provide time step and time of next observation
     END SUBROUTINE PDAF_assimilate_seek
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_seek_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_seek_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_seik(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_prepoststep, U_prodRinvA, U_init_obsvar, flag)
       INTEGER, INTENT(out) :: flag   ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &      ! Initialize dimension of observation vector
            U_obs_op, &            ! Observation operator
            U_init_obsvar, &       ! Initialize mean observation error variance
            U_init_obs, &          ! Initialize observation vector
            U_prepoststep, &       ! User supplied pre/poststep routine
            U_prodRinvA            ! Provide product R^-1 A
     END SUBROUTINE PDAF_put_state_seik
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_seik_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_seik_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_seik(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_prodRinvA, &
          U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &   ! Routine to distribute a state vector
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_obs_op, &             ! Observation operator
            U_init_obs, &           ! Initialize observation vector
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_prodRinvA, &          ! Provide product R^-1 HV
            U_next_observation      ! Provide time step and time of next observation
     END SUBROUTINE PDAF_assimilate_seik
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_seik_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_seik_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_enkf(U_collect_state, U_init_dim_obs, U_obs_op,  &
          U_init_obs, U_prepoststep, U_add_obs_err, U_init_obs_covar, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_obs_op, &             ! Observation operator
            U_init_obs, &           ! Initialize observation vector
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_add_obs_err, &        ! Add obs error covariance R to HPH in EnKF
            U_init_obs_covar        ! Initialize obs. error cov. matrix R in EnKF
     END SUBROUTINE PDAF_put_state_enkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_enkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_enkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_enkf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_add_obs_error, &
          U_init_obs_covar, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_distribute_state, &     ! Routine to distribute a state vector
            U_init_dim_obs, &      ! Initialize dimension of observation vector
            U_obs_op, &            ! Observation operator
            U_init_obs_covar, &    ! Initialize obs. error cov. matrix R in EnKF
            U_init_obs, &          ! Initialize observation vector
            U_prepoststep, &       ! User supplied pre/poststep routine
            U_add_obs_error, &     ! Add obs error covariance R to HPH in EnKF
            U_next_observation     ! Provide time step and time of next observation
     END SUBROUTINE PDAF_assimilate_enkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_enkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_enkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_lseik(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
          U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
          U_init_obsvar, U_init_obsvar_l, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_obs_op, &             ! Observation operator
            U_init_n_domains_p, &   ! Provide number of local analysis domains
            U_init_dim_l, &         ! Init state dimension for local ana. domain
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_init_dim_obs_l, &     ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs, &           ! Initialize PE-local observation vector
            U_init_obs_l, &         ! Init. observation vector on local analysis domain
            U_init_obsvar, &        ! Initialize mean observation error variance
            U_init_obsvar_l, &      ! Initialize local mean observation error variance
            U_g2l_state, &          ! Get state on local ana. domain from full state
            U_l2g_state, &          ! Init full state from state on local analysis domain
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prepoststep           ! User supplied pre/poststep routine
     END SUBROUTINE PDAF_put_state_lseik
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_lseik_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_lseik_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_lseik(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
          U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
          U_g2l_state, U_l2g_state, U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
          U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &      ! Routine to distribute a state vector
            U_obs_op, &             ! Observation operator
            U_init_n_domains_p, &   ! Provide number of local analysis domains
            U_init_dim_l, &         ! Init state dimension for local ana. domain
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_init_dim_obs_l, &     ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs, &           ! Initialize PE-local observation vector
            U_init_obs_l, &         ! Init. observation vector on local analysis domain
            U_init_obsvar, &        ! Initialize mean observation error variance
            U_init_obsvar_l, &      ! Initialize local mean observation error variance
            U_g2l_state, &          ! Get state on local ana. domain from full state
            U_l2g_state, &          ! Init full state from state on local analysis domain
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_next_observation      ! Provide time step and time of next observation
     END SUBROUTINE PDAF_assimilate_lseik
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_lseik_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_lseik_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_etkf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_prepoststep, U_prodRinvA, U_init_obsvar, flag)
       INTEGER, INTENT(out) :: flag   ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &      ! Initialize dimension of observation vector
            U_obs_op, &            ! Observation operator
            U_init_obsvar, &       ! Initialize mean observation error variance
            U_init_obs, &          ! Initialize observation vector
            U_prepoststep, &       ! User supplied pre/poststep routine
            U_prodRinvA            ! Provide product R^-1 A
     END SUBROUTINE PDAF_put_state_etkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_etkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_etkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_etkf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_prodRinvA, &
          U_init_obsvar, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &   ! Routine to distribute a state vector
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_obs_op, &             ! Observation operator
            U_init_obsvar, &        ! Initialize mean observation error variance
            U_init_obs, &           ! Initialize observation vector
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_prodRinvA, &          ! Provide product R^-1 HV
            U_next_observation      ! Provide time step and time of next observation
     END SUBROUTINE PDAF_assimilate_etkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_etkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_etkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_letkf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
          U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
          U_init_obsvar, U_init_obsvar_l, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_obs_op, &             ! Observation operator
            U_init_n_domains_p, &   ! Provide number of local analysis domains
            U_init_dim_l, &         ! Init state dimension for local ana. domain
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_init_dim_obs_l, &     ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs, &           ! Initialize PE-local observation vector
            U_init_obs_l, &         ! Init. observation vector on local analysis domain
            U_init_obsvar, &        ! Initialize mean observation error variance
            U_init_obsvar_l, &      ! Initialize local mean observation error variance
            U_g2l_state, &          ! Get state on local ana. domain from full state
            U_l2g_state, &          ! Init full state from state on local analysis domain
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prepoststep           ! User supplied pre/poststep routine
     END SUBROUTINE PDAF_put_state_letkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_letkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_letkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_letkf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
          U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
          U_g2l_state, U_l2g_state, U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
          U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &      ! Routine to distribute a state vector
            U_obs_op, &             ! Observation operator
            U_init_n_domains_p, &   ! Provide number of local analysis domains
            U_init_dim_l, &         ! Init state dimension for local ana. domain
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_init_dim_obs_l, &     ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs, &           ! Initialize PE-local observation vector
            U_init_obs_l, &         ! Init. observation vector on local analysis domain
            U_init_obsvar, &        ! Initialize mean observation error variance
            U_init_obsvar_l, &      ! Initialize local mean observation error variance
            U_g2l_state, &          ! Get state on local ana. domain from full state
            U_l2g_state, &          ! Init full state from state on local analysis domain
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_next_observation      ! Provide time step and time of next observation
     END SUBROUTINE PDAF_assimilate_letkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_letkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_letkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_estkf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_prepoststep, U_prodRinvA, U_init_obsvar, flag)
       INTEGER, INTENT(out) :: flag   ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &      ! Initialize dimension of observation vector
            U_obs_op, &            ! Observation operator
            U_init_obsvar, &       ! Initialize mean observation error variance
            U_init_obs, &          ! Initialize observation vector
            U_prepoststep, &       ! User supplied pre/poststep routine
            U_prodRinvA            ! Provide product R^-1 A
     END SUBROUTINE PDAF_put_state_estkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_estkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_estkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_estkf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_prodRinvA, &
          U_init_obsvar, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &   ! Routine to distribute a state vector
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_init_obsvar, &        ! Initialize mean observation error variance
            U_obs_op, &             ! Observation operator
            U_init_obs, &           ! Initialize observation vector
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_prodRinvA, &          ! Provide product R^-1 HV
            U_next_observation      ! Provide time step and time of next observation
     END SUBROUTINE PDAF_assimilate_estkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_estkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_estkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_lestkf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
          U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
          U_init_obsvar, U_init_obsvar_l, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_obs_op, &             ! Observation operator
            U_init_n_domains_p, &   ! Provide number of local analysis domains
            U_init_dim_l, &         ! Init state dimension for local ana. domain
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_init_dim_obs_l, &     ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs, &           ! Initialize PE-local observation vector
            U_init_obs_l, &         ! Init. observation vector on local analysis domain
            U_init_obsvar, &        ! Initialize mean observation error variance
            U_init_obsvar_l, &      ! Initialize local mean observation error variance
            U_g2l_state, &          ! Get state on local ana. domain from full state
            U_l2g_state, &          ! Init full state from state on local analysis domain
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prepoststep           ! User supplied pre/poststep routine
     END SUBROUTINE PDAF_put_state_lestkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_lestkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_lestkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_lestkf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
          U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
          U_g2l_state, U_l2g_state, U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
          U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &      ! Routine to distribute a state vector
            U_obs_op, &             ! Observation operator
            U_init_n_domains_p, &   ! Provide number of local analysis domains
            U_init_dim_l, &         ! Init state dimension for local ana. domain
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_init_dim_obs_l, &     ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs, &           ! Initialize PE-local observation vector
            U_init_obs_l, &         ! Init. observation vector on local analysis domain
            U_init_obsvar, &        ! Initialize mean observation error variance
            U_init_obsvar_l, &      ! Initialize local mean observation error variance
            U_g2l_state, &          ! Get state on local ana. domain from full state
            U_l2g_state, &          ! Init full state from state on local analysis domain
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_next_observation      ! Provide time step and time of next observation
     END SUBROUTINE PDAF_assimilate_lestkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_lestkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_lestkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_lenkf(U_collect_state, U_init_dim_obs, U_obs_op,  &
          U_init_obs, U_prepoststep, U_localize, U_add_obs_err, U_init_obs_covar, &
          flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_obs_op, &             ! Observation operator
            U_init_obs, &           ! Initialize observation vector
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_localize, &           ! Apply localization to HP and HPH^T
            U_add_obs_err, &        ! Add obs error covariance R to HPH in EnKF
            U_init_obs_covar        ! Initialize obs. error cov. matrix R in EnKF
     END SUBROUTINE PDAF_put_state_lenkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_lenkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_lenkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_lenkf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_localize, &
          U_add_obs_error, U_init_obs_covar, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_distribute_state, &     ! Routine to distribute a state vector
            U_init_dim_obs, &      ! Initialize dimension of observation vector
            U_obs_op, &            ! Observation operator
            U_init_obs_covar, &    ! Initialize obs. error cov. matrix R in EnKF
            U_init_obs, &          ! Initialize observation vector
            U_prepoststep, &       ! User supplied pre/poststep routine
            U_localize, &          ! Apply localization to HP and HPH^T
            U_add_obs_error, &     ! Add obs error covariance R to HPH in EnKF
            U_next_observation     ! Provide time step and time of next observation
     END SUBROUTINE PDAF_assimilate_lenkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_lenkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_lenkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_netf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_prepoststep, U_likelihood, flag)
       INTEGER, INTENT(out) :: flag   ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &      ! Initialize dimension of observation vector
            U_obs_op, &            ! Observation operator
            U_init_obs, &          ! Initialize observation vector
            U_prepoststep, &       ! User supplied pre/poststep routine
            U_likelihood           ! Compute observation likelihood for an ensemble member
     END SUBROUTINE PDAF_put_state_netf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_netf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_netf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_netf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, &
          U_likelihood, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &   ! Routine to distribute a state vector
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_obs_op, &             ! Observation operator
            U_init_obs, &           ! Initialize observation vector
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_likelihood, &         ! Compute observation likelihood for an ensemble member
            U_next_observation      ! Provide time step and time of next observation
     END SUBROUTINE PDAF_assimilate_netf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_netf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_netf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_lnetf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_init_obs_l, U_prepoststep, U_likelihood_l, U_init_n_domains_p, &
          U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
          outflag)
       INTEGER, INTENT(out) :: outflag  ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_obs_op, &             ! Observation operator
            U_init_n_domains_p, &   ! Provide number of local analysis domains
            U_init_dim_l, &         ! Init state dimension for local ana. domain
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_init_dim_obs_l, &     ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs, &           ! Initialize PE-local observation vector
            U_init_obs_l, &         ! Init. observation vector on local analysis domain
            U_g2l_state, &          ! Get state on local ana. domain from full state
            U_l2g_state, &          ! Init full state from state on local analysis domain
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_likelihood_l, &       ! Compute observation likelihood for an ensemble member
            U_prepoststep           ! User supplied pre/poststep routine
     END SUBROUTINE PDAF_put_state_lnetf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_lnetf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_lnetf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_lnetf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs_l, U_prepoststep, &
          U_likelihood_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
          U_g2l_state, U_l2g_state, U_g2l_obs, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &      ! Routine to distribute a state vector
            U_obs_op, &             ! Observation operator
            U_init_n_domains_p, &   ! Provide number of local analysis domains
            U_init_dim_l, &         ! Init state dimension for local ana. domain
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_init_dim_obs_l, &     ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs_l, &         ! Init. observation vector on local analysis domain
            U_g2l_state, &          ! Get state on local ana. domain from full state
            U_l2g_state, &          ! Init full state from state on local analysis domain
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_likelihood_l, &       ! Compute observation likelihood for an ensemble member
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_next_observation      ! Provide time step and time of next observation
     END SUBROUTINE PDAF_assimilate_lnetf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_lnetf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_lnetf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_lknetf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_prodRinvA_hyb_l, &
          U_init_n_domains_p, &
          U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
          U_init_obsvar, U_init_obsvar_l, U_likelihood_l, U_likelihood_hyb_l, outflag)
       INTEGER, INTENT(out) :: outflag  ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_obs_op, &             ! Observation operator
            U_init_n_domains_p, &   ! Provide number of local analysis domains
            U_init_dim_l, &         ! Init state dimension for local ana. domain
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_init_dim_obs_l, &     ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs, &           ! Initialize PE-local observation vector
            U_init_obs_l, &         ! Init. observation vector on local analysis domain
            U_init_obsvar, &        ! Initialize mean observation error variance
            U_init_obsvar_l, &      ! Initialize local mean observation error variance
            U_g2l_state, &          ! Get state on local ana. domain from full state
            U_l2g_state, &          ! Init full state from state on local analysis domain
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prodRinvA_hyb_l, &    ! Provide product R^-1 A on local analysis domain with hybrid weight
            U_likelihood_l, &       ! Compute likelihood
            U_likelihood_hyb_l, &   ! Compute likelihood with hybrid weight
            U_prepoststep           ! User supplied pre/poststep routine
     END SUBROUTINE PDAF_put_state_lknetf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_lknetf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_lknetf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_lknetf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
          U_prodRinvA_l, U_prodRinvA_hyb_l, U_init_n_domains_p, U_init_dim_l, &
          U_init_dim_obs_l, &
          U_g2l_state, U_l2g_state, U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
          U_likelihood_l, U_likelihood_hyb_l, &
          U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_obs_op, &             ! Observation operator
            U_init_n_domains_p, &   ! Provide number of local analysis domains
            U_init_dim_l, &         ! Init state dimension for local ana. domain
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_init_dim_obs_l, &     ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs, &           ! Initialize PE-local observation vector
            U_init_obs_l, &         ! Init. observation vector on local analysis domain
            U_init_obsvar, &        ! Initialize mean observation error variance
            U_init_obsvar_l, &      ! Initialize local mean observation error variance
            U_g2l_state, &          ! Get state on local ana. domain from full state
            U_l2g_state, &          ! Init full state from state on local analysis domain
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prodRinvA_hyb_l, &    ! Provide product R^-1 A on local analysis domain with hybrid weight
            U_likelihood_l, &       ! Compute likelihood
            U_likelihood_hyb_l, &   ! Compute likelihood with hybrid weight
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_next_observation, &   ! Routine to provide time step, time and dimension
                                    !   of next observation
            U_distribute_state      ! Routine to distribute a state vector
     END SUBROUTINE PDAF_assimilate_lknetf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_lknetf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_lknetf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_generate_obs(U_collect_state, U_init_dim_obs_f, U_obs_op_f, &
          U_get_obs_f, U_init_obserr_f, U_prepoststep, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_init_dim_obs_f, &        ! Initialize dimension of observation vector
            U_obs_op_f, &              ! Observation operator
            U_get_obs_f, &             ! Provide observation vector to user
            U_init_obserr_f, &         ! Initialize vector of observation errors
            U_prepoststep              ! User supplied pre/poststep routine
     END SUBROUTINE PDAF_put_state_generate_obs
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_generate_obs_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_generate_obs_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_generate_obs(U_collect_state, U_distribute_state, &
          U_init_dim_obs_f, U_obs_op_f, U_get_obs_f, U_init_obserr_f, U_prepoststep, &
          U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &      ! Routine to distribute a state vector
            U_init_dim_obs_f, &        ! Initialize dimension of observation vector
            U_obs_op_f, &              ! Observation operator
            U_get_obs_f, &             ! Provide observation vector to user
            U_init_obserr_f, &         ! Initialize vector of observation error standard deviations
            U_prepoststep, &           ! User supplied pre/poststep routine
            U_next_observation         ! Provide time step and time of next observation
     END SUBROUTINE PDAF_generate_obs
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_generate_obs_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_generate_obs_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_pf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_prepoststep, U_likelihood, flag)
       INTEGER, INTENT(out) :: flag   ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &      ! Initialize dimension of observation vector
            U_obs_op, &            ! Observation operator
            U_init_obs, &          ! Initialize observation vector
            U_prepoststep, &       ! User supplied pre/poststep routine
            U_likelihood           ! Compute observation likelihood for an ensemble member
     END SUBROUTINE PDAF_put_state_pf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_pf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_pf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_pf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, &
          U_likelihood, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &   ! Routine to distribute a state vector
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_obs_op, &             ! Observation operator
            U_init_obs, &           ! Initialize observation vector
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_likelihood, &         ! Compute observation likelihood for an ensemble member
            U_next_observation      ! Provide time step and time of next observation
     END SUBROUTINE PDAF_assimilate_pf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_pf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_pf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_prepost_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_prepost_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_prepost(U_collect_state, U_prepoststep, flag)
       INTEGER, INTENT(out) :: flag   ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_prepoststep             ! User supplied pre/poststep routine
     END SUBROUTINE PDAF_put_state_prepost
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_prepost(U_collect_state, U_distribute_state, &
          U_prepoststep, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &      ! Routine to distribute a state vector
            U_prepoststep, &           ! User supplied pre/poststep routine
            U_next_observation         ! Provide time step and time of next observation
     END SUBROUTINE PDAF_assimilate_prepost
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_prepost_si(flag)
       INTEGER, INTENT(inout) :: flag  ! Status flag
     END SUBROUTINE PDAF_assimilate_prepost_si
  END INTERFACE

! 3D-Var

  INTERFACE
     SUBROUTINE PDAF_put_state_3dvar(U_collect_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
          U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
          U_prepoststep, outflag)
       INTEGER, INTENT(out) :: outflag  ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &      ! Initialize dimension of observation vector
            U_obs_op, &            ! Observation operator
            U_init_obsvar, &       ! Initialize mean observation error variance
            U_init_obs, &          ! Initialize observation vector
            U_prepoststep, &       ! User supplied pre/poststep routine
            U_prodRinvA, &         ! Provide product R^-1 A
            U_cvt, &               ! Apply control vector transform matrix to control vector
            U_cvt_adj, &           ! Apply adjoint control vector transform matrix
            U_obs_op_lin, &        ! Linearized observation operator
            U_obs_op_adj           ! Adjoint observation operator
     END SUBROUTINE PDAF_put_state_3dvar
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_en3dvar_estkf(U_collect_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
          U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
          U_init_obsvar, U_prepoststep, outflag)
       INTEGER, INTENT(out) :: outflag  ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &      ! Initialize dimension of observation vector
            U_obs_op, &            ! Observation operator
            U_init_obs, &          ! Initialize observation vector
            U_prepoststep, &       ! User supplied pre/poststep routine
            U_prodRinvA, &         ! Provide product R^-1 A
            U_cvt_ens, &           ! Apply control vector transform matrix (ensemble)
            U_cvt_adj_ens, &       ! Apply adjoint control vector transform matrix (ensemble var)
            U_obs_op_lin, &        ! Linearized observation operator
            U_obs_op_adj           ! Adjoint observation operator
       EXTERNAL :: U_init_obsvar   ! Initialize mean observation error variance
     END SUBROUTINE PDAF_put_state_en3dvar_estkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_en3dvar_lestkf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_prodRinvA, &
          U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
          U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
          U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, &
          U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
          U_prepoststep, outflag)
       INTEGER, INTENT(out) :: outflag  ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &      ! Initialize dimension of observation vector
            U_obs_op, &            ! Observation operator
            U_init_obs, &          ! Initialize observation vector
            U_prepoststep, &       ! User supplied pre/poststep routine
            U_prodRinvA, &         ! Provide product R^-1 A
            U_cvt_ens, &           ! Apply control vector transform matrix (ensemble)
            U_cvt_adj_ens, &       ! Apply adjoint control vector transform matrix (ensemble var)
            U_obs_op_lin, &        ! Linearized observation operator
            U_obs_op_adj           ! Adjoint observation operator
       EXTERNAL :: U_obs_op_f, &    ! Observation operator
            U_init_n_domains_p, &   ! Provide number of local analysis domains
            U_init_dim_l, &         ! Init state dimension for local ana. domain
            U_init_dim_obs_f, &     ! Initialize dimension of observation vector
            U_init_dim_obs_l, &     ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs_f, &         ! Initialize PE-local observation vector
            U_init_obs_l, &         ! Init. observation vector on local analysis domain
            U_init_obsvar, &        ! Initialize mean observation error variance
            U_init_obsvar_l, &      ! Initialize local mean observation error variance
            U_g2l_state, &          ! Get state on local ana. domain from full state
            U_l2g_state, &          ! Init full state from state on local analysis domain
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l           ! Provide product R^-1 A on local analysis domain
     END SUBROUTINE PDAF_put_state_en3dvar_lestkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_hyb3dvar_estkf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_prodRinvA, &
          U_cvt, U_cvt_adj, U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
          U_init_obsvar, U_prepoststep, outflag)
       INTEGER, INTENT(out) :: outflag  ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &      ! Initialize dimension of observation vector
            U_obs_op, &            ! Observation operator
            U_init_obs, &          ! Initialize observation vector
            U_prepoststep, &       ! User supplied pre/poststep routine
            U_prodRinvA, &         ! Provide product R^-1 A
            U_cvt_ens, &           ! Apply control vector transform matrix (ensemble)
            U_cvt_adj_ens, &       ! Apply adjoint control vector transform matrix (ensemble var)
            U_cvt, &               ! Apply control vector transform matrix to control vector
            U_cvt_adj, &           ! Apply adjoint control vector transform matrix
            U_obs_op_lin, &        ! Linearized observation operator
            U_obs_op_adj           ! Adjoint observation operator
       EXTERNAL :: U_init_obsvar   ! Initialize mean observation error variance
     END SUBROUTINE PDAF_put_state_hyb3dvar_estkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_hyb3dvar_lestkf(U_collect_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
          U_cvt_ens, U_cvt_adj_ens, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
          U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
          U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, &
          U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
          U_prepoststep, outflag)
       INTEGER, INTENT(out) :: outflag  ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &      ! Initialize dimension of observation vector
            U_obs_op, &            ! Observation operator
            U_init_obs, &          ! Initialize observation vector
            U_prepoststep, &       ! User supplied pre/poststep routine
            U_prodRinvA, &         ! Provide product R^-1 A
            U_cvt_ens, &           ! Apply control vector transform matrix (ensemble)
            U_cvt_adj_ens, &       ! Apply adjoint control vector transform matrix (ensemble var)
            U_cvt, &               ! Apply control vector transform matrix to control vector
            U_cvt_adj, &           ! Apply adjoint control vector transform matrix
            U_obs_op_lin, &        ! Linearized observation operator
            U_obs_op_adj           ! Adjoint observation operator
       EXTERNAL :: U_obs_op_f, &   ! Observation operator
            U_init_n_domains_p, &  ! Provide number of local analysis domains
            U_init_dim_l, &        ! Init state dimension for local ana. domain
            U_init_dim_obs_f, &    ! Initialize dimension of observation vector
            U_init_dim_obs_l, &    ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs_f, &        ! Initialize PE-local observation vector
            U_init_obs_l, &        ! Init. observation vector on local analysis domain
            U_init_obsvar, &       ! Initialize mean observation error variance
            U_init_obsvar_l, &     ! Initialize local mean observation error variance
            U_g2l_state, &         ! Get state on local ana. domain from full state
            U_l2g_state, &         ! Init full state from state on local analysis domain
            U_g2l_obs, &           ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l          ! Provide product R^-1 A on local analysis domain
     END SUBROUTINE PDAF_put_state_hyb3dvar_lestkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_3dvar(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
          U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
          U_prepoststep, U_next_observation, outflag)
       INTEGER, INTENT(out) :: outflag  ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &         ! Initialize dimension of observation vector
            U_obs_op, &               ! Observation operator
            U_init_obsvar, &          ! Initialize mean observation error variance
            U_init_obs, &             ! Initialize observation vector
            U_prepoststep, &          ! User supplied pre/poststep routine
            U_prodRinvA, &            ! Provide product R^-1 A
            U_next_observation, &     ! Routine to provide time step, time and dimension
                                 !   of next observation
            U_distribute_state, &     ! Routine to distribute a state vector
            U_cvt, &                  ! Apply control vector transform matrix to control vector
            U_cvt_adj, &              ! Apply adjoint control vector transform matrix
            U_obs_op_lin, &           ! Linearized observation operator
            U_obs_op_adj              ! Adjoint observation operator
     END SUBROUTINE PDAF_assimilate_3dvar
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_en3dvar_estkf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
          U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
          U_init_obsvar, U_prepoststep, U_next_observation, outflag)
       INTEGER, INTENT(out) :: outflag  ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &         ! Initialize dimension of observation vector
            U_obs_op, &               ! Observation operator
            U_init_obsvar, &          ! Initialize mean observation error variance
            U_init_obs, &             ! Initialize observation vector
            U_prepoststep, &          ! User supplied pre/poststep routine
            U_prodRinvA, &            ! Provide product R^-1 A
            U_next_observation, &     ! Routine to provide time step, time and dimension
                                 !   of next observation
            U_distribute_state, &     ! Routine to distribute a state vector
            U_cvt_ens, &              ! Apply control vector transform matrix (ensemble)
            U_cvt_adj_ens, &          ! Apply adjoint control vector transform matrix (ensemble var)
            U_cvt, &                  ! Apply control vector transform matrix to control vector
            U_cvt_adj, &              ! Apply adjoint control vector transform matrix
            U_obs_op_lin, &           ! Linearized observation operator
            U_obs_op_adj              ! Adjoint observation operator
     END SUBROUTINE PDAF_assimilate_en3dvar_estkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_en3dvar_lestkf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
          U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
          U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
          U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, &
          U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
          U_prepoststep, U_next_observation, outflag)
       INTEGER, INTENT(inout) :: outflag  ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &         ! Initialize dimension of observation vector
            U_obs_op, &               ! Observation operator
            U_init_obsvar, &          ! Initialize mean observation error variance
            U_init_obs, &             ! Initialize observation vector
            U_prepoststep, &          ! User supplied pre/poststep routine
            U_prodRinvA, &            ! Provide product R^-1 A
            U_next_observation, &     ! Routine to provide time step, time and dimension
                                 !   of next observation
            U_distribute_state, &     ! Routine to distribute a state vector
            U_cvt_ens, &              ! Apply control vector transform matrix (ensemble)
            U_cvt_adj_ens, &          ! Apply adjoint control vector transform matrix (ensemble var)
            U_cvt, &                  ! Apply control vector transform matrix to control vector
            U_cvt_adj, &              ! Apply adjoint control vector transform matrix
            U_obs_op_lin, &           ! Linearized observation operator
            U_obs_op_adj              ! Adjoint observation operator
       EXTERNAL :: U_obs_op_f, &      ! Observation operator
            U_init_n_domains_p, &     ! Provide number of local analysis domains
            U_init_dim_l, &           ! Init state dimension for local ana. domain
            U_init_dim_obs_f, &       ! Initialize dimension of observation vector
            U_init_dim_obs_l, &       ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs_f, &           ! Initialize PE-local observation vector
            U_init_obs_l, &           ! Init. observation vector on local analysis domain
            U_init_obsvar_l, &        ! Initialize local mean observation error variance
            U_g2l_state, &            ! Get state on local ana. domain from full state
            U_l2g_state, &            ! Init full state from state on local analysis domain
            U_g2l_obs, &              ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l             ! Provide product R^-1 A on local analysis domain
     END SUBROUTINE PDAF_assimilate_en3dvar_lestkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_hyb3dvar_estkf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
          U_cvt_ens, U_cvt_adj_ens, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
          U_init_obsvar, U_prepoststep, U_next_observation, outflag)
       INTEGER, INTENT(out) :: outflag  ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &         ! Initialize dimension of observation vector
            U_obs_op, &               ! Observation operator
            U_init_obsvar, &          ! Initialize mean observation error variance
            U_init_obs, &             ! Initialize observation vector
            U_prepoststep, &          ! User supplied pre/poststep routine
            U_prodRinvA, &            ! Provide product R^-1 A
            U_next_observation, &     ! Routine to provide time step, time and dimension
                                 !   of next observation
            U_distribute_state, &     ! Routine to distribute a state vector
            U_cvt_ens, &              ! Apply control vector transform matrix (ensemble)
            U_cvt_adj_ens, &          ! Apply adjoint control vector transform matrix (ensemble var)
            U_cvt, &                  ! Apply control vector transform matrix to control vector
            U_cvt_adj, &              ! Apply adjoint control vector transform matrix
            U_obs_op_lin, &           ! Linearized observation operator
            U_obs_op_adj              ! Adjoint observation operator
     END SUBROUTINE PDAF_assimilate_hyb3dvar_estkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_hyb3dvar_lestkf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
          U_cvt_ens, U_cvt_adj_ens, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
          U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
          U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, &
          U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
          U_prepoststep, U_next_observation, outflag)
       INTEGER, INTENT(out) :: outflag  ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &         ! Initialize dimension of observation vector
            U_obs_op, &               ! Observation operator
            U_init_obsvar, &          ! Initialize mean observation error variance
            U_init_obs, &             ! Initialize observation vector
            U_prepoststep, &          ! User supplied pre/poststep routine
            U_prodRinvA, &            ! Provide product R^-1 A
            U_next_observation, &     ! Routine to provide time step, time and dimension
                                 !   of next observation
            U_distribute_state, &     ! Routine to distribute a state vector
            U_cvt_ens, &              ! Apply control vector transform matrix (ensemble)
            U_cvt_adj_ens, &          ! Apply adjoint control vector transform matrix (ensemble var)
            U_cvt, &                  ! Apply control vector transform matrix to control vector
            U_cvt_adj, &              ! Apply adjoint control vector transform matrix
            U_obs_op_lin, &           ! Linearized observation operator
            U_obs_op_adj              ! Adjoint observation operator
       EXTERNAL :: U_obs_op_f, &      ! Observation operator
            U_init_n_domains_p, &     ! Provide number of local analysis domains
            U_init_dim_l, &           ! Init state dimension for local ana. domain
            U_init_dim_obs_f, &       ! Initialize dimension of observation vector
            U_init_dim_obs_l, &       ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs_f, &           ! Initialize PE-local observation vector
            U_init_obs_l, &           ! Init. observation vector on local analysis domain
            U_init_obsvar_l, &        ! Initialize local mean observation error variance
            U_g2l_state, &            ! Get state on local ana. domain from full state
            U_l2g_state, &            ! Init full state from state on local analysis domain
            U_g2l_obs, &              ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l             ! Provide product R^-1 A on local analysis domain
     END SUBROUTINE PDAF_assimilate_hyb3dvar_lestkf
  END INTERFACE

END MODULE PDAF_assim_interfaces
