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
!$Id$
!BOP
!
! !ROUTINE: PDAF_interfaces_module --- Interface definitions for PDAF
!
! !INTERFACE:
MODULE PDAF_interfaces_module

! !DESCRIPTION:
! Module providing interface definition of the PDAF routines that
! are called from the model code.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2012-05 - Lars Nerger - Initial code
! Later revisions - see svn log
!EOP
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAFlocal_interfaces

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
     SUBROUTINE PDAF_print_info(printtype)
  INTEGER, INTENT(in) :: printtype    ! Type of screen output
     END SUBROUTINE PDAF_print_info
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_deallocate()
     END SUBROUTINE PDAF_deallocate
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
          U_init_obs_l, U_prepoststep, U_likelihood_l, U_init_n_domains_p, &
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

! Other routines

  INTERFACE
     SUBROUTINE PDAF_gather_dim_obs_f(dim_obs_p, dim_obs_f)
       INTEGER, INTENT(in)  :: dim_obs_p    ! PE-local observation dimension
       INTEGER, INTENT(out) :: dim_obs_f    ! Full observation dimension
     END SUBROUTINE PDAF_gather_dim_obs_f
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_seik_TtimesA(rank, dim_col, A, B)
       INTEGER, INTENT(in) :: rank         ! Rank of initial covariance matrix
       INTEGER, INTENT(in) :: dim_col            ! Number of columns in A and B
       REAL, INTENT(in)    :: A(rank, dim_col)   ! Input matrix
       REAL, INTENT(out)   :: B(rank+1, dim_col) ! Output matrix (TA)
     END SUBROUTINE PDAF_seik_TtimesA
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_etkf_Tleft(dim_ens, dim, A)
       INTEGER, INTENT(in) :: dim_ens      ! Rank of initial covariance matrix
       INTEGER, INTENT(in) :: dim               ! Number of columns in A and B
       REAL, INTENT(inout) :: A(dim_ens, dim)   ! Input/output matrix
     END SUBROUTINE PDAF_etkf_Tleft
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_estkf_OmegaA(rank, dim_col, A, B)
       INTEGER, INTENT(in) :: rank         ! Rank of initial covariance matrix
       INTEGER, INTENT(in) :: dim_col            ! Number of columns in A and B
       REAL, INTENT(in)    :: A(rank, dim_col)   ! Input matrix
       REAL, INTENT(out)   :: B(rank+1, dim_col) ! Output matrix (TA)
     END SUBROUTINE PDAF_estkf_OmegaA
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_enkf_Omega(seed, r, dim_ens, Omega, norm, &
          otype, screen)
       INTEGER, INTENT(in) :: seed(4)  ! Seed for random number generation
       INTEGER, INTENT(in) :: r        ! Approximated rank of covar matrix
       INTEGER, INTENT(in) :: dim_ens  ! Ensemble size
       REAL, INTENT(inout) :: Omega(dim_ens,r)  ! Random matrix
       REAL, INTENT(inout) :: norm     ! Norm for ensemble transformation
       INTEGER, INTENT(in) :: otype    ! Type of Omega
       INTEGER, INTENT(in) :: screen    ! Verbosity flag
     END SUBROUTINE PDAF_enkf_Omega
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_seik_Omega(rank, Omega, Omegatype, screen)
       INTEGER, INTENT(in) :: rank      ! Approximated rank of covar matrix
       REAL, INTENT(inout) :: Omega(rank+1, rank) ! Matrix Omega
       INTEGER, INTENT(in) :: Omegatype ! Select type of Omega
       INTEGER, INTENT(in) :: screen    ! Verbosity flag
     END SUBROUTINE PDAF_seik_Omega
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_get_memberid(memberid)
       INTEGER,INTENT(inout) :: memberid ! Index in the local ensemble
     END SUBROUTINE PDAF_get_memberid
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_get_obsmemberid(memberid)
       INTEGER,INTENT(inout) :: memberid ! Index in the local observed ensemble
     END SUBROUTINE PDAF_get_obsmemberid
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_local_weight(wtype, rtype, cradius, sradius, distance, &
          nrows, ncols, A, var_obs, weight, verbose)
       INTEGER, INTENT(in) :: wtype      ! Type of weight function
       INTEGER, INTENT(in) :: rtype      ! Type of regulated weighting
       REAL, INTENT(in)    :: cradius    ! Cut-off radius
       REAL, INTENT(in)    :: sradius    ! Support radius 
       REAL, INTENT(in)    :: distance   ! Distance to observation
       INTEGER, INTENT(in) :: nrows      ! Number of rows in matrix A
       INTEGER, INTENT(in) :: ncols      ! Number of columns in matrix A
       REAL, INTENT(in) :: A(nrows, ncols) ! Input matrix
       REAL, INTENT(in)    :: var_obs    ! Observation variance
       REAL, INTENT(out)   :: weight     ! Weights
       INTEGER, INTENT(in) :: verbose    ! Verbosity flag
     END SUBROUTINE PDAF_local_weight
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_incremental(steps, U_dist_stateinc)
       INTEGER, INTENT(in) :: steps ! Time steps over which increment is distributed
       EXTERNAL :: U_dist_stateinc  ! Add state increment during integration
     END SUBROUTINE PDAF_incremental
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_add_increment(dim_p, state_p)
       INTEGER,INTENT(in) :: dim_p          ! State dimension
       REAL,INTENT(inout) :: state_p(dim_p) ! State vector
     END SUBROUTINE PDAF_add_increment
  END INTERFACE

  INTERFACE 
     SUBROUTINE PDAF_get_smootherens(sens_point, maxlag, status)
       REAL, POINTER, INTENT(out) :: sens_point(:,:,:)  ! Pointer to smoother array
       INTEGER, INTENT(out)       :: maxlag  ! Number of past timesteps processed in sens
       INTEGER, INTENT(out)       :: status  ! Status flag 
     END SUBROUTINE PDAF_get_smootherens
  END INTERFACE

  INTERFACE 
     SUBROUTINE PDAF_set_smootherens(sens_point, maxlag, status)
       REAL, POINTER, INTENT(out) :: sens_point(:,:,:)  ! Pointer to smoother array
       INTEGER, INTENT(in)        :: maxlag  ! Number of past timesteps processed in sens
       INTEGER, INTENT(out)       :: status  ! Status flag 
     END SUBROUTINE PDAF_set_smootherens
  END INTERFACE

  INTERFACE 
     SUBROUTINE PDAF_get_ensstats(skew_ptr, kurt_ptr, status)
       REAL, POINTER, INTENT(out) :: skew_ptr(:)  ! Pointer to skewness array
       REAL, POINTER, INTENT(out) :: kurt_ptr(:)  ! Pointer to kurtosis array
       INTEGER, INTENT(out)       :: status  ! Status flag 
     END SUBROUTINE PDAF_get_ensstats
  END INTERFACE

  INTERFACE 
     SUBROUTINE PDAF_reset_forget(forget_in)
       REAL, INTENT(in) :: forget_in    ! New value of forgetting factor
     END SUBROUTINE PDAF_reset_forget
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

  INTERFACE
     SUBROUTINE PDAF_set_ens_pointer(ens_point, status)
       REAL, POINTER, INTENT(out) :: ens_point(:,:)  ! Pointer to smoother array
       INTEGER, INTENT(out)       :: status  ! Status flag
     END SUBROUTINE PDAF_set_ens_pointer
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_SampleEns(dim, dim_ens, modes, svals, state, ens, verbose, flag)
       INTEGER, INTENT(in) :: dim                   ! Size of state vector
       INTEGER, INTENT(in) :: dim_ens               ! Size of ensemble
       REAL, INTENT(inout) :: modes(dim, dim_ens-1) ! Array of EOF modes
       REAL, INTENT(in)    :: svals(dim_ens-1)      ! Vector of singular values
       REAL, INTENT(inout) :: state(dim)            ! PE-local model state
       REAL, INTENT(out)   :: ens(dim, dim_ens)     ! State ensemble
       INTEGER, INTENT(in) :: verbose               ! Verbosity flag
       INTEGER, INTENT(inout) :: flag               ! Status flag
     END SUBROUTINE PDAF_SampleEns
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_eofcovar(dim_state, nstates, nfields, dim_fields, offsets, &
          remove_mstate, do_mv, states, stddev, svals, svec, meanstate, verbose, status)
       INTEGER, INTENT(in) :: dim_state           ! Dimension of state vector
       INTEGER, INTENT(in) :: nstates             ! Number of state vectors
       INTEGER, INTENT(in) :: nfields             ! Number of fields in state vector
       INTEGER, INTENT(in) :: dim_fields(nfields) ! Size of each field
       INTEGER, INTENT(in) :: offsets(nfields)    ! Start position of each field
       INTEGER, INTENT(in) :: do_mv               ! 1: Do multivariate scaling; 0: no scaling
                     ! nfields, dim_fields and offsets are only used if do_mv=1
       INTEGER, INTENT(in) :: remove_mstate       ! 1: subtract mean state from states
                                                  ! before computing EOFs; 0: don't remove
       REAL, INTENT(inout)  :: states(dim_state, nstates) ! State perturbations
       REAL, INTENT(out) :: stddev(nfields)       ! Standard deviation of field variability
                     ! Without multivariate scaling (do_mv=0), it is stddev = 1.0
       REAL, INTENT(out) :: svals(nstates)        ! Singular values divided by sqrt(nstates-1)
       REAL, INTENT(out) :: svec(dim_state, nstates)   ! Singular vectors
       REAL, INTENT(inout) :: meanstate(dim_state)     ! Mean state (only changed if remove_mstate=1)
       INTEGER, INTENT(in) :: verbose             ! Verbosity flag
       INTEGER, INTENT(out) :: status             ! Status flag
     END SUBROUTINE PDAF_eofcovar
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_diag_histogram(ncall, dim, dim_ens, element, &
          state, ens, hist, delta, status)
       INTEGER, INTENT(in) :: ncall              ! Number of calls to routine
       INTEGER, INTENT(in) :: dim                ! State dimension
       INTEGER, INTENT(in) :: dim_ens            ! Ensemble size
       INTEGER, INTENT(in) :: element            ! Element of vector used for histogram
          ! If element=0, all elements are used
       REAL, INTENT(in)   :: state(dim)          ! State vector
       REAL, INTENT(in)   :: ens(dim, dim_ens)   ! State ensemble
       INTEGER, INTENT(inout) :: hist(dim_ens+1) ! Histogram about the state
       REAL, INTENT(out)     :: delta            ! deviation measure from flat histogram
       INTEGER, INTENT(out)   :: status          ! Status flag (0=success)
     END SUBROUTINE PDAF_diag_histogram
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_diag_ensstats(dim, dim_ens, element, &
          state, ens, skewness, kurtosis, status)
       INTEGER, INTENT(in) :: dim               ! PE-local state dimension
       INTEGER, INTENT(in) :: dim_ens           ! Ensemble size
       INTEGER, INTENT(in) :: element           ! ID of element to be used
          ! If element=0, mean values over all elements are computed
       REAL, INTENT(in)    :: state(dim)        ! State vector
       REAL, INTENT(in)    :: ens(dim, dim_ens) ! State ensemble
       REAL, INTENT(out)   :: skewness          ! Skewness of ensemble
       REAL, INTENT(out)   :: kurtosis          ! Kurtosis of ensemble
       INTEGER, INTENT(out) :: status           ! Status flag (0=success)
     END SUBROUTINE PDAF_diag_ensstats
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_diag_effsample(dim_sample, weights, effSample)
       INTEGER, INTENT(in)  :: dim_sample         ! Sample size
       REAL, INTENT(in)    :: weights(dim_sample) ! weights of the samples
       REAL, INTENT(out)   :: effsample           ! effecfive sample size
     END SUBROUTINE PDAF_diag_effsample
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_diag_crps(dim_p, dim_ens, element, oens, obs, &
          CRPS, reli, pot_CRPS, uncert, status)!
       IMPLICIT NONE
       INTEGER, INTENT(in) :: dim_p                !< PE-local state dimension
       INTEGER, INTENT(in) :: dim_ens              !< Ensemble size
       INTEGER, INTENT(in) :: element              !< index of element in full state vector
       !< If element=0, mean values over dim_p grid points/cases are computed
       REAL, INTENT(in)    :: oens(dim_p, dim_ens) !< State ensemble
       REAL, INTENT(in)    :: obs(dim_p)           !< Observation / truth
       REAL, INTENT(out)   :: CRPS                 !< CRPS
       REAL, INTENT(out)   :: reli                 !< Reliability
       REAL, INTENT(out)   :: pot_CRPS             !< potential CRPS
       REAL, INTENT(out)   :: uncert               !< uncertainty
       INTEGER, INTENT(out) :: status              !< Status flag (0=success)
     END SUBROUTINE PDAF_diag_crps
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_diag_crps_mpi(dim_p, dim_ens, element, oens, obs, &
          COMM_filter, mype_filter, npes_filter, &
          CRPS, reli, pot_CRPS, uncert, status)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: dim_p                !< PE-local state dimension
       INTEGER, INTENT(in) :: dim_ens              !< Ensemble size
       INTEGER, INTENT(in) :: element              !< index of element in full state vector
       !< If element=0, mean values over dim_p grid points/cases are computed
       INTEGER, INTENT(in) :: COMM_filter          !< MPI communicator for filter
       INTEGER, INTENT(in) :: mype_filter          !< rank of MPI communicator
       INTEGER, INTENT(in) :: npes_filter          !< size of MPI communicator
       REAL, INTENT(in)    :: oens(dim_p, dim_ens) !< State ensemble
       REAL, INTENT(in)    :: obs(dim_p)           !< Observation / truth
       REAL, INTENT(out)   :: CRPS                 !< CRPS
       REAL, INTENT(out)   :: reli                 !< Reliability
       REAL, INTENT(out)   :: pot_CRPS             !< potential CRPS
       REAL, INTENT(out)   :: uncert               !< uncertainty
       INTEGER, INTENT(out) :: status              !< Status flag (0=success)
     END SUBROUTINE PDAF_diag_crps_mpi
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_diag_CRPS_nompi(dim, dim_ens, element, oens, obs, &
          CRPS, reli, resol, uncert, status)!
       IMPLICIT NONE
       INTEGER, INTENT(in) :: dim                !< PE-local state dimension
       INTEGER, INTENT(in) :: dim_ens            !< Ensemble size
       INTEGER, INTENT(in) :: element            !< ID of element to be used
       !< If element=0, mean values over all elements are computed
       REAL, INTENT(in)    :: oens(dim, dim_ens) !< State ensemble
       REAL, INTENT(in)    :: obs(dim)           !< State ensemble
       REAL, INTENT(out)   :: CRPS               !< CRPS
       REAL, INTENT(out)   :: reli               !< Reliability
       REAL, INTENT(out)   :: resol              !< resolution
       REAL, INTENT(out)   :: uncert             !< uncertainty
       INTEGER, INTENT(out) :: status            !< Status flag (0=success)
     END SUBROUTINE PDAF_diag_CRPS_nompi
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_gather_obs_f(obs_p, obs_f, status)
       USE PDAF_mod_filtermpi, &
            ONLY: dimobs_p, dimobs_f
       IMPLICIT NONE
       REAL, INTENT(in)  :: obs_p(dimobs_p)  ! PE-local vector
       REAL, INTENT(out) :: obs_f(dimobs_f)  ! Full gathered vector
       INTEGER, INTENT(out) :: status   ! Status flag: 
                                        ! (0) no error
                                        ! (1) when PDAF_gather_dim_obs_f not executed before
     END SUBROUTINE PDAF_gather_obs_f
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_gather_obs_f2(coords_p, coords_f, nrows, status)
       USE PDAF_mod_filtermpi, &
            ONLY: dimobs_p, dimobs_f
       IMPLICIT NONE
       INTEGER, INTENT(in) :: nrows     ! Number of rows in array
       REAL, INTENT(in)  :: coords_p(nrows, dimobs_p)  ! PE-local array
       REAL, INTENT(out) :: coords_f(nrows, dimobs_f)  ! Full gathered array
       INTEGER, INTENT(out) :: status   ! Status flag: 
                                        ! (0) no error
                                        ! (1) when PDAF_gather dim_obs_f not executed before
     END SUBROUTINE PDAF_gather_obs_f2
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_get_assim_flag(did_assim)
       IMPLICIT NONE
       INTEGER,INTENT(out) :: did_assim    ! Flag: (1) for assimilation; (0) else
     END SUBROUTINE PDAF_get_assim_flag
  END INTERFACE


  INTERFACE
     SUBROUTINE PDAF_get_localfilter(lfilter)
       INTEGER, INTENT(out) :: lfilter   ! Whether the filter is domain-localized
     END SUBROUTINE PDAF_get_localfilter
  END INTERFACE

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

  INTERFACE 
     SUBROUTINE PDAF_set_debug_flag(debugval)
       INTEGER, INTENT(in)        :: debugval  ! Value of debugging flag; print debug information for >0
     END SUBROUTINE PDAF_set_debug_flag
  END INTERFACE

  INTERFACE 
     SUBROUTINE PDAF_set_offline_mode(screen)
       INTEGER, INTENT(in)        :: screen    ! Verbosity flag
     END SUBROUTINE PDAF_set_offline_mode
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_g2l(dim_p, dim_l, idx_l_in_p, state_p, state_l)
       INTEGER, INTENT(in) :: dim_p          !< PE-local full state dimension
       INTEGER, INTENT(in) :: dim_l          !< Local state dimension
       INTEGER, INTENT(in) :: idx_l_in_p(dim_l)     !< Index array for projection
       REAL, INTENT(in)    :: state_p(dim_p) !< PE-local full state vector 
       REAL, INTENT(out)   :: state_l(dim_l) !< State vector on local analysis domain
     END SUBROUTINE PDAF_g2l
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_l2g(dim_p, dim_l, idx_l_in_p, state_p, state_l)
       INTEGER, INTENT(in) :: dim_p          !< PE-local full state dimension
       INTEGER, INTENT(in) :: dim_l          !< Local state dimension
       INTEGER, INTENT(in) :: idx_l_in_p(dim_l)     !< Index array for projection
       REAL, INTENT(inout) :: state_p(dim_p) !< PE-local full state vector 
       REAL, INTENT(in)    :: state_l(dim_l) !< State vector on local analysis domain
     END SUBROUTINE PDAF_l2g
  END INTERFACE


! OMI INTERFACES ---------------------

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

END MODULE PDAF_interfaces_module
