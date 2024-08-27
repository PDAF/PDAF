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
! !ROUTINE: PDAFlocal_interfaces --- Interface definitions for PDAFlocal
!
! !INTERFACE:
MODULE PDAFlocal_interfaces

! !DESCRIPTION:
! Module providing interface definition of the PDAFlocal routines that
! are called from the model code.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2024-08 - Lars Nerger - Initial code
! Later revisions - see svn log
!EOP
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  INTERFACE
     SUBROUTINE PDAFlocal_put_state_lseik(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
          U_init_dim_l, U_init_dim_obs_l, U_g2l_obs, &
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
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prepoststep           ! User supplied pre/poststep routine
     END SUBROUTINE PDAFlocal_put_state_lseik
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_put_state_lseik_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAFlocal_put_state_lseik_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_assimilate_lseik(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
          U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
          U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
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
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_next_observation      ! Provide time step and time of next observation
     END SUBROUTINE PDAFlocal_assimilate_lseik
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_assimilate_lseik_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAFlocal_assimilate_lseik_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_put_state_letkf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
          U_init_dim_l, U_init_dim_obs_l, U_g2l_obs, &
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
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prepoststep           ! User supplied pre/poststep routine
     END SUBROUTINE PDAFlocal_put_state_letkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_put_state_letkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAFlocal_put_state_letkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_assimilate_letkf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
          U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
          U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
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
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_next_observation      ! Provide time step and time of next observation
     END SUBROUTINE PDAFlocal_assimilate_letkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_assimilate_letkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAFlocal_assimilate_letkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_put_state_lestkf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
          U_init_dim_l, U_init_dim_obs_l, U_g2l_obs, &
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
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prepoststep           ! User supplied pre/poststep routine
     END SUBROUTINE PDAFlocal_put_state_lestkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_put_state_lestkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAFlocal_put_state_lestkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_assimilate_lestkf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
          U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
          U_g2l_obs, U_init_obsvar, U_init_obsvar_l, U_next_observation, flag)
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
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_next_observation      ! Provide time step and time of next observation
     END SUBROUTINE PDAFlocal_assimilate_lestkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_assimilate_lestkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAFlocal_assimilate_lestkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_put_state_lnetf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs_l, U_prepoststep, U_likelihood_l, U_init_n_domains_p, &
          U_init_dim_l, U_init_dim_obs_l, U_g2l_obs, &
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
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_likelihood_l, &       ! Compute observation likelihood for an ensemble member
            U_prepoststep           ! User supplied pre/poststep routine
     END SUBROUTINE PDAFlocal_put_state_lnetf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_put_state_lnetf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAFlocal_put_state_lnetf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_assimilate_lnetf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs_l, U_prepoststep, &
          U_likelihood_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
          U_g2l_obs, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &      ! Routine to distribute a state vector
            U_obs_op, &             ! Observation operator
            U_init_n_domains_p, &   ! Provide number of local analysis domains
            U_init_dim_l, &         ! Init state dimension for local ana. domain
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_init_dim_obs_l, &     ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs_l, &         ! Init. observation vector on local analysis domain
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_likelihood_l, &       ! Compute observation likelihood for an ensemble member
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_next_observation      ! Provide time step and time of next observation
     END SUBROUTINE PDAFlocal_assimilate_lnetf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_assimilate_lnetf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAFlocal_assimilate_lnetf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_put_state_lknetf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_prodRinvA_hyb_l, &
          U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_obs, &
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
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prodRinvA_hyb_l, &    ! Provide product R^-1 A on local analysis domain with hybrid weight
            U_likelihood_l, &       ! Compute likelihood
            U_likelihood_hyb_l, &   ! Compute likelihood with hybrid weight
            U_prepoststep           ! User supplied pre/poststep routine
     END SUBROUTINE PDAFlocal_put_state_lknetf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_put_state_lknetf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAFlocal_put_state_lknetf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_assimilate_lknetf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
          U_prodRinvA_l, U_prodRinvA_hyb_l, U_init_n_domains_p, U_init_dim_l, &
          U_init_dim_obs_l, U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
          U_likelihood_l, U_likelihood_hyb_l, U_next_observation, flag)
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
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prodRinvA_hyb_l, &    ! Provide product R^-1 A on local analysis domain with hybrid weight
            U_likelihood_l, &       ! Compute likelihood
            U_likelihood_hyb_l, &   ! Compute likelihood with hybrid weight
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_next_observation, &   ! Routine to provide time step, time and dimension
                                    !   of next observation
            U_distribute_state      ! Routine to distribute a state vector
     END SUBROUTINE PDAFlocal_assimilate_lknetf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_assimilate_lknetf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAFlocal_assimilate_lknetf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_put_state_en3dvar_lestkf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_prodRinvA, &
          U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
          U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
          U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_obs, &
          U_init_obsvar, U_init_obsvar_l, U_prepoststep, outflag)
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
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l           ! Provide product R^-1 A on local analysis domain
     END SUBROUTINE PDAFlocal_put_state_en3dvar_lestkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_assimilate_en3dvar_lestkf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
          U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
          U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
          U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
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
            U_g2l_obs, &              ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l             ! Provide product R^-1 A on local analysis domain
     END SUBROUTINE PDAFlocal_assimilate_en3dvar_lestkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_put_state_hyb3dvar_lestkf(U_collect_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
          U_cvt_ens, U_cvt_adj_ens, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
          U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
          U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_obs, &
          U_init_obsvar, U_init_obsvar_l, U_prepoststep, outflag)
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
            U_g2l_obs, &           ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l          ! Provide product R^-1 A on local analysis domain
     END SUBROUTINE PDAFlocal_put_state_hyb3dvar_lestkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_assimilate_hyb3dvar_lestkf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
          U_cvt_ens, U_cvt_adj_ens, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
          U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
          U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_obs, &
          U_init_obsvar, U_init_obsvar_l, U_prepoststep, U_next_observation, outflag)
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
            U_g2l_obs, &              ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l             ! Provide product R^-1 A on local analysis domain
     END SUBROUTINE PDAFlocal_assimilate_hyb3dvar_lestkf
  END INTERFACE

! PDAFlocal INTERFACES ---------------------

  INTERFACE
     SUBROUTINE PDAFlocalomi_put_state(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_prepoststep, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
          flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_obs_op, &                ! Observation operator
            U_init_n_domains_p, &      ! Provide number of local analysis domains
            U_init_dim_l, &            ! Init state dimension for local ana. domain
            U_init_dim_obs, &          ! Initialize dimension of observation vector
            U_init_dim_obs_l, &        ! Initialize dim. of obs. vector for local ana. domain
            U_prepoststep              ! User supplied pre/poststep routine
     END SUBROUTINE PDAFlocalomi_put_state
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_put_state_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFlocalomi_put_state_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_assimilate(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_prepoststep, U_init_n_domains_p, U_init_dim_l, &
          U_l2g_state, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &      ! Routine to distribute a state vector
            U_obs_op, &                ! Observation operator
            U_init_n_domains_p, &      ! Provide number of local analysis domains
            U_init_dim_l, &            ! Init state dimension for local ana. domain
            U_init_dim_obs, &          ! Initialize dimension of observation vector
            U_init_dim_obs_l, &        ! Initialize dim. of obs. vector for local ana. domain
            U_prepoststep, &           ! User supplied pre/poststep routine
            U_next_observation         ! Provide time step and time of next observation
     END SUBROUTINE PDAFlocalomi_assimilate
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_assimilate_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFlocalomi_assimilate_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_put_state_nondiagR(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_prepoststep, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_prodRinvA_l, &
          flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_obs_op, &                ! Observation operator
            U_init_n_domains_p, &      ! Provide number of local analysis domains
            U_init_dim_l, &            ! Init state dimension for local ana. domain
            U_init_dim_obs, &          ! Initialize dimension of observation vector
            U_init_dim_obs_l, &        ! Initialize dim. of obs. vector for local ana. domain
            U_prodRinvA_l, &           ! Provide product R^-1 A on local analysis domain
            U_prepoststep              ! User supplied pre/poststep routine
     END SUBROUTINE PDAFlocalomi_put_state_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_put_state_nondiagR_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFlocalomi_put_state_nondiagR_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_assimilate_nondiagR(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_prepoststep, U_init_n_domains_p, U_init_dim_l, &
          U_init_dim_obs_l, U_prodRinvA_l, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &      ! Routine to distribute a state vector
            U_obs_op, &                ! Observation operator
            U_init_n_domains_p, &      ! Provide number of local analysis domains
            U_init_dim_l, &            ! Init state dimension for local ana. domain
            U_init_dim_obs, &          ! Initialize dimension of observation vector
            U_init_dim_obs_l, &        ! Initialize dim. of obs. vector for local ana. domain
            U_prodRinvA_l, &           ! Provide product R^-1 A on local analysis domain
            U_prepoststep, &           ! User supplied pre/poststep routine
            U_next_observation         ! Provide time step and time of next observation
     END SUBROUTINE PDAFlocalomi_assimilate_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_assimilate_nondiagR_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFlocalomi_assimilate_nondiagR_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_put_state_lnetf_nondiagR(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_prepoststep, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_likelihood_l, &
          flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_obs_op, &                ! Observation operator
            U_init_n_domains_p, &      ! Provide number of local analysis domains
            U_init_dim_l, &            ! Init state dimension for local ana. domain
            U_init_dim_obs, &          ! Initialize dimension of observation vector
            U_init_dim_obs_l, &        ! Initialize dim. of obs. vector for local ana. domain
            U_likelihood_l, &          ! Compute likelihood and apply localization
            U_prepoststep              ! User supplied pre/poststep routine
     END SUBROUTINE PDAFlocalomi_put_state_lnetf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_put_state_lnetf_nondiagR_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFlocalomi_put_state_lnetf_nondiagR_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_assimilate_lnetf_nondiagR(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_prepoststep, U_init_n_domains_p, U_init_dim_l, &
          U_init_dim_obs_l, U_likelihood_l, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &      ! Routine to distribute a state vector
            U_obs_op, &                ! Observation operator
            U_init_n_domains_p, &      ! Provide number of local analysis domains
            U_init_dim_l, &            ! Init state dimension for local ana. domain
            U_init_dim_obs, &          ! Initialize dimension of observation vector
            U_init_dim_obs_l, &        ! Initialize dim. of obs. vector for local ana. domain
            U_likelihood_l, &          ! Compute likelihood and apply localization
            U_prepoststep, &           ! User supplied pre/poststep routine
            U_next_observation         ! Provide time step and time of next observation
     END SUBROUTINE PDAFlocalomi_assimilate_lnetf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_assimilate_lnetf_nondiagR_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFlocalomi_assimilate_lnetf_nondiagR_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_put_state_lknetf_nondiagR(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_prepoststep, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
          U_prodRinvA_l, U_prodRinvA_hyb_l, U_likelihood_l, U_likelihood_hyb_l, flag)
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
            U_prepoststep              ! User supplied pre/poststep routine
     END SUBROUTINE PDAFlocalomi_put_state_lknetf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_put_state_lknetf_nondiagR_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFlocalomi_put_state_lknetf_nondiagR_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_assimilate_lknetf_nondiagR(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_prepoststep, U_init_n_domains_p, U_init_dim_l, &
          U_init_dim_obs_l, U_prodRinvA_l, U_prodRinvA_hyb_l, U_likelihood_l, U_likelihood_hyb_l, &
          U_next_observation, flag)
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
            U_prepoststep, &           ! User supplied pre/poststep routine
            U_next_observation         ! Provide time step and time of next observation
     END SUBROUTINE PDAFlocalomi_assimilate_lknetf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_assimilate_lknetf_nondiagR_si(flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
     END SUBROUTINE PDAFlocalomi_assimilate_lknetf_nondiagR_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_put_state_en3dvar_lestkf(collect_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          prepoststep_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag    ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
       prepoststep_pdaf                     ! User supplied pre/poststep routine
       EXTERNAL :: cvt_ens_pdaf, &          ! Apply control vector transform matrix to control vector
            cvt_adj_ens_pdaf, &             ! Apply adjoint control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
       EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
            init_dim_l_pdaf, &              ! Init state dimension for local ana. domain
            init_dim_obs_f_pdaf, &          ! Initialize dimension of full observation vector
            obs_op_f_pdaf, &                ! Full observation operator
            init_dim_obs_l_pdaf             ! Initialize local dimimension of obs. vector
     END SUBROUTINE PDAFlocalomi_put_state_en3dvar_lestkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_assimilate_en3dvar_lestkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          prepoststep_pdaf, next_observation_pdaf, outflag)
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
            init_dim_obs_f_pdaf, &          ! Initialize dimension of full observation vector
            obs_op_f_pdaf, &                ! Full observation operator
            init_dim_obs_l_pdaf             ! Initialize local dimimension of obs. vector
     END SUBROUTINE PDAFlocalomi_assimilate_en3dvar_lestkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_put_state_en3dvar_lestkf_nondiagR(collect_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, prodRinvA_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          prepoststep_pdaf, outflag)
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
            init_dim_obs_f_pdaf, &          ! Initialize dimension of full observation vector
            obs_op_f_pdaf, &                ! Full observation operator
            init_dim_obs_l_pdaf, &          ! Initialize local dimimension of obs. vector
            prodRinvA_l_pdaf                ! Provide product R^-1 A with localization
     END SUBROUTINE PDAFlocalomi_put_state_en3dvar_lestkf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_assimilate_en3dvar_lestkf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, prodRinvA_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          prepoststep_pdaf, next_observation_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag    ! Status flag
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
            init_dim_obs_f_pdaf, &          ! Initialize dimension of full observation vector
            obs_op_f_pdaf, &                ! Full observation operator
            init_dim_obs_l_pdaf, &          ! Initialize local dimimension of obs. vector
            prodRinvA_l_pdaf                ! Provide product R^-1 A with localization
     END SUBROUTINE PDAFlocalomi_assimilate_en3dvar_lestkf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_put_state_hyb3dvar_lestkf(collect_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          prepoststep_pdaf, outflag)
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
            init_dim_obs_f_pdaf, &          ! Initialize dimension of full observation vector
            obs_op_f_pdaf, &                ! Full observation operator
            init_dim_obs_l_pdaf             ! Initialize local dimimension of obs. vector
     END SUBROUTINE PDAFlocalomi_put_state_hyb3dvar_lestkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_assimilate_hyb3dvar_lestkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          prepoststep_pdaf, next_observation_pdaf, outflag)
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
            init_dim_obs_f_pdaf, &          ! Initialize dimension of full observation vector
            obs_op_f_pdaf, &                ! Full observation operator
            init_dim_obs_l_pdaf             ! Initialize local dimimension of obs. vector
     END SUBROUTINE PDAFlocalomi_assimilate_hyb3dvar_lestkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_put_state_hyb3dvar_lestkf_nondiagR(collect_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, prodRinvA_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          prepoststep_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag    ! Status flag
       EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
       prepoststep_pdaf                     ! User supplied pre/poststep routine
       EXTERNAL :: cvt_ens_pdaf, &          ! Apply control vector transform matrix to control vector
            cvt_adj_ens_pdaf, &             ! Apply adjoint control vector transform matrix
            cvt_pdaf, &                     ! Apply control vector transform matrix to control vector
            cvt_adj_pdaf, &                 ! Apply adjoint control vector transform matrix
            obs_op_lin_pdaf, &              ! Linearized observation operator
            obs_op_adj_pdaf                 ! Adjoint observation operator
       EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
            init_dim_l_pdaf, &              ! Init state dimension for local ana. domain
            prodRinvA_pdaf, &               ! Provide product R^-1 A
            init_dim_obs_f_pdaf, &          ! Initialize dimension of full observation vector
            obs_op_f_pdaf, &                ! Full observation operator
            init_dim_obs_l_pdaf, &          ! Initialize local dimimension of obs. vector
            prodRinvA_l_pdaf                ! Provide product R^-1 A with localization
     END SUBROUTINE PDAFlocalomi_put_state_hyb3dvar_lestkf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocalomi_assimilate_hyb3dvar_lestkf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, prodRinvA_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          prepoststep_pdaf, next_observation_pdaf, outflag)
       INTEGER, INTENT(inout) :: outflag    ! Status flag
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
            init_dim_obs_f_pdaf, &          ! Initialize dimension of full observation vector
            obs_op_f_pdaf, &                ! Full observation operator
            init_dim_obs_l_pdaf, &          ! Initialize local dimimension of obs. vector
            prodRinvA_l_pdaf                ! Provide product R^-1 A with localization
     END SUBROUTINE PDAFlocalomi_assimilate_hyb3dvar_lestkf_nondiagR
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_set_indices(dim_l, map)
       INTEGER, INTENT(in) :: dim_l          !< Dimension of local state vector
       INTEGER, INTENT(in) :: map(dim_l)     !< Index array for mapping
     END SUBROUTINE PDAFlocal_set_indices
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_set_increment_weights(dim_l, weights)
       INTEGER, INTENT(in) :: dim_l          !< Dimension of local state vector
       REAL, INTENT(in) :: weights(dim_l)    !< Weights array
     END SUBROUTINE PDAFlocal_set_increment_weights
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAFlocal_clear_increment_weights()
     END SUBROUTINE PDAFlocal_clear_increment_weights
  END INTERFACE

END MODULE PDAFlocal_interfaces
