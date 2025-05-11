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
!!
!! __Revision history:__
!! * 2025-04 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
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

  ABSTRACT INTERFACE 
     SUBROUTINE collect_cb(dim_p, state_p)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
       REAL, INTENT(inout) :: state_p(dim_p)  !< local state vector
     END SUBROUTINE collect_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE distribute_cb(dim_p, state_p)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
       REAL, INTENT(inout) :: state_p(dim_p)  !< PE-local state vector
     END SUBROUTINE distribute_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE init_dim_obs_cb(step, dim_obs)
       IMPLICIT NONE
       INTEGER, INTENT(in)  :: step     !< Current time step
       INTEGER, INTENT(out) :: dim_obs  !< Dimension of full observation vector
     END SUBROUTINE init_dim_obs_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE obs_op_cb(step, dim_p, dim_obs, state_p, ostate)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: step                 !< Current time step
       INTEGER, INTENT(in) :: dim_p                !< PE-local state dimension
       INTEGER, INTENT(in) :: dim_obs              !< Dimension of full observed state
       REAL, INTENT(in)    :: state_p(dim_p)       !< PE-local model state
       REAL, INTENT(inout) :: ostate(dim_obs)      !< PE-local full observed state
     END SUBROUTINE obs_op_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE init_n_domains_cb(step, n_domains_p)
       IMPLICIT NONE
       INTEGER, INTENT(in)  :: step        !< Current time step
       INTEGER, INTENT(out) :: n_domains_p !< PE-local number of analysis domains
     END SUBROUTINE init_n_domains_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE init_dim_l_cb(step, domain_p, dim_l)
       IMPLICIT NONE
       INTEGER, INTENT(in)  :: step     !< Current time step
       INTEGER, INTENT(in)  :: domain_p !< Current local analysis domain
       INTEGER, INTENT(out) :: dim_l    !< Local state dimension
     END SUBROUTINE init_dim_l_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE init_dim_obs_l_cb(domain_p, step, dim_obs, dim_obs_l)
       IMPLICIT NONE
       INTEGER, INTENT(in)  :: domain_p   !< Index of current local analysis domain
       INTEGER, INTENT(in)  :: step       !< Current time step
       INTEGER, INTENT(in)  :: dim_obs    !< Full dimension of observation vector
       INTEGER, INTENT(out) :: dim_obs_l  !< Local dimension of observation vector
     END SUBROUTINE init_dim_obs_l_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE prepost_cb(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
          state_p, Uinv, ens_p, flag)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: step        !< Current time step (negative for call after forecast)
       INTEGER, INTENT(in) :: dim_p       !< PE-local state dimension
       INTEGER, INTENT(in) :: dim_ens     !< Size of state ensemble
       INTEGER, INTENT(in) :: dim_ens_p   !< PE-local size of ensemble
       INTEGER, INTENT(in) :: dim_obs_p   !< PE-local dimension of observation vector
       REAL, INTENT(inout) :: state_p(dim_p) !< PE-local forecast/analysis state
       REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) !< Inverse of matrix U
       REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      !< PE-local state ensemble
       INTEGER, INTENT(in) :: flag        !< PDAF status flag
     END SUBROUTINE prepost_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE next_obs_cb(stepnow, nsteps, doexit, time)
       IMPLICIT NONE
       INTEGER, INTENT(in)  :: stepnow  !< Number of the current time step
       INTEGER, INTENT(out) :: nsteps   !< Number of time steps until next obs
       INTEGER, INTENT(out) :: doexit   !< Whether to exit forecasting (1 for exit)
       REAL, INTENT(out)    :: time     !< Current model (physical) time
     END SUBROUTINE next_obs_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE g2l_state_cb(step, domain_p, dim_p, state_p, dim_l, state_l)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: step           !< Current time step
       INTEGER, INTENT(in) :: domain_p       !< Current local analysis domain
       INTEGER, INTENT(in) :: dim_p          !< PE-local full state dimension
       INTEGER, INTENT(in) :: dim_l          !< Local state dimension
       REAL, INTENT(in)    :: state_p(dim_p) !< PE-local full state vector 
       REAL, INTENT(out)   :: state_l(dim_l) !< State vector on local analysis domain
     END SUBROUTINE g2l_state_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE l2g_state_cb(step, domain_p, dim_l, state_l, dim_p, state_p)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: step           !< Current time step
       INTEGER, INTENT(in) :: domain_p       !< Current local analysis domain
       INTEGER, INTENT(in) :: dim_l          !< Local state dimension
       INTEGER, INTENT(in) :: dim_p          !< PE-local full state dimension
       REAL, INTENT(in)    :: state_l(dim_l) !< State vector on local analysis domain
       REAL, INTENT(inout) :: state_p(dim_p) !< PE-local full state vector 
     END SUBROUTINE l2g_state_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE localize_cb(dim_p, dim_obs, HP_p, HPH)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: dim_p                 !< Process-local state dimension
       INTEGER, INTENT(in) :: dim_obs               !< Number of observations
       REAL, INTENT(inout) :: HP_p(dim_obs, dim_p)  !< Process-local part of matrix HP
       REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH
     END SUBROUTINE localize_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE get_obs_cb(step, dim_obs, observation)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: step                 !< Current time step
       INTEGER, INTENT(in) :: dim_obs              !< Dimension of obs. vector
       REAL, INTENT(out)   :: observation(dim_obs) !< Observation vector
     END SUBROUTINE get_obs_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE localize_serial_cb(iobs, dim_p, dim_obs, HP_p, HXY_p)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: iobs           !< Index of current observation
       INTEGER, INTENT(in) :: dim_p          !< Process-local state dimension
       INTEGER, INTENT(in) :: dim_obs        !< Number of observations
       REAL, INTENT(inout) :: HP_p(dim_p)    !< Process-local part of matrix HP for observation iobs
       REAL, INTENT(inout) :: HXY_p(dim_obs) !< Process-local part of matrix HX(HX_all) for full observations
     END SUBROUTINE localize_serial_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE obs_op_lin_cb(step, dim_p, dim_obs, state_p, ostate)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: step                 !< Current time step
       INTEGER, INTENT(in) :: dim_p                !< PE-local state dimension
       INTEGER, INTENT(in) :: dim_obs              !< Dimension of full observed state
       REAL, INTENT(in)    :: state_p(dim_p)       !< PE-local model state
       REAL, INTENT(inout) :: ostate(dim_obs)      !< PE-local full observed state
     END SUBROUTINE obs_op_lin_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE obs_op_adj_cb(step, dim_p, dim_obs, ostate, state_p)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: step                 !< Current time step
       INTEGER, INTENT(in) :: dim_p                !< PE-local state dimension
       INTEGER, INTENT(in) :: dim_obs              !< Dimension of full observed state
       REAL, INTENT(in)    :: ostate(dim_obs)      !< PE-local full observed state
       REAL, INTENT(inout) :: state_p(dim_p)       !< PE-local model state
     END SUBROUTINE obs_op_adj_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE cvt_cb(iter, dim_p, dim_cvec, v_p, Vv_p)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: iter          !< Iteration of optimization
       INTEGER, INTENT(in) :: dim_p         !< PE-local observation dimension
       INTEGER, INTENT(in) :: dim_cvec      !< Dimension of control vector
       REAL, INTENT(in)    :: v_p(dim_cvec) !< PE-local control vector
       REAL, INTENT(inout) :: Vv_p(dim_p)   !< PE-local result vector
     END SUBROUTINE cvt_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE cvt_adj_cb(iter, dim_p, dim_cvec, Vv_p, v_p)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: iter          !< Iteration of optimization
       INTEGER, INTENT(in) :: dim_p         !< PE-local observation dimension
       INTEGER, INTENT(in) :: dim_cvec      !< Dimension of control vector
       REAL, INTENT(in)    :: Vv_p(dim_p)   !< PE-local input vector
       REAL, INTENT(inout) :: v_p(dim_cvec) !< PE-local result vector
     END SUBROUTINE cvt_adj_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE cvt_ens_cb(iter, dim_p, dim_ens, dim_cvec_ens, ens_p, &
          v_p, Vv_p)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: iter               !< Iteration of optimization
       INTEGER, INTENT(in) :: dim_p              !< PE-local dimension of state
       INTEGER, INTENT(in) :: dim_ens            !< Ensemble size
       INTEGER, INTENT(in) :: dim_cvec_ens       !< Dimension of control vector
       REAL, INTENT(in) :: ens_p(dim_p, dim_ens) !< PE-local ensemble
       REAL, INTENT(in) :: v_p(dim_cvec_ens)     !< PE-local control vector
       REAL, INTENT(inout) :: Vv_p(dim_p)        !< PE-local state increment
     END SUBROUTINE cvt_ens_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE cvt_adj_ens_cb(iter, dim_p, dim_ens, dim_cvec_ens, ens_p, &
          Vv_p, v_p)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: iter               !< Iteration of optimization
       INTEGER, INTENT(in) :: dim_p              !< PE-local dimension of state
       INTEGER, INTENT(in) :: dim_ens            !< Ensemble size
       INTEGER, INTENT(in) :: dim_cvec_ens       !< Number of columns in HV_p
       REAL, INTENT(in) :: ens_p(dim_p, dim_ens) !< PE-local ensemble
       REAL, INTENT(in)    :: Vv_p(dim_p)        !< PE-local input vector
       REAL, INTENT(inout) :: v_p(dim_cvec_ens)  !< PE-local result vector
     END SUBROUTINE cvt_adj_ens_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE prodRinvA_l_cb(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: domain_p          !< Index of current local analysis domain
       INTEGER, INTENT(in) :: step              !< Current time step
       INTEGER, INTENT(in) :: dim_obs_l         !< Dimension of local observation vector
       INTEGER, INTENT(in) :: rank              !< Rank of initial covariance matrix
       REAL, INTENT(in)    :: obs_l(dim_obs_l)  !< Local vector of observations
       REAL, INTENT(inout) :: A_l(dim_obs_l, rank) !< Input matrix
       REAL, INTENT(out)   :: C_l(dim_obs_l, rank) !< Output matrix
     END SUBROUTINE prodRinvA_l_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE likelihood_l_cb(domain_p, step, dim_obs_l, obs_l, resid_l, lhood_l)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: domain_p           ! Current local analysis domain
       INTEGER, INTENT(in) :: step               !< Current time step
       INTEGER, INTENT(in) :: dim_obs_l          !< PE-local dimension of obs. vector
       REAL, INTENT(in)    :: obs_l(dim_obs_l)   !< PE-local vector of observations
       REAL, INTENT(inout) :: resid_l(dim_obs_l) !< Input vector of residuum
       REAL, INTENT(out)   :: lhood_l            !< Output vector - log likelihood
     END SUBROUTINE likelihood_l_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE prodRinvA_hyb_l_cb(domain_p, step, dim_obs_l, rank, obs_l, alpha, A_l, C_l)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: domain_p          !< Index of current local analysis domain
       INTEGER, INTENT(in) :: step              !< Current time step
       INTEGER, INTENT(in) :: dim_obs_l         !< Dimension of local observation vector
       INTEGER, INTENT(in) :: rank              !< Rank of initial covariance matrix
       REAL, INTENT(in)    :: obs_l(dim_obs_l)  !< Local vector of observations
       REAL, INTENT(in)    :: alpha             !< Hybrid weight
       REAL, INTENT(inout) :: A_l(dim_obs_l, rank) !< Input matrix
       REAL, INTENT(out)   :: C_l(dim_obs_l, rank) !< Output matrix
     END SUBROUTINE prodRinvA_hyb_l_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE likelihood_hyb_l_cb(domain_p, step, dim_obs_l, obs_l, resid_l, alpha, lhood_l)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: domain_p           !< Current local analysis domain
       INTEGER, INTENT(in) :: step               !< Current time step
       INTEGER, INTENT(in) :: dim_obs_l          !< PE-local dimension of obs. vector
       REAL, INTENT(in)    :: obs_l(dim_obs_l)   !< PE-local vector of observations
       REAL, INTENT(inout) :: resid_l(dim_obs_l) !< Input vector of residuum
       REAL, INTENT(in)    :: alpha              !< Hybrid weight
       REAL, INTENT(out)   :: lhood_l            !< Output vector - log likelihood
     END SUBROUTINE likelihood_hyb_l_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE prodRinvA_cb(step, dim_obs_p, ncol, obs_p, A_p, C_p)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: step              !< Current time step
       INTEGER, INTENT(in) :: dim_obs_p         !< Dimension of PE-local observation vector
       INTEGER, INTENT(in) :: ncol              !< Number of columns in A_p and C_p
       REAL, INTENT(in)    :: obs_p(dim_obs_p)  !< PE-local vector of observations
       REAL, INTENT(in)    :: A_p(dim_obs_p, ncol) !< Input matrix
       REAL, INTENT(out)   :: C_p(dim_obs_p, ncol) !< Output matrix
     END SUBROUTINE prodRinvA_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE likelihood_cb(step, dim_obs, obs, resid, lhood)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: step             !< Current time step
       INTEGER, INTENT(in) :: dim_obs          !< PE-local dimension of obs. vector
       REAL, INTENT(in)    :: obs(dim_obs)     !< PE-local vector of observations
       REAL, INTENT(in)    :: resid(dim_obs)   !< Input vector of residuum
       REAL, INTENT(out)   :: lhood            !< Output vector - log likelihood
     END SUBROUTINE likelihood_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE add_obs_error_cb(step, dim_obs_p, C_p)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: step              !< Current time step
       INTEGER, INTENT(in) :: dim_obs_p         !< Dimension of PE-local observation vector
       REAL, INTENT(inout) :: C_p(dim_obs_p,dim_obs_p) ! Matrix to which R is added
     END SUBROUTINE add_obs_error_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE init_obscovar_cb(step, dim_obs, dim_obs_p, covar, m_state_p, &
          isdiag)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: step                 !< Current time step
       INTEGER, INTENT(in) :: dim_obs              !< Dimension of observation vector
       INTEGER, INTENT(in) :: dim_obs_p            !< PE-local dimension of obs. vector
       REAL, INTENT(out) :: covar(dim_obs,dim_obs) !< Observation error covar. matrix
       REAL, INTENT(in) :: m_state_p(dim_obs_p)    !< Observation vector
       LOGICAL, INTENT(out) :: isdiag              !< Whether matrix R is diagonal
     END SUBROUTINE init_obscovar_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE init_obserr_cb(step, dim_obs_f, obs_f, obserr_f)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: step                !< Current time step
       INTEGER, INTENT(in) :: dim_obs_f           !< Full dimension of observation vector
       REAL, INTENT(in)    :: obs_f(dim_obs_f)    !< Full observation vector
       REAL, INTENT(out)   :: obserr_f(dim_obs_f) !< Full observation error stddev
     END SUBROUTINE init_obserr_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE init_obs_cb(step, dim_obs_f, observation_f)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: step        !< Current time step
       INTEGER, INTENT(in) :: dim_obs_f   !< Dimension of full observation vector
       REAL, INTENT(out)   :: observation_f(dim_obs_f) !< Full observation vector
     END SUBROUTINE init_obs_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE init_obs_l_cb(domain_p, step, dim_obs_l, observation_l)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: domain_p   !< Index of current local analysis domain index
       INTEGER, INTENT(in) :: step       !< Current time step
       INTEGER, INTENT(in) :: dim_obs_l  !< Local dimension of observation vector
       REAL, INTENT(out)   :: observation_l(dim_obs_l) !< Local observation vector
     END SUBROUTINE init_obs_l_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE init_obsvar_cb(step, dim_obs_p, obs_p, meanvar)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: step          !< Current time step
       INTEGER, INTENT(in) :: dim_obs_p     !< PE-local dimension of observation vector
       REAL, INTENT(in) :: obs_p(dim_obs_p) !< PE-local observation vector
       REAL, INTENT(out)   :: meanvar       !< Mean observation error variance
     END SUBROUTINE init_obsvar_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE init_obsvar_l_cb(domain_p, step, dim_obs_l, obs_l, meanvar_l)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: domain_p      !< Index of current local analysis domain
       INTEGER, INTENT(in) :: step          !< Current time step
       INTEGER, INTENT(in) :: dim_obs_l     !< Local dimension of observation vector
       REAL, INTENT(in) :: obs_l(dim_obs_l) !< Local observation vector
       REAL, INTENT(out)   :: meanvar_l     !< Mean local observation error variance
     END SUBROUTINE init_obsvar_l_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE g2l_obs_cb(domain_p, step, dim_obs_f, dim_obs_l, ostate_f, ostate_l)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: domain_p   !< Index of current local analysis domain
       INTEGER, INTENT(in) :: step       !< Current time step
       INTEGER, INTENT(in) :: dim_obs_f  !< Dimension of full PE-local observation vector
       INTEGER, INTENT(in) :: dim_obs_l  !< Dimension of local observation vector
       REAL, INTENT(in)    :: ostate_f(dim_obs_f)   !< Full PE-local obs.ervation vector
       REAL, INTENT(out)   :: ostate_l(dim_obs_l)   !< Observation vector on local domain
     END SUBROUTINE g2l_obs_cb
  END INTERFACE

  ABSTRACT INTERFACE 
     SUBROUTINE init_obsvars_cb(step, dim_obs_f, var_f)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: step           !< Current time step
       INTEGER, INTENT(in) :: dim_obs_f      !< Dimension of full observation vector
       REAL, INTENT(out) :: var_f(dim_obs_f) !< vector of observation error variances
     END SUBROUTINE init_obsvars_cb
  END INTERFACE


END MODULE PDAF_cb_procedures
