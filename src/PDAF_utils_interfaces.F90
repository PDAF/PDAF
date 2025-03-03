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
!! * 2025-02 - Lars Nerger - Split from PDAF_interfaces_module for better coding overview
!! * Other revisions - see repository log
MODULE PDAF_utils_interfaces

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
     SUBROUTINE PDAF_gather_dim_obs_f(dim_obs_p, dim_obs_f)
       INTEGER, INTENT(in)  :: dim_obs_p    ! PE-local observation dimension
       INTEGER, INTENT(out) :: dim_obs_f    ! Full observation dimension
     END SUBROUTINE PDAF_gather_dim_obs_f
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
     SUBROUTINE PDAF_reset_forget(forget_in)
       REAL, INTENT(in) :: forget_in    ! New value of forgetting factor
     END SUBROUTINE PDAF_reset_forget
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_sampleens(dim, dim_ens, modes, svals, state, ens, verbose, flag)
       INTEGER, INTENT(in) :: dim                   ! Size of state vector
       INTEGER, INTENT(in) :: dim_ens               ! Size of ensemble
       REAL, INTENT(inout) :: modes(dim, dim_ens-1) ! Array of EOF modes
       REAL, INTENT(in)    :: svals(dim_ens-1)      ! Vector of singular values
       REAL, INTENT(inout) :: state(dim)            ! PE-local model state
       REAL, INTENT(out)   :: ens(dim, dim_ens)     ! State ensemble
       INTEGER, INTENT(in) :: verbose               ! Verbosity flag
       INTEGER, INTENT(inout) :: flag               ! Status flag
     END SUBROUTINE PDAF_sampleens
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_mvnormalize(mode, dim_state, dim_field, offset, &
          ncol, states, stddev, status)
       INTEGER, INTENT(in) :: mode       ! Mode: (1) normalize, (2) re-scale
       INTEGER, INTENT(in) :: dim_state  ! Dimension of state vector
       INTEGER, INTENT(in) :: dim_field  ! Dimension of a field in state vector
       INTEGER, INTENT(in) :: offset     ! Offset of field in state vector
       INTEGER, intent(in) :: ncol       ! Number of columns in array states
       REAL, INTENT(inout) :: states(dim_state, ncol)  ! State vector array
       REAL, INTENT(inout) :: stddev     ! Standard deviation of field
                                         ! stddev is output for mode=1 and input for mode=2
       INTEGER, INTENT(out) :: status    ! Status flag (0=success)
     END SUBROUTINE PDAF_mvnormalize
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

! PDAF_diag

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

! PDAF_gather

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

! PDAF_g2l/l2g

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

! PDAF_get

  INTERFACE
     SUBROUTINE PDAF_get_assim_flag(did_assim)
       IMPLICIT NONE
       INTEGER,INTENT(out) :: did_assim    ! Flag: (1) for assimilation; (0) else
     END SUBROUTINE PDAF_get_assim_flag
  END INTERFACE

  INTERFACE 
     SUBROUTINE PDAF_get_ensstats(skew_ptr, kurt_ptr, status)
       REAL, POINTER, INTENT(out) :: skew_ptr(:)  ! Pointer to skewness array
       REAL, POINTER, INTENT(out) :: kurt_ptr(:)  ! Pointer to kurtosis array
       INTEGER, INTENT(out)       :: status  ! Status flag 
     END SUBROUTINE PDAF_get_ensstats
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_get_localfilter(lfilter)
       INTEGER, INTENT(out) :: lfilter   ! Whether the filter is domain-localized
     END SUBROUTINE PDAF_get_localfilter
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
     SUBROUTINE PDAF_get_smootherens(sens_point, maxlag, status)
       REAL, POINTER, INTENT(out) :: sens_point(:,:,:)  ! Pointer to smoother array
       INTEGER, INTENT(out)       :: maxlag  ! Number of past timesteps processed in sens
       INTEGER, INTENT(out)       :: status  ! Status flag 
     END SUBROUTINE PDAF_get_smootherens
  END INTERFACE

! PDAF_set

  INTERFACE
     SUBROUTINE PDAF_set_comm_pdaf(in_COMM_pdaf)
       INTEGER,INTENT(in) :: in_COMM_pdaf    !< MPI communicator for PDAF
     END SUBROUTINE PDAF_set_comm_pdaf
  END INTERFACE

  INTERFACE 
     SUBROUTINE PDAF_set_debug_flag(debugval)
       INTEGER, INTENT(in)        :: debugval  ! Value of debugging flag; print debug information for >0
     END SUBROUTINE PDAF_set_debug_flag
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_set_ens_pointer(ens_point, status)
       REAL, POINTER, INTENT(out) :: ens_point(:,:)  ! Pointer to smoother array
       INTEGER, INTENT(out)       :: status          ! Status flag
     END SUBROUTINE PDAF_set_ens_pointer
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_set_iparam(id, value, flag)
       INTEGER, INTENT(in) :: id       !< Index of parameter
       INTEGER, INTENT(in) :: value    !< Parameter value
       INTEGER, INTENT(out) :: flag    !< Status flag: 0 for no error
     END SUBROUTINE PDAF_set_iparam
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_set_memberid(memberid)
       INTEGER,INTENT(inout) :: memberid    !< Index in the local ensemble
     END SUBROUTINE PDAF_set_memberid
  END INTERFACE

  INTERFACE 
     SUBROUTINE PDAF_set_offline_mode(screen)
       INTEGER, INTENT(in)        :: screen    ! Verbosity flag
     END SUBROUTINE PDAF_set_offline_mode
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_set_rparam(id, value, flag)
       INTEGER, INTENT(in)  :: id       !< Index of parameter
       REAL, INTENT(in)     :: value    !< Parameter value
       INTEGER, INTENT(out) :: flag     !< Status flag: 0 for no error
     END SUBROUTINE PDAF_set_rparam
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_set_seedset(seedset_in)
       INTEGER,INTENT(in) :: seedset_in    !< Seedset index (1-20)
     END SUBROUTINE PDAF_set_seedset
  END INTERFACE

  INTERFACE 
     SUBROUTINE PDAF_set_smootherens(sens_point, maxlag, status)
       REAL, POINTER, INTENT(out) :: sens_point(:,:,:)  ! Pointer to smoother array
       INTEGER, INTENT(in)        :: maxlag  ! Number of past timesteps processed in sens
       INTEGER, INTENT(out)       :: status  ! Status flag 
     END SUBROUTINE PDAF_set_smootherens
  END INTERFACE

END MODULE PDAF_utils_interfaces
