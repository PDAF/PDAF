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
!> Module providing shared variables for ensemble framework
!!
!! This module provides variables shared between the
!! subroutines of PDAF.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2003-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE PDAF_mod_filter
  
  IMPLICIT NONE
  SAVE

  ! *** Dimensions ***
  INTEGER :: dim_ens              !< Ensemble size 
  INTEGER :: dim_eof              !< Rank (number of columns of ens in SEEK)
  INTEGER :: dim_p                !< State dimension for PE-local domain
  INTEGER :: dim_bias_p=0         !< Dimension of bias vector
  INTEGER :: dim_lag=0            !< Number of past time instances considered for smoother
  INTEGER :: local_dim_ens        !< Local ensemble sizes (including state forecast)

  ! *** Variables for time stepping ***
  INTEGER :: step                 !< Current time step
  INTEGER :: step_obs             !< Time step of next observation
  INTEGER :: nsteps               !< Number of time steps to perform
  INTEGER :: cnt_steps            !< Number of time steps in current forecast phase

  ! *** Status and output variables ***
  INTEGER :: flag=0               !< PDAF status flag
  INTEGER :: screen=0             !< Control verbosity of filter routines
                                  !< (0) quiet; (1) normal output; (2); plus timings; (3) debug output
  INTEGER :: debug=0              !< Debugging flag: print debug information if >0

  ! *** Variables controlling ensemble forecasts ***
  LOGICAL :: offline_mode=.false. !< Whether to use PDAF offline mode
  INTEGER :: firsttime=1          !< Are the filter routines called for the first time?
  INTEGER :: initevol=1           !< Initialize a new forecast phase?
  INTEGER :: member=1             !< Which member of sub-ensemble to evolve
  INTEGER :: member_get=1         !< Which member of sub-ensemble to evolve (used in PDAF_get_state)
  INTEGER :: end_forecast         !< Whether to exit the forecasting
  INTEGER :: assim_flag=0         !< (1) if assimilation was done at this time stepn in PDAF_assimilate, (0) if not
  INTEGER :: incremental=0        !< Whether to perform incremental updating

  ! *** Specification of type and subtype of DA method ***
  INTEGER :: type_filter          !< Type of Filter
                                  !< (0) SEEK  (Pham et al., 1998a)
                                  !< (1) SEIK  (Pham et al., 1998b)
                                  !< (2) EnKF  (Evensen, 1994)
                                  !< (3) LSEIK (Nerger et al., 2007)
                                  !< (4) ETKF  (Bishop et al., 2001) 
                                  !<     (ETKF uses symmetric square roots like LETKF)
                                  !< (5) LETKF (Hunt et al., 2007)
  INTEGER :: subtype_filter       !< Sub-type of Filter
                                  !< Subtype of SEEK: 
                                  !<     (0) Evolve with finite difference approx to TLM
                                  !<     (1) Scaled modes, unit U
                                  !<     (2) Fixed basis (V); variable U matrix
                                  !<     (3) Fixed covar matrix (V,U kept static)
                                  !<     (5) PDAF offline mode
                                  !< Subtype of SEIK: 
                                  !<     (0) Usual SEIK with mean forecast, new formulation;
                                  !<     (1) Usual SEIK with mean forecast, old formulation;
                                  !<     (2) Fixed error space basis
                                  !<     (3) Fixed state covariance matrix
                                  !<     (4) SEIK with ensemble transformation (like ETKF)
                                  !<     (5) PDAF offline mode
                                  !< Subtype of EnKF forecast and update step: 
                                  !<     (0) Mean forecast & representer analysis for large dim_obs;
                                  !<     (1) Mean forecast & representer analysis for small dim_obs;
                                  !<     (5) PDAF offline mode
                                  !< Subtype of LSEIK: 
                                  !<     (0) Mean forecast;
                                  !<     (2) Fixed error space basis
                                  !<     (3) Fixed state covariance matrix
                                  !<     (4) LSEIK with ensemble transformation (like LETKF)
                                  !<     (5) PDAF offline mode
                                  !< Subtype of ETKF:
                                  !<     (0) ETKF using T-matrix like SEIK
                                  !<     (1) ETKF following Hunt et al. (2007)
                                  !<       There are no fixed basis/covariance cases, as
                                  !<       these are equivalent to SEIK subtypes 2/3
                                  !<     (5) PDAF offline mode
                                  !< Subtype of LETKF:
                                  !<     (0) ETKF using T-matrix like SEIK
                                  !<       There are no fixed basis/covariance cases, as
                                  !<       these are equivalent to LSEIK subtypes 2/3
                                  !<     (5) PDAF offline mode


  ! *** Control variables for DA method ***
  LOGICAL :: ensemblefilter=.true.          !< Whether the chosen filter is ensemble-based
  INTEGER :: localfilter=0                  !< Whether the chosen filter is domain-localized (1: yes)
  INTEGER :: globalobs=0                    !< Whether the chosen filter needs global observations (1: yes)
  CHARACTER(len=10) :: filterstr            !< String defining the filter type
  INTEGER :: cnt_maxlag=0                   !< Smoother: Count maximum number of past time instances
  LOGICAL :: inloop=.false.                 !< Whether the program is in the local analysis loop
  LOGICAL :: use_PDAF_assim = .false.       !< Whether we use PDAF_assimilate

  ! *** Information variables for filter ***
  INTEGER :: member_save = 1                !< Store member index for query with PDAF_get_memberid
  INTEGER :: obs_member = 0                 !< Ensemble member when calling the observation operator routine

  ! *** Filter fields ***
  REAL, ALLOCATABLE :: state(:)             !< PE-local model state
  REAL, ALLOCATABLE :: state_inc(:)         !< PE-local analysis increment for inc. updating
  REAL, ALLOCATABLE :: Ainv(:,:)            !< Transform matrix or matrix of eigenvalues
  REAL, ALLOCATABLE :: bias(:)              !< Model bias vector
  REAL, TARGET, ALLOCATABLE :: ens(:,:)     !< Ensemble matrix
                                            !<   or matrix of eigenvectors from EOF computation
  REAL, TARGET, ALLOCATABLE :: sens(:,:,:)  !< Ensemble matrix holding past times for smoothing
  REAL, TARGET, ALLOCATABLE :: skewness(:)  !< Skewness of ensemble for each local domain
  REAL, TARGET, ALLOCATABLE :: kurtosis(:)  !< Kurtosis of ensemble for each local domain

!$OMP THREADPRIVATE(cnt_maxlag, obs_member, debug)

END MODULE PDAF_mod_filter
