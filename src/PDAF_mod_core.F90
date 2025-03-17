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
!! subroutines of PDAF. These variables are internal to
!! PDAF, and should not directly be access by the user.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2003-06 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAF_mod_core
  
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

  ! *** Specification of type and subtype of DA method ***
  INTEGER :: type_filter          !< Type of Filter
  INTEGER :: subtype_filter       !< Sub-type of Filter

  ! *** Control variables for DA method ***
  LOGICAL :: ensemblefilter=.true.          !< Whether the chosen filter is ensemble-based
  INTEGER :: localfilter=0                  !< Whether the chosen filter is domain-localized (1: yes)
  INTEGER :: covarloc=0                     !< Whether the chosen filter uses covariance localization (1: yes)
  CHARACTER(len=10) :: filterstr            !< String defining the filter type
  INTEGER :: cnt_maxlag=0                   !< Smoother: Count maximum number of past time instances
  LOGICAL :: inloop=.false.                 !< Whether the program is in the local analysis loop
  LOGICAL :: use_PDAF_assim=.false.         !< Whether we use PDAF_assimilate
  INTEGER :: seedset=1                      !< Seed set for PDAF_generate_rndmat; can be set with PDAF_set_seedset
  LOGICAL :: new_seedset=.FALSE.            !< Whether the seetset was reset by PDAF_set_seedset

  ! *** Information variables for ensemble loop operations ***
  INTEGER :: member_save = 1                !< Store member index for query with PDAF_get_memberid
  INTEGER :: obs_member = 0                 !< Ensemble member when calling the observation operator routine

  ! *** Filter fields ***
  REAL, ALLOCATABLE :: state(:)             !< PE-local model state
  REAL, ALLOCATABLE :: Ainv(:,:)            !< Transform matrix or matrix of eigenvalues
  REAL, ALLOCATABLE :: bias(:)              !< Model bias vector
  REAL, TARGET, ALLOCATABLE :: ens(:,:)     !< Ensemble matrix
  REAL, TARGET, ALLOCATABLE :: sens(:,:,:)  !< Ensemble matrix holding past times for smoothing
  REAL, TARGET, ALLOCATABLE :: skewness(:)  !< Skewness of ensemble for each local domain
  REAL, TARGET, ALLOCATABLE :: kurtosis(:)  !< Kurtosis of ensemble for each local domain

!$OMP THREADPRIVATE(cnt_maxlag, obs_member, debug)

END MODULE PDAF_mod_core
