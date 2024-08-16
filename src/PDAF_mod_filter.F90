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
! !MODULE:
MODULE PDAF_mod_filter
  
! !DESCRIPTION:
! This module provides variables shared between the
! subroutines of PDAF.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2003-06 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  SAVE

! !PUBLIC DATA MEMBERS:
  INTEGER :: dim_eof       ! Rank (number of columns of eofV in SEEK)
  INTEGER :: dim_ens       ! Ensemble size 
  INTEGER :: rank          ! Rank of initial covariance matrix
  INTEGER :: dim_p         ! State dimension for PE-local domain
  INTEGER :: dim_bias_p=0  ! Dimension of bias vector
  REAL    :: forget        ! Forgetting factor
  LOGICAL :: offline_mode=.false.   ! Wether to use PDAF offline mode
  INTEGER :: type_filter   ! Type of Filter
                           ! (0) SEEK  (Pham et al., 1998a)
                           ! (1) SEIK  (Pham et al., 1998b)
                           ! (2) EnKF  (Evensen, 1994)
                           ! (3) LSEIK (Nerger et al., 2007)
                           ! (4) ETKF  (Bishop et al., 2001) 
                           !     (ETKF uses symmetric square roots like LETKF)
                           ! (5) LETKF (Hunt et al., 2007)
  INTEGER :: subtype_filter  ! Sub-type of Filter
                   ! Subtype of SEEK: 
                   !     (0) Evolve with finite difference approx to TLM
                   !     (1) Scaled modes, unit U
                   !     (2) Fixed basis (V); variable U matrix
                   !     (3) Fixed covar matrix (V,U kept static)
                   !     (5) PDAF offline mode
                   ! Subtype of SEIK: 
                   !     (0) Usual SEIK with mean forecast, new formulation;
                   !     (1) Usual SEIK with mean forecast, old formulation;
                   !     (2) Fixed error space basis
                   !     (3) Fixed state covariance matrix
                   !     (4) SEIK with ensemble transformation (like ETKF)
                   !     (5) PDAF offline mode
                   ! Subtype of EnKF forecast and update step: 
                   !     (0) Mean forecast & representer analysis for large dim_obs;
                   !     (1) Mean forecast & representer analysis for small dim_obs;
                   !     (5) PDAF offline mode
                   ! Subtype of LSEIK: 
                   !     (0) Mean forecast;
                   !     (2) Fixed error space basis
                   !     (3) Fixed state covariance matrix
                   !     (4) LSEIK with ensemble transformation (like LETKF)
                   !     (5) PDAF offline mode
                   ! Subtype of ETKF:
                   !     (0) ETKF using T-matrix like SEIK
                   !     (1) ETKF following Hunt et al. (2007)
                   !       There are no fixed basis/covariance cases, as
                   !       these are equivalent to SEIK subtypes 2/3
                   !     (5) PDAF offline mode
                   ! Subtype of LETKF:
                   !     (0) ETKF using T-matrix like SEIK
                   !       There are no fixed basis/covariance cases, as
                   !       these are equivalent to LSEIK subtypes 2/3
                   !     (5) PDAF offline mode
  INTEGER :: type_trans=0  ! Type of ensemble transformation
                           ! For SEIK/LSEIK:
                           ! (0) use deterministic Omega
                           ! (1) use random orthonormal Omega orthogonal to (1,...,1)^T
                           ! (2) use product of (0) with random orthonomal matrix with
                           !     eigenvector (1,...,1)^T
                           ! For ETKF/LETKF:
                           ! (0) use deterministic symmetric transformation
                           ! (2) use product of (0) with random orthonomal matrix with
                           !     eigenvector (1,...,1)^T
  INTEGER :: step          ! Current time step
  INTEGER :: step_obs      ! Time step of next observation
  INTEGER :: dim_obs       ! Dimension of next observation
  INTEGER :: screen=0      ! Control verbosity of filter routines
                   ! (0) quiet; (1) normal output; (2); plus timings; (3) debug output
  INTEGER :: debug=0       ! Debugging flag: print debug information if >0
  INTEGER :: incremental=0 ! Whether to perform incremental updating
  INTEGER :: type_forget=0 ! Type of forgetting factor
                           ! (0): fixed; (1) global adaptive; (2) local adaptive
  INTEGER :: type_sqrt=0   ! Type of sqrt of U in SEIK/LSEIK-trans or A in ESTKF/LESTKF
                           ! (0): symmetric sqrt; (1): Cholesky decomposition
                           ! In SEIK/LSEIK the default is 1
  INTEGER :: dim_lag = 0   ! Number of past time instances considered for smoother

  ! SEEK
  INTEGER :: int_rediag=1  ! Interval for perform rediagonalization (SEEK)
  REAL    :: epsilon=0.1   ! Epsilon for approximated TLM evolution

  ! EnKF/LEnKF
  INTEGER :: rank_ana_enkf ! Rank to be considered for inversion of HPH
                           !   in analysis of EnKF
  ! NETF and PF
  INTEGER :: type_winf=0   ! Type of weights inflation for NETF
                           ! (0): none; (1) inflate for N_eff/N > limit_winf
  REAL :: limit_winf = 0.0 ! Limit to weights inflation

  ! LKNETF
  INTEGER :: type_hyb = 0  ! Type of hybrid weight: (2) adaptive
  REAL :: hyb_g = 1.0      ! Hybrid weight for state in LKNEF (1.0 for LETKF; 0.0 for LNETF)
  REAL :: hyb_k = 50.0     ! Hybrid weight norm for using skewness and kurtosis
  LOGICAL :: store_rndmat = .false.  ! Whether to recompute or store the random matrix

  ! PF
  INTEGER :: restype = 1     ! Resampling type for particle filters
                             ! (1) probabilistic resampling, (2) stochastic universal resampling
                             ! (3) residual resampling
  INTEGER :: noise_type = 0  ! Type of perturbing noise in PF
                             ! (1) constant variance, (2) amplitude relative to ensemble std.
  REAL :: pf_noise_amp = 0.0 ! Amplitudy of noise in PF
  
  ! Variational
  INTEGER :: type_opt = 0     ! Type of minimizer for 3DVar
                              ! (0) LBFGS, (1) CG+, (-1) steepest descent
  INTEGER :: dim_cvec = 0     ! Size of control vector (fixed part)
  INTEGER :: dim_cvec_ens = 0 ! Size of control vector (ensemble part)
  REAL :: beta_3dvar = 0.5    ! Hybrid weight for hybrid 3D-Var
  INTEGER :: m_lbfgs_var=5        ! Parameter 'm' of LBFGS
  INTEGER :: method_cgplus_var=2  ! Parameter 'method' of CG+
  INTEGER :: irest_cgplus_var=1   ! Parameter 'irest' of CG+
  INTEGER :: maxiter_cg_var=200   ! Parameter 'maxiter' of CG
  REAL :: eps_cg_var = 1.0e-6     ! Parameter 'EPS' of  CG
  REAL :: eps_cgplus_var = 1.0e-5  ! Parameter 'EPS' of CG+
  REAL :: pgtol_lbfgs_var=1.0e-5  ! Parameter 'pgtol' of LBFGS
  REAL :: factr_lbfgs_var=1.0e7   ! Parameter 'factr' of LBFGS

  ! *** Control variables for filter ***
  INTEGER :: firsttime = 1  ! Are the filter routines called for the first time?
  INTEGER :: initevol = 1   ! Initialize a new forecast phase?
  INTEGER :: member = 1     ! Which member of sub-ensemble to evolve
  INTEGER :: member_get = 1 ! Which member of sub-ensemble to evolve (used in PDAF_get_state)
  INTEGER :: member_save = 1 ! Store member index for quewry with PDAF_get_memberid
  INTEGER :: nsteps         ! Number of time steps to perform
  INTEGER :: cnt_steps      ! Number of time steps in current forecast phase
  INTEGER :: end_forecast   ! Whether to exit the forecasting
  INTEGER :: local_dim_ens  ! Local ensemble sizes (including state forecast)
  INTEGER :: flag = 0       ! Status flag
  INTEGER :: cnt_maxlag = 0 ! Count maximum number of past time instances for smoother
  INTEGER :: Nm1vsN=1       ! Flag which definition of P ist used in SEIK
  INTEGER :: obs_member=0   ! Ensemble member when calling the observation operator routine
  LOGICAL :: observe_ens=.false.  ! Whether (F) to apply H to ensemble mean to compute residual
                            ! or (T) apply H to X, compute mean of HX and then residual
  INTEGER :: assim_flag=0   ! (1) if assimilation done at this time step, (0) if not
  ! (0): Factor N^-1; (1): Factor (N-1)^-1 - Recommended is 1 for 
  ! a real ensemble filter, 0 is for compatibility with older PDAF versions
  LOGICAL :: ensemblefilter ! Whether the chosen filter is ensemble-based
  INTEGER :: localfilter = 0 ! Whether the chosen filter is domain-localized (1: yes)
  INTEGER :: globalobs = 0  ! Whether the chosen filter needs global observations (1: yes)
  CHARACTER(len=10) :: filterstr   ! String defining the filter type
  REAL    :: forget_l       ! Forgetting factor in local analysis loop
  LOGICAL :: inloop=.false. ! Whether the program is in the local analysis loop
  LOGICAL :: use_PDAF_assim = .false. ! Whether we use PDAF_assimilate

  ! *** Filter fields ***
  REAL, ALLOCATABLE :: state(:)     ! PE-local model state
  REAL, ALLOCATABLE :: state_inc(:) ! PE-local analysis increment for inc. updating
  REAL, ALLOCATABLE :: eofU(:,:)    ! Matrix of eigenvalues from EOF computation
  REAL, TARGET, ALLOCATABLE :: eofV(:,:)    ! Ensemble matrix
                                            !   or matrix of eigenvectors from EOF computation
  REAL, TARGET, ALLOCATABLE :: sens(:,:,:)  ! Ensemble matrix holding past times for smoothing
  REAL, TARGET, ALLOCATABLE :: skewness(:)  ! Skewness of ensemble for each local domain
  REAL, TARGET, ALLOCATABLE :: kurtosis(:)  ! Kurtosis of ensemble for each local domain
  REAL, ALLOCATABLE :: bias(:)      ! Model bias vector
!EOP

!$OMP THREADPRIVATE(cnt_maxlag, obs_member, forget_l, debug)

END MODULE PDAF_mod_filter
