!$Id: mod_assimilation.F90 1866 2017-12-21 09:05:27Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_assimilation

! !DESCRIPTION:
! This module provides variables needed for the 
! assimilation within the routines of the dummy model.
! For simplicity, all assimilation-related variables
! are stored here, even if they are only used in
! the main program for the filter initialization.
! Most variables can be specified as a command line 
! argument.
!
! Implementation for the 2D offline example
! without parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  SAVE
!EOP

! *** Variables specific for offline tutorial example ***

  INTEGER :: nx, ny              ! Size of 2D grid

! *** Model- and data specific variables ***

  INTEGER :: dim_state           ! Global model state dimension
  INTEGER :: dim_state_p         ! Model state dimension for PE-local domain
  INTEGER, ALLOCATABLE :: local_dims(:)  ! Array for local state dimensions

  INTEGER :: dim_obs_p                    ! Process-local number of observations
  REAL, ALLOCATABLE    :: obs_p(:)        ! Vector holding process-local observations
  INTEGER, ALLOCATABLE :: obs_index_p(:)  ! Vector holding state-vector indices of observations
  REAL, ALLOCATABLE    :: obs_f(:)        ! Vector holding full vector of observations
  REAL, ALLOCATABLE :: coords_obs_f(:,:)  ! Array for full observation coordinates


! *** Below are the generic variables used for configuring PDAF ***
! *** Their values are set in init_PDAF_offline                 ***

! !PUBLIC MEMBER FUNCTIONS:
! ! Settings for time stepping - available as command line options
  LOGICAL :: model_error   ! Control application of model error
  REAL    :: model_err_amp ! Amplitude for model error

! ! Settings for observations - available as command line options
  INTEGER :: delt_obs      ! time step interval between assimilation steps
  LOGICAL :: twin_experiment   !< Whether to run an twin experiment with synthetic observations
  REAL    :: rms_obs       ! RMS error size for observation generation
  INTEGER :: dim_obs       ! Number of observations

! ! General control of PDAF - available as command line options
  INTEGER :: screen       ! Control verbosity of PDAF
                          ! (0) no outputs, (1) progess info, (2) add timings
                          ! (3) debugging output
  INTEGER :: dim_ens      ! Size of ensemble for SEIK/LSEIK/EnKF/ETKF
                          ! Number of EOFs to be used for SEEK
  INTEGER :: filtertype   ! Select filter algorithm:
                          !   SEEK (0), SEIK (1), EnKF (2), LSEIK (3), ETKF (4)
                          !   LETKF (5), ESTKF (6), LESTKF (7), NETF (9), LNETF (10)
  INTEGER :: subtype      ! Subtype of filter algorithm
                          !   SEEK: 
                          !     (0) evolve normalized modes
                          !     (1) evolve scaled modes with unit U
                          !     (2) fixed basis (V); variable U matrix
                          !     (3) fixed covar matrix (V,U kept static)
                          !   SEIK:
                          !     (0) ensemble forecast; new formulation
                          !     (1) ensemble forecast; old formulation
                          !     (2) fixed error space basis
                          !     (3) fixed state covariance matrix
                          !     (4) SEIK with ensemble transformation
                          !   EnKF:
                          !     (0) analysis for large observation dimension
                          !     (1) analysis for small observation dimension
                          !   LSEIK:
                          !     (0) ensemble forecast;
                          !     (2) fixed error space basis
                          !     (3) fixed state covariance matrix
                          !     (4) LSEIK with ensemble transformation
                          !   ETKF:
                          !     (0) ETKF using T-matrix like SEIK
                          !     (1) ETKF following Hunt et al. (2007)
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to SEIK subtypes 2/3
                          !   LETKF:
                          !     (0) LETKF using T-matrix like SEIK
                          !     (1) LETKF following Hunt et al. (2007)
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to LSEIK subtypes 2/3
                          !   ESTKF:
                          !     (0) standard ESTKF 
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to SEIK subtypes 2/3
                          !   LESTKF:
                          !     (0) standard LESTKF 
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to LSEIK subtypes 2/3
                          !   NETF:
                          !     (0) standard NETF 
                          !   LNETF:
                          !     (0) standard LNETF 
  INTEGER :: incremental  ! Perform incremental updating in LSEIK
  INTEGER :: dim_lag      ! Number of time instances for smoother

! ! Filter settings - available as command line options
!    ! General
  INTEGER :: type_forget  ! Type of forgetting factor
  REAL    :: forget       ! Forgetting factor for filter analysis
  INTEGER :: dim_bias     ! dimension of bias vector
!    ! SEEK
  INTEGER :: int_rediag   ! Interval to perform re-diagonalization in SEEK
  REAL    :: epsilon      ! Epsilon for gradient approx. in SEEK forecast
!    ! ENKF
  INTEGER :: rank_ana_enkf  ! Rank to be considered for inversion of HPH
!    ! SEIK/ETKF/ESTKF/LSEIK/LETKF/LESTKF
  INTEGER :: type_trans    ! Type of ensemble transformation
                           ! SEIK/LSEIK:
                           ! (0) use deterministic omega
                           ! (1) use random orthonormal omega orthogonal to (1,...,1)^T
                           ! (2) use product of (0) with random orthonormal matrix with
                           !     eigenvector (1,...,1)^T
                           ! ETKF/LETKF with subtype=4:
                           ! (0) use deterministic symmetric transformation
                           ! (2) use product of (0) with random orthonormal matrix with
                           !     eigenvector (1,...,1)^T
                           ! ESTKF/LESTKF:
                           ! (0) use deterministic omega
                           ! (1) use random orthonormal omega orthogonal to (1,...,1)^T
                           ! (2) use product of (0) with random orthonormal matrix with
                           !     eigenvector (1,...,1)^T
                           ! NETF/LNETF:
                           ! (0) use random orthonormal transformation orthogonal to (1,...,1)^T
                           ! (1) use identity transformation
!    ! LSEIK/LETKF/LESTKF
  REAL    :: cradius       ! Cut-off radius for local observation domain
  INTEGER :: locweight     ! Type of localizing weighting of observations
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRADIUS
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
  REAL    :: sradius       ! Support radius for 5th order polynomial
                           !   or radius for 1/e for exponential weighting
!    ! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
  INTEGER :: type_sqrt     ! Type of the transform matrix square-root 
                    !   (0) symmetric square root, (1) Cholesky decomposition
!    ! NETF/LNETF/PF
  INTEGER :: type_winf     !< Set weights inflation: 
                           !<   (0) no weights inflation
                           !<   (1) use N_eff/N>limit_winf
  REAL    :: limit_winf    !< Limit for weights inflation: N_eff/N>limit_winf
!    ! hybrid LKNETF
  INTEGER :: type_hyb      !< * Type of hybrid weight:
                           !<   (0) use fixed hybrid weight hyb_gamma
                           !<   (1) use gamma_lin: (1 - N_eff/N_e)*hyb_gamma
                           !<   (2) use gamma_alpha: hybrid weight from N_eff/N>=hyb_gamma
                           !<   (3) use gamma_ska: 1 - min(s,k)/sqrt(hyb_kappa) with N_eff/N>=hyb_gamma
                           !<   (4) use gamma_sklin: 1 - min(s,k)/sqrt(hyb_kappa) >= 1-N_eff/N>=hyb_gamma
  REAL    :: hyb_gamma     !< Hybrid filter weight for state (1.0: LETKF, 0.0 LNETF)
  REAL    :: hyb_kappa     !< Hybrid norm for using skewness and kurtosis
!    ! Particle filter
  INTEGER :: pf_res_type   !< * Resampling type for PF
                           !<   (1) probabilistic resampling
                           !<   (2) stochastic universal resampling
                           !<   (3) residual resampling        
  INTEGER :: pf_noise_type !< * Resampling type for PF
                           !<   (0) no perturbations, (1) constant stddev, 
                           !<   (2) amplitude of stddev relative of ensemble variance
  REAL :: pf_noise_amp     !< Noise amplitude (>=0.0, only used if pf_noise_type>0)

!    ! 3D-Var
  INTEGER :: type_opt      !< * Type of minimizer for 3DVar
                           !<   (1) LBFGS (default)
                           !<   (2) CG+
                           !<   (3) plain CG
                           !<   (12) CG+ parallelized
                           !<   (13) plain CG parallelized
  INTEGER :: dim_cvec = 0  !< Size of control vector (parameterized part; for subtypes 0,1)
  INTEGER :: dim_cvec_ens = 0   !< Size of control vector (ensemble part; for subtypes 1,2)
  INTEGER :: mcols_cvec_ens = 1 !< Multiplication factor for number of columns for ensemble control vector
  REAL :: beta_3dvar = 0.5 !< Hybrid weight for hybrid 3D-Var
  INTEGER :: solver_iparam1 = 2 !< Solver specific parameter
                                !<  LBFGS: parameter m (default=5)
                                !<       Number of corrections used in limited memory matrix; 3<=m<=20
                                !<  CG+: parameter method (default=2)
                                !<       (1) Fletcher-Reeves, (2) Polak-Ribiere, (3) positive Polak-Ribiere
                                !<  CG: maximum number of iterations (default=200)
  INTEGER :: solver_iparam2 = 1 !< Solver specific parameter
                                !<  LBFGS: - not used - 
                                !<  CG+: parameter irest (default=1)
                                !<       (0) no restarts; (n>0) restart every n steps
                                !<  CG: - not used -
  REAL :: solver_rparam1 = 1.0e-6 !< Solver specific parameter
                                !<  LBFGS: limit for stopping iterations 'pgtol' (default=1.0e-5)
                                !<  CG+: convergence parameter 'eps' (default=1.0e-5)
                                !<  CG: conpergence parameter 'eps' (default=1.0e-6)
  REAL :: solver_rparam2 = 1.0e+7 !< Solver specific parameter
                                !<  LBFGS: tolerance in termination test 'factr' (default=1.0e+7) 
                                !<  CG+: - not used -
                                !<  CG: - not used -

!    ! File output - available as a command line option
  CHARACTER(len=110) :: filename  ! file name for assimilation output

!    ! Other variables - _NOT_ available as command line options!
  INTEGER :: covartype     ! For SEIK: Definition of ensemble covar matrix
                           ! (0): Factor (r+1)^-1 (or N^-1)
                           ! (1): Factor r^-1 (or (N-1)^-1) - real ensemble covar.
                           ! This setting is only for the model part; The definition
                           ! of P has also to be specified in PDAF_filter_init.
                           ! Only for upward-compatibility of PDAF!
  REAL    :: time          ! model time
  REAL :: coords_l(2)      ! Coordinates of local analysis domain
  INTEGER, ALLOCATABLE :: id_lstate_in_pstate(:) ! Indices of local state vector in PE-local global state vector
  INTEGER, ALLOCATABLE :: id_lobs_in_fobs(:)  ! Indices of local observations in full obs. vector
  REAL, ALLOCATABLE    :: distance_l(:)   ! Vector holding distances of local observations

!$OMP THREADPRIVATE(coords_l, id_lstate_in_pstate, id_lobs_in_fobs, distance_l)

END MODULE mod_assimilation
