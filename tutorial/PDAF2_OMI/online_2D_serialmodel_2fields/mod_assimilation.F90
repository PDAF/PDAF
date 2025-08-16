!>  Module for assimilation variables
!!
!! This module provides variables needed for the 
!! assimilation within the routines of the dummy model.
!! For simplicity, all assimilation-related variables
!! are stored here, even if they are only used in
!! the main program for the filter initialization.
!! Most variables can be specified as a command line 
!! argument.
!!
!! Implementation for the 2D online example
!! with or without parallelization.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE mod_assimilation

  IMPLICIT NONE
  SAVE


! *** Variables specific for online tutorial example ***

  INTEGER :: ensgroup=1    !< Type of initial ensemble

! *** Variables specific for model setup ***

  REAL :: coords_l(2)      !< Coordinates of local analysis domain

! *** Variables to handle multiple fields in the state vector ***

  INTEGER :: n_fields      !< number of fields in state vector
  INTEGER, ALLOCATABLE :: off_fields(:) !< Offsets of fields in state vector
  INTEGER, ALLOCATABLE :: dim_fields(:) !< Dimension of fields in state vector

  ! Declare Fortran type holding the indices of model fields in the state vector
  ! This can be extended to any number of fields - it severs to give each field a name
  TYPE field_ids
     INTEGER :: fieldA 
     INTEGER :: fieldB
  END TYPE field_ids

  ! Type variable holding field IDs in state vector
  TYPE(field_ids) :: id

!$OMP THREADPRIVATE(coords_l)


! -----------------------------------------------------------------
! *** Below are the generic variables used for configuring PDAF ***
! *** Their values are set in init_PDAF                         ***

! Settings for state vector size
  INTEGER :: dim_state     !< Global model state dimension
  INTEGER :: dim_state_p   !< Model state dimension for PE-local domain

! Settings for time stepping - available as command line options
  LOGICAL :: model_error   !< Control application of model error
  REAL    :: model_err_amp !< Amplitude for model error

! Settings for observations - available as command line options
  INTEGER :: delt_obs      !< time step interval between assimilation steps
  LOGICAL :: twin_experiment  !< Whether to run an twin experiment with synthetic observations

! General control of PDAF - available as command line options
  INTEGER :: screen       !< Control verbosity of PDAF
                          !< * (0) no outputs
                          !< * (1) progress info
                          !< * (2) add timings
                          !< * (3) debugging output
  INTEGER :: dim_ens      !< Size of ensemble
  INTEGER :: filtertype   !< Select filter algorithm:
                          !<   * SEIK (1), EnKF (2), LSEIK (3), ETKF (4)
                          !<   LETKF (5), ESTKF (6), LESTKF (7), NETF (9), LNETF (10)
                          !<   LKNETF (11), PF (12), GENOBS (100), 3DVAR (200)
  INTEGER :: subtype      !< Subtype of filter algorithm
                          !<   * SEIK:
                          !<     (0) ensemble forecast; new formulation
                          !<     (1) ensemble forecast; old formulation
                          !<     (2) fixed error space basis
                          !<     (3) fixed state covariance matrix
                          !<     (4) SEIK with ensemble transformation
                          !<   * EnKF:
                          !<     (0) analysis for large observation dimension
                          !<     (1) analysis for small observation dimension
                          !<   * LSEIK:
                          !<     (0) ensemble forecast;
                          !<     (2) fixed error space basis
                          !<     (3) fixed state covariance matrix
                          !<     (4) LSEIK with ensemble transformation
                          !<   * ETKF:
                          !<     (0) ETKF using T-matrix like SEIK
                          !<     (1) ETKF following Hunt et al. (2007)
                          !<       There are no fixed basis/covariance cases, as
                          !<       these are equivalent to SEIK subtypes 2/3
                          !<   * LETKF:
                          !<     (0) LETKF using T-matrix like SEIK
                          !<     (1) LETKF following Hunt et al. (2007)
                          !<       There are no fixed basis/covariance cases, as
                          !<       these are equivalent to LSEIK subtypes 2/3
                          !<   * ESTKF:
                          !<     (0) standard ESTKF 
                          !<       There are no fixed basis/covariance cases, as
                          !<       these are equivalent to SEIK subtypes 2/3
                          !<   * LESTKF:
                          !<     (0) standard LESTKF 
                          !<       There are no fixed basis/covariance cases, as
                          !<       these are equivalent to LSEIK subtypes 2/3
                          !<   * NETF:
                          !<     (0) standard NETF 
                          !<   * LNETF:
                          !<     (0) standard LNETF 
                          !<   * LKNETF:
                          !<     (0) HNK: 2-step LKNETF with NETF before LETKF
                          !<     (1) HKN: 2-step LKNETF with LETKF before NETF
                          !<     (4) HSync: LKNETF synchronous
                          !<   * PF:
                          !<     (0) standard PF 
                          !<   * 3D-Var:
                          !<     (0) parameterized 3D-Var
                          !<     (1) 3D Ensemble Var using LESTKF for ensemble update
                          !<     (4) 3D Ensemble Var using ESTKF for ensemble update
                          !<     (6) hybrid 3D-Var using LESTKF for ensemble update
                          !<     (7) hybrid 3D-Var using ESTKF for ensemble update
  INTEGER :: incremental  !< SEIK/LSEIK: (1) Perform incremental updating
  INTEGER :: dim_lag      !< Number of time instances for smoother

! Filter settings - available as command line options
!    ! General
  INTEGER :: type_forget  !< Type of forgetting factor
                          !<  SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
                          !<   (0) fixed
                          !<   (1) global adaptive
                          !<   (2) local adaptive for LSEIK/LETKF/LESTKF
                          !<  NETF/LNETF/PF
                          !<   (0) apply inflation on forecast ensemble
                          !<   (2) apply inflation on analysis ensemble
  REAL    :: forget       !< Forgetting factor for filter analysis
  INTEGER :: dim_bias     !< dimension of bias vector
!    ! All localized filters
  REAL    :: cradius       !< Cut-off radius for local observation domain
  INTEGER :: locweight     !< * Type of localizing weighting of observations
                           !<   (0) constant weight of 1
                           !<   (1) exponentially decreasing with SRADIUS
                           !<   (2) use 5th-order polynomial
                           !<   (3) regulated localization of R with mean error variance
                           !<   (4) regulated localization of R with single-point error variance
  REAL    :: sradius       !< Support radius for 5th order polynomial
                           !<   or radius for 1/e for exponential weighting
!    ! ENKF
  INTEGER :: rank_ana_enkf !< Rank to be considered for inversion of HPH in analysis of EnKF
                           !<  (0) for analysis w/o eigendecomposition
!    ! SEIK/ETKF/ESTKF/LSEIK/LETKF/LESTKF/NETF/LNETF/LKNETF
  INTEGER :: type_trans    !< Type of ensemble transformation 
                           !< * SEIK/LSEIK: 
                           !< (0) use deterministic omega
                           !< (1) use random orthonormal omega orthogonal to (1,...,1)^T
                           !< (2) use product of (0) with random orthonormal matrix with
                           !<     eigenvector (1,...,1)^T 
                           !< * ETKF/LETKF with subtype=4: 
                           !< (0) use deterministic symmetric transformation
                           !< (2) use product of (0) with random orthonormal matrix with
                           !<     eigenvector (1,...,1)^T 
                           !< * ESTKF/LESTKF:
                           !< (0) use deterministic omega
                           !< (1) use random orthonormal omega orthogonal to (1,...,1)^T
                           !< (2) use product of (0) with random orthonormal matrix with
                           !<     eigenvector (1,...,1)^T
                           !< * NETF/LNETF:
                           !< (0) use random orthonormal transformation orthogonal to (1,...,1)^T
                           !< (1) use identity transformation
                           !< * LKNETF:
                           !< (0) use random orthonormal transformation orthogonal to (1,...,1)^T
                           !< (1) use identity transformation
!    ! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
  INTEGER :: type_sqrt     !< * Type of the transform matrix square-root 
                           !<   (0) symmetric square root
                           !<   (1) Cholesky decomposition
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

!    ! Other variables - _NOT_ available as command line options!
  REAL    :: time               !< model time

END MODULE mod_assimilation
