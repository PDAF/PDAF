!$Id$
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
!! Implementation for the 2D offline example
!! with or without parallelization.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE mod_assimilation

  IMPLICIT NONE
  SAVE

! *** Variables specific for offline tutorial example ***

  INTEGER :: nx, ny                     !< Size of 2D grid
  INTEGER, ALLOCATABLE :: local_dims(:) !< Array for local state dimensions

  ! *** Specific for the 3D-Var tutorial cases
  REAL, ALLOCATABLE    :: Vmat_p(:,:)            !< square-root of P for 3D-Var
  REAL, ALLOCATABLE    :: Vmat_ens_p(:,:)        !< square-root of P for ensemble 3D-Var

! *** Model- and data specific variables ***

  INTEGER :: dim_state     !< Global model state dimension
  INTEGER :: dim_state_p   !< Model state dimension for PE-local domain


! *** Below are the generic variables used for configuring PDAF ***
! *** Their values are set in init_PDAF_offline                 ***

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
  INTEGER :: dim_ens      !< Size of ensemble for SEIK/LSEIK/EnKF/ETKF \n
                          !< Number of EOFs to be used for SEEK
  INTEGER :: filtertype   !< Select filter algorithm:
                          !<   * SEEK (0), SEIK (1), EnKF (2), LSEIK (3), ETKF (4)
                          !<   LETKF (5), ESTKF (6), LESTKF (7), NETF (9), LNETF (10)
                          !<   PF (12), GENOBS (100), 3DVAR (200)
  INTEGER :: subtype      !< Subtype of filter algorithm
                          !<   * SEEK: 
                          !<     (0) evolve normalized modes
                          !<     (1) evolve scaled modes with unit U
                          !<     (2) fixed basis (V); variable U matrix
                          !<     (3) fixed covar matrix (V,U kept static)
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
                          !<   * PF:
                          !<     (0) standard PF 
                          !<   * 3D-Var:
                          !<     (0) parameterized 3D-Var
                          !<     (1) 3D Ensemble Var using LESTKF for ensemble update
                          !<     (4) 3D Ensemble Var using ESTKF for ensemble update
                          !<     (6) hybrid 3D-Var using LESTKF for ensemble update
                          !<     (7) hybrid 3D-Var using ESTKF for ensemble update
  INTEGER :: incremental  !< Perform incremental updating in LSEIK
  INTEGER :: dim_lag      !< Number of time instances for smoother

! Filter settings - available as command line options
!    ! General
  INTEGER :: type_forget  !< Type of forgetting factor
  REAL    :: forget       !< Forgetting factor for filter analysis
  INTEGER :: dim_bias     !< dimension of bias vector
!    ! SEEK
  INTEGER :: int_rediag   !< Interval to perform re-diagonalization in SEEK
  REAL    :: epsilon      !< Epsilon for gradient approx. in SEEK forecast
!    ! ENKF
  INTEGER :: rank_analysis_enkf  !< Rank to be considered for inversion of HPH
!    ! SEIK/ETKF/ESTKF/LSEIK/LETKF/LESTKF
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
!    ! LSEIK/LETKF/LESTKF
  REAL    :: local_range   !< Range for local observation domain
  INTEGER :: locweight     !< Type of localizing weighting of observations
                    !<   * (0) constant weight of 1
                    !<   * (1) exponentially decreasing with SRANGE
                    !<   * (2) use 5th-order polynomial
                    !<   * (3) regulated localization of R with mean error variance
                    !<   * (4) regulated localization of R with single-point error variance
  REAL    :: srange        !< Support range for 5th order polynomial
                           !<   or radius for 1/e for exponential weighting
!    ! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
  INTEGER :: type_sqrt     !< Type of the transform matrix square-root 
                    !<   * (0) symmetric square root
                    !<   * (1) Cholesky decomposition
!    ! 3D-Var
  INTEGER :: type_3dvar !< Type of 3D-Var method
                    !<   * (0) Parameterized 3D-Var
                    !<   * (1) Ensemble 3D-Var using LETKF for ensemble transformation
                    !<   * (4) Ensemble 3D-Var using global ETKF for ensemble transformation
                    !<   * (6) Hybrid 3D-Var using LETKF for ensemble transformation
                    !<   * (7) Hybrid 3D-Var using global ETKF for ensemble transformation
  INTEGER :: type_opt      !< Type of minimizer for 3DVar
                    !<   * (1) LBFGS (default)
                    !<   * (2) CG+
                    !<   * (3) plain CG
                    !<   * (12) CG+ parallelized
                    !<   * (13) plain CG parallelized
  INTEGER :: dim_cvec = 0  !< Size of control vector (parameterized part; for subtypes 0,1)
  INTEGER :: dim_cvec_ens = 0   !< Size of control vector (ensemble part; for subtypes 1,2)
  INTEGER :: mcols_cvec_ens = 1 !< Multiplication factor for number of columns for ensemble control vector
  REAL :: beta_3dvar = 0.5 !< Hybrid weight for hybrid 3D-Var
!    ! NETF/LNETF
  INTEGER :: type_winf     ! Set weights inflation: (1) activate
  REAL    :: limit_winf    ! Limit for weights inflation: N_eff/N>limit_winf
!    ! Particle filter
  INTEGER :: pf_res_type   ! Resampling type for PF
                           ! (1) probabilistic resampling
                           ! (2) stochastic universal resampling
                           ! (3) residual resampling        
  INTEGER :: pf_noise_type    ! Resampling type for PF
                           ! (0) no perturbations, (1) constant stddev, 
                           ! (2) amplitude of stddev relative of ensemble variance
  REAL :: pf_noise_amp     ! Noise amplitude (>=0.0, only used if pf_noise_type>0)

!    ! File output - available as a command line option
  CHARACTER(len=110) :: filename  !< file name for assimilation output

!    ! Other variables - _NOT_ available as command line options!
  INTEGER :: covartype     !< For SEIK: Definition of ensemble covar matrix
                           !<   * (0): Factor (r+1)^-1 (or N^-1)
                           !<   * (1): Factor r^-1 (or (N-1)^-1) - real ensemble covar.
                           !< This setting is only for the model part; The definition
                           !< of P has also to be specified in PDAF_filter_init.
                           !< Only for upward-compatibility of PDAF!
  REAL    :: time          !< model time

  REAL :: coords_l(2)      !< Coordinates of local analysis domain
  INTEGER, ALLOCATABLE :: id_lstate_in_pstate(:) !< Indices of local state vector in PE-local global state vector

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

!$OMP THREADPRIVATE(coords_l, id_lstate_in_pstate)

END MODULE mod_assimilation
