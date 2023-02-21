!$Id$
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
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  SAVE
!EOP

! *** Model- and data specific variables ***

  INTEGER :: dim_state           ! Global model state dimension
  INTEGER :: dim_state_p         ! Model state dimension for PE-local domain

  INTEGER :: dim_obs_p                    ! Process-local number of observations
  REAL, ALLOCATABLE    :: obs_p(:)        ! Vector holding process-local observations
  INTEGER, ALLOCATABLE :: obs_index_p(:)  ! Vector holding state-vector indices of observations
  REAL, ALLOCATABLE    :: obs_f(:)        ! Vector holding full vector of observations
  REAL, ALLOCATABLE :: coords_obs_f(:,:)  ! Array for full observation coordinates


! *** Below are the generic variables used for configuring PDAF ***
! *** Their values are set in init_PDAF                         ***

! !PUBLIC MEMBER FUNCTIONS:
! ! Settings for time stepping - available as command line options
  LOGICAL :: model_error   ! Control application of model error
  REAL    :: model_err_amp ! Amplitude for model error

! ! Settings for observations - available as command line options
  INTEGER :: delt_obs      ! time step interval between assimilation steps
  REAL    :: rms_obs       ! RMS error size for observation generation
  INTEGER :: dim_obs       ! Number of observations

! ! General control of PDAF - available as command line options
  INTEGER :: screen       ! Control verbosity of PDAF
                          ! (0) no outputs, (1) progess info, (2) add timings
                          ! (3) debugging output
  INTEGER :: dim_ens      ! Size of ensemble for SEIK/LSEIK/EnKF/ETKF
                          ! Number of EOFs to be used for SEEK
  INTEGER :: filtertype   ! Select filter algorithm:
                          !   SEEK (0), SEIK (1), EnKF (2), LSEIK (3), ETKF (4), LETKF (5)
                          !   ESTKF (6), LESTKF (7), LEnKF (8), NETF (9), LNETF (10), PF (12)
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
                          !   LEnKF:
                          !     (0) Standard form of EnKF with covariance localization
                          !   NETF:
                          !     (0) standard NETF 
                          !   LNETF:
                          !     (0) standard LNETF 
  INTEGER :: incremental  ! Perform incremental updating in LSEIK
  INTEGER :: dim_lag      ! Number of time instances for smoother

! ! Filter settings - available as command line options
!    ! General
  INTEGER :: type_forget  ! Type of forgetting factor
                          !   (0) fixed
                          !   (1) global adaptive
                          !   (2) local adaptive for LSEIK/LETKF/LESTKF
  REAL    :: forget       ! Forgetting factor for filter analysis
  INTEGER :: dim_bias     ! dimension of bias vector
!    ! SEEK
  INTEGER :: int_rediag   ! Interval to perform re-diagonalization in SEEK
  REAL    :: epsilon      ! Epsilon for gradient approx. in SEEK forecast
!    ! ENKF
  INTEGER :: rank_analysis_enkf  ! Rank to be considered for inversion of HPH
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
                    ! For LESTKF, LETKF, and LSEIK
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRADIUS
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
                    ! For LEnKF
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRADIUS
                    !   (2) 5th-order polynomial weight function
  REAL    :: sradius        ! Support radius for 5th order polynomial
                           !   or radius for 1/e for exponential weighting
!    ! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
  INTEGER :: type_sqrt     ! Type of the transform matrix square-root 
                    !   (0) symmetric square root, (1) Cholesky decomposition
!    ! NETF/LNETF
  INTEGER :: type_winf     ! Set weights inflation: (1) activate
  REAL    :: limit_winf    ! Limit for weights inflation: N_eff/N>limit_winf
!    ! hybrid LKNETF
  INTEGER :: type_hyb      ! Type of hybrid weight:
                    !   (0) use fixed hybrid weight hyb_gamma
                    !   (1) use gamma_lin: (1 - N_eff/N_e)*hyb_gamma
                    !   (2) use gamma_alpha: hybrid weight from N_eff/N>=hyb_gamma
                    !   (3) use gamma_ska: 1 - min(s,k)/sqrt(hyb_kappa) with N_eff/N>=hyb_gamma
                    !   (4) use gamma_sklin: 1 - min(s,k)/sqrt(hyb_kappa) >= 1-N_eff/N>=hyb_gamma
  REAL    :: hyb_gamma     ! Hybrid filter weight for state (1.0: LETKF, 0.0 LNETF)
  REAL    :: hyb_kappa     ! Hybrid norm for using skewness and kurtosis
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
  CHARACTER(len=110) :: filename  ! file name for assimilation output

!    ! Other variables - _NOT_ available as command line options!
  REAL    :: time          ! model time
  REAL :: coords_l(2)      ! Coordinates of local analysis domain
  INTEGER, ALLOCATABLE :: id_lstate_in_pstate(:) ! Indices of local state vector in PE-local global state vector
  INTEGER, ALLOCATABLE :: id_lobs_in_fobs(:)  ! Indices of local observations in full obs. vector
  REAL, ALLOCATABLE    :: distance_l(:)   ! Vector holding distances of local observations

!$OMP THREADPRIVATE(coords_l, id_lstate_in_pstate, id_lobs_in_fobs, distance_l)

END MODULE mod_assimilation
