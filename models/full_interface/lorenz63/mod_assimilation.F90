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

! !PUBLIC MEMBER FUNCTIONS:
! ! Settings for time stepping - available as command line options
  LOGICAL :: model_error   ! Control application of model error
  REAL    :: model_err_amp ! Amplitude for model error

! ! Settings for observations - available as command line options
  INTEGER :: delt_obs      ! time step interval between assimilation steps
  REAL    :: rms_obs       ! RMS error size for observation generation
  LOGICAL :: use_obs_mask  ! Whether to use a mask for observation gaps
  LOGICAL :: use_maskfile  ! Whether to read mask from a file
  INTEGER :: numobs        ! If not read from file, use this number of obs. (1 to numobs)
  INTEGER :: dx_obs        ! Of not read from file, use this grid point distance of obs.
  INTEGER :: obs_err_type  ! Type of observation error: (0) Gaussian, (1) double-exponential
  LOGICAL :: twin_experiment  ! Wether to run an twin experiment with synthetic observations

! ! General control of PDAF - available as command line options
  INTEGER :: screen       ! Control verbosity of PDAF
                          ! (0) no outputs, (1) progess info, (2) add timings
                          ! (3) debugging output
  INTEGER :: dim_ens      ! Size of ensemble for SEIK/LSEIK/EnKF/ETKF
                          ! Number of EOFs to be used for SEEK
  INTEGER :: filtertype   ! Select filter algorithm:
                          ! (0) SEEK, (1) SEIK, (2) EnKF, (3) LSEIK, 
                          ! (4) ETKF, (5) LETKF, (6) ESTKF, (7) LESTKF
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
                          !     (0) ETKF using T-matrix like SEIK
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
  INTEGER :: dim_lag      ! Number of time instances for smoother

! ! Filter settings - available as command line options
!    ! General
  INTEGER :: type_forget  ! Type of forgetting factor
  REAL    :: forget       ! Forgetting factor for filter analysis
  INTEGER :: dim_bias     ! dimension of bias vector
  CHARACTER(len=3) :: type_ensinit ! Type of ensemble initialization:
                          ! 'eof': Initialize by 2nd-order exact sampling from EOFs
                          ! 'rnd': Initialize by random sampling from state trajectory
  INTEGER :: seedset      ! Select set of seeds for random numbers (only for 'rnd')
  INTEGER :: init_dt      ! Time step interval used for random ensemble initialization
  INTEGER :: init_maxtime ! Maximum time step to pick from in random ensemble init
  INTEGER :: stepnull_means=0 ! Steps at which computation of time mean error is started
!    ! SEEK
  INTEGER :: int_rediag   ! Interval to perform re-diagonalization in SEEK
  REAL    :: epsilon      ! Epsilon for gradient approx. in SEEK forecast
!    ! ENKF
  INTEGER :: rank_analysis_enkf  ! Rank to be considered for inversion of HPH
!    ! SEIK/ETKF/ESTKF/LSEIK/LETKF/LESTKF
  INTEGER :: type_trans    ! Type of ensemble transformation
                           ! SEIK/LSEIK and ESTKF/LESTKF:
                           ! (0) use deterministic omega
                           ! (1) use random orthonormal omega orthogonal to (1,...,1)^T
                           ! (2) use product of (0) with random orthonormal matrix with
                           !     eigenvector (1,...,1)^T
                           ! ETKF/LETKF:
                           ! (0) use deterministic symmetric transformation
                           ! (2) use product of (0) with random orthonormal matrix with
                           !     eigenvector (1,...,1)^T
!    ! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
  INTEGER :: type_sqrt     ! Type of the transform matrix square-root 
                    !   (0) symmetric square root, (1) Cholesky decomposition
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

!    ! File names - available as a command line option
  CHARACTER(len=110) :: file_ini  ! netcdf file holding distributed initial
                                  ! state and covariance matrix
  CHARACTER(len=110) :: file_obs  ! netcdf file holding observations
  CHARACTER(len=110) :: file_obs_mask  ! ASCII file holding observation mask
  CHARACTER(len=110) :: file_syntobs   ! netcdf file holding synthetic observations

!    ! Other variables - _NOT_ available as command line options!
  INTEGER :: covartype     ! For SEIK: Definition of ensemble covar matrix
                           ! (0): Factor (r+1)^-1 (or N^-1)
                           ! (1): Factor r^-1 (or (N-1)^-1) - real ensemble covar.
                           ! This setting is only for the model part; The definition
                           ! of P has also to be specified in PDAF_filter_init.
                           ! Only for upward-compatibility of PDAF!
  REAL    :: time          ! model time
  INTEGER :: obsfile_laststep(1) ! Last time step in observation file
  INTEGER :: delt_obs_file     ! Observation interval in input file
  LOGICAL :: have_obs          ! Flag whether we consider observations
                               ! at next possible analysis time
  INTEGER :: fileid_state      ! Netcdf ID of the file holding true state trajectory
  REAL, ALLOCATABLE :: state_true(:)    ! Array holding true state
  REAL, ALLOCATABLE :: observation_g(:) ! For local filter: global observation vector
  INTEGER, ALLOCATABLE :: obs_mask(:) ! Mask array for observation availability
  INTEGER, ALLOCATABLE :: obsindx(:)  ! Index array for observations
  INTEGER, ALLOCATABLE :: obsindx_l(:) ! Index array for local observations
             ! This array contains the index of local observation in the 
             ! reordered gappy global observation vector
!EOP


END MODULE mod_assimilation
