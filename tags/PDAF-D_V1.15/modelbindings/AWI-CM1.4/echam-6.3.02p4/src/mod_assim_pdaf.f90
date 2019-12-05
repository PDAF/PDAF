!$Id: mod_assim_pdaf.f90 2135 2019-11-22 18:56:29Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_assim_pdaf

! !DESCRIPTION:
! This module provides variables needed for the 
! assimilation within the routines of the dummy model.
! For simplicity, all assimilation-related variables
! are stored here, even if they are only used in
! the main program for the filter initialization.
! Most variables can be specified as a command line 
! argument.
!
! Implementation for the 2D online example
! with parallelization.
!
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mo_kind_pdaf, ONLY: dp

  IMPLICIT NONE
  SAVE
!EOP


! *** Below are the generic variables used for configuring PDAF ***
! *** Their values are set in init_PDAF                         ***

! !PUBLIC MEMBER FUNCTIONS:

! ! Settings for time stepping - available as namelist read-in
  INTEGER :: step_null = 0 ! initial time step of assimilation

! ! Settings for observations - available as command line options
  INTEGER :: delt_obs_ocn  ! time step interval between assimilation steps - Ocean
  INTEGER :: delt_obs_atm  ! time step interval between assimilation steps - Atmosphere
  REAL(dp):: rms_obs       ! RMS error size for observation generation  --- SST
  REAL(dp):: rms_obs_T     ! RMS error size for observation generation  --- T
  REAL(dp):: rms_obs_S     ! RMS error size for observation generation  --- S
  INTEGER :: dim_obs       ! Number of observations
  REAL(dp):: bias_obs      ! Assumed observation bias
  REAL    :: peak_obs_error  ! Peak value used to define the observation error
  INTEGER :: obs_err_type  ! Observation errors: (0) for Gaussian (1) for double-exponential
  INTEGER :: proffiles_o   ! (0) don't generate  them; 
                           ! (1) generate distributed profile files
                           ! (2) generate global profile file
  LOGICAL :: assim_o_sst   ! Whether to assimilate SST data
  LOGICAL :: assim_o_en4_t ! Whether to assimilate temperature profiles
  LOGICAL :: assim_o_en4_s ! Whether to assimilate salinity profiles
  INTEGER :: writeprofile  ! (0) No profile is written (1) Write Temperature profile 
                           ! (2) Write Salinity profile (3) Write both Temperature and salinity profiles
  LOGICAL :: sst_exclude_ice ! Whether to exclude SST observations at grid points with ice
  REAL    :: sst_exclude_diff ! Limit difference beyond which observations are excluded (0.0 to deactivate)
  REAL    :: prof_exclude_diff ! Limit difference beyond which observations are excluded (0.0 to deactivate)
  LOGICAL :: twin_experiment = .false.   ! Whether to perform a twin experiment with synthetic observations
  INTEGER :: dim_obs_max   ! Expect max. number of observations for synthetic obs.

! ! General control of PDAF - available as command line options
  INTEGER :: screen       ! Control verbosity of PDAF
                          ! (0) no outputs, (1) progess info, (2) add timings
                          ! (3) debugging output
  INTEGER :: dim_ens      ! Size of ensemble for SEIK/LSEIK/EnKF/ETKF
                          ! Number of EOFs to be used for SEEK
  INTEGER :: filtertype   ! Select filter algorithm:
                          ! SEEK (0), SEIK (1), EnKF (2), LSEIK (3), ETKF (4), LETKF (5)
  INTEGER :: subtype      ! Subtype of filter algorithm
                          !   SEIK:
                          !     (0) ensemble forecast; new formulation
                          !     (1) ensemble forecast; old formulation
                          !     (2) fixed error space basis
                          !     (3) fixed state covariance matrix
                          !     (4) SEIK with ensemble transformation
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
  INTEGER :: incremental  ! Perform incremental updating in LSEIK
  INTEGER :: dim_lag      ! Number of time instances for smoother
  INTEGER :: DA_couple_type=0 ! (0) for weakly-coupled, (1) for strongly-coupled assimilation

! ! Filter settings - available as command line options
!    ! General
  INTEGER :: type_forget  ! Type of forgetting factor
  REAL(dp) :: forget      ! Forgetting factor for filter analysis
  INTEGER :: dim_bias     ! dimension of bias vector
  REAL(dp) :: varscale=1.0 ! Scaling factor for initial ensemble variance
!    ! SEIK/ETKF/LSEIK/ETKFS
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

!    ! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
  INTEGER :: type_sqrt     ! Type of the transform matrix square-root 
                           !   (0) symmetric square root, (1) Cholesky decomposition
!    ! Localization - LSEIK/LETKF/LESTKF
  REAL(dp) :: local_range  ! Range for local observation domain
  INTEGER :: locweight     ! Type of localizing weighting of observations
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRANGE
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
  REAL(dp) :: srange        ! Support range for 5th order polynomial
                           !   or radius for 1/e for exponential weighting

!    ! Specific for ECHAM
  INTEGER :: dim_state              ! Global size of model state
  INTEGER :: dim_state_p            ! PE-local size of model state
  INTEGER :: istep_asml             ! Time step at end of an forecat phase
  LOGICAL :: flag_final=.false.     ! Whether the current is the final analysis step
  REAL(dp), ALLOCATABLE :: all_ssh_p(:)          ! Complete surface SSH field
  REAL(dp), ALLOCATABLE :: obs_ssh(:)            ! Global observed surface SSH field
  INTEGER, ALLOCATABLE :: local_dimobs(:)        ! PE-local observation dimensions
  INTEGER, ALLOCATABLE :: id_observed_ssh(:)     ! IDs of observed surface nodes

!    ! Specific for local filters
  INTEGER, ALLOCATABLE :: index_local_domain(:)  ! Node indices for local state vector
  INTEGER, ALLOCATABLE :: local_obs_nod2d(:)     ! Global node indices of local obs. domain
  REAL(dp), ALLOCATABLE :: ocoord_n2d(:,:)       ! Coordinates of full observation vector
!    ! Variables for adaptive localization radius
  REAL(dp), ALLOCATABLE :: eff_dim_obs(:)        ! Effective observation dimension
  REAL(dp), ALLOCATABLE :: loc_radius(:)         ! Effective observation dimension
  INTEGER :: loctype       ! Type of localization
                    !   (0) Fixed radius defined by local_range
                    !   (1) Variable radius for constant effective observation dimension
  REAL(dp) :: loc_ratio    ! Choose local_range so the effective observation dim. is loc_ratio times dim_ens


!    ! File output and input - available as as namelist read-in
  LOGICAL :: read_inistate = .false.            ! Whether to read initial state from separate file
  CHARACTER(len=100) :: path_obs_sst  = '.'     ! Path to SST observations
  CHARACTER(len=110) :: file_sst_prefix  = ''   ! file name prefix for SST observations 
  CHARACTER(len=110) :: file_sst_suffix  = '.nc'! file name suffix for SST observations 
  CHARACTER(len=100) :: path_obs_prof  = '.'     ! Path to profile observations
  CHARACTER(len=110) :: file_prof_prefix  = ''   ! file name prefix for profile observations 
  CHARACTER(len=110) :: file_prof_suffix  = '.nc'! file name suffix for profile observations   
  CHARACTER(len=100) :: path_init = '.'      ! Path to initialization files
  CHARACTER(len=110) :: file_init = 'covar_' ! netcdf file holding distributed initial
                                    ! state and covariance matrix (added is _XX.nc)
  CHARACTER(len=110) :: file_inistate = 'state_ini_' ! netcdf file holding distributed initial
                                    ! state (added is _XX.nc)
  CHARACTER(len=110) :: file_syntobs = 'syntobs.nc' ! File name for synthetic observations

!    ! Other variables - _NOT_ available as command line options!
  REAL(dp)  :: time                 ! model time
  INTEGER, ALLOCATABLE :: offset(:) ! Offsets of fields in state vector

END MODULE mod_assim_pdaf
