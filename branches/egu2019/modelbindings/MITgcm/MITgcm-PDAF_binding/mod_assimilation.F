C     $Id$
C!BOP
!
C     !MODULE:
      MODULE mod_assimilation

C     !DESCRIPTION:
C     This module provides variables needed for the 
C     assimilation within the routines of the dummy model.
C     For simplicity, all assimilation-related variables
C     are stored here, even if they are only used in
C     the main program for the filter initialization.
C     Most variables can be specified as a command line 
C     argument.
C
C     Implementation for the 2D online example
C     with parallelization.
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !USES:
      IMPLICIT NONE
      SAVE
C!EOP

C     *** Variables specific for MITgcm ***

      INTEGER :: myIter
      INTEGER :: myThid
      REAL*8  :: myTime

C     *** Variables specific for state and observation vectors ***

      INTEGER :: dim_state           ! Global model state dimension
      INTEGER :: dim_state_p         ! Model state dimension for PE-local domain
      INTEGER, ALLOCATABLE :: off_fields(:)  ! Offsets of fields in state vector
      INTEGER, ALLOCATABLE :: dim_fields(:)  ! Field dimensions in state vector
      REAL*8, ALLOCATABLE  :: obs(:)            ! Vector holding all observations
      INTEGER, ALLOCATABLE :: obs_index(:)      ! Vector holding state-vector indices of observations
      INTEGER, ALLOCATABLE :: obs_index_l(:)    ! Vector holding local state-vector indices of observations
      REAL*8, ALLOCATABLE :: obs_coords(:,:)    ! Array for observation coordinates
      REAL*8, ALLOCATABLE :: distance_l(:)      ! Vector of distances of local observations
      INTEGER :: domain_index(2)                ! i,j indices of local analysis domain
      REAL*8 :: coords_l(2)                     ! Coordinates of analysis domain
      INTEGER, ALLOCATABLE :: local_dims_obs(:) ! Array for process-local observation dimensions


C     *** Below are the generic variables used for configuring PDAF ***
C     *** Their values are set in init_PDAF                         ***

C     !PUBLIC MEMBER FUNCTIONS:
C     Settings for time stepping - available as command line options
      LOGICAL :: model_error   ! Control application of model error
      REAL*8  :: model_err_amp ! Amplitude for model error

C     Settings for observations - available as command line options
      INTEGER :: nsteps_between_obs ! time step interval between assimilation steps
      REAL*8  :: rms_obs       ! RMS error size for observation generation
      INTEGER :: dim_obs       ! Number of observations

C     General control of PDAF - available as command line options
      INTEGER :: screen       ! Control verbosity of PDAF
                          ! (0) no outputs, (1) progess info, (2) add timings
                          ! (3) debugging output
      INTEGER :: dim_ens      ! Size of ensemble for SEIK/LSEIK/EnKF/ETKF
                          ! Number of EOFs to be used for SEEK
      INTEGER :: filtertype   ! Select filter algorithm:
                          ! SEEK (0), SEIK (1), EnKF (2), LSEIK (3), ETKF (4), LETKF (5)
      INTEGER :: subtype      ! Subtype of filter algorithm
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
                          !     (0) ESTKF using T-matrix 
                          !   LESTKF:
                          !     (0) LESTKF using T-matrix 
      INTEGER :: incremental  ! Perform incremental updating in LSEIK
      INTEGER :: dim_lag      ! Number of time instances for smoother

C     Filter settings - available as command line options
C     General
      INTEGER :: type_forget  ! Type of forgetting factor
      REAL*8  :: forget       ! Forgetting factor for filter analysis
      INTEGER :: dim_bias     ! dimension of bias vector
C     SEIK/ETKF/LSEIK/ETKFS
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
C     LSEIK/LETKF
      REAL*8  :: local_range   ! Range for local observation domain
      INTEGER :: locweight     ! Type of localizing weighting of observations
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRANGE
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
      REAL*8  :: srange        ! Support range for 5th order polynomial
                               ! or radius for 1/e for exponential weighting
C     SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
      INTEGER :: type_sqrt     ! Type of the transform matrix square-root 
                    !   (0) symmetric square root, (1) Cholesky decomposition

C     File output - available as a command line option
      LOGICAL :: write_ens   ! Whether to write full ensemble into files

C     Other variables - _NOT_ available as command line options!
      INTEGER :: iter_null   ! Initial time step
      REAL*8  :: time        ! model time
      INTEGER :: nfields     ! Number of modelfields in state vector

      END MODULE mod_assimilation
