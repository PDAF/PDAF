C     $Id$
C!BOP
C
C     !DESCRIPTION:
C     This module provides variables needed for the
C     assimilation within the routines of the dummy model.
C     For simplicity, all assimilation-related variables
C     are stored here, even if they are only used in
C     the main program for the filter initialization.
C     Most variables can be specified as a command line
C     argument.
C
C     Implementation for MITgcm without parallelization.
C
C     !REVISION HISTORY:
C     2018-01 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !MODULE:
      MODULE mod_assimilation

C     !USES:
      IMPLICIT NONE
      SAVE
C!EOP

C     *** Variables specific for MITgcm ***

      INTEGER :: myIter       ! Iteration counter
      INTEGER :: myThid       ! Thread ID
      REAL(kind=8) :: myTime  ! Time

C     *** Variables specific for state and observation vectors ***

      INTEGER :: dim_state                   ! Global model state dimension
      INTEGER :: dim_state_p                 ! Model state dimension for PE-local domain
      INTEGER :: iter_null                   ! Initial time step
      REAL(kind=8) :: time                   ! model time
      INTEGER :: n_fields                    ! number of fields in state vector
      INTEGER, ALLOCATABLE :: off_fields(:)  ! Offsets of fields in state vector
      INTEGER, ALLOCATABLE :: dim_fields(:)  ! Dimension of fields in state vector
      REAL(kind=8) :: coords_l(2)            ! Coordinates of analysis domain
      INTEGER, ALLOCATABLE :: id_lstate_in_pstate(:) ! Indices of local state vector in PE-local global state vector

C     Declare Fortran type holding the indices of model fields in the state vector
C     This can be extended to any number of fields - it severs to give each field a name
      TYPE field_ids
        INTEGER :: etan 
        INTEGER :: theta
        INTEGER :: salt
        INTEGER :: uvel
        INTEGER :: vvel
      END TYPE field_ids

C     Type variable holding field IDs in state vector
      TYPE(field_ids) :: id


C     *** Below are the generic variables used for configuring PDAF ***
C     *** Their values are set in init_PDAF                         ***

C     !PUBLIC MEMBER FUNCTIONS:
C     Settings for time stepping - available as command line options
      LOGICAL :: model_error         ! Control application of model error
      REAL(kind=8)  :: model_err_amp ! Amplitude for model error

C     Settings for observations - available as command line options
      INTEGER :: nsteps_between_obs  ! time step interval between assimilation steps

C     General control of PDAF - available as command line options
      INTEGER :: screen       ! Control verbosity of PDAF
                          ! (0) no outputs, (1) progess info, (2) add timings
                          ! (3) debugging output
      INTEGER :: dim_ens      ! Size of ensemble for SEIK/LSEIK/EnKF/ETKF
                          ! Number of EOFs to be used for SEEK
      INTEGER :: filtertype   ! Select filter algorithm:
                          !   SEIK (1), EnKF (2), LSEIK (3), ETKF (4), LETKF (5)
                          !   ESTKF (6), LESTKF (7), LEnKF (8), NETF (9), LNETF (10), PF (12)
      INTEGER :: subtype      ! Subtype of filter algorithm
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
                          !   PF:
                          !     (0) standard PF with resampling
      INTEGER :: incremental  ! Perform incremental updating in LSEIK
      INTEGER :: dim_lag      ! Number of time instances for smoother

C     Filter settings - available as command line options
C     General
      INTEGER :: type_forget  ! Type of forgetting factor
                          !   (0) fixed
                          !   (1) global adaptive
                          !   (2) local adaptive for LSEIK/LETKF/LESTKF
      REAL(kind=8) :: forget  ! Forgetting factor for filter analysis
      INTEGER :: dim_bias     ! dimension of bias vector
C     ENKF
      INTEGER :: rank_analysis_enkf ! Rank to be considered for inversion of HPH
C     SEIK/ETKF/ESTKF/LSEIK/LETKF/LESTKF
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
C     LSEIK/LETKF/LESTKF/LNETF
      REAL(kind=8)  :: cradius  ! Cut-off radius for local observation domain
      INTEGER :: locweight      ! Type of localizing weighting of observations
                    ! For LESTKF, LETKF, LSEIK, LNETF
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRADIUS
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
                    ! For LEnKF
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRADIUS
                    !   (2) 5th-order polynomial weight function
      REAL(kind=8)  :: sradius  ! Support radius for 5th order polynomial
                                ! or radius for 1/e for exponential weighting
C     SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
      INTEGER :: type_sqrt     ! Type of the transform matrix square-root
                    !   (0) symmetric square root, (1) Cholesky decomposition

C     File output - available as a command line option
      LOGICAL :: write_ens   ! Whether to write full ensemble into files

      END MODULE mod_assimilation
