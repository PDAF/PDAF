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
! Implementation for the 2D online example
! with parallelization.
!
! !REVISION HISTORY:
! 2013-09 - Lars Nerger - Initial code
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
  INTEGER, ALLOCATABLE :: obs_index_l(:)  ! Vector holding local state-vector indices of observations
  REAL, ALLOCATABLE    :: distance_l(:)   ! Vector holding distances of local observations


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
  INTEGER :: ensgroup     ! Type of initial ensemble

! ! Filter settings - available as command line options
!    ! General
  REAL    :: forget       ! Forgetting factor for filter analysis
!    ! LSEIK/LETKF/LESTKF
  REAL    :: local_range   ! Range for local observation domain
  INTEGER :: locweight     ! Type of localizing weighting of observations
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRANGE
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
  REAL    :: srange        ! Support range for 5th order polynomial
                           !   or radius for 1/e for exponential weighting

!    ! File output - available as a command line option
  CHARACTER(len=110) :: filename  ! file name for assimilation output

!    ! Other variables - _NOT_ available as command line options!
  REAL    :: time          ! model time

END MODULE mod_assimilation
