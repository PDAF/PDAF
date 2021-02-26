!$Id$
!>  Module holding variables for assimilation
!!
!! This module provides variables needed for the 
!! assimilation to be shared within the different
!! user-supplied routines for PDAF.
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
MODULE mod_assim_pdaf

  USE mod_assim_atm_pdaf, ONLY: dp

  IMPLICIT NONE
  SAVE

! *** Below are the generic variables used for configuring PDAF ***
! *** Their values are set in init_PDAF                         ***

! *** Observation-specific variables are kepts in the different ***
! *** OMI observation modules                                   ***

! Settings for time stepping - available as namelist read-in
  INTEGER :: step_null = 0       ! initial time step of assimilation
  LOGICAL :: restart = .false.   ! Whether to restart the data assimilation from previous run

! General settings for observations - available as namelist read-in
  INTEGER :: use_global_obs=1 ! Whether to use global full obs, of full obs limited to process domains

! General control of PDAF - available as namelist read-in
  INTEGER :: screen       ! Control verbosity of PDAF
                          ! (0) no outputs, (1) progess info, (2) add timings
                          ! (3) debugging output
  INTEGER :: dim_ens      ! Size of ensemble
  INTEGER :: filtertype   ! Select filter algorithm:
                          !   SEIK (1), EnKF (2), LSEIK (3), ETKF (4), LETKF (5), 
                          !   ESTKF (6), LESTKF (7), LEnKF (8), NETF (9), LNETF (10), PF (12)
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
  INTEGER :: incremental  ! Perform incremental updating
  INTEGER :: dim_lag      ! Number of time instances for smoother
  INTEGER :: DA_couple_type=0 ! (0) for weakly-coupled, (1) for strongly-coupled assimilation

! *** Filter settings - available as namelist read-in
! General
  INTEGER :: type_forget  ! Type of forgetting factor
  REAL(dp) :: forget      ! Forgetting factor for filter analysis
  INTEGER :: dim_bias     ! dimension of bias vector
  REAL(dp) :: varscale=1.0 ! Scaling factor for initial ensemble variance
! SEIK/ETKF/LSEIK/LETKF
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
! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
  INTEGER :: type_sqrt     ! Type of the transform matrix square-root 
                           !   (0) symmetric square root, (1) Cholesky decomposition
! NETF/LNETF
  INTEGER :: type_winf     ! Set weights inflation: (1) activate
  REAL    :: limit_winf    ! Limit for weights inflation: N_eff/N>limit_winf
! Localization - LSEIK/LETKF/LESTKF/LNETF
  INTEGER :: locweight     ! Type of localizing weighting of observations
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRADIUS
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance

! File output and input - available as as namelist read-in
  LOGICAL :: read_inistate = .false.            ! Whether to read initial state from separate file
  CHARACTER(len=100) :: path_init = '.'         ! Path to initialization files
  CHARACTER(len=110) :: file_init = 'covar_'    ! netcdf file holding distributed initial
                                                ! state and covariance matrix (added is _XX.nc)
  CHARACTER(len=110) :: file_inistate = 'state_ini_' ! netcdf file holding distributed initial
                                                ! state (added is _XX.nc)

! Settings for state
  INTEGER :: dim_state              ! Global size of model state
  INTEGER :: dim_state_p            ! PE-local size of model state
  INTEGER :: istep_asml             ! Time step at end of an forecast phase
  LOGICAL :: flag_final=.false.     ! Whether the current is the final analysis step

! Variables for adaptive localization radius
  REAL(dp), ALLOCATABLE :: eff_dim_obs(:)        ! Effective observation dimension
  INTEGER :: loctype       ! Type of localization
                           !   (0) Fixed radius defined by lradius
                           !   (1) Variable radius for constant effective observation dimension
  REAL(dp) :: loc_ratio    ! Choose lradius so the effective observation dim. is loc_ratio times dim_ens

! Other variables - _NOT_ available in the namelist
  REAL(dp)  :: time                       ! model time
  INTEGER :: n_fields                     ! Number of model fields in state vector
  INTEGER, ALLOCATABLE :: dim_fields_p(:) ! Dimension of fields in process-local state vector
  INTEGER, ALLOCATABLE :: dim_fields(:)   ! Dimension of fields in global state vector
  INTEGER, ALLOCATABLE :: off_fields_p(:) ! Offsets of fields in state vector
  REAL :: coords_l(2)                     ! Coordinates of local analysis domain
  INTEGER, ALLOCATABLE :: id_lstate_in_pstate(:) ! Indices of local state vector in PE-local global state vector
  REAL, PARAMETER :: pi=3.141592653589793

END MODULE mod_assim_pdaf
