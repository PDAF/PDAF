!$Id$
!> PDAF-OMI template observation module 
!!
!! This module handles operations for one data type (called 'module-type' below).
!! 
!! Module type here: simulated airt temperature observations (taken from model)
!!
!! The subroutines in this module are for the particular handling of
!! a single observation type.
!! The routines are called by the different call-back routines of PDAF.
!! Most of the routines are generic so that in practice only 2 routines
!! need to be adapted for a particular data type. These are the routines
!! for the initialization of the observation information (init_dim_obs)
!! and for the observation operator (obs_op).
!!
!! The module and the routines are named according to the observation type.
!! This allows to distinguish the observation type and the routines in this
!! module from other observation types.
!!
!! The module uses two derived data type (obs_f and obs_l), which contain
!! all information about the full and local observations. Only variables
!! of the type obs_f need to be initialized in this module. The variables
!! in the type obs_l are initilized by the generic routines from PDAFomi.
!!
!!
!! These 3 routines need to be adapted for the particular observation type:
!! * init_dim_obs_TYPE \n
!!           Count number of process-local and full observations; 
!!           initialize vector of observations and their inverse variances;
!!           initialize coordinate array and index array for indices of
!!           observed elements of the state vector.
!! * obs_op_TYPE \n
!!           observation operator to get full observation vector of this type. Here
!!           one has to choose a proper observation operator or implement one.
!! * deallocate_obs_TYPE \n
!!           Deallocate observation type
!!
!! In addition, there are two optional routine, which are required if filters 
!! with localization are used:
!! * init_dim_obs_l_TYPE \n
!!           Only required if domain-localized filters (e.g. LESTKF, LETKF) are used:
!!           Count number of local observations of module-type according to
!!           their coordinates (distance from local analysis domain). Initialize
!!           module-internal distances and index arrays.
!! * localize_covar_TYPE \n
!!           Only required if the localized EnKF is used:
!!           Apply covariance localization in the LEnKF.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE obs_airt_pdafomi

  USE mod_parallel_pdaf, &
       ONLY: mype_filter_echam    ! Rank of filter process
  USE PDAFomi, &
       ONLY: obs_f, obs_l         ! Declaration of observation data types
 
  IMPLICIT NONE
  SAVE

  ! Variables which are inputs to the module (usually set in init_pdaf)
  LOGICAL :: assim_a_airt      !< Whether to assimilate this data type
  REAL    :: rms_obs_airt      !< Observation error standard deviation (for constant errors)

  REAL    :: lradius_airt      ! Local observation radius for airT
  REAL    :: sradius_airt      ! Support radius for 5th order polynomial or radius for 1/e for exponential weighting

  CHARACTER(len=100) :: path_obs_airt  = ''      ! Path to observation files
  CHARACTER(len=110) :: file_airt_prefix  = ''   ! file name prefix for observations 
  CHARACTER(len=110) :: file_airt_suffix  = '.nc'! file name suffix for observations 
  CHARACTER(len=110) :: file_syntobs_airt = 'syntobs_airt.nc' ! File name for synthetic observations


! ***********************************************************************
! *** The following two data types are used in PDAFomi                ***
! *** They are declared in PDAFomi and only listed here for reference ***
! ***********************************************************************

! Data type to define the full observations by internally shared variables of the module
!   TYPE obs_f
!           Mandatory variables to be set in INIT_DIM_OBS
!      INTEGER :: doassim                   ! Whether to assimilate this observation type
!      INTEGER :: disttype                  ! Type of distance computation to use for localization
!                                           ! (0) Cartesian, (1) Cartesian periodic
!                                           ! (2) simplified geographic, (3) geographic haversine function
!      INTEGER :: ncoord                    ! Number of coordinates use for distance computation
!      INTEGER, ALLOCATABLE :: id_obs_p(:,:) ! Indices of observed field in state vector (process-local)
!           
!           Optional variables - they can be set in INIT_DIM_OBS
!      REAL, ALLOCATABLE :: icoeff_p(:,:)   ! Interpolation coefficients for obs. operator
!      REAL, ALLOCATABLE :: domainsize(:)   ! Size of domain for periodicity (<=0 for no periodicity)
!
!           Variables with predefined values - they can be changed in INIT_DIM_OBS
!      INTEGER :: obs_err_type=0            ! Type of observation error: (0) Gauss, (1) Laplace
!      INTEGER :: use_global_obs=1          ! Whether to use (1) global full obs. 
!                                           ! or (0) obs. restricted to those relevant for a process domain
!
!           The following variables are set in the routine PDAFomi_gather_obs
!      INTEGER :: dim_obs_p                 ! number of PE-local observations
!      INTEGER :: dim_obs_f                 ! number of full observations
!      INTEGER :: dim_obs_g                 ! global number of observations
!      INTEGER :: off_obs_f                 ! Offset of this observation in overall full obs. vector
!      INTEGER :: off_obs_g                 ! Offset of this observation in overall global obs. vector
!      INTEGER :: obsid                     ! Index of observation over all assimilated observations
!      REAL, ALLOCATABLE :: obs_f(:)        ! Full observed field
!      REAL, ALLOCATABLE :: ocoord_f(:,:)   ! Coordinates of full observation vector
!      REAL, ALLOCATABLE :: ivar_obs_f(:)   ! Inverse variance of full observations
!      INTEGER, ALLOCATABLE :: id_obs_f_lim(:) ! Indices of domain-relevant full obs. in global vector of obs.
!                                           ! (only if full obs. are restricted to process domain))
!   END TYPE obs_f

! Data type to define the local observations by internally shared variables of the module
!   TYPE obs_l
!      INTEGER :: dim_obs_l                 ! number of local observations
!      INTEGER :: off_obs_l                 ! Offset of this observation in overall local obs. vector
!      INTEGER, ALLOCATABLE :: id_obs_l(:)  ! Indices of local observations in full obs. vector 
!      REAL, ALLOCATABLE :: distance_l(:)   ! Distances of local observations
!      REAL, ALLOCATABLE :: ivar_obs_l(:)   ! Inverse variance of local observations
!      INTEGER :: locweight                 ! Specify localization function
!      REAL :: lradius                      ! localization radius
!      REAL :: sradius                      ! support radius for localization function
!   END TYPE obs_l
! ***********************************************************************

! Declare instances of observation data types used here
! We use generic names here, but one could renamed the variables
  TYPE(obs_f), TARGET, PUBLIC :: thisobs      ! full observation
  TYPE(obs_l), TARGET, PUBLIC :: thisobs_l    ! local observation

!$OMP THREADPRIVATE(thisobs_l)


!-------------------------------------------------------------------------------

CONTAINS

!> Initialize information on the module-type observation
!!
!! The routine is called by each filter process.
!! at the beginning of the analysis step before 
!! the loop through all local analysis domains.
!! 
!! It has to count the number of observations of the
!! observation type handled in this module according
!! to the current time step for all observations 
!! required for the analyses in the loop over all local 
!! analysis domains on the PE-local state domain.
!!
!! The following four variables have to be initialized in this routine
!! * thisobs\%doassim     - Whether to assimilate this type of observations
!! * thisobs\%disttype    - type of distance computation for localization with this observaton
!! * thisobs\%ncoord      - number of coordinates used for distance computation
!! * thisobs\%id_obs_p    - index of module-type observation in PE-local state vector
!!
!! Optional is the use of
!! * thisobs\%icoeff_p    - Interpolation coefficients for obs. operator (only if interpolation is used)
!! * thisobs\%domainsize  - Size of domain for periodicity for disttype=1 (<0 for no periodicity)
!! * thisobs\%obs_err_type - Type of observation errors for particle filter and NETF (default: 0=Gaussian)
!! * thisobs\%use_global obs - Whether to use global observations or restrict the observations to the relevant ones
!!                          (default: 1=use global full observations)
!!
!! Further variables are set when the routine PDAFomi_gather_obs is called.
!!
  SUBROUTINE init_dim_obs_airt(step, dim_obs)

    USE PDAFomi, ONLY: PDAFomi_gather_obs
    USE mod_assim_pdaf, ONLY: filtertype, pi
    USE mo_geoloc, ONLY: philat_2d, philon_2d
    USE mo_memory_g1a,    ONLY: tm1
    USE mo_decomposition, ONLY: dc=>local_decomposition

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)    :: step       !< Current time step
    INTEGER, INTENT(inout) :: dim_obs    !< Dimension of full observation vector

! *** Local variables ***
    INTEGER :: i, j                      ! Counters
    INTEGER :: k, jk, jrow, jl           ! Counters
    INTEGER :: nproma 
    INTEGER :: dim_obs_p                 ! Number of process-local observations
    REAL, ALLOCATABLE :: obs_p(:)        ! Process-local observation vector
    REAL, ALLOCATABLE :: ivar_obs_p(:)   ! Process-local inverse observation error variance
    REAL, ALLOCATABLE :: ocoord_p(:,:)   ! Process-local observation coordinates


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    IF (mype_filter_echam==0) &
         WRITE (*,'(8x,a)') 'Assimilate observations - OBS_airt'

    ! Store whether to assimilate this observation type (used in routines below)
    IF (assim_a_airt) thisobs%doassim = 1

    ! Specify type of distance computation
    thisobs%disttype = 2   ! 2=Geographic

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2


! **********************************
! *** Read PE-local observations ***
! **********************************

    ! For testing purpose, generate the virtual observations using the bottom
    ! layer values with perturbations
    ! thus we don't read anything here
    

! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

    ! *** Count valid observations that lie within the process sub-domain ***

    dim_obs_p = dc%nglat * dc%nglon

    ! *** Initialize vector of observations on the process sub-domain ***

    ALLOCATE(obs_p(dim_obs_p))

    k = 1
    jk = dc%nlev
    DO jrow = 1, dc%ngpblks
       IF ( jrow == dc%ngpblks ) THEN
          nproma = dc%npromz
       ELSE
          nproma = dc%nproma
       END IF
       DO jl = 1, nproma
          obs_p(k) = tm1(jl,jk,jrow)
          k = k + 1
       END DO
    END DO

    obs_p = obs_p + 1.0


    ! *** Initialize coordinate array of observations on the process sub-domain ***

    ALLOCATE(ocoord_p(2, dim_obs_p))

    k = 1
    DO jrow = 1, dc%ngpblks

       IF ( jrow == dc%ngpblks ) THEN
          nproma = dc% npromz
       ELSE
          nproma = dc% nproma
       END IF

       DO jl = 1, nproma
          ocoord_p(1, k) = philon_2d(jl,jrow) * pi / 180.0
          ocoord_p(2, k) = philat_2d(jl,jrow) * pi / 180.0
          k = k + 1

       END DO
    END DO


    ! *** Initialize process local index array                         ***
    ! *** This array holds the information which elements of the state ***
    ! *** vector are used in the observation optation.                 ***
    ! *** It has a many rows as required for the observation operator, ***
    ! *** i.e. 1 if observations are at grid points; >1 if             ***
    ! *** interpolation is required                                    ***

    ALLOCATE(thisobs%id_obs_p(1, dim_obs_p))

    DO i = 1, dc%nglat * dc%nglon
       thisobs%id_obs_p(1, i)=i
    END DO


! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

    ALLOCATE(ivar_obs_p(dim_obs_p))
    
    ! Set inverse observation error variance
    DO i = 1, dim_obs_p
       ivar_obs_p(i) = 1.0 / (rms_obs_airt*rms_obs_airt)
    END DO


! ****************************************
! *** Gather global observation arrays ***
! ****************************************

    CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
         thisobs%ncoord, lradius_airt, dim_obs)


! ********************
! *** Finishing up ***
! ********************

    ! Deallocate all local arrays
    DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

    ! Arrays in THISOBS have to be deallocated after the analysis step
    ! by a call to deallocate_obs() in prepoststep_pdaf.

  END SUBROUTINE init_dim_obs_airt



!-------------------------------------------------------------------------------
!> Implementation of observation operator 
!!
!! This routine applies the full observation operator
!! for the type of observations handled in this module.
!!
!! One can choose a proper observation operator from
!! PDAFOMI_OBS_OP or add one to that module or 
!! implement another observation operator here.
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE obs_op_airt(dim_p, dim_obs, state_p, ostate)

    USE PDAFomi, &
         ONLY: PDAFomi_obs_op_gridpoint

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
    REAL, INTENT(in)    :: state_p(dim_p)        !< PE-local model state
    REAL, INTENT(inout) :: ostate(dim_obs)       !< Full observed state


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

    IF (thisobs%doassim==1) THEN
       CALL PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)
    END IF

  END SUBROUTINE obs_op_airt



!-------------------------------------------------------------------------------
!> Initialize local information on the module-type observation
!!
!! The routine is called during the loop over all local
!! analysis domains. It has to initialize the information
!! about local observations of the module type. It returns
!! number of local observations of the module type for the
!! current local analysis domain in DIM_OBS_L and the full
!! and local offsets of the observation in the overall
!! observation vector.
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type and  local analysis domain.
!!
  SUBROUTINE init_dim_obs_l_airt(domain_p, step, dim_obs, dim_obs_l)

    ! Include PDAFomi function
    USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l

    ! Include localization radius and local coordinates
    USE mod_assim_pdaf, ONLY: coords_l, locweight

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)  :: domain_p     !< Index of current local analysis domain
    INTEGER, INTENT(in)  :: step         !< Current time step
    INTEGER, INTENT(in)  :: dim_obs      !< Full dimension of observation vector
    INTEGER, INTENT(inout) :: dim_obs_l  !< Local dimension of observation vector


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    IF (thisobs%doassim==1) THEN
       CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
            locweight, lradius_airt, lradius_airt, dim_obs_l)
    END IF

  END SUBROUTINE init_dim_obs_l_airt



!-------------------------------------------------------------------------------
!> Deallocate observation type
!!
!! This routine calls the routine PDAFomi_deallocate_obs
!!
  SUBROUTINE deallocate_obs_airt()

    ! Include PDAFomi function
    USE PDAFomi, ONLY: PDAFomi_deallocate_obs


! *************************************
! *** Deallocate observation arrays ***
! *************************************

    CALL PDAFomi_deallocate_obs(thisobs)

  END SUBROUTINE deallocate_obs_airt

END MODULE obs_airt_pdafomi
