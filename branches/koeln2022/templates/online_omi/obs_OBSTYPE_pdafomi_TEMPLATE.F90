!$Id$
!> PDAF-OMI template observation module 
!!
!! This module handles operations for one data type (called 'module-type' below).
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
!! The module uses two derived data types (obs_f and obs_l), which contain
!! all information about the full and local observations. Only variables
!! of the type obs_f need to be initialized in this module. The variables
!! in the type obs_l are initilized by the generic routines from PDAFomi.
!!
!!
!! **Using this template:**
!!   To be able to distinguish the observation type and the routines in this module,
!!   we recommend to rename the module according to the observation module-type.
!!   Further,we recommend to replace 'TYPE' in the routine names according to the
!!   type of the observation so that they can be identified when calling them from 
!!   the call-back routines.
!!
!!
!! These 2 routines need to be adapted for the particular observation type:
!! * init_dim_obs_OBSTYPE \n
!!           Count number of process-local and full observations; 
!!           initialize vector of observations and their inverse variances;
!!           initialize coordinate array and index array for indices of
!!           observed elements of the state vector.
!! * obs_op_OBSTYPE \n
!!           observation operator to get full observation vector of this type. Here
!!           one has to choose a proper observation operator or implement one.
!!
!! In addition, there are two optional routine, which are required if filters 
!! with localization are used:
!! * init_dim_obs_l_OBSTYPE \n
!!           Only required if domain-localized filters (e.g. LESTKF, LETKF) are used:
!!           Count number of local observations of module-type according to
!!           their coordinates (distance from local analysis domain). Initialize
!!           module-internal distances and index arrays.
!! * localize_covar_OBSTYPE \n
!!           Only required if the localized EnKF is used:
!!           Apply covariance localization in the LEnKF.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE obs_OBSTYPE_pdafomi

  USE mod_parallel_pdaf, &
       ONLY: mype_filter    ! Rank of filter process
  USE PDAFomi, &
       ONLY: obs_f, obs_l   ! Declaration of observation data types
 
  IMPLICIT NONE
  SAVE

  ! Variables which are inputs to the module (usually set in init_pdaf)
  LOGICAL :: assim_OBSTYPE        !< Whether to assimilate this data type
  REAL    :: rms_obs_OBSTYPE      !< Observation error standard deviation (for constant errors)

  ! One can declare further variables, e.g. for file names which can
  ! be use-included in init_pdaf() and initialized there.


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
!! **Adapting the template**
!! In this routine the variables listed above have to be initialized. One
!! can include modules from the model with 'use', e.g. for mesh information.
!! Alternatively one could include these as subroutine arguments
!!
  SUBROUTINE init_dim_obs_OBSTYPE(step, dim_obs)

    USE PDAFomi, &
         ONLY: PDAFomi_gather_obs
    USE mod_assimilation, &
         ONLY: filtertype, local_range

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)    :: step       !< Current time step
    INTEGER, INTENT(inout) :: dim_obs    !< Dimension of full observation vector

! *** Local variables ***
    INTEGER :: i, j                      ! Counters
    INTEGER :: dim_obs_p                 ! Number of process-local observations
    REAL, ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
    REAL, ALLOCATABLE :: ivar_obs_p(:)   ! PE-local inverse observation error variance
    REAL, ALLOCATABLE :: ocoord_p(:,:)   ! PE-local observation coordinates 


    ! Template reminder - delete when implementing functionality
    WRITE (*,*) 'TEMPLATE init_OBSTYPE_pdafomi_TEMPLATE.F90: Initialize observations'

! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    IF (mype_filter==0) &
         WRITE (*,'(8x,a)') 'Assimilate observations - OBS_OBSTYPE'

    ! Store whether to assimilate this observation type (used in routines below)
    IF (assim_OBSTYPE) thisobs%doassim = 1

    ! Specify type of distance computation
    thisobs%disttype = 0   ! 0=Cartesian

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2


! +++++ This is a dummy implementation for a single observation
! +++++ Its only purpose is to let the program run with segmentation fault

    IF (mype_filter==0) &
         WRITE (*,'(8x,a)') 'Dummy implementation of a single observation'

    dim_obs_p = 1
    ALLOCATE(obs_p(dim_obs_p))
    obs_p = 1.0
    ALLOCATE(ocoord_p(thisobs%ncoord, dim_obs_p))
    ocoord_p = 1.0
    ALLOCATE(thisobs%id_obs_p(1 , dim_obs_p))
    thisobs%id_obs_p = 1
    ALLOCATE(ivar_obs_p(dim_obs_p))
    ivar_obs_p = 1.0
! +++++ END OF DUMMY OBSERVATION


! **********************************
! *** Read PE-local observations ***
! **********************************

  ! read observation values and their coordinates
  ! also read observation error information if available


! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

    ! *** Count valid observations that lie within the process sub-domain ***

!    dim_obs_p = ...
    

    ! *** Initialize vector of observations on the process sub-domain ***

!    ALLOCATE(obs_p(dim_obs_p))

!    obs_p = ....


    ! *** Initialize coordinate array of observations on the process sub-domain ***

  ! The coordinates are only used in case of the local filters or to compute
  ! interpolation coefficients (see below). In any case, the array must be
  ! allocated because it's an argument of PDAFomi_gather_obs called below.

!    ALLOCATE(ocoord_p(thisobs%ncoord, dim_obs_p))

!    ocoord_p = ....

    ! *** Initialize process local index array                         ***
    ! *** This array holds the information which elements of the state ***
    ! *** vector are used in the observation operator.                 ***
    ! *** It has a many rows as required for the observation operator, ***
    ! *** i.e. 1 if observations are at grid points; >1 if             ***
    ! *** interpolation is required                                    ***

  ! The initialization is done locally for each process sub-domain and later
  ! used in the observation operator. 
  ! Examples:
  ! 1. If the observations are model fields located at grid points, one should
  !   initialize the index array thisobs%id_obs_p with one row so that it contains 
  !   the indices of the observed field values in the process-local state vector
  !   (state_p). Then one can use the observation operator OBS_OP_GRIDPOINT 
  !   provided by the module PDAFomi.
  ! 2. If the observations are the average of model fields located at grid points,
  !   one should initialize the index array thisobs%id_obs_p with as many rows as 
  !   values to be averaged. Each column of the arrays then contains the indices of
  !   the elements of the process-local state vector that have to be averaged. With
  !   this index array one can then use the observation operator OBS_OP_GRIDAVG
  !   provided by the module PDAFomi.
  ! 3. If model values need to be interpolated to the observation location
  !   one should initialize the index array thisobs%id_obs_p with as many rows as 
  !   values are required in the interpolationto be averaged. Each column of the 
  !   array then contains the indices of elements of the process-local state vector 
  !   that are used in the interpolation.
  ! Below, you need to replace NROWS by the number of required rows
  ! Note: This array is only used in the observation operator routines. If you
  !   replace the observation operator routines provided by PDAF-OMI by a custom
  !   observation operator, you might not need this array.

!    ALLOCATE(thisobs%id_obs_p( NROWS , dim_obs_p))

!    thisobs%id_obs_p = ...


! **********************************************************************
! *** Initialize interpolation coefficients for observation operator ***
! **********************************************************************

  ! This initialization is only required if an observation operator
  ! with interpolation is used. The coefficients should be determined
  ! here instead of the observation operator, because the operator is
  ! called for each ensemble member while init_dim_obs is only called
  ! once.

  ! Allocate array of interpolation coefficients. As thisobs%id_obs_p, the number
  ! of rows corresponds to the number of grid points using the the interpolation

!    ALLOCATE(thisobs%icoeff_p( NROWS , dim_obs_p))

  ! Ensure that the order of the coefficients is consistent with the
  ! indexing in thisobs%id_obs_p. Further ensure that the order is consistent
  ! with the assumptions used in the observation operator.

!    thisobs%icoeff_p = ...


! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

!    ALLOCATE(ivar_obs_p(dim_obs_p))
    
!    ivar_obs_p = ...


! ****************************************
! *** Gather global observation arrays ***
! ****************************************

    ! NOTE FOR DIM_OBS_P=0
    ! For the call to PDAFomi_gather_obs_f, obs_p, ivar_obs_p, ocoord_p,
    ! and thisobs%id_obs_p need to be allocated. Thus, if dim_obs_p=0 can 
    ! happen in your application you should explicitly handle this case.
    ! You can introduce an IF block in the initializations above:
    !  IF dim_obs_p>0 THEN
    !     regular allocation and initialization of obs_p, ivar_obs_p, ocoord_p
    !  ELSE
    !     allocate obs_p, ivar_obs_p, ocoord_p, thisobs%id_obs_p with size=1
    !  ENDIF


    ! This routine is generic for the case that only the observations, 
    ! inverse variances and observation coordinates are gathered

    CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
         thisobs%ncoord, local_range, dim_obs)


! *********************************************************
! *** For twin experiment: Read synthetic observations  ***
! *********************************************************

!   IF (twin_experiment .AND. filtertype/=100) THEN
!      CALL read_syn_obs(file_syntobs_OBSTYPE, dim_obs, thisobs%obs_f, 0, 1-mype_filter)
!   END IF


! ********************
! *** Finishing up ***
! ********************

    ! Deallocate all local arrays
!    DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

    ! Arrays in THISOBS have to be deallocated after the analysis step
    ! by a call to deallocate_obs() in prepoststep_pdaf.

  END SUBROUTINE init_dim_obs_OBSTYPE



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
  SUBROUTINE obs_op_OBSTYPE(dim_p, dim_obs, state_p, ostate)

    USE PDAFomi, &
         ONLY: PDAFomi_obs_op_gridpoint

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
    REAL, INTENT(in)    :: state_p(dim_p)        !< PE-local model state
    REAL, INTENT(inout) :: ostate(dim_obs)       !< Full observed state


    ! Template reminder - delete when implementing functionality
    WRITE (*,*) 'TEMPLATE init_OBSTYPE_pdafomi_TEMPLATE.F90: Apply observation operator'

! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

    IF (thisobs%doassim==1) THEN
    
       !+++  Choose suitable observation operator from the
       !+++  module PDAFomi_obs_op or implement your own

       ! Example: Observation operator for observed grid point values
       CALL PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)

    END IF

  END SUBROUTINE obs_op_OBSTYPE



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
  SUBROUTINE init_dim_obs_l_OBSTYPE(domain_p, step, dim_obs, dim_obs_l)

    ! Include PDAFomi function
    USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l

    ! Include localization radius and local coordinates
    ! one can also set observation-specific values for the localization.
    USE mod_assimilation, &   
         ONLY: coords_l, local_range, locweight, srange

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)  :: domain_p     !< Index of current local analysis domain
    INTEGER, INTENT(in)  :: step         !< Current time step
    INTEGER, INTENT(in)  :: dim_obs      !< Full dimension of observation vector
    INTEGER, INTENT(inout) :: dim_obs_l  !< Local dimension of observation vector


    ! Template reminder - delete when implementing functionality
    WRITE (*,*) 'TEMPLATE init_OBSTYPE_pdafomi_TEMPLATE.F90: Initialize local observations'

! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    ! Here one has to specify the coordinates of the local analysis domain
    ! (coords_l) and the localization variables, which can be different for
    ! each observation type and can be made dependent on the index DOMAIN_P.
    ! coords_l should be set in the call-back routine init_dim_l.

    CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
         locweight, local_range, srange, dim_obs_l)

  END SUBROUTINE init_dim_obs_l_OBSTYPE



!-------------------------------------------------------------------------------
!> Perform covariance localization for local EnKF on the module-type observation
!!
!! The routine is called in the analysis step of the localized
!! EnKF. It has to apply localization to the two matrices
!! HP and HPH of the analysis step for the module-type
!! observation.
!!
!! This routine calls the routine PDAFomi_localize_covar
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type.
!!
  SUBROUTINE localize_covar_OBSTYPE(dim_p, dim_obs, HP_p, HPH, coords_p)

    ! Include PDAFomi function
    USE PDAFomi, ONLY: PDAFomi_localize_covar

    ! Include localization radius and local coordinates
    USE mod_assimilation, &   
         ONLY: local_range, locweight, srange

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< Dimension of observation vector
    REAL, INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
    REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH
    REAL, INTENT(in)    :: coords_p(:,:)         !< Coordinates of state vector elements


    ! Template reminder - delete when implementing functionality
    WRITE (*,*) 'TEMPLATE init_OBSTYPE_pdafomi_TEMPLATE.F90: Apply covariance localization'

! *************************************
! *** Apply covariance localization ***
! *************************************

    ! Here one has to specify the three localization variables
    ! which can be different for each observation type.

    CALL PDAFomi_localize_covar(thisobs, dim_p, locweight, local_range, srange, &
         coords_p, HP_p, HPH)

  END SUBROUTINE localize_covar_OBSTYPE

END MODULE obs_OBSTYPE_pdafomi
