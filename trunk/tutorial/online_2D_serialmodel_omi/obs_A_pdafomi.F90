!$Id$
!> PDAF-OMI observation module for type A observations
!!
!! This module handles operations for one data type (called 'module-type' below):
!! TYPE = A
!!
!! __Observation type A:__
!! The observation type A in this tutorial is the full set of observations except
!! for three observations at the locations (8,5), (12,15), and (4,30) which are
!! removed and only used for observation type B.
!!
!! The different subroutines in this module are for the particular handling of
!! a single observation type.
!! The routines are called by the different call-back routines of PDAF or the
!! interface routines in interface_pdafomi. Most of the routines are generic,
!! so that in practice only 2 routines need to be adapted for a particular
!! data type. There are the routines for the initialization of the observation
!! information (init_dim_obs_f) and for the observation operator (obs_op_f).
!!
!! The module and the routines are named according to the observation type.
!! This allows to distinguish the observation type and the routines in this
!! module from other observation types.
!!
!! These 2 routines usually need to be adapted for the particular observation type:
!! * init_dim_obs_f_TYPE \n
!!           Count number of process-local and full observations; 
!!           initialize vector of observations and their inverse variances;
!!           initialize coordinate array and index array for indices of
!!           observed elements of the state vector.
!! * obs_op_f_TYPE \n
!!           observation operator to get full observation vector of this type. Here
!!           one has to choose a proper observation operator or implement one.
!!
!! The following routines are usually generic. Usually one does not need to modify
!! them, except for the name of the subroutine, which should indicate the 
!! observation type. These routines are part of the module to be able to access
!! module-internal variables and to name them according to the data type.
!!
!! * deallocate_obs_TYPE \n
!!           Deallocate observation arrays after the analysis step. The routine
!!           is mainly generic, but might also deallocate some arrays that are
!!           specific to the module-type observation.
!!
!! These 2 generic routine perform operations for the full observation vector: 
!! * init_obs_f_TYPE \n
!!           Fill the provided full vector of observations with values for
!!           the observation of this type starting from provided offset
!! * init_obsvar_TYPE \n
!!           Compute the mean observation error variance. This is only used with
!!           an adaptive forgetting factor.
!!
!! These 5 routines perform operations for localization:
!! * init_dim_obs_l_TYPE \n
!!           Count number of local observations of module-type according to
!!           their coordinates (distance from local analysis domain). Initialize
!!           module-internal distances and index arrays. Further count offsets
!!           of this observation in full and local observation vectors
!! * init_obs_l_TYPE \n
!!           Initialize module-type part of local observation vector and local
!!           vector of inverse observation error variances
!! * g2l_obs_TYPE \n
!!           Initialize the modules-specific part of the overall local
!!           observation vector from a full observation vector
!! * prodRinvA_l_TYPE \n
!!           Compute the product of the inverse of the local observation error
!!           covariance matrix with the matrix of locally observed ensemble 
!!           perturbations. In addition a localizing weighting can be applied.
!!           The product is computed for the part corresponding to the 
!!           module-type observation.
!! * init_obsvar_l_TYPE \n
!!           Compute the mean observation error variance for local observations. 
!!           This is only used with a local adaptive forgetting factor.
!! 
!! For global square-root filters, this generic routine performs operations for
!! the process-local observation vector
!! * prodRinvA_TYPE \n
!!        Multiply an intermediate matrix of the global filter analysis
!!        with the inverse of the observation error covariance matrix
!! 
!! For the stochastic EnKF, these routines are used for global observations
!! * add_obs_error_TYPE \n
!!        Add the observation error covariance matrix to some other matrix
!! * init_obscovar_TYPE \n
!!        Initialize the full observation error covariance matrix
!! * localize_covar_TYPE \n
!!        Apply covariance localization in localized EnKF
!! 
!! For the particle filter and NETF, and for the localizated NETF these routines
!! are used
!! * likelihood_TYPE \n
!!        Compute global likelihood
!! * likelihood_l_TYPE \n
!!        Compute local likelihood
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE obs_A_pdafomi

  USE mod_parallel_pdaf, &
       ONLY: mype_filter    ! Rank of filter process
  USE PDAFomi_obs_f, &
       ONLY: obs_f          ! Declaration of data type obs_f
  USE PDAFomi_obs_l, &
       ONLY: obs_l          ! Declaration of data type obs_l
 
  IMPLICIT NONE
  SAVE

  ! Variables which are inputs to the module (usually set in init_pdaf)
  LOGICAL :: assim_A        !< Whether to assimilate this data type
  REAL    :: rms_obs_A      !< Observation error standard deviation (for constant errors)

  ! One can declare further variables, e.g. for file names which can
  ! be use-included in init_pdaf() and initialized there.


! *** The following two data types are used inside this module. ***
! *** They are declared in PDAFomi_obs_f and PDAFomi_obs_l and  ***
! *** only listed here for reference

! Data type to define the full observations by internally shared variables of the module
!   TYPE obs_f
!           Mandatory variables to be set in init_dim_obs_f
!      INTEGER :: doassim                   ! Whether to assimilate this observation type
!      INTEGER :: disttype                  ! Type of distance computation to use for localization
!      INTEGER :: ncoord                    ! Number of coordinates use for distance computation
!      LOGICAL :: use_global_obs=.true.     ! Whether to use (T) global full obs. 
!                                           ! or (F) obs. restricted to those relevant for a process domain
!      INTEGER, ALLOCATABLE :: id_obs_p(:,:) ! indices of observed field in state vector
!           
!           Optional variables - they can be set in init_dim_obs_f
!      LOGICAL :: localfilter=.true.        ! Whether a localized filter is used
!      REAL, ALLOCATABLE :: icoeff_p(:,:)   ! Interpolation coefficients for obs. operator
!      REAL, ALLOCATABLE :: domainsize(:)   ! Size of domain for periodicity (<=0 for no periodicity)
!      INTEGER :: obs_err_type=0            ! Type of observation error: (0) Gauss, (1) Laplace
!
!           The following variables are set in the routine PDAFomi_gather_obs_f
!      INTEGER :: dim_obs_p                 ! number of PE-local observations
!      INTEGER :: dim_obs_f                 ! number of full observations
!      REAL, ALLOCATABLE :: obs_f(:)        ! Full observed field
!      REAL, ALLOCATABLE :: ocoord_f(:,:)   ! Coordinates of full observation vector
!      REAL, ALLOCATABLE :: ivar_obs_f(:)   ! Inverse variance of full observations
!      INTEGER :: dim_obs_g                 ! global number of observations 
!                                           ! (only if full obs. are restricted to process domain))
!      INTEGER, ALLOCATABLE :: id_obs_f_lim(:) ! Indices of domain-relevant full obs. in global vector of obs.
!                                           ! (only if full obs. are restricted to process domain))
!
!           Mandatory variable to be set in obs_op_f
!      INTEGER :: off_obs_f                 ! Offset of this observation in overall full obs. vector
!   END TYPE obs_f

! Data type to define the local observations by internally shared variables of the module
!   TYPE obs_l
!      INTEGER :: dim_obs_l                 ! number of local observations
!      INTEGER :: off_obs_l                 ! Offset of this observation in overall local obs. vector
!      INTEGER, ALLOCATABLE :: id_obs_l(:)  ! Indices of local observations in full obs. vector 
!      REAL, ALLOCATABLE :: distance_l(:)   ! Distances of local observations
!      REAL, ALLOCATABLE :: ivar_obs_l(:)   ! Inverse variance of local observations
!   END TYPE obs_l

! Declare instances of observation data types used here
  type(obs_f), private :: thisobs      ! full observation
  type(obs_l), private :: thisobs_l    ! local observation

!$OMP THREADPRIVATE(thisobs_l)


!-------------------------------------------------------------------------------

CONTAINS

!> Set full dimension of observations of module-type
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
!! * thisobs\%localfilter - Whether we use a localized filter 
!!                          (default: .true.; only relevant if the model uses domain decomposition)
!! * thisobs\%domainsize  - Size of domain for periodicity for disttype=1 (<0 for no periodicity)
!! * thisobs\%obs_err_type - Type of observation errors for particle filter and NETF
!! * thisobs\%use_global obs - Whether to use global observations or restrict the observations to the relevant ones
!!                          (default: .true.: use global full observations)
!!
!! The following variables are set in the routine gather_obs_f
!! * thisobs\%dim_obs_p   - PE-local number of module-type observations
!! * thisobs\%dim_obs_f   - full number of module-type observations
!! * thisobs\%obs_f       - full vector of module-type observations
!! * thisobs\%ocoord_f    - coordinates of observations in OBS_MOD_F
!! * thisobs\%ivar_obs_f  - full vector of inverse obs. error variances of module-type
!! * thisobs\%dim_obs_g   - Number of global observations (only if if use_global_obs=.false)
!! * thisobs\%id_obs_f_lim - Ids of full observations in global observations (if use_global_obs=.false)
!!
  SUBROUTINE init_dim_obs_f_A(step, dim_obs_f)

    USE PDAFomi_obs_f, &
         ONLY: PDAFomi_gather_obs_f
    USE mod_model, &
         ONLY: nx, ny
    USE mod_assimilation, &
         ONLY: filtertype, local_range

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)    :: step       !< Current time step
    INTEGER, INTENT(inout) :: dim_obs_f  !< Dimension of full observation vector

! *** Local variables ***
    INTEGER :: i, j                      ! Counters
    INTEGER :: cnt, cnt0                 ! Counters
    INTEGER :: off_nx                    ! Offset of local grid in global domain in x-direction
    INTEGER :: dim_obs_p                 ! Number of process-local observations
    INTEGER :: status                    ! Status flag
    REAL, ALLOCATABLE :: obs_field(:,:)  ! Observation field read from file
    REAL, ALLOCATABLE :: obs_p(:)        ! PE-local observed SST field
    REAL, ALLOCATABLE :: ivar_obs_p(:)   ! PE-local inverse observation error variance
    REAL, ALLOCATABLE :: ocoord_p(:,:)   ! PE-local coordinates of observed SST field
    CHARACTER(len=2) :: stepstr          ! String for time step


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    IF (mype_filter==0) &
         WRITE (*,'(8x,a)') 'Assimilate observations - obs type A'

    ! Store whether to assimilate this observation type (used in routines below)
    IF (assim_A) thisobs%doassim = 1

    ! Specify type of distance computation
    thisobs%disttype = 0   ! 0=Cartesian

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2

    ! Whether we use a local filter (default: .true.) 
    IF (filtertype==6) thisobs%localfilter = .false.  ! (for ESTKF)


! **********************************
! *** Read PE-local observations ***
! **********************************

    ! Read observation field from file
    ALLOCATE(obs_field(ny, nx))

    IF (step<10) THEN
       WRITE (stepstr, '(i1)') step
    ELSE
       WRITE (stepstr, '(i2)') step
    END IF

    OPEN (12, file='../inputs_online/obs_step'//TRIM(stepstr)//'.txt', status='old')
    DO i = 1, ny
       READ (12, *) obs_field(i, :)
    END DO
    CLOSE (12)

    ! Make the observations at (8,5), (12,15) and (4,30) invalid
    ! They will be used in observation type B
    obs_field(8, 5) = -1000.0 
    obs_field(12, 15) = -1000.0 
    obs_field(4, 30) = -1000.0 


! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

    ! *** Count valid observations that lie within the process sub-domain ***

    cnt = 0
    DO j = 1, nx
       DO i= 1, ny
          IF (obs_field(i,j) > -999.0) cnt = cnt + 1
       END DO
    END DO
    dim_obs_p = cnt
    dim_obs_f = cnt

    IF (mype_filter==0) &
         WRITE (*,'(8x, a, i6)') '--- number of full observations', dim_obs_f


    ! *** Initialize vector of observations on the process sub-domain ***
    ! *** Initialize coordinate array of observations on the process sub-domain ***

    ! Allocate process-local observation arrays
    ALLOCATE(obs_p(dim_obs_p))
    ALLOCATE(ivar_obs_p(dim_obs_p))
    ALLOCATE(ocoord_p(2, dim_obs_p))

    ! Allocate process-local index array
    ! This array has a many rows as required for the observation operator
    ! 1 if observations are at grid points; >1 if interpolation is required
    ALLOCATE(thisobs%id_obs_p(1, dim_obs_p))

    cnt = 0
    cnt0 = 0
    DO j = 1, nx
       DO i= 1, ny
          cnt0 = cnt0 + 1
          IF (obs_field(i,j) > -999.0) THEN
             cnt = cnt + 1
             thisobs%id_obs_p(1, cnt) = cnt0
             obs_p(cnt) = obs_field(i, j)
             ocoord_p(1, cnt) = REAL(j)
             ocoord_p(2, cnt) = REAL(i)
          END IF
       END DO
    END DO


! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

    ! *** Set inverse observation error variances ***

    ivar_obs_p(:) = 1.0 / (rms_obs_A*rms_obs_A)


! ****************************************
! *** Gather global observation arrays ***
! ****************************************

    CALL PDAFomi_gather_obs_f(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
         thisobs%ncoord, local_range, dim_obs_f)


! *********************************************************
! *** For twin experiment: Read synthetic observations  ***
! *********************************************************

!     IF (twin_experiment .AND. filtertype/=11) THEN
!        CALL read_syn_obs(file_syntobs_TYPE, dim_obs_f, thisobs%obs_f, 0, 1-mype_filter)
!     END IF


! ********************
! *** Finishing up ***
! ********************

    ! Deallocate all local arrays
    DEALLOCATE(obs_field)
    DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

    ! Arrays in THISOBS have to be deallocated after the analysis step
    ! by a call to deallocate_obs() in prepoststep_pdaf.

  END SUBROUTINE init_dim_obs_f_A



!-------------------------------------------------------------------------------
!> Implementation of observation operator 
!!
!! This routine applies the full observation operator
!! for the type of observations handled in this module
!! It has to append the observations to obsstate_f from
!! position OFFSET_OBS+1. For the return value OFFSET_OBS
!! has to be incremented by the number of added observations.
!!
!! One can choose a proper observation operator from
!! PDAFOMI_OBS_OP or add one to that module or 
!! implement another observation operator here.
!!
!! The order of the calls to this routine for different modules
!! is important because it influences the offset of the 
!! module-type observation in the overall full observation vector.
!!
!! Outputs for within the module are:
!! * thisobs\%off_obs_f - Offset of full module-type observation in overall full obs. vector
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE obs_op_f_A(dim_p, dim_obs_f, state_p, obsstate_f, offset_obs)

    USE PDAFomi_obs_op, &
         ONLY: PDAFomi_obs_op_f_gridpoint

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs_f             !< Dimension of full observed state (all observed fields)
    REAL, INTENT(in)    :: state_p(dim_p)        !< PE-local model state
    REAL, INTENT(inout) :: obsstate_f(dim_obs_f) !< Full observed state
    INTEGER, INTENT(inout) :: offset_obs         !< input: offset of module-type observations in obsstate_f
                                                 !< output: input + number of added observations


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

    IF (thisobs%doassim==1) THEN

       ! Store offset
       thisobs%off_obs_f = offset_obs

       ! observation operator for observed grid point values
       CALL PDAFomi_obs_op_f_gridpoint(thisobs%localfilter, dim_p, dim_obs_f, thisobs%dim_obs_p, &
            thisobs%dim_obs_f, thisobs%id_obs_p, state_p, obsstate_f, offset_obs)
    END IF

  END SUBROUTINE obs_op_f_A




!-------------------------------------------------------------------------------
!++++++      THE FOLLOWING ROUTINES SHOULD BE USABLE WITHOUT CHANGES       +++++
!-------------------------------------------------------------------------------

!> Deallocate observation arrays
!!
!! This routine is called after the analysis step
!! (usually in prepoststep) to deallocate observation
!! arrays before going into the next forecast phase.
!!
!! One only needs to adapt this routine if one introduces
!! additional module arrays in init_dim_obs_f.
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE deallocate_obs_A()

    USE PDAFomi_obs_f, &
         ONLY: PDAFomi_deallocate_obs

    IMPLICIT NONE

    ! Deallocate arrays in full observation type
    CALL PDAFomi_deallocate_obs(thisobs)

  END SUBROUTINE deallocate_obs_A




!-------------------------------------------------------------------------------
!> Initialize full vector of observations
!!
!! This routine initializes the part of the full vector of
!! observations for the current observation type.
!! It has to fill the observations to obsstate_f from
!! position OFFSET_OBS+1. For the return value OFFSET_OBS
!! has to be incremented by the number of added observations.
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE init_obs_f_A(dim_obs_f, obsstate_f, offset_obs)

    USE PDAFomi_obs_f, &
         ONLY: PDAFomi_init_obs_f

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_obs_f             !< Dimension of full observed state (all observed fields)
    REAL, INTENT(inout) :: obsstate_f(dim_obs_f) !< Full observation vector
    INTEGER, INTENT(inout) :: offset_obs         !< input: offset of module-type observations in obsstate_f
                                                 !< output: input + number of added observations


! ******************************************
! *** Initialize full observation vector ***
! ******************************************

    IF (thisobs%doassim == 1) THEN
       CALL PDAFomi_init_obs_f(thisobs, dim_obs_f, obsstate_f, offset_obs)
    END IF

  END SUBROUTINE init_obs_f_A



!-------------------------------------------------------------------------------
!> Compute mean observation error variance
!!
!! This routine will only be called, if the adaptive
!! forgetting factor feature is used. Please note that
!! this is an experimental feature.
!!
!! The routine is called in global filters (like ESTKF)
!! during the analysis or in local filters (e.g. LESTKF)
!! before the loop over local analysis domains 
!! by the routine PDAF_set_forget that estimates an 
!! adaptive forgetting factor.  The routine has to 
!! initialize the mean observation error variance.  
!! For global filters this should be the global mean,
!! while for local filters it should be the mean for the
!! PE-local  sub-domain. (init_obsvar_l_A is the 
!! localized variant for local filters)
!!
!! The implemented functionality is generic. There 
!! should be no changes required as long as the 
!! observation error covariance matrix is diagonal.
!!
!! If the observation counter is zero the computation
!! of the mean variance is initialized. The output is 
!! always the mean variance. If the observation counter
!! is >0 first the variance sum is computed by 
!! multiplying with the observation counter.
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE init_obsvar_A(meanvar, cnt_obs)

    USE PDAFomi_obs_f, &
         ONLY: PDAFomi_init_obsvar_f

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: cnt_obs      !< Observation counter
    REAL, INTENT(inout) :: meanvar         !< Mean variance


! *****************************
! *** Compute mean variance ***
! *****************************

    IF (thisobs%doassim == 1) THEN
       CALL PDAFomi_init_obsvar_f(thisobs, meanvar, cnt_obs)
    END IF

  END SUBROUTINE init_obsvar_A



!-------------------------------------------------------------------------------
!> Set dimension of local obs. vector current type
!!
!! The routine is called during the loop over all local
!! analysis domains. It has to initialize the number of
!! local observations of the module type for the current
!! local analysis domain which is returned in DIM_OBS_L.
!! The routine further stores this number in the module-
!! internal variable thisobs_l%dim_obs_l for later use
!! within the module.
!!
!! The routine further initialized the internal array
!! THISOBS_L\%ID_OBS_L, which stores the indices of the local
!! observations in the full observation vector of the module
!! type. In addition THISOBS_L%DISTANCE_L is initialied,
!! which stores the distances of the local observations
!! from the local analysis domain
!! The offset of the current observation type in the full
!! observation vector is given by OFFSET_OBS_F.
!! Likewise the offset of the current observation type
!! in the local observation vector is given by OFFSET_OBS_L.
!! For their return values the added number of full and
!! local observations has to be added.
!!
!! The implemented functionality using the routine
!! INIT_DIM_OBS_L is generic. There should be no changes
!! required for other observation types.
!!
!! Outputs for within the module are:
!! * thisobs_l\%dim_obs_l  - Local number of module-type observations
!! * thisobs_l\%off_obs_l  - Offset of local module-type observation in overall local obs. vector
!! * thisobs_l\%id_obs_l   - Index module-type local observation in module-type full obs. vector
!! * thisobs_l\%distance_l - Distance of observation from local analysis domain
!!
!! The routine is called by each filter process.
!!
  SUBROUTINE init_dim_obs_l_A(coords_l, lradius, dim_obs_l, offset_obs_l, offset_obs_f)

  USE PDAFomi_obs_l, &
       ONLY: PDAFomi_init_dim_obs_l

    IMPLICIT NONE

! *** Arguments ***
    REAL, INTENT(in) :: coords_l(:)        !< Coordinates of local analysis domain
    REAL, INTENT(in) :: lradius            !< Localization radius
    INTEGER, INTENT(out) :: dim_obs_l      !< Local number of module-type observations
    INTEGER, INTENT(inout) :: offset_obs_l !< input: Offset of module-type obs. in local obs. vector;
                                           !< output: input + number of added observations
    INTEGER, INTENT(inout) :: offset_obs_f !< input: Offset of module-type obs. in full obs. vector;
                                           !< output: input + number of added observations


! ********************************************************************
! *** Initialize local observation dimension and local obs. arrays ***
! ********************************************************************

    IF (thisobs%doassim == 1) THEN

       ! *** Check offset in full observation vector ***
       IF (offset_obs_f /= thisobs%off_obs_f) THEN
          WRITE (*,*) 'ERROR: INCONSISTENT ORDER of observation calls in OBS_OP_F and INIT_DIM_OBS_L! SST'
       END IF

       ! Initialize
       CALL PDAFomi_init_dim_obs_l(thisobs, thisobs_l, coords_l, lradius, dim_obs_l, &
            offset_obs_l, offset_obs_f)
    ELSE
       dim_obs_l = 0
    END IF

  END SUBROUTINE init_dim_obs_l_A




!-------------------------------------------------------------------------------
!> Initialize local observations and inverse variances
!!
!! This routine is called during the loop over
!! all local analysis domains. It has to initialize
!! the local vector of observations for the current
!! local analysis domain and the corresponding vector
!! of inverse observation variances. 
!!
!! The implemented functionality using the routine
!! INIT_OBS_L is generic. There should be no changes
!! required for other observation types as long as
!! the observation error covariance matrix is diagonal.
!!
!! Outputs for within the module are:
!! * thisobs_l\%ivar_obs_l  - Local inverse observation variances
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE init_obs_l_A(dim_obs_l, obs_l)

  USE PDAFomi_obs_l, &
       ONLY: PDAFomi_init_obs_l

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_obs_l         !< Local dimension of observation vector
    REAL, INTENT(inout) :: obs_l(dim_obs_l)  !< Local observation vector


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

    IF (thisobs%doassim == 1) THEN
       CALL PDAFomi_init_obs_l(dim_obs_l, thisobs_l, thisobs, obs_l)
    END IF

  END SUBROUTINE init_obs_l_A



!-------------------------------------------------------------------------------
!> Restrict an observation vector to local analysis domain
!!
!! This routine is called during the loop over
!! all local analysis domains. It has to initialize
!! the local vector of observations for the current
!! local analysis domain and the corresponding vector
!! of inverse observation variances. Further,
!! OFFSET_OBS_L is the offset of the observation of the 
!! module type in the local state vector holding all
!! observation type. The routine has to add the number
!! of module-type observations to it for the return value.
!!
!! The implemented functionality using the routine
!! G2L_OBS is generic. There should be no changes
!! required for other observation types as long as
!! the observation error covariance matrix is diagonal.
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE g2l_obs_A(dim_obs_l, dim_obs_f, obs_f, obs_l)

  USE PDAFomi_obs_l, &
       ONLY: PDAFomi_g2l_obs

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_obs_l         !< Local dimension of observation vector
    INTEGER, INTENT(in) :: dim_obs_f         !< Full dimension of observation vector
    REAL, INTENT(in) :: obs_f(dim_obs_f)     !< Full observation vector
    REAL, INTENT(inout) :: obs_l(dim_obs_l)  !< Local observation vector


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

    IF (thisobs%doassim == 1) THEN
       CALL PDAFomi_g2l_obs(dim_obs_l, thisobs_l%dim_obs_l, thisobs%dim_obs_f, thisobs_l%id_obs_l, &
            obs_f(thisobs%off_obs_f+1:thisobs%off_obs_f+thisobs%dim_obs_f), &
            thisobs_l%off_obs_l, obs_l)
    END IF

  END SUBROUTINE g2l_obs_A



!-------------------------------------------------------------------------------
!> Compute product of inverse of R with some matrix and apply localization
!!
!! The routine is called during the analysis step
!! on each local analysis domain. It has to 
!! compute the product of the inverse of the local
!! observation error covariance matrix with
!! the matrix of locally observed ensemble 
!! perturbations.
!!
!! Next to computing the product, a localizing 
!! weighting ('obervation localization) can be applied
!! to matrix A.
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE prodRinvA_l_A(verbose, dim_obs_l, ncol, locweight, lradius, &
       sradius, A_l, C_l)

  USE PDAFomi_obs_l, &
       ONLY: PDAFomi_prodRinvA_l

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: verbose           !< Verbosity flag
    INTEGER, INTENT(in) :: dim_obs_l         !< Local dimension of observation vector
    INTEGER, INTENT(in) :: ncol              !< Rank of initial covariance matrix
    INTEGER, INTENT(in) :: locweight         !< Localization weight type
    REAL, INTENT(in)    :: lradius           !< localization radius
    REAL, INTENT(in)    :: sradius           !< support radius for weight functions
    REAL, INTENT(inout) :: A_l(:, :)         !< Input matrix
    REAL, INTENT(out)   :: C_l(:, :)         !< Output matrix

! *** Local variable ***
    INTEGER :: idummy   ! Dummy variable to prevent compiler warning

    idummy = dim_obs_l


! ***********************
! *** Compute product ***
! ***********************

    IF (thisobs%doassim == 1) THEN
       CALL PDAFomi_prodRinvA_l(verbose, thisobs_l%dim_obs_l, ncol, &
            locweight, lradius, sradius, &
            thisobs_l%ivar_obs_l, thisobs_l%distance_l, &
            A_l(thisobs_l%off_obs_l+1 : thisobs_l%off_obs_l+thisobs_l%dim_obs_l, :), &
            C_l(thisobs_l%off_obs_l+1 : thisobs_l%off_obs_l+thisobs_l%dim_obs_l, :))
    END IF

  END SUBROUTINE prodRinvA_l_A



!-------------------------------------------------------------------------------
!> Compute local mean observation error variance
!!
!! This routine will only be called, if the local 
!! adaptive forgetting factor feature is used.
!!
!! The routine is called in the loop over all
!! local analysis domains during each analysis
!! by the routine PDAF_set_forget_local that 
!! estimates a local adaptive forgetting factor.
!! The routine has to initialize the mean observation 
!! error variance for the current local analysis 
!! domain.  (See init_obsvar_A for a global variant)
!!
!! The implemented functionality is generic. There
!! should be no changes required as long as the
!! observation error covariance matrix is diagonal.
!!
!! If the observation counter is zero the computation
!! of the mean variance is initialized. The output is 
!! always the mean variance. If the observation counter
!! is >0 first the variance sum is computed by 
!! multiplying with the observation counter.
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE init_obsvar_l_A(meanvar_l, cnt_obs_l)

    USE PDAFomi_obs_l, &
         ONLY: PDAFomi_init_obsvar_l

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: cnt_obs_l      !< Observation counter
    REAL, INTENT(inout) :: meanvar_l         !< Mean variance


! ***********************************
! *** Compute local mean variance ***
! ***********************************

    IF (thisobs%doassim == 1) THEN
       CALL PDAFomi_init_obsvar_l(thisobs_l, meanvar_l, cnt_obs_l)
    END IF

  END SUBROUTINE init_obsvar_l_A



!-------------------------------------------------------------------------------
!> Compute product of inverse of R with some matrix
!!
!! The routine is called during the analysis step
!! of the global square-root filters. It has to 
!! compute the product of the inverse of the
!! process-local observation error covariance matrix
!! with the matrix of process-local observed ensemble 
!! perturbations.
!!
!! This routine assumes a diagonal observation error
!! covariance matrix, but allows for varying observation
!! error variances.
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE prodRinvA_A(ncol, A_p, C_p)

  USE PDAFomi_obs_f, &
         ONLY: PDAFomi_prodRinvA

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: ncol        !< Rank of initial covariance matrix
    REAL, INTENT(in   ) :: A_p(:, :)   !< Input matrix
    REAL, INTENT(out)   :: C_p(:, :)   !< Output matrix


! ***********************
! *** Compute product ***
! ***********************

    IF (thisobs%doassim == 1) THEN

       ! Check process-local observation dimension

       IF (thisobs%dim_obs_p /= thisobs%dim_obs_f) THEN
          ! This error usually happens when thisobs%localfilter=.true.
          IF (thisobs%localfilter) THEN
             WRITE (*,*) 'ERROR: INCONSISTENT value for DIM_OBS_P because localfilter=true.'
          ELSE
             WRITE (*,*) 'ERROR: INCONSISTENT value for DIM_OBS_P'
          END IF
       END IF

       ! compute product
       CALL PDAFomi_prodRinvA(thisobs%dim_obs_f, ncol, thisobs%ivar_obs_f, &
            A_p(thisobs%off_obs_f+1 : thisobs%off_obs_f+thisobs%dim_obs_f, :), &
            C_p(thisobs%off_obs_f+1 : thisobs%off_obs_f+thisobs%dim_obs_f, :))
    END IF

  END SUBROUTINE prodRinvA_A



!-------------------------------------------------------------------------------
!> Add observation error to some matrix
!!
!! The routine is called during the analysis step
!! of the stochastic EnKF. It it provided with a
!! matrix in observation space and has to add the 
!! observation error covariance matrix.
!!
!! This routine assumes a diagonal observation error
!! covariance matrix, but allows for varying observation
!! error variances.
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE add_obs_error_A(dim_obs, C_p)

    USE PDAFomi_obs_f, &
         ONLY: PDAFomi_add_obs_error

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_obs    !< Number of observations
    REAL, INTENT(inout) :: C_p(:, :)  !< Matrix to which R is added


! *************************************************
! *** Check process-local observation dimension ***
! *************************************************

    IF (thisobs%doassim == 1) THEN

       IF (thisobs%dim_obs_p /= thisobs%dim_obs_f) THEN
          ! This error usually happens when thisobs%localfilter=.true.
          IF (thisobs%localfilter) THEN
             WRITE (*,*) 'ERROR: INCONSISTENT value for DIM_OBS_P because localfilter=true.'
          ELSE
             WRITE (*,*) 'ERROR: INCONSISTENT value for DIM_OBS_P'
          END IF
       END IF


! ***********************************************
! *** Add observation error covariance matrix ***
! ***********************************************

       CALL PDAFomi_add_obs_error(thisobs%dim_obs_f, thisobs%ivar_obs_f, C_p, thisobs%off_obs_f)

    END IF

  END SUBROUTINE add_obs_error_A



!-------------------------------------------------------------------------------
!> Initialize global observation error covariance matrix
!!
!! The routine is called during the analysis
!! step when an ensemble of observations is
!! generated by PDAF_enkf_obs_ensemble. 
!! It has to initialize the global observation 
!! error covariance matrix.
!!
!! This routine assumes a diagonal observation error
!! covariance matrix, but allows for varying observation
!! error variances.
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE init_obscovar_A(dim_obs, covar, isdiag)

    USE PDAFomi_obs_f, &
         ONLY: PDAFomi_init_obscovar

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_obs      !< Number of observations
    REAL, INTENT(out) :: covar(:, :)    !< covariance matrix R
    LOGICAL, INTENT(out) :: isdiag      !< Whether matrix R is diagonal


! ***********************************************
! *** Add observation error covariance matrix ***
! ***********************************************

    IF (thisobs%doassim == 1) THEN
       CALL PDAFomi_init_obscovar(thisobs%dim_obs_f, thisobs%ivar_obs_f, thisobs%off_obs_f, &
            covar, isdiag)
    END IF

  END SUBROUTINE init_obscovar_A



!-------------------------------------------------------------------------------
!> Compute likelihood for an ensemble member
!!
!! The routine is called during the analysis step
!! of the NETF or a particle filter.
!! It has to compute the likelihood of the
!! ensemble according to the difference from the
!! observation (residual) and the error distribution
!! of the observations.
!!
!! In general this routine is similar to the routine
!! prodRinvA used for ensemble square root Kalman
!! filters. As an addition to this routine, we here have
!! to evaluate the likelihood weight according the
!! assumed observation error statistics.
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE likelihood_A(dim_obs, obs, resid, lhood)

    USE PDAFomi_obs_f, &
         ONLY: PDAFomi_likelihood

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_obs        !< Number of observations
    REAL, INTENT(in)    :: obs(dim_obs)   !< PE-local vector of observations
    REAL, INTENT(in)    :: resid(dim_obs) !< Input vector of residuum
    REAL, INTENT(inout)   :: lhood        !< Output vector - log likelihood


! **************************
! *** Compute likelihood ***
! **************************

    IF (thisobs%doassim == 1) THEN
       CALL PDAFomi_likelihood(dim_obs, obs, resid, thisobs%ivar_obs_f, lhood, thisobs%obs_err_type)
    END IF

  END SUBROUTINE likelihood_A



!-------------------------------------------------------------------------------
!> Compute local likelihood for an ensemble member
!!
!! The routine is called during the analysis step
!! of the localized NETF.
!! It has to compute the likelihood of the
!! ensemble according to the difference from the
!! observation (residual) and the error distribution
!! of the observations.
!!
!! In addition, a localizing weighting of the 
!! inverse of R by expotential decrease or a 5-th order 
!! polynomial of compact support can be applied. This is 
!! defined by the variables 'locweight', 'local_range, 
!! 'local_range2' and 'srange' in the main program.
!!
!! In general this routine is similar to the routine
!! prodRinvA_l used for ensemble square root Kalman
!! filters. As an addition to this routine, we here have
!! to evaluate the likelihood weight according the
!! assumed observation error statistics.
!!
  SUBROUTINE likelihood_l_A(verbose, dim_obs_l, obs_l, resid_l, &
       locweight, lradius, sradius, lhood_l)

    USE PDAFomi_obs_l, &
         ONLY: PDAFomi_likelihood_l

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: verbose            !< Verbosity flag
    INTEGER, INTENT(in) :: dim_obs_l          !< Number of observations
    REAL, INTENT(in)    :: obs_l(dim_obs_l)   !< PE-local vector of observations
    REAL, INTENT(inout) :: resid_l(dim_obs_l) !< Input vector of residuum
    INTEGER, INTENT(in) :: locweight          !< Localization weight type
    REAL, INTENT(in)    :: lradius            !< localization radius
    REAL, INTENT(in)    :: sradius            !< support radius for weight functions
    REAL, INTENT(inout)   :: lhood_l            !< Output vector - log likelihood


! **************************
! *** Compute likelihood ***
! **************************

    IF (thisobs%doassim == 1) THEN
       CALL PDAFomi_likelihood_l(verbose, dim_obs_l, obs_l, resid_l, locweight, lradius,  &
            sradius, thisobs_l%ivar_obs_l, thisobs_l%distance_l, lhood_l, thisobs%obs_err_type)
    END IF

  END SUBROUTINE likelihood_l_A



!-------------------------------------------------------------------------------
!> Apply covariance localization in LenKF
!!
!! This routine applies a localization matrix B
!! to the matrices HP and HPH^T of the localized EnKF.
!!
  SUBROUTINE localize_covar_A(verbose, dim, dim_obs, &
       locweight, lradius, sradius, coords, HP, HPH, offset_obs)

    USE PDAFomi_obs_l, &
         ONLY: PDAFomi_localize_covar

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: verbose        !< Verbosity flag
    INTEGER, INTENT(in) :: dim            !< State dimension
    INTEGER, INTENT(in) :: dim_obs        !< Number of observations (one or all obs. types)
    INTEGER, INTENT(in) :: locweight      !< Localization weight type
    REAL, INTENT(in)    :: lradius        !< localization radius
    REAL, INTENT(in)    :: sradius        !< support radius for weight functions
    REAL, INTENT(in)    :: coords(:,:)    !< Coordinates of state vector elements
    REAL, INTENT(inout) :: HP(dim_obs, dim)      !< Matrix HP
    REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH
    INTEGER, INTENT(inout) :: offset_obs         !< input: offset of module-type observations in obsstate_f
                                                 !< output: input + number of added observations


! **************************
! *** Compute likelihood ***
! **************************

    IF (thisobs%doassim == 1) THEN
       CALL PDAFomi_localize_covar(verbose, thisobs, dim, &
            locweight, lradius, sradius, coords, HP, HPH, offset_obs)
    END IF

  END SUBROUTINE localize_covar_A

END MODULE obs_A_pdafomi
