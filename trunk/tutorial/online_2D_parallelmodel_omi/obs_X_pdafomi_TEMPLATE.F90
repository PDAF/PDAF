!$Id$
!BOP
!
! !MODULE:
MODULE obs_X_pdafomi_TEMPLATE
!
! !DESCRIPTION:
! This module handles operations for one data type (called 'module-type' below).
!
! There are 10 subroutines for the particular handling of a single data type:
!
! These 3 routines usually need to be adapted for the particular observation type:
! init_dim_obs_f_TYPE 
!           Count number of process-local and full observations; 
!           initialize vector of observations and their inverse variances;
!           initialize coordinate array and index array for indices of
!           observed elements of the state vector.
! obs_op_f_TYPE 
!           observation operator to get full observation vector of this type. Here
!           one has to choose a proper observation operator or implement one.
! deallocate_obs_TYPE 
!           Deallocate observation arrays after the analysis step. The routine
!           is mainly generic, but might also deallocate some additional arrays
!           that are specific to the module-type observation.
!
! The following routines are usually generic so that one does not need to modify
! them, except for the name of the subroutine, which should indicate the 
! observation type. These routines are part of the module to be able to access
! module-internal variables and to name them according to the data type.
!
! These 2 generic routine perform operations for the full observation vector: 
! init_obs_f_TYPE
!           Fill the provided full vector of observations with values for
!           the observation of this type starting from provided offset
! init_obsvar_TYPE
!           Compute the mean observation error variance. This is only used with
!           an adaptive forgetting factor.
!
! These 5 routines perform operations for localization:
! init_dim_obs_l_TYPE 
!           Count number of local observations of module-type according to
!           their coordinates (distance from local analysis domain). Initialize
!           module-internal distances and index arrays. Further count offsets
!           of this observation in full and local observation vectors
! init_obs_l_TYPE 
!           Initialize module-type part of local observation vector and local
!           vector of inverse observation error variances
! g2l_obs_TYPE
!           Initialize the modules-specific part of the overall local
!           observation vector from a full observation vector
! prodRinvA_l_TYPE
!           Compute the product of the inverse of the local observation error
!           covariance matrix with the matrix of locally observed ensemble 
!           perturbations. In addition a localizing weighting can be applied.
!           The product is computed for the part corresponding to the 
!           module-type observation.
! init_obsvar_l_TYPE
!           Compute the mean observation error variance for local observations. 
!           This is only used with a local adaptive forgetting factor.
!
! The routines are called by the different call-back routines of PDAF. To be
! able to distinguish the observation type and the routines in this module,
! we recommend to rename the module according to the observation module-type.
! Further,we recommend to replace 'TYPE' in the routine names according to the
! type of the observation so that they can be identified when 
! 
! !REVISION HISTORY:
! 2019-06 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, &
       ONLY: mype_filter
  USE PDAFomi_obs_f, &
       ONLY: obs_f        ! Declaration of data type obs_f
  USE PDAFomi_obs_l, &
       ONLY: obs_l        ! Declaration of data type obs_l
 
  IMPLICIT NONE
  SAVE

  ! Variables which are inputs to the module (usually set in init_pdaf)
  LOGICAL :: assim_TYPE        ! Whether to assimilate this data type
  REAL    :: rms_obs_TYPE      ! Observation error standard deviation (for constant errors)

  ! One can declare further variables, e.g. for file names which can
  ! be use-included in init_pdaf() and initialized there.


! *** The following two data types are used inside this module. ***
! *** They are declared in PDAFomi_obs_f and PDAFomi_obs_l and  ***
! *** only listed here for reference

! Data type to define the full observations by internally shared variables of the module
!   type obs_f
!      INTEGER :: dim_obs_p                 ! number of PE-local observations
!      INTEGER :: dim_obs_f                 ! number of full observations
!      INTEGER :: off_obs_f                 ! Offset of this observation in overall full obs. vector
!      INTEGER, ALLOCATABLE :: id_obs_p(:,:) ! indices of observed field in state vector
!      REAL, ALLOCATABLE :: icoeff_p(:,:)   ! Interpolation coefficients for obs. operator
!      REAL, ALLOCATABLE :: obs_f(:)        ! Full observed field
!      REAL, ALLOCATABLE :: ocoord_f(:,:)   ! Coordinates of full observation vector
!      REAL, ALLOCATABLE :: ivar_obs_f(:)   ! Inverse variance of full observations
!      INTEGER :: disttype                  ! Type of distance computation to use for localization
!      INTEGER :: ncoord                    ! Number of coordinates use for distance computation
!   end type obs_f

! Data type to define the local observations by internally shared variables of the module
!   type obs_l
!      INTEGER :: dim_obs_l                 ! number of local observations
!      INTEGER :: off_obs_l                 ! Offset of this observation in overall local obs. vector
!      INTEGER, ALLOCATABLE :: id_obs_l(:)  ! Indices of local observations in full obs. vector 
!      REAL, ALLOCATABLE :: distance_l(:)   ! Distances of local observations
!      REAL, ALLOCATABLE :: ivar_obs_l(:)   ! Inverse variance of local observations
!   end type obs_l

! Declare instances of observation data types used here
  type(obs_f), private :: thisobs   
  type(obs_l), private :: thisobs_l

!$OMP THREADPRIVATE(thisobs_l)

! EOP
!-------------------------------------------------------------------------------

CONTAINS

!BOP
!
! !ROUTINE: init_dim_obs_f_TYPE --- Set full dimension of observations of module-type
!
! !INTERFACE:
  SUBROUTINE init_dim_obs_f_TYPE(step, dim_obs_f)

! !DESCRIPTION:
! The routine is called by each filter process.
! at the beginning of the analysis step before 
! the loop through all local analysis domains.
! 
! It has to count the number of observations of the
! observation type handled in this module according
! to the current time step for all observations 
! required for the analyses in the loop over all local 
! analysis domains on the PE-local state domain.
!
! Outputs for within the module are:
! thisobs%dim_obs_p  - PE-local number of module-type observations
! thisobs%dim_obs_f  - full number of module-type observations
! thisobs%id_obs_p   - index of module-type observation in PE-local state vector
! thisobs%obs_f      - full vector of module-type observations
! thisobs%ocoord_f   - coordinates of observations in OBS_MOD_F
! thisobs%ivar_obs_f - full vector of inverse obs. error variances of module-type
! thisobs%disttype   - type of distance computation for localization with this observaton
! thisobs%ncoord     - number of coordinates used for distance computation
!
! Optional is the use of
! thisobs%icoeff_p   - Interpolation coefficients for obs. operator (only if interpolation is used)
!
! !REVISION HISTORY:
! 2019-06 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Here, you can use-include modules from the model, e.g. for mesh information.
! Alternatively you could include these as subroutine arguments

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in)    :: step       ! Current time step
    INTEGER, INTENT(inout) :: dim_obs_f  ! Dimension of full observation vector
!EOP

! Local variables
    INTEGER :: i, j                      ! Counters
    INTEGER :: dim_obs_p                 ! Number of process-local observations
    INTEGER :: status                    ! Status flag
    REAL, ALLOCATABLE :: obs_p(:)        ! PE-local observed SST field
    REAL, ALLOCATABLE :: ivar_obs_p(:)   ! PE-local inverse observation error variance
    REAL, ALLOCATABLE :: ocoord_p(:,:)   ! PE-local coordinates of observed SST field


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    IF (mype_filter==0) &
         WRITE (*,'(8x,a)') 'Assimilate observations - obs type A'

    ! Specify type of distance computation
    thisobs%disttype = 0   ! 0=Cartesian

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2


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

!    IF (ALLOCATED(obs_p)) DEALLOCATE(obs_p)
!    ALLOCATE(obs_p(dim_obs_p))

!    obs_p = ....
!    ivar_obs_p = ....


  ! *** Initialize coordinate array of observations on the process sub-domain ***

!    IF (ALLOCATED(ocoord_p)) DEALLOCATE(ocoord_p)
!    ALLOCATE(ocoord_p(dim_obs_p))


  ! *** Initialize process local index array                                ***
  ! *** This array holds the information which elements of the state vector ***
  ! *** are used in the observation optation. This array has a many rows as ***
  ! *** required for the observation operator, i.e. 1 if observations are   ***
  ! *** at grid points; >1 if interpolation is required                     ***

  ! The initialization is done locally for each process sub-domain and later
  ! used in the observation operator. 
  ! Examples:
  ! 1. If the observations are model fields located at grid points, one should
  !   initialize the index array thisobs%id_obs_p with one row so that it contains 
  !   the indices of the observed field values in the process-local state vector
  !   (state_p). Then one can use the observation operator OBS_OP_F_GRIDPOINT 
  !   provided by MOD_OBS_OP_GENERAL_PDAF.
  ! 2. If the observations are the average of model fields located at grid points,
  !   one should initialize the index array thisobs%id_obs_p with as many rows as 
  !   values to be averaged. Each column of the arrays then contains the indices of
  !   the elements of the process-local state vector that have to be averaged. With
  !   this index array one can then use the observation operator OBS_OP_F_GRIDAVG
  !   provided by MOD_OBS_OP_GENERAL_PDAF.
  ! 3. If model values need to be interpolated to the observation location
  !   one should initialize the index array thisobs%id_obs_p with as many rows as 
  !   values are required in the interpolationto be averaged. Each column of the 
  !   array then contains the indices of elements of the process-local state vector 
  !   that are used in the interpolation.
  ! Below, you need to replace NROWS by the number of required rows

!    IF (ALLOCATED(thisobs%id_obs_p)) DEALLOCATE(thisobs%id_obs_p)
!    ALLOCATE(thisobs%id_obs_p( NROWS , dim_obs_p))

!    thisobs%id_obs_p = ...


! **********************************************************************
! *** Initialize interpolation coefficients for observation operator ***
! **********************************************************************

  ! This initialization is only required if an observation operator
  ! with interpolation is used. The coefficients should be determined
  ! here instead of the observation operator, because the operator is
  ! called for each ensemble member while init_dim_obs_f is only called
  ! once.

  ! Allocate array of interpolation coefficients. As thisobs%id_obs_p, the number
  ! of rows corresponds to the number of grid points using the the interpolation

!    IF (ALLOCATED(thisobs%icoeff_p)) DEALLOCATE(thisobs%icoeff_p)
!    ALLOCATE(thisobs%icoeff_p(4, dim_obs_p))

  ! Ensure that the order of the coefficients is consistent with the
  ! indexing in thisobs%id_obs_p. Further ensure that the order is consistent
  ! with the assumptions used in the observation operator.

!    thisobs%icoeff_p = ...


! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

!    IF (ALLOCATED(ivar_obs_p)) DEALLOCATE(ivar_obs_p)
!    ALLOCATE(ivar_obs_p(dim_obs_p))
    
!    ivar_obs_p = ...


! ****************************************
! *** Gather global observation arrays ***
! ****************************************

  ! This part should be the same for all observations 
  ! You need to replace NROWS by the number of rows in the coordinate array

  ! *** Initialize global dimension of observation vector ***
  CALL PDAF_gather_dim_obs_f(dim_obs_p, dim_obs_f)

  ! *** Gather full observation vector and corresponding coordinates ***

  ! Allocate full observation arrays
  ! The arrays are deallocated in deallocate_obs in this module
  ALLOCATE(thisobs%obs_f(dim_obs_f))
  ALLOCATE(thisobs%ivar_obs_f(dim_obs_f))
  ALLOCATE(thisobs%ocoord_f(NROWS, dim_obs_f))

  CALL PDAF_gather_obs_f_flex(dim_obs_p, dim_obs_f, obs_p, thisobs%obs_f, status)
  CALL PDAF_gather_obs_f_flex(dim_obs_p, dim_obs_f, ivar_obs_p, thisobs%ivar_obs_f, status)
  CALL PDAF_gather_obs_f2_flex(dim_obs_p, dim_obs_f, ocoord_p, thisobs%ocoord_f, NROWS, status)


! *********************************************************
! *** For twin experiment: Read synthetic observations  ***
! *********************************************************

!   IF (twin_experiment .AND. filtertype/=11) THEN
!      CALL read_syn_obs(file_syntobs_TYPE, dim_obs_f, thisobs%obs_f, 0, 1-mype_filter)
!   END IF


! ********************
! *** Finishing up ***
! ********************

    ! store full and PE-local observation dimension in module variables
    thisobs%dim_obs_p = dim_obs_p
    thisobs%dim_obs_f = dim_obs_f

    ! Clean up arrays
!    DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

  END SUBROUTINE init_dim_obs_f_TYPE



!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: obs_op_f_TYPE --- Implementation of observation operator 
!
! !INTERFACE:
  SUBROUTINE obs_op_f_TYPE(dim_p, dim_obs_f, state_p, obsstate_f, offset_obs)

! !DESCRIPTION:
! This routine applies the full observation operator
! for the type of observations handled in this module
! It has to append the observations to obsstate_f from
! position OFFSET_OBS+1. For the return value OFFSET_OBS
! has to be incremented by the number of added observations.
!
! One can choose a proper observation operator from
! PDAFOMI_OBS_OP or add one to that module or 
! implement another observation operator here.
!
! The order of the calls to this routine for different modules
! is important because it influences the offset of the 
! module-type observation in the overall full observation vector.
!
! Outputs for within the module are:
! thisobs%off_obs_f    - Offset of full module-type observation in overall full obs. vector
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2019-06 - Lars Nerger - Initial code from restructuring observation routines
! Later revisions - see svn log
!
! !USES:
    USE PDAFomi_obs_op, &
         ONLY: obs_op_f_gridpoint

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: dim_p                 ! PE-local dimension of state
    INTEGER, INTENT(in) :: dim_obs_f             ! Dimension of full observed state (all observed fields)
    REAL, INTENT(in)    :: state_p(dim_p)        ! PE-local model state
    REAL, INTENT(inout) :: obsstate_f(dim_obs_f) ! Full observed state
    INTEGER, INTENT(inout) :: offset_obs         ! input: offset of module-type observations in obsstate_f
                                                 ! output: input + number of added observations
!EOP


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

    ! Store offset
    thisobs%off_obs_f = offset_obs

    
    ! Choose suitable observation operator from the module PDAFomi_obs_op
    ! or implement your own

    ! observation operator for observed grid point values
    CALL obs_op_f_gridpoint(dim_p, dim_obs_f, thisobs%dim_obs_p, thisobs%dim_obs_f, &
         thisobs%id_obs_p, state_p, obsstate_f, offset_obs)

  END SUBROUTINE obs_op_f_TYPE



!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: deallocate_obs_TYPE --- Deallocate observation errors
!
! !INTERFACE:
  SUBROUTINE deallocate_obs_TYPE()

! !DESCRIPTION:
! This routine is called after the analysis step
! (usually in prepoststep) to deallocate observation
! arrays before going into the next forecast phase.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2019-06 - Lars Nerger - Initial code from restructuring observation routines
! Later revisions - see svn log
!
! !USES:
    USE PDAFomi_obs_f, &
         ONLY: deallocate_obs

    IMPLICIT NONE

    ! Deallocate arrays in full observation type
    CALL deallocate_obs(thisobs)

  END SUBROUTINE deallocate_obs_TYPE




!-------------------------------------------------------------------------------
!++++++      THE FOLLOWING ROUTINES SHOULD BE USABLE WITHOUT CHANGES       +++++
!-------------------------------------------------------------------------------

!BOP
!
! !ROUTINE: init_obs_f_TYPE --- Initialize full vector of observations
!
! !INTERFACE:
  SUBROUTINE init_obs_f_TYPE(dim_obs_f, obsstate_f, offset_obs)

! !DESCRIPTION:
! This routine initializes the part of the full vector of
! observations for the current observation type.
! It has to fill the observations to obsstate_f from
! position OFFSET_OBS+1. For the return value OFFSET_OBS
! has to be incremented by the number of added observations.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2019-09 - Lars Nerger - Initial code from restructuring observation routines
! Later revisions - see svn log
!
! !USES:
    USE PDAFomi_obs_f, &
         ONLY: init_obs_f

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: dim_obs_f             ! Dimension of full observed state (all observed fields)
    REAL, INTENT(inout) :: obsstate_f(dim_obs_f) ! Full observation vector
    INTEGER, INTENT(inout) :: offset_obs         ! input: offset of module-type observations in obsstate_f
                                                 ! output: input + number of added observations
!EOP


! ******************************************
! *** Initialize full observation vector ***
! ******************************************

    CALL init_obs_f(thisobs, dim_obs_f, obsstate_f, offset_obs)

  END SUBROUTINE init_obs_f_TYPE



!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_obsvar_TYPE --- Compute mean observation error variance
!
! !INTERFACE:
  SUBROUTINE init_obsvar_TYPE(meanvar, cnt_obs)

! !DESCRIPTION:
! This routine will only be called, if the adaptive
! forgetting factor feature is used. Please note that
! this is an experimental feature.
!
! The routine is called in global filters (like ESTKF)
! during the analysis or in local filters (e.g. LESTKF)
! before the loop over local analysis domains 
! by the routine PDAF\_set\_forget that estimates an 
! adaptive forgetting factor.  The routine has to 
! initialize the mean observation error variance.  
! For global filters this should be the global mean,
! while for local filters it should be the mean for the
! PE-local  sub-domain. (init_obsvar_l_TYPE is the 
! localized variant for local filters)
!
! The implemented functionality is generic. There 
! should be no changes required as long as the 
! observation error covariance matrix is diagonal.
!
! If the observation counter is zero the computation
! of the mean variance is initialized. The output is 
! always the mean variance. If the observation counter
! is >0 first the variance sum is computed by 
! multiplying with the observation counter.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2019-09 - Lars Nerger - Initial code from restructuring observation routines
! Later revisions - see svn log
!
! !USES:
    USE PDAFomi_obs_f, &
         ONLY: init_obsvar_f

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(inout) :: cnt_obs      ! Observation counter
    REAL, INTENT(inout) :: meanvar         ! Mean variance
!EOP


! *****************************
! *** Compute mean variance ***
! *****************************

    CALL init_obsvar_f(thisobs, meanvar, cnt_obs)

  END SUBROUTINE init_obsvar_TYPE



!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_dim_obs_l_TYPE --- Set dimension of local obs. vector current type
!
! !INTERFACE:
  SUBROUTINE init_dim_obs_l_TYPE(coords_l, lradius, dim_obs_l, offset_obs_l, offset_obs_f)

! !DESCRIPTION:
! The routine is called during the loop over all local
! analysis domains. It has to initialize the number of
! local observations of the module type for the current
! local analysis domain which is returned in DIM_OBS_L.
! The routine further stores this number in the module-
! internal variable thisobs_l%dim_obs_l for later use
! within the module.
!
! The routine further initializes the internal array
! THISOBS_L%ID_OBS_L, which stores the indices of the local
! observations in the full observation vector of the module
! type. In addition THISOBS_L%DISTANCE_L is initialied,
! which stores the distances of the local observations
! from the local analysis domain
! The offset of the current observation type in the full
! observation vector is given by OFFSET_OBS_F.
! Likewise the offset of the current observation type
! in the local observation vector is given by OFFSET_OBS_L.
! For their return values the added number of full and
! local observations has to be added.
!
! The implemented functionality using the routine
! INIT_DIM_OBS_L is generic. There should be no changes
! required for other observation types.
!
! Outputs for within the module are:
! thisobs_l%dim_obs_l  - Local number of module-type observations
! thisobs_l%off_obs_l  - Offset of local module-type observation in overall local obs. vector
! thisobs_l%id_obs_l   - Index module-type local observation in module-type full obs. vector
! thisobs_l%distance_l - Distance of observation from local analysis domain
!
! The routine is called by each filter process.
!
! !REVISION HISTORY:
! 2019-06 - Lars Nerger - Initial code from restructuring observation routines
! Later revisions - see svn log
!
! !USES:
  USE PDAFomi_obs_l, &
       ONLY: init_dim_obs_l

    IMPLICIT NONE

! !ARGUMENTS:
    REAL, INTENT(in) :: coords_l(:)        ! Coordinates of local analysis domain
    REAL, INTENT(in) :: lradius            ! Localization radius
    INTEGER, INTENT(out) :: dim_obs_l      ! Local number of module-type observations
    INTEGER, INTENT(inout) :: offset_obs_l ! input: Offset of module-type obs. in local obs. vector
                                           ! output: input + number of added observations
    INTEGER, INTENT(inout) :: offset_obs_f ! input: Offset of module-type obs. in full obs. vector
                                           ! output: input + number of added observations
!EOP


! ***********************************************
! *** Check offset in full observation vector ***
! ***********************************************

    IF (offset_obs_f /= thisobs%off_obs_f) THEN
       WRITE (*,*) 'ERROR: INCONSISTENT ORDER of observation calls in OBS_OP_F and INIT_DIM_OBS_L!'
    END IF


! ********************************************************************
! *** Initialize local observation dimension and local obs. arrays ***
! ********************************************************************

    CALL init_dim_obs_l(thisobs, thisobs_l, coords_l, lradius, dim_obs_l, &
         offset_obs_l, offset_obs_f)

  END SUBROUTINE init_dim_obs_l_TYPE




!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_obs_l_TYPE --- Initialize local observations and inverse variances
!
! !INTERFACE:
  SUBROUTINE init_obs_l_TYPE(dim_obs_l, obs_l)

! !DESCRIPTION:
! This routine is called during the loop over
! all local analysis domains. It has to initialize
! the local vector of observations for the current
! local analysis domain and the corresponding vector
! of inverse observation variances. 
!
! The implemented functionality using the routine
! INIT_OBS_L is generic. There should be no changes
! required for other observation types as long as
! the observation error covariance matrix is diagonal.
!
! Outputs for within the module are:
! thisobs_l%ivar_obs_l  - Local inverse observation variances
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2019-06 - Lars Nerger - Initial code from restructuring observation routines
! Later revisions - see svn log
!
! !USES:
  USE PDAFomi_obs_l, &
       ONLY: init_obs_l

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: dim_obs_l         ! Local dimension of observation vector
    REAL, INTENT(inout) :: obs_l(dim_obs_l)  ! Local observation vector
!EOP


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

    CALL init_obs_l(dim_obs_l, thisobs_l, thisobs, obs_l)

  END SUBROUTINE init_obs_l_TYPE



!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: g2l_obs_TYPE --- Restrict an obs. vector to local analysis domain
!
! !INTERFACE:
  SUBROUTINE g2l_obs_TYPE(dim_obs_l, dim_obs_f, obs_f, obs_l)

! !DESCRIPTION:
! This routine is called during the loop over
! all local analysis domains. It has to initialize
! the local vector of observations for the current
! local analysis domain and the corresponding vector
! of inverse observation variances. Further,
! OFFSET_OBS_L is the offset of the observation of the 
! module type in the local state vector holding all
! observation type. The routine has to add the number
! of module-type observations to it for the return value.
!
! The implemented functionality using the routine
! G2L_OBS is generic. There should be no changes
! required for other observation types as long as
! the observation error covariance matrix is diagonal.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2019-06 - Lars Nerger - Initial code from restructuring observation routines
! Later revisions - see svn log
!
! !USES:
  USE PDAFomi_obs_l, &
       ONLY: g2l_obs

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: dim_obs_l         ! Local dimension of observation vector
    INTEGER, INTENT(in) :: dim_obs_f         ! Full dimension of observation vector
    REAL, INTENT(in) :: obs_f(dim_obs_f)     ! Full observation vector
    REAL, INTENT(inout) :: obs_l(dim_obs_l)  ! Local observation vector
!EOP


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

    ! Initialize local observations
    CALL g2l_obs(dim_obs_l, thisobs_l%dim_obs_l, thisobs%dim_obs_f, thisobs_l%id_obs_l, &
         obs_f(thisobs%off_obs_f+1:thisobs%off_obs_f+thisobs%dim_obs_f), &
         thisobs_l%off_obs_l, obs_l)

  END SUBROUTINE g2l_obs_TYPE



!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: prodRinvA_l_TYPE --- Restrict an obs. vector to local analysis domain
!
! !INTERFACE:
  SUBROUTINE prodRinvA_l_TYPE(verbose, dim_obs_l, ncol, locweight, lradius, &
       sradius, A_l, C_l)

! !DESCRIPTION:
! This routine is called during the loop over
! all local analysis domains. It has to initialize
! the local vector of observations for the current
! local analysis domain and the corresponding vector
! of inverse observation variances. Further,
! OFFSET_OBS_L is the offset of the observation of the 
! module type in the local state vector holding all
! observation type. The routine has to add the number
! of module-type observations to it for the return value.
!
! The implemented functionality using the routine
! INIT_OBS_L is generic. There should be no changes
! required for other observation types as long as
! the observation error covariance matrix is diagonal.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2019-06 - Lars Nerger - Initial code from restructuring observation routines
! Later revisions - see svn log
!
! !USES:
  USE PDAFomi_obs_l, &
       ONLY: prodRinvA_l

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: verbose           ! Verbosity flag
    INTEGER, INTENT(in) :: dim_obs_l         ! Local dimension of observation vector
    INTEGER, INTENT(in) :: ncol              ! Rank of initial covariance matrix
    INTEGER, INTENT(in) :: locweight         ! Localization weight type
    REAL, INTENT(in)    :: lradius           ! localization radius
    REAL, INTENT(in)    :: sradius           ! support radius for weight functions
    REAL, INTENT(inout) :: A_l(:, :)         ! Input matrix
    REAL, INTENT(out)   :: C_l(:, :)         ! Output matrix
!EOP

! Local variable
    INTEGER :: idummy   ! Dummy variable to present compiler warning

    idummy = dim_obs_l


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

    ! Initialize local observations
    CALL prodRinvA_l(verbose, thisobs_l%dim_obs_l, ncol, locweight, lradius, sradius, &
         thisobs_l%ivar_obs_l, thisobs_l%distance_l, &
         A_l(thisobs_l%off_obs_l+1 : thisobs_l%off_obs_l+thisobs_l%dim_obs_l, :), &
         C_l(thisobs_l%off_obs_l+1 : thisobs_l%off_obs_l+thisobs_l%dim_obs_l, :))

  END SUBROUTINE prodRinvA_l_TYPE



!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_obsvar_l_TYPE --- Compute local mean observation error variance
!
! !INTERFACE:
  SUBROUTINE init_obsvar_l_TYPE(meanvar_l, cnt_obs_l)

! !DESCRIPTION:
! This routine will only be called, if the local 
! adaptive forgetting factor feature is used.
!
! The routine is called in the loop over all
! local analysis domains during each analysis
! by the routine PDAF\_set\_forget\_local that 
! estimates a local adaptive forgetting factor.
! The routine has to initialize the mean observation 
! error variance for the current local analysis 
! domain.  (See init_obsvar_TYPE for a global variant)
!
! The implemented functionality is generic. There
! should be no changes required as long as the
! observation error covariance matrix is diagonal.
!
! If the observation counter is zero the computation
! of the mean variance is initialized. The output is 
! always the mean variance. If the observation counter
! is >0 first the variance sum is computed by 
! multiplying with the observation counter.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2019-09 - Lars Nerger - Initial code from restructuring observation routines
! Later revisions - see svn log
!
! !USES:
    USE PDAFomi_obs_l, &
         ONLY: init_obsvar_l

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(inout) :: cnt_obs_l      ! Observation counter
    REAL, INTENT(inout) :: meanvar_l         ! Mean variance
!EOP


! ***********************************
! *** Compute local mean variance ***
! ***********************************

    CALL init_obsvar_l(thisobs_l, meanvar_l, cnt_obs_l)

  END SUBROUTINE init_obsvar_l_TYPE

END MODULE obs_X_pdafomi_TEMPLATE
