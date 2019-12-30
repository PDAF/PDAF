!$Id: obs_C_pdafomi.F90 251 2019-11-19 08:43:39Z lnerger $
!> PDAF-OMI observation module for type B observations
!!
!! This module handles operations for one data type (called 'module-type' below):
!! TYPE = C
!!
!! __Observation type C:__
!! The observation type C in this tutorial is a set of observations that are 
!! placed in between the grid points. It demonstrates the use of the observation
!! operator for bi-linear interpolation.
!!
!! There are 10 subroutines for the particular handling of a single data type.
!! The routines are called by the different call-back routines of PDAF or the
!! interface routines in interface_pdafomi. 
!!
!! The module and the routines are named according to the observation type.
!! This allows to distinguish the observation type and the routines in this
!! module from other observation types.
!!
!! These 3 routines usually need to be adapted for the particular observation type:
!! * init_dim_obs_f_TYPE \n
!!           Count number of process-local and full observations; 
!!           initialize vector of observations and their inverse variances;
!!           initialize coordinate array and index array for indices of
!!           observed elements of the state vector.
!! * obs_op_f_TYPE \n
!!           observation operator to get full observation vector of this type. Here
!!           one has to choose a proper observation operator or implement one.
!! * deallocate_obs_TYPE \n
!!           Deallocate observation arrays after the analysis step. The routine
!!           is mainly generic, but might also deallocate some arrays that are
!!           specific to the module-type observation.
!!
!! The following routines are usually generic so that one does not need to modify
!! them, except for the name of the subroutine, which should indicate the 
!! observation type. These routines are part of the module to be able to access
!! module-internal variables and to name them according to the data type.
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
!! __Revision history:__
!! * 2019-12 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE obs_C_pdafomi

  USE mod_parallel_pdaf, &
       ONLY: mype_filter    ! Rank of filter process
  USE PDAFomi_obs_f, &
       ONLY: obs_f          ! Declaration of data type obs_f
  USE PDAFomi_obs_l, &
       ONLY: obs_l          ! Declaration of data type obs_l
 
  IMPLICIT NONE
  SAVE

  ! Variables which are inputs to the module (usually set in init_pdaf)
  LOGICAL :: assim_C        !< Whether to assimilate this data type
  REAL    :: rms_obs_C      !< Observation error standard deviation (for constant errors)

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
  type(obs_f), private :: thisobs      ! full observation
  type(obs_l), private :: thisobs_l    ! local observation

!$OMP THREADPRIVATE(thisobs_l)

! EOP
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
!! Outputs for within the module are:
!! * thisobs\%dim_obs_p  - PE-local number of module-type observations
!! * thisobs\%dim_obs_f  - full number of module-type observations
!! * thisobs\%id_obs_p   - index of module-type observation in PE-local state vector
!! * thisobs\%obs_f      - full vector of module-type observations
!! * thisobs\%ocoord_f   - coordinates of observations in OBS_MOD_F
!! * thisobs\%ivar_obs_f - full vector of inverse obs. error variances of module-type
!! * thisobs\%disttype   - type of distance computation for localization with this observaton
!! * thisobs\%ncoord     - number of coordinates used for distance computation
!!
!! Optional is the use of
!! * thisobs\%icoeff_p   - Interpolation coefficients for obs. operator (only if interpolation is used)
!!
  SUBROUTINE init_dim_obs_f_C(step, dim_obs_f)

    USE mod_model, &
         ONLY: ny
    USE PDAFomi_obs_op, &
         ONLY: get_interp_coeff_lin

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)    :: step      !< Current time step
    INTEGER, INTENT(inout) :: dim_obs_f !< Dimension of full observation vector

! *** Local variables ***
    INTEGER :: i                        ! Counters
    INTEGER :: nobs                     ! Number of observations in file
    INTEGER :: dim_obs_p                ! number of process-local observations
    REAL, ALLOCATABLE :: obs_list(:,:)  ! List of observations field read from file
    CHARACTER(len=2) :: stepstr         ! String for time step
    REAL :: gcoords(4,2)                ! Grid point coordinated to compute interpolation coeffs


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    IF (mype_filter==0) &
         WRITE (*,'(8x,a)') 'Assimilate observations - obs type C: interpolated observations'

    ! Specify type of distance computation
    thisobs%disttype = 0   ! 0=Cartesian

    ! Number of coordinates used for distance computation
    thisobs%ncoord = 2


! **********************************
! *** Read PE-local observations ***
! **********************************

    ! Open file
    IF (step<10) THEN
       WRITE (stepstr, '(i1)') step
    ELSE
       WRITE (stepstr, '(i2)') step
    END IF
    OPEN (12, file='../inputs_online/iobs_step'//TRIM(stepstr)//'.txt', status='old')

    ! Read number of observations
    READ (12, *) nobs

    ! Read observations and coordinates
    ALLOCATE(obs_list(nobs, 3))
    DO i = 1, nobs
       READ (12, *) obs_list(i, :)
    END DO
    CLOSE (12)


! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

    ! *** Count valid observations that lie within the process sub-domain ***
    dim_obs_p = nobs
    dim_obs_f = nobs

    IF (mype_filter==0) &
         WRITE (*,'(8x, a, i6)') '--- number of full observations', dim_obs_f


    ! *** Initialize vector of observations on the process sub-domain ***
    ! *** Initialize coordinate array of observations on the process sub-domain ***

    ! Allocate full observation arrays
    ! The arrays are deallocated in deallocate_obs in this module
    ALLOCATE(thisobs%obs_f(dim_obs_f))
    ALLOCATE(thisobs%ivar_obs_f(dim_obs_f))
    ALLOCATE(thisobs%ocoord_f(2, dim_obs_f))

    ! Allocate process-local index array
    ! This array has a many rows as required for the observation operator
    ! 1 if observations are at grid points; >1 if interpolation is required
    IF (ALLOCATED(thisobs%id_obs_p)) DEALLOCATE(thisobs%id_obs_p)
    ALLOCATE(thisobs%id_obs_p(4, dim_obs_p))

    DO i= 1, dim_obs_p

       ! Observation value and coordinates
       thisobs%obs_f(i) = obs_list(i,1)
       thisobs%ocoord_f(1, i) = obs_list(i,2)
       thisobs%ocoord_f(2, i) = obs_list(i,3)
     
       ! State vector indices of 4 grid points in which box the observation lies
       ! Note: These indices have to be consistent with the coordinates used to
       ! compute the interpolation coefficients (see below)
       thisobs%id_obs_p(1, i) = (FLOOR(obs_list(i,2))-1)*ny + FLOOR(obs_list(i,3))
       thisobs%id_obs_p(2, i) = (FLOOR(obs_list(i,2)))*ny + FLOOR(obs_list(i,3))
       thisobs%id_obs_p(3, i) = thisobs%id_obs_p(1, i) + 1
       thisobs%id_obs_p(4, i) = thisobs%id_obs_p(2, i) + 1

    END DO

    DEALLOCATE(obs_list)


! **********************************************************************
! *** Initialize interpolation coefficients for observation operator ***
! **********************************************************************

    ! Allocate array of interpolation coefficients. As ID_OBS_P, the number
    ! of rows corresponds to the number of grid points using the the interpolation
    IF (ALLOCATED(thisobs%icoeff_p)) DEALLOCATE(thisobs%icoeff_p)
    ALLOCATE(thisobs%icoeff_p(4, dim_obs_p))

    DO i= 1, dim_obs_p
       ! Determine coordinates of grid points around observation

       ! Note: The computation of the coefficients assumes that the
       ! grid points 1 and 2 (likewise 3 and 4) differ only in their
       ! first coordinate, and grid points 1 and 3 (likewise 2 and 4)
       ! differ only in the second coordinate. The setup has to be 
       ! consistent with thisobs%id_obs_p initialized above.
       gcoords(1,1) = REAL(FLOOR(thisobs%ocoord_f(1, i)))
       gcoords(1,2) = REAL(FLOOR(thisobs%ocoord_f(2, i)))
       gcoords(2,1) = gcoords(1,1) + 1.0
       gcoords(2,2) = gcoords(1,2)
       gcoords(3,1) = gcoords(1,1)
       gcoords(3,2) = gcoords(1,2) + 1.0
       gcoords(4,1) = gcoords(1,1) + 1.0
       gcoords(4,2) = gcoords(1,2) + 1.0

       ! Compute interpolation coefficients
       CALL get_interp_coeff_lin(4, 2, gcoords, thisobs%ocoord_f(:, i), thisobs%icoeff_p(:,i))

    END DO


! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

    ! *** Set inverse observation error variances ***

    thisobs%ivar_obs_f(:) = 1.0 / (rms_obs_C*rms_obs_C)


! ****************************************
! *** Gather global observation arrays ***
! ****************************************

  ! Nothing to be done here as the model is serial


! *********************************************************
! *** For twin experiment: Read synthetic observations  ***
! *********************************************************

!     IF (twin_experiment .AND. filtertype/=11) THEN
!        CALL read_syn_obs(file_syntobs_TYPE, dim_obs_f, thisobs%obs_f, 0, 1-mype_filter)
!     END IF


! ********************
! *** Finishing up ***
! ********************

    ! store full and PE-local observation dimension in module variables
    thisobs%dim_obs_p = dim_obs_p
    thisobs%dim_obs_f = dim_obs_f

  END SUBROUTINE init_dim_obs_f_C



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
  SUBROUTINE obs_op_f_C(dim_p, dim_obs_f, state_p, obsstate_f, offset_obs)

    USE PDAFomi_obs_op, &
         ONLY: obs_op_f_interp_lin

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

    ! Store offset
    thisobs%off_obs_f = offset_obs

    ! observation operator for bi-linear interpolation
    CALL obs_op_f_interp_lin(dim_p, dim_obs_f, thisobs%dim_obs_p, thisobs%dim_obs_f, 4, &
         thisobs%id_obs_p, thisobs%icoeff_p, state_p, obsstate_f, offset_obs)

  END SUBROUTINE obs_op_f_C



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
  SUBROUTINE deallocate_obs_C()

    USE PDAFomi_obs_f, &
         ONLY: deallocate_obs

    IMPLICIT NONE

    ! Deallocate arrays in full observation type
    CALL deallocate_obs(thisobs)

  END SUBROUTINE deallocate_obs_C




!-------------------------------------------------------------------------------
!++++++      THE FOLLOWING ROUTINES SHOULD BE USABLE WITHOUT CHANGES       +++++
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
  SUBROUTINE init_obs_f_C(dim_obs_f, obsstate_f, offset_obs)

    USE PDAFomi_obs_f, &
         ONLY: init_obs_f

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_obs_f             !< Dimension of full observed state (all observed fields)
    REAL, INTENT(inout) :: obsstate_f(dim_obs_f) !< Full observation vector
    INTEGER, INTENT(inout) :: offset_obs         !< input: offset of module-type observations in obsstate_f
                                                 !< output: input + number of added observations


! ******************************************
! *** Initialize full observation vector ***
! ******************************************

    CALL init_obs_f(thisobs, dim_obs_f, obsstate_f, offset_obs)

  END SUBROUTINE init_obs_f_C



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
!! PE-local  sub-domain. (init_obsvar_l_TYPE is the 
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
  SUBROUTINE init_obsvar_C(meanvar, cnt_obs)

    USE PDAFomi_obs_f, &
         ONLY: init_obsvar_f

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: cnt_obs      !< Observation counter
    REAL, INTENT(inout) :: meanvar         !< Mean variance


! *****************************
! *** Compute mean variance ***
! *****************************

    CALL init_obsvar_f(thisobs, meanvar, cnt_obs)

  END SUBROUTINE init_obsvar_C



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
  SUBROUTINE init_dim_obs_l_C(coords_l, lradius, dim_obs_l, offset_obs_l, offset_obs_f)

  USE PDAFomi_obs_l, &
       ONLY: init_dim_obs_l

    IMPLICIT NONE

! *** Arguments ***
    REAL, INTENT(in) :: coords_l(:)        !< Coordinates of local analysis domain
    REAL, INTENT(in) :: lradius            !< Localization radius
    INTEGER, INTENT(out) :: dim_obs_l      !< Local number of module-type observations
    INTEGER, INTENT(inout) :: offset_obs_l !< input: Offset of module-type obs. in local obs. vector;
                                           !< output: input + number of added observations
    INTEGER, INTENT(inout) :: offset_obs_f !< input: Offset of module-type obs. in full obs. vector;
                                           !< output: input + number of added observations


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

  END SUBROUTINE init_dim_obs_l_C




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
  SUBROUTINE init_obs_l_C(dim_obs_l, obs_l)

  USE PDAFomi_obs_l, &
       ONLY: init_obs_l

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_obs_l         !< Local dimension of observation vector
    REAL, INTENT(inout) :: obs_l(dim_obs_l)  !< Local observation vector


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

    CALL init_obs_l(dim_obs_l, thisobs_l, thisobs, obs_l)

  END SUBROUTINE init_obs_l_C



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
  SUBROUTINE g2l_obs_C(dim_obs_l, dim_obs_f, obs_f, obs_l)

  USE PDAFomi_obs_l, &
       ONLY: g2l_obs

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_obs_l         !< Local dimension of observation vector
    INTEGER, INTENT(in) :: dim_obs_f         !< Full dimension of observation vector
    REAL, INTENT(in) :: obs_f(dim_obs_f)     !< Full observation vector
    REAL, INTENT(inout) :: obs_l(dim_obs_l)  !< Local observation vector


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

     CALL g2l_obs(dim_obs_l, thisobs_l%dim_obs_l, thisobs%dim_obs_f, thisobs_l%id_obs_l, &
         obs_f(thisobs%off_obs_f+1:thisobs%off_obs_f+thisobs%dim_obs_f), &
         thisobs_l%off_obs_l, obs_l)

  END SUBROUTINE g2l_obs_C



!-------------------------------------------------------------------------------
!> Compute product of inverse of R with some matrix
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
  SUBROUTINE prodRinvA_l_C(verbose, dim_obs_l, ncol, locweight, lradius, &
       sradius, A_l, C_l)

  USE PDAFomi_obs_l, &
       ONLY: prodRinvA_l

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
    INTEGER :: idummy   ! Dummy variable to present compiler warning

    idummy = dim_obs_l


! ***********************
! *** Compute product ***
! ***********************

    CALL prodRinvA_l(verbose, thisobs_l%dim_obs_l, ncol, locweight, lradius, sradius, &
         thisobs_l%ivar_obs_l, thisobs_l%distance_l, &
         A_l(thisobs_l%off_obs_l+1 : thisobs_l%off_obs_l+thisobs_l%dim_obs_l, :), &
         C_l(thisobs_l%off_obs_l+1 : thisobs_l%off_obs_l+thisobs_l%dim_obs_l, :))

  END SUBROUTINE prodRinvA_l_C



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
!! domain.  (See init_obsvar_TYPE for a global variant)
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
  SUBROUTINE init_obsvar_l_C(meanvar_l, cnt_obs_l)

    USE PDAFomi_obs_l, &
         ONLY: init_obsvar_l

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: cnt_obs_l      !< Observation counter
    REAL, INTENT(inout) :: meanvar_l         !< Mean variance


! ***********************************
! *** Compute local mean variance ***
! ***********************************

    CALL init_obsvar_l(thisobs_l, meanvar_l, cnt_obs_l)

  END SUBROUTINE init_obsvar_l_C

END MODULE obs_C_pdafomi
