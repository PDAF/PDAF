!$Id$
!> PDAF-OMI observation module for type B observations
!!
!! This module handles operations for one data type (called 'module-type' below):
!! TYPE = B
!!
!! __Observation type B:__
!! The observation type B in this tutorial are the only the observations at
!! the locations (8,5), (12,15), and (4,30). 
!!
!! The subroutines in this module are for the particular handling of
!! a single observation type.
!! The routines are called by the different call-back routines of PDAF.
!! Most of the routines are generic so that in practice only 2 routines
!! need to be adapted for a particular data type. These are the routines
!! for the initialization of the observation information (init_dim_obs_f)
!! and for the observation operator (obs_op_f).
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
!! These 2 routines need to be adapted for the particular observation type:
!! * init_dim_obs_f_TYPE \n
!!           Count number of process-local and full observations; 
!!           initialize vector of observations and their inverse variances;
!!           initialize coordinate array and index array for indices of
!!           observed elements of the state vector.
!! * obs_op_f_TYPE \n
!!           observation operator to get full observation vector of this type. Here
!!           one has to choose a proper observation operator or implement one.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE obs_B_pdafomi

  USE mod_parallel_pdaf, &
       ONLY: mype_filter    ! Rank of filter process
  USE PDAFomi, &
       ONLY: obs_f, obs_l   ! Declaration of observation data types
 
  IMPLICIT NONE
  SAVE

  ! Variables which are inputs to the module (usually set in init_pdaf)
  LOGICAL :: assim_B        !< Whether to assimilate this data type
  REAL    :: rms_obs_B      !< Observation error standard deviation (for constant errors)

  ! One can declare further variables, e.g. for file names which can
  ! be use-included in init_pdaf() and initialized there.


! ***********************************************************************
! *** The following two data types are used in PDAFomi                ***
! *** They are declared in PDAFomi and only listed here for reference ***

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
! ***********************************************************************

! Declare instances of observation data types used here
! We use generic names here, but one could renamed the variables
  type(obs_f), target, public :: thisobs      ! full observation
  TYPE(obs_l), target, PUBLIC :: thisobs_l    ! local observation

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
  SUBROUTINE init_dim_obs_f_B(step, dim_obs_f)

    USE PDAFomi, &
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
    REAL, ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
    REAL, ALLOCATABLE :: ivar_obs_p(:)   ! PE-local inverse observation error variance
    REAL, ALLOCATABLE :: ocoord_p(:,:)   ! PE-local observation coordinates 
    CHARACTER(len=2) :: stepstr          ! String for time step
    REAL :: obs_tmp(3)                   ! Temporary storage of observation values


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    IF (mype_filter==0) &
         WRITE (*,'(8x,a)') 'Assimilate observations - obs type B'

    ! Store whether to assimilate this observation type (used in routines below)
    IF (assim_B) thisobs%doassim = 1

    ! Specify type of distance computation
    thisobs%disttype = 0   ! 0=Cartesian

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2


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

    ! Just keep observations at grid points (8,5), (12,15), (4,30)
    ! The other observations are used as observation type A
    obs_tmp(1) = obs_field(8,5)
    obs_tmp(2) = obs_field(12,15)
    obs_tmp(3) = obs_field(4,30)
    obs_field = -1000.0
    obs_field(8,5) = obs_tmp(1)
    obs_field(12,15) = obs_tmp(2)
    obs_field(4,30) = obs_tmp(3)


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

    ivar_obs_p(:) = 1.0 / (rms_obs_B*rms_obs_B)


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

  END SUBROUTINE init_dim_obs_f_B



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
  SUBROUTINE obs_op_f_B(dim_p, dim_obs_f, state_p, obsstate_f, offset_obs)

    USE PDAFomi, &
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
       CALL PDAFomi_obs_op_f_gridpoint(dim_p, dim_obs_f, thisobs%dim_obs_p, &
            thisobs%dim_obs_f, thisobs%id_obs_p, state_p, obsstate_f, offset_obs)
    END IF

  END SUBROUTINE obs_op_f_B

END MODULE obs_B_pdafomi
