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
MODULE obs_C_pdafomi

  USE mod_parallel_pdaf, &
       ONLY: mype_filter    ! Rank of filter process
  USE PDAFomi, &
       ONLY: obs_f, obs_l   ! Declaration of observation data types
 
  IMPLICIT NONE
  SAVE

  ! Variables which are inputs to the module (usually set in init_pdaf)
  LOGICAL :: assim_C        !< Whether to assimilate this data type
  REAL    :: rms_obs_C      !< Observation error standard deviation (for constant errors)

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
  type(obs_f), public :: thisobs      ! full observation
  type(obs_l), public :: thisobs_l    ! local observation

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
  SUBROUTINE init_dim_obs_f_C(step, dim_obs_f)

    USE PDAFomi, &
         ONLY: PDAFomi_gather_obs_f, &
         PDAFomi_get_interp_coeff_lin
    USE mod_model, &
         ONLY: ny
    USE mod_assimilation, &
         ONLY: filtertype, local_range

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)    :: step       !< Current time step
    INTEGER, INTENT(inout) :: dim_obs_f  !< Dimension of full observation vector

! *** Local variables ***
    INTEGER :: i                         ! Counters
    INTEGER :: nobs                      ! Number of observations in file
    INTEGER :: dim_obs_p                 ! number of process-local observations
    REAL, ALLOCATABLE :: obs_list(:,:)   ! List of observations field read from file
    REAL, ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
    REAL, ALLOCATABLE :: ivar_obs_p(:)   ! PE-local inverse observation error variance
    REAL, ALLOCATABLE :: ocoord_p(:,:)   ! PE-local observation coordinates 
    CHARACTER(len=2) :: stepstr          ! String for time step
    REAL :: gcoords(4,2)                 ! Grid point coordinated to compute interpolation coeffs


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    IF (mype_filter==0) &
         WRITE (*,'(8x,a)') 'Assimilate observations - obs type C: interpolated observations'

    ! Store whether to assimilate this observation type (used in routines below)
    IF (assim_C) thisobs%doassim = 1

    ! Specify type of distance computation
    thisobs%disttype = 0   ! 0=Cartesian

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
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
    ALLOCATE(obs_p(dim_obs_p))
    ALLOCATE(ivar_obs_p(dim_obs_p))
    ALLOCATE(ocoord_p(2, dim_obs_p))

    ! Allocate process-local index array
    ! This array has a many rows as required for the observation operator
    ! 1 if observations are at grid points; >1 if interpolation is required
    ALLOCATE(thisobs%id_obs_p(4, dim_obs_p))

    DO i= 1, dim_obs_p

       ! Observation value and coordinates
       obs_p(i) = obs_list(i,1)
       ocoord_p(1, i) = obs_list(i,2)
       ocoord_p(2, i) = obs_list(i,3)
     
       ! State vector indices of 4 grid points in which box the observation lies
       ! Note: These indices have to be consistent with the coordinates used to
       ! compute the interpolation coefficients (see below)
       thisobs%id_obs_p(1, i) = (FLOOR(obs_list(i,2))-1)*ny + FLOOR(obs_list(i,3))
       thisobs%id_obs_p(2, i) = (FLOOR(obs_list(i,2)))*ny + FLOOR(obs_list(i,3))
       thisobs%id_obs_p(3, i) = thisobs%id_obs_p(1, i) + 1
       thisobs%id_obs_p(4, i) = thisobs%id_obs_p(2, i) + 1

    END DO


! **********************************************************************
! *** Initialize interpolation coefficients for observation operator ***
! **********************************************************************

    ! Allocate array of interpolation coefficients. As ID_OBS_P, the number
    ! of rows corresponds to the number of grid points using the the interpolation
    ALLOCATE(thisobs%icoeff_p(4, dim_obs_p))

    DO i= 1, dim_obs_p
       ! Determine coordinates of grid points around observation

       ! Note: The computation of the coefficients assumes that the
       ! grid points 1 and 2 (likewise 3 and 4) differ only in their
       ! first coordinate, and grid points 1 and 3 (likewise 2 and 4)
       ! differ only in the second coordinate. The setup has to be 
       ! consistent with thisobs%id_obs_p initialized above.
       gcoords(1,1) = REAL(FLOOR(ocoord_p(1, i)))
       gcoords(1,2) = REAL(FLOOR(ocoord_p(2, i)))
       gcoords(2,1) = gcoords(1,1) + 1.0
!        gcoords(2,2) = gcoords(1,2)
        gcoords(3,1) = gcoords(1,1)
!        gcoords(3,2) = gcoords(1,2) + 1.0
!        gcoords(4,1) = gcoords(1,1) + 1.0
!        gcoords(4,2) = gcoords(1,2) + 1.0

       ! Compute interpolation coefficients
       CALL PDAFomi_get_interp_coeff_lin(4, 2, gcoords, ocoord_p(:, i), thisobs%icoeff_p(:,i))

    END DO


! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

    ! *** Set inverse observation error variances ***

    ivar_obs_p(:) = 1.0 / (rms_obs_C*rms_obs_C)


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
    DEALLOCATE(obs_list)
    DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

    ! Arrays in THISOBS have to be deallocated after the analysis step
    ! by a call to deallocate_obs() in prepoststep_pdaf.

  END SUBROUTINE init_dim_obs_f_C



!-------------------------------------------------------------------------------
!> Implementation of observation operator 
!!
!! This routine applies the full observation operator
!! for the type of observations handled in this module
!! It has to append the observations to ostate_f from
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
  SUBROUTINE obs_op_f_C(dim_p, dim_obs_f, state_p, ostate_f, offset_obs)

    USE PDAFomi, &
         ONLY: PDAFomi_obs_op_f_interp_lin

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs_f             !< Dimension of full observed state (all observed fields)
    REAL, INTENT(in)    :: state_p(dim_p)        !< PE-local model state
    REAL, INTENT(inout) :: ostate_f(dim_obs_f)   !< Full observed state
    INTEGER, INTENT(inout) :: offset_obs         !< input: offset of module-type observations in ostate_f
                                                 !< output: input + number of added observations


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

    IF (thisobs%doassim==1) THEN
       ! observation operator for bi-linear interpolation
       CALL PDAFomi_obs_op_f_interp_lin(thisobs, 4, state_p, ostate_f, offset_obs)
    END IF

  END SUBROUTINE obs_op_f_C

END MODULE obs_C_pdafomi
