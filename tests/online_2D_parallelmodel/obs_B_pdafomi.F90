!> PDAF-OMI observation module for type B observations
!!
!! This module handles operations for one data type (called 'module-type' below):
!! OBSTYPE = B
!!
!! __Observation type B:__
!! The observation type B in this tutorial are 6 observations at specified 
!! model grid points.
!!
!! The subroutines in this module are for the particular handling of
!! a single observation type.
!! The routines are called by the different call-back routines of PDAF
!! usually by callback_obs_pdafomi.F90
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
!! In addition, there are two optional routines, which are required if filters 
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


! *********************************************************
! *** Data type obs_f defines the full observations by  ***
! *** internally shared variables of the module         ***
! *********************************************************

! Relevant variables that can be modified by the user:
!   TYPE obs_f
!      ---- Mandatory variables to be set in INIT_DIM_OBS ----
!      INTEGER :: doassim                    ! Whether to assimilate this observation type
!      INTEGER :: disttype                   ! Type of distance computation to use for localization
!                                            ! (0) Cartesian, (1) Cartesian periodic
!                                            ! (2) simplified geographic, (3) geographic haversine function
!      INTEGER :: ncoord                     ! Number of coordinates use for distance computation
!      INTEGER, ALLOCATABLE :: id_obs_p(:,:) ! Indices of observed field in state vector (process-local)
!
!      ---- Optional variables - they can be set in INIT_DIM_OBS ----
!      REAL, ALLOCATABLE :: icoeff_p(:,:)   ! Interpolation coefficients for obs. operator
!      REAL, ALLOCATABLE :: domainsize(:)   ! Size of domain for periodicity (<=0 for no periodicity)
!
!      ---- Variables with predefined values - they can be changed in INIT_DIM_OBS  ----
!      INTEGER :: obs_err_type=0            ! Type of observation error: (0) Gauss, (1) Laplace
!      INTEGER :: use_global_obs=1          ! Whether to use (1) global full obs. 
!                                           ! or (0) obs. restricted to those relevant for a process domain
!      REAL :: inno_omit=0.0                ! Omit obs. if squared innovation larger this factor times
!                                           !     observation variance
!      REAL :: inno_omit_ivar=1.0e-12       ! Value of inverse variance to omit observation
!   END TYPE obs_f

! Data type obs_l defines the local observations by internally shared variables of the module

! ***********************************************************************

! Declare instances of observation data types used here
! We use generic names here, but one could rename the variables
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
!! * thisobs\%id_obs_p    - array with indices of module-type observation in process-local state vector
!!
!! Optional is the use of
!! * thisobs\%icoeff_p    - Interpolation coefficients for obs. operator (only if interpolation is used)
!! * thisobs\%domainsize  - Size of domain for periodicity for disttype=1 (<0 for no periodicity)
!! * thisobs\%obs_err_type - Type of observation errors for particle filter and NETF (default: 0=Gaussian)
!! * thisobs\%use_global obs - Whether to use global observations or restrict the observations to the relevant ones
!!                          (default: 1=use global full observations)
!! * thisobs\%inno_omit   - Omit obs. if squared innovation larger this factor times observation variance
!!                          (default: 0.0, omission is active if >0) 
!! * thisobs\%inno_omit_ivar - Value of inverse variance to omit observation
!!                          (default: 1.0e-12, change this if this value is not small compared to actual obs. error)
!!
!! Further variables are set when the routine PDAFomi_gather_obs is called.
!!
  SUBROUTINE init_dim_obs_B(step, dim_obs)

    USE PDAF, &
         ONLY: PDAF_local_type
    USE PDAFomi, &
         ONLY: PDAFomi_gather_obs, PDAFomi_set_localize_covar
    USE mod_assimilation, &
         ONLY: filtertype, dim_state_p, locweight, cradius, sradius
    USE mod_model, &
         ONLY: nx, ny, nx_p, ndim, coords_p

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)    :: step       !< Current time step
    INTEGER, INTENT(inout) :: dim_obs    !< Dimension of full observation vector

! *** Local variables ***
    INTEGER :: i, j                      ! Counters
    INTEGER :: cnt_p, cnt0_p             ! Counters
    INTEGER :: off_nx                    ! Offset of local grid in global domain in x-direction
    INTEGER :: dim_obs_p                 ! Number of process-local observations
    INTEGER :: localtype                 ! Localization type index (2 or 3 for covariance localization)
    REAL, ALLOCATABLE :: obs_field(:,:)  ! Observation field read from file
    REAL, ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
    REAL, ALLOCATABLE :: ivar_obs_p(:)   ! PE-local inverse observation error variance
    REAL, ALLOCATABLE :: ocoord_p(:,:)   ! PE-local observation coordinates 
    CHARACTER(len=2) :: stepstr          ! String for time step


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

    IF (nx==36) THEN
       OPEN (12, file='../inputs_online.18x36/obsB_step'//TRIM(stepstr)//'.txt', status='old')
    ELSEIF (nx==256) THEN
       OPEN (12, file='../inputs_online.128x256/obsB_step'//TRIM(stepstr)//'.txt', status='old')
    ELSE
       OPEN (12, file='../inputs_online.512x2048obsB_step'//TRIM(stepstr)//'.txt', status='old')
    END IF
    DO i = 1, ny
       READ (12, *) obs_field(i, :)
    END DO
    CLOSE (12)


! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

    ! *** Count valid observations that lie within the process sub-domain ***

    ! Get offset of local domain in global domain in x-direction
    off_nx = 0
    DO i = 1, mype_filter
       off_nx = off_nx + nx_p
    END DO

    ! Count process-local observations
    cnt_p = 0
    DO j = 1 + off_nx, nx_p + off_nx
       DO i= 1, ny
          IF (obs_field(i,j) > -999.0) cnt_p = cnt_p + 1
       END DO
    END DO

    ! Set number of local observations
    dim_obs_p = cnt_p


    ! *** Initialize vector of observations on the process sub-domain ***
    ! *** Initialize coordinate array of observations on the process sub-domain ***

    haveobs: IF (dim_obs_p > 0) THEN

       ! Allocate process-local observation arrays
       ALLOCATE(obs_p(dim_obs_p))
       ALLOCATE(ivar_obs_p(dim_obs_p))
       ALLOCATE(ocoord_p(2, dim_obs_p))

       ! Allocate process-local index array
       ! This array has a many rows as required for the observation operator
       ! 1 if observations are at grid points; >1 if interpolation is required
       ALLOCATE(thisobs%id_obs_p(1, dim_obs_p))

       cnt_p = 0
       cnt0_p = 0
       DO j = 1 + off_nx, nx_p + off_nx
          DO i= 1, ny
             cnt0_p = cnt0_p + 1
             IF (obs_field(i,j) > -999.0) THEN
                cnt_p = cnt_p + 1
                thisobs%id_obs_p(1, cnt_p) = cnt0_p
                obs_p(cnt_p) = obs_field(i, j)
                ocoord_p(1, cnt_p) = REAL(j)
                ocoord_p(2, cnt_p) = REAL(i)
             END IF
          END DO
       END DO


! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

       ! *** Set inverse observation error variances ***

       ivar_obs_p(:) = 1.0 / (rms_obs_B*rms_obs_B)

    ELSE haveobs

       ! *** For dim_obs_p=0 allocate arrays with minimum size

       ALLOCATE(obs_p(1))
       ALLOCATE(ivar_obs_p(1))
       ALLOCATE(ocoord_p(2, 1))
       ALLOCATE(thisobs%id_obs_p(1, 1))

    END IF haveobs


! ****************************************
! *** Gather global observation arrays ***
! ****************************************

    CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
         thisobs%ncoord, cradius, dim_obs)


! ************************************************************
! *** Provide localization information for LEnKF and ENSRF ***
! ************************************************************

    IF (PDAF_local_type() > 1) THEN
       CALL PDAFomi_set_localize_covar(thisobs, dim_state_p, ndim, coords_p, &
            locweight, cradius, sradius)
    END IF


! *********************************************************
! *** For twin experiment: Read synthetic observations  ***
! *********************************************************

!     IF (twin_experiment .AND. filtertype/=100) THEN
!        CALL read_syn_obs(file_syntobs_OBSTYPE, dim_obs, thisobs%obs_f, 0, 1-mype_filter)
!     END IF


! ********************
! *** Finishing up ***
! ********************

    ! Deallocate all local arrays
    DEALLOCATE(obs_field)
    DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

    ! Arrays in THISOBS have to be deallocated after the analysis step
    ! by a call to deallocate_obs() in prepoststep_pdaf.

  END SUBROUTINE init_dim_obs_B



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
  SUBROUTINE obs_op_B(dim_p, dim_obs, state_p, ostate)

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

    ! observation operator for observed grid point values
    CALL PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)

  END SUBROUTINE obs_op_B



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
  SUBROUTINE init_dim_obs_l_B(domain_p, step, dim_obs, dim_obs_l)

    ! Include PDAFomi function
    USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l

    ! Include localization radius and local coordinates
    USE mod_assimilation, &   
         ONLY: coords_l, cradius, locweight, sradius

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)  :: domain_p     !< Index of current local analysis domain
    INTEGER, INTENT(in)  :: step         !< Current time step
    INTEGER, INTENT(in)  :: dim_obs      !< Full dimension of observation vector
    INTEGER, INTENT(inout) :: dim_obs_l  !< Local dimension of observation vector


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
         locweight, cradius, sradius, dim_obs_l)

  END SUBROUTINE init_dim_obs_l_B



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
  SUBROUTINE localize_covar_B(dim_p, dim_obs, HP_p, HPH, coords_p)

    ! Include PDAFomi function
    USE PDAFomi, ONLY: PDAFomi_localize_covar

    ! Include localization radius and local coordinates
    USE mod_assimilation, &   
         ONLY: cradius, locweight, sradius

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< Dimension of observation vector
    REAL, INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
    REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH
    REAL, INTENT(in)    :: coords_p(:,:)         !< Coordinates of state vector elements


! *************************************
! *** Apply covariance localization ***
! *************************************

    CALL PDAFomi_localize_covar(thisobs, dim_p, locweight, cradius, sradius, &
         coords_p, HP_p, HPH)

  END SUBROUTINE localize_covar_B



!-------------------------------------------------------------------------------
!> Perform covariance localization for ENSRF on the module-type observation
!!
!! The routine is called in the analysis step of the ENSRF
!! with serial observation processing. It has to apply localization
!! to the two vectors HP and HXY for the single observation with
!! index iobs
!!
!! This routine calls the routine PDAFomi_localize_covar_serial
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type.
!!
  SUBROUTINE localize_covar_serial_B(iobs, dim_p, dim_obs, HP_p, HXY_p, coords_p)

    ! Include PDAFomi function
    USE PDAFomi, ONLY: PDAFomi_localize_covar_serial

    ! Include localization radius and local coordinates
    USE mod_assimilation, &   
         ONLY: cradius, locweight, sradius

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: iobs                  !< Index of current observation
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< Dimension of observation vector
    REAL, INTENT(inout) :: HP_p(dim_p)           !< Process-local part of matrix HP for observation iobs
    REAL, INTENT(inout) :: HXY_p(dim_obs)        !< Process-local part of matrix HX(HX_all) for full observations
    REAL, INTENT(in)    :: coords_p(:,:)         !< Coordinates of state vector elements


! *************************************
! *** Apply covariance localization ***
! *************************************

   CALL PDAFomi_localize_covar_serial(thisobs, iobs, dim_p, dim_obs, locweight, &
        cradius, sradius, coords_p, HP_p, HXY_p)

  END SUBROUTINE localize_covar_serial_B

END MODULE obs_B_pdafomi
