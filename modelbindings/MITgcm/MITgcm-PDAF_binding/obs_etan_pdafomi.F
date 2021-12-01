C     $Id: obs_etan_pdafomi.F90 561 2020-11-21 09:59:46Z lnerger $
C     !> PDAF-OMI observation module for MITgcm

C     ! This module handles operations for one data type (called 'module-type' below).
C     ! Here this type are simulated measurements of sea surface height (etan in MITgcm).
C     !  
C     ! The subroutines in this module are for the particular handling of
C     ! a single observation type.
C     ! The routines are called by the different call-back routines of PDAF.
C     ! Most of the routines are generic so that in practice only 2 routines
C     ! need to be adapted for a particular data type. These are the routines
C     ! for the initialization of the observation information (init_dim_obs)
C     ! and for the observation operator (obs_op).
C     !
C     ! The module and the routines are named according to the observation type.
C     ! This allows to distinguish the observation type and the routines in this
C     ! module from other observation types.
C     !
C     ! The module uses two derived data types (obs_f and obs_l), which contain
C     ! all information about the full and local observations. Only variables
C     ! of the type obs_f need to be initialized in this module. The variables
C     ! in the type obs_l are initilized by the generic routines from PDAFomi.
C     !
C     !
C     ! These 2 routines need to be adapted for the particular observation type:
C     ! * init_dim_obs_TYPE \n
C     !           Count number of process-local and full observations;
C     !           initialize vector of observations and their inverse variances;
C     !           initialize coordinate array and index array for indices of
C     !           observed elements of the state vector.
C     ! * obs_op_TYPE \n
C     !           observation operator to get full observation vector of this type. Here
C     !           one has to choose a proper observation operator or implement one.
C     !
C     ! In addition, there are two optional routine, which are required if filters
C     ! with localization are used:
C     ! * init_dim_obs_l_TYPE \n
C     !           Only required if domain-localized filters (e.g. LESTKF, LETKF) are used:
C     !           Count number of local observations of module-type according to
C     !           their coordinates (distance from local analysis domain). Initialize
C     !           module-internal distances and index arrays.
C     ! * localize_covar_TYPE \n
C     !           Only required if the localized EnKF is used:
C     !           Apply covariance localization in the LEnKF.
C     !
C     ! __Revision history:__
C     ! * 2021-03 - Nabir Mamnun - Initial code
C     ! * Later revisions - see repository log

#include "CPP_OPTIONS.h"

      MODULE obs_etan_pdafomi

      USE mod_parallel_pdaf, ONLY: mype_filter    ! Rank of filter process
      USE PDAFomi, ONLY: obs_f, obs_l   ! Declaration of observation data types

      IMPLICIT NONE
      SAVE

C     ! Variables which are inputs to the module (usually set in init_pdaf)
      LOGICAL       :: assim_etan        !< Whether to assimilate this data type
      REAL(kind=8)  :: rms_obs_etan      !< Observation error standard deviation (for constant errors)

C     ! One can declare further variables, e.g. for file names which can
C     ! be use-included in init_pdaf() and initialized there.

C     !***********************************************************************
C     !*** The following two data types are used in PDAFomi                ***
C     !*** They are declared in PDAFomi and only listed here for reference ***
C     !***********************************************************************

C     !Data type to define the full observations by internally shared variables of the module
C     !  TYPE obs_f
C     !          Mandatory variables to be set in INIT_DIM_OBS
C     !     INTEGER :: doassim                  ! Whether to assimilate this observation type
C     !     INTEGER :: disttype                  ! Type of distance computation to use for localization
C     !                                          ! (0) Cartesian, (1) Cartesian periodic
C     !                                          ! (2) simplified geographic, (3) geographic haversine function
C     !     INTEGER :: ncoord                    ! Number of coordinates use for distance computation
C     !     INTEGER, ALLOCATABLE :: id_obs_p(:,:) ! Indices of observed field in state vector (process-local)

C     !          Optional variables - they can be set in INIT_DIM_OBS
C     !     REAL(kind=8), ALLOCATABLE :: icoeff_p(:,:)   ! Interpolation coefficients for obs. operator
C     !     REAL(kind=8), ALLOCATABLE :: domainsize(:)   ! Size of domain for periodicity (<=0 for no periodicity)

C     !          Variables with predefined values - they can be changed in INIT_DIM_OBS
C     !     INTEGER :: obs_err_type=0            ! Type of observation error: (0) Gauss, (1) Laplace
C     !     INTEGER :: use_global_obs=1          ! Whether to use (1) global full obs.
C     !                                          ! or (0) obs. restricted to those relevant for a process domain

C     !          The following variables are set in the routine PDAFomi_gather_obs
C     !     INTEGER :: dim_obs_p                 ! number of PE-local observations
C     !     INTEGER :: dim_obs_f                 ! number of full observations
C     !     INTEGER :: dim_obs_g                 ! global number of observations
C     !     INTEGER :: off_obs_f                 ! Offset of this observation in overall full obs. vector
C     !     INTEGER :: off_obs_g                 ! Offset of this observation in overall global obs. vector
C     !     INTEGER :: obsid                     ! Index of observation over all assimilated observations
C     !     REAL(kind=8), ALLOCATABLE :: obs_f(:)        ! Full observed field
C     !     REAL(kind=8), ALLOCATABLE :: ocoord_f(:,:)   ! Coordinates of full observation vector
C     !     REAL(kind=8), ALLOCATABLE :: ivar_obs_f(:)   ! Inverse variance of full observations
C     !     INTEGER, ALLOCATABLE :: id_obs_f_lim(:) ! Indices of domain-relevant full obs. in global vector of obs.
C     !                                          ! (only if full obs. are restricted to process domain))
C     !  END TYPE obs_f

C     !Data type to define the local observations by internally shared variables of the module
C     !  TYPE obs_l
C     !     INTEGER :: dim_obs_l                 ! number of local observations
C     !     INTEGER :: off_obs_l                 ! Offset of this observation in overall local obs. vector
C     !     INTEGER, ALLOCATABLE :: id_obs_l(:)  ! Indices of local observations in full obs. vector
C     !     REAL(kind=8), ALLOCATABLE :: distance_l(:)   ! Distances of local observations
C     !     REAL(kind=8), ALLOCATABLE :: ivar_obs_l(:)   ! Inverse variance of local observations
C     !     INTEGER :: locweight                 ! Specify localization function
C     !     REAL(kind=8) :: lradius                      ! localization radius
C     !     REAL(kind=8) :: sradius                      ! support radius for localization function
C     ! END TYPE obs_l
C     !***********************************************************************

C     !Declare instances of observation data types used here
C     !We use generic names here, but one could renamed the variables
      TYPE(obs_f), TARGET, PUBLIC :: thisobs      ! full observation
      TYPE(obs_l), TARGET, PUBLIC :: thisobs_l    ! local observation

C$OMP THREADPRIVATE(thisobs_l)

C     !-------------------------------------------------------------------------------

      CONTAINS

C     ! > Initialize information on the module-type observation
C     !
C     ! The routine is called by each filter process.
C     ! at the beginning of the analysis step before
C     ! the loop through all local analysis domains.

C     ! It has to count the number of observations of the
C     ! observation type handled in this module according
C     ! to the current time step for all observations
C     ! required for the analyses in the loop over all local
C     ! analysis domains on the PE-local state domain.

C     ! The following four variables have to be initialized in this routine
C     ! * thisobs\%doassim     - Whether to assimilate this type of observations
C     ! * thisobs\%disttype    - type of distance computation for localization with this observaton
C     ! * thisobs\%ncoord      - number of coordinates used for distance computation
C     ! * thisobs\%id_obs_p    - index of module-type observation in PE-local state vector

C     ! Optional is the use of
C     ! * thisobs\%icoeff_p    - Interpolation coefficients for obs. operator (only if interpolation is used)
C     ! * thisobs\%domainsize  - Size of domain for periodicity for disttype=1 (<0 for no periodicity)
C     ! * thisobs\%obs_err_type - Type of observation errors for particle filter and NETF (default: 0=Gaussian)
C     ! * thisobs\%use_global obs - Whether to use global observations or restrict the observations to the relevant ones
C     !                          (default: 1=use global full observations)

C     ! Further variables are set when the routine PDAFomi_gather_obs is called.
C     !
      SUBROUTINE init_dim_obs_etan(step, dim_obs)


      USE PDAFomi, ONLY: PDAFomi_gather_obs
      USE mod_assimilation,
     &     ONLY: filtertype, local_range, id, off_fields

      IMPLICIT NONE

#include "SIZE.h"
#include "GRID.h"

C     !*** Arguments ***
      INTEGER, INTENT(in)    :: step       !< Current time step
      INTEGER, INTENT(inout) :: dim_obs    !< Dimension of full observation vector

C     !*** Local variables ***
      INTEGER :: i, j                              ! Counters
      INTEGER :: cnt, cnt0                         ! Counters
      INTEGER :: dim_obs_p                         ! Number of process-local observations
      REAL(kind=8), ALLOCATABLE :: obs_field(:,:)  ! Observation field read from file
      REAL(kind=8), ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
      REAL(kind=8), ALLOCATABLE :: ivar_obs_p(:)   ! PE-local inverse observation error variance
      REAL(kind=8), ALLOCATABLE :: ocoord_p(:,:)   ! PE-local observation coordinates


C     !*********************************************
C     !*** Initialize full observation dimension ***
C     !*********************************************

      IF (mype_filter==0)
     &  WRITE (*,'(8x,a)') 'Assimilate observations - obs type: etan'

C     ! Store whether to assimilate this observation type (used in routines below)
      IF (assim_etan) thisobs%doassim = 1

C     ! Specify type of distance computation
      thisobs%disttype = 0   ! 0=Cartesian

C     ! Number of coordinates used for distance computation
C     ! The distance compution starts from the first row
      thisobs%ncoord = 2


C     !**********************************
C     !*** Read PE-local observations ***
C     !**********************************

C     ! Read observation field from file
      ALLOCATE(obs_field(sNy, sNx))

      ! Here, we just set a value
      obs_field = 0.5d0


C     !***********************************************************
C     !*** Count available observations for the process domain ***
C     !*** and initialize index and coordinate arrays.         ***
C     !***********************************************************

C     ! *** Count valid observations that lie within the process sub-domain ***

      cnt = 0
      DO j = 1, sNx
        DO i= 1, sNy
          IF (obs_field(i,j) > -999.0) cnt = cnt + 1
        END DO
      END DO
      dim_obs_p = cnt
      dim_obs = cnt

      IF (mype_filter==0) WRITE (*,'(8x, a, i6)')
     &      '--- number of full observations', dim_obs


C     ! *** Initialize vector of observations on the process sub-domain ***
C     ! *** Initialize coordinate array of observations on the process sub-domain ***

C     ! Allocate process-local observation arrays
      ALLOCATE(obs_p(dim_obs_p))
      ALLOCATE(ivar_obs_p(dim_obs_p))
      ALLOCATE(ocoord_p(2, dim_obs_p))

C     ! Allocate process-local index array
C     ! This array has a many rows as required for the observation operator
C     ! 1 if observations are at grid points; >1 if interpolation is required
      ALLOCATE(thisobs%id_obs_p(1, dim_obs_p))

        cnt = 0
        cnt0 = 0
        DO j = 1, sNx
          DO i= 1, sNy
            cnt0 = cnt0 + 1
            IF (obs_field(i,j) > -999.0) THEN
              cnt = cnt + 1
              thisobs%id_obs_p(1, cnt) = cnt0 + off_fields(id%etan)
              obs_p(cnt) = obs_field(i, j)
              ocoord_p(1, cnt) = xC(i, j, 1, 1)
              ocoord_p(2, cnt) = yC(i, j, 1, 1)
            END IF
          END DO
        END DO


C     !****************************************************************
C     !*** Define observation errors for process-local observations ***
C     !****************************************************************

C    ! *** Set inverse observation error variances ***

      ivar_obs_p(:) = 1.0 / (rms_obs_etan*rms_obs_etan)


C     !****************************************
C     !*** Gather global observation arrays ***
C     !****************************************

      CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p,
     &     ocoord_p, thisobs%ncoord, local_range, dim_obs)


C     !*********************************************************
C     !*** For twin experiment: Read synthetic observations  ***
C     !*********************************************************

C     !    IF (twin_experiment .AND. filtertype/=100) THEN
C     !       CALL read_syn_obs(file_syntobs_TYPE, dim_obs, thisobs%obs_f, 0, 1-mype_filter)
C     !    END IF


C     !********************
C     !*** Finishing up ***
C     !********************

C     ! Deallocate all local arrays
      DEALLOCATE(obs_field)
      DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

      END



C     -------------------------------------------------------------------------------
C     > Implementation of observation operator
C     !
C     ! This routine applies the full observation operator
C     ! for the type of observations handled in this module.
C     !
C     ! One can choose a proper observation operator from
C     ! PDAFOMI_OBS_OP or add one to that module or
C     ! implement another observation operator here.
C     !
C     ! The routine is called by all filter processes.
C     !
      SUBROUTINE obs_op_etan(dim_p, dim_obs, state_p, ostate)

      USE PDAFomi, ONLY: PDAFomi_obs_op_gridpoint

      IMPLICIT NONE

C     !*** Arguments ***
      INTEGER, INTENT(in) :: dim_p     !< PE-local state dimension
      INTEGER, INTENT(in) :: dim_obs   !< Dimension of full observed state (all observed fields)
      REAL(kind=8), INTENT(in)    :: state_p(dim_p)   !< PE-local model state
      REAL(kind=8), INTENT(inout) :: ostate(dim_obs)  !< Full observed state

C     !******************************************************
C     !*** Apply observation operator H on a state vector ***
C     !******************************************************

      IF (thisobs%doassim==1) THEN
C       ! observation operator for observed grid point values
        CALL PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)
      END IF

      END


C     !-------------------------------------------------------------------------------
C     !> Initialize local information on the module-type observation
C     !
C     ! The routine is called during the loop over all local
C     ! analysis domains. It has to initialize the information
C     ! about local observations of the module type. It returns
C     ! number of local observations of the module type for the
C     ! current local analysis domain in DIM_OBS_L and the full
C     ! and local offsets of the observation in the overall
C     ! observation vector.
C     !
C     ! This routine calls the routine PDAFomi_init_dim_obs_l
C     ! for each observation type. The call allows to specify a
C     ! different localization radius and localization functions
C     ! for each observation type and  local analysis domain.
C     !
      SUBROUTINE init_dim_obs_l_etan(domain_p, step, dim_obs,
     &    dim_obs_l)

C     ! Include PDAFomi function
      USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l
C    ! Include localization radius and local coordinates
      USE mod_assimilation, ONLY: coords_l, local_range,
     &    locweight, srange

      IMPLICIT NONE

C     !*** Arguments ***
      INTEGER, INTENT(in)  :: domain_p     !< Index of current local analysis domain
      INTEGER, INTENT(in)  :: step         !< Current time step
      INTEGER, INTENT(in)  :: dim_obs      !< Full dimension of observation vector
      INTEGER, INTENT(inout) :: dim_obs_l  !< Local dimension of observation vector


C     !**********************************************
C     !*** Initialize local observation dimension ***
C     !**********************************************

      CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l,
     &      locweight, local_range, srange, dim_obs_l)

      END


C     !-------------------------------------------------------------------------------
C     !> Perform covariance localization for local EnKF on the module-type observation
C     !
C     ! The routine is called in the analysis step of the localized
C     ! EnKF. It has to apply localization to the two matrices
C     ! HP and HPH of the analysis step for the module-type
C     ! observation.
C     !
C     ! This routine calls the routine PDAFomi_localize_covar
C     ! for each observation type. The call allows to specify a
C     ! different localization radius and localization functions
C     ! for each observation type.
C     !
      SUBROUTINE localize_covar_etan(dim_p,dim_obs,HP_p,HPH,coords_p)

C     ! Include PDAFomi function
      USE PDAFomi, ONLY: PDAFomi_localize_covar
C     ! Include localization radius and local coordinates
      USE mod_assimilation, ONLY: local_range, locweight, srange

      IMPLICIT NONE

C     ! *** Arguments ***
      INTEGER, INTENT(in) :: dim_p                         !< PE-local state dimension
      INTEGER, INTENT(in) :: dim_obs                       !< Dimension of observation vector
      REAL(kind=8), INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
      REAL(kind=8), INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH
      REAL(kind=8), INTENT(in)    :: coords_p(:,:)         !< Coordinates of state vector elements


C     ! *************************************
C     ! *** Apply covariance localization ***
C     ! *************************************

      CALL PDAFomi_localize_covar(thisobs, dim_p, locweight,
     &    local_range, srange, coords_p, HP_p, HPH)


      END

      END MODULE obs_etan_pdafomi
