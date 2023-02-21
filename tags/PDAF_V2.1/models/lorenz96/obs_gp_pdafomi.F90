!$Id: obs_gp_pdafomi.F90 488 2020-06-07 11:00:55Z lnerger $
!> PDAF-OMI observation module 
!!
!! This module handles operations for one data type (called 'module-type' below).
!!
!! __Observation type GP:__
!! The observations are located at model grid points. The state can be fully
!! or only partially observed.
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
!! These 2 routines need to be adapted for the particular observation type:
!! * init_dim_obs_TYPE \n
!!           Count number of process-local and full observations; 
!!           initialize vector of observations and their inverse variances;
!!           initialize coordinate array and index array for indices of
!!           observed elements of the state vector.
!! * obs_op_TYPE \n
!!           observation operator to get full observation vector of this type. Here
!!           one has to choose a proper observation operator or implement one.
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
MODULE obs_gp_pdafomi

  USE mod_parallel, &
       ONLY: mype_filter    ! Rank of filter process
  USE PDAFomi, &
       ONLY: obs_f, obs_l   ! Declaration of observation data types
 
  IMPLICIT NONE
  SAVE

  ! Variables which are inputs to the module (usually set in init_pdaf)
  LOGICAL :: assim_gp=.TRUE.        !< Whether to assimilate this data type
  REAL    :: rms_obs                !< Observation error standard deviation (for constant errors)

  ! One can declare further variables, e.g. for file names which can
  ! be use-included in init_pdaf() and initialized there.
  CHARACTER(len=110) :: file_obs  ! netcdf file holding observations
  CHARACTER(len=110) :: file_obs_mask  ! ASCII file holding observation mask
  CHARACTER(len=110) :: file_syntobs   ! netcdf file holding synthetic observations
  INTEGER :: obsfile_laststep(1)  ! Last time step in observation file
  INTEGER :: delt_obs_file     ! Observation interval in input file
  LOGICAL :: have_obs          ! Flag whether we consider observations
                               ! at next possible analysis time
  LOGICAL :: use_obs_mask  ! Whether to use a mask for observation gaps
  LOGICAL :: use_maskfile  ! Whether to read mask from a file
  INTEGER :: numobs        ! If not read from file, use this number of obs. (1 to numobs)
  INTEGER :: dx_obs        ! Of not read from file, use this grid point distance of obs.
  INTEGER :: obs_err_type  ! Type of observation error: (0) Gaussian, (1) double-exponential
  INTEGER, ALLOCATABLE :: obs_mask(:) ! Mask array for observation availability


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
  SUBROUTINE init_dim_obs_gp(step, dim_obs)

    USE netcdf
    USE PDAFomi, &
         ONLY: PDAFomi_gather_obs
    USE mod_assimilation, &
         ONLY: twin_experiment, filtertype, cradius
    USE mod_model, &
         ONLY: dim_state, step_null

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)    :: step       !< Current time step
    INTEGER, INTENT(inout) :: dim_obs    !< Dimension of full observation vector

! *** Local variables ***
    INTEGER :: i, j, s                   ! Counters
    INTEGER :: dim_obs_p                 ! Number of process-local observations
    INTEGER :: stat(50)                  ! Array for status flag
    INTEGER :: fileid                    ! Id of netcdf file
    INTEGER :: id_step                   ! ID for step information
    INTEGER :: id_obs                    ! File ID for observation
    INTEGER :: nsteps_file               ! Number of time steps in file
    INTEGER :: obs_step1and2(2)          ! Array holding first and second time step index
    INTEGER :: pos(2)                    ! Position index for writing
    INTEGER :: cnt(2)                    ! Count index for writing
    REAL, ALLOCATABLE :: obs_p(:)        ! PE-local observed SST field
    REAL, ALLOCATABLE :: ivar_obs_p(:)   ! PE-local inverse observation error variance
    REAL, ALLOCATABLE :: ocoord_p(:,:)   ! PE-local coordinates of observed SST field
    REAL, ALLOCATABLE :: obs_g(:)        ! Global observation vector
    INTEGER, ALLOCATABLE :: obsindx(:)   ! Index array for observations


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    IF (mype_filter==0) &
         WRITE (*,'(8x,a)') 'Assimilate observations - OBS_GP'

    ! Store whether to assimilate this observation type (used in routines below)
    IF (assim_gp) thisobs%doassim = 1

    ! Specify type of distance computation
    thisobs%disttype = 1   ! 1=Cartesian with periodicity

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 1

    ! Define size of domain for periodicity
    ALLOCATE(thisobs%domainsize(thisobs%ncoord))
    thisobs%domainsize(1) = REAL(dim_state)


! **********************************
! *** Read PE-local observations ***
! **********************************

    IF (.NOT.(filtertype==100 .OR. twin_experiment)) THEN
      ! *** If we don't generate observations with PDAF or run the twin experiment ***
      ! *** we use the observations generated using tools/generate_obs.F90         ***

      ! Retrieve information on observations from file
       s = 1
       stat(s) = NF90_OPEN(TRIM(file_obs), NF90_NOWRITE, fileid)

      ! Read number of time steps in file
       s = s + 1
       stat(s) = NF90_INQ_DIMID(fileid, 'timesteps', id_step)
       s = s + 1
       stat(s) = NF90_Inquire_dimension(fileid, id_step, len=nsteps_file)

      ! Read time step information
       s = s + 1
       stat(s) = NF90_INQ_VARID(fileid, 'step', id_step)

       pos(1) = 1
       cnt(1) = 2
       s = s + 1
       stat(s) = NF90_GET_VAR(fileid, id_step, obs_step1and2, start=pos(1:1), count=cnt(1:1))

       pos(1) = nsteps_file
       cnt(1) = 1
       s = s + 1
       stat(s) = NF90_GET_VAR(fileid, id_step, obsfile_laststep, start=pos(1:1), count=cnt(1:1))

       DO i = 1,  s
          IF (stat(i) /= NF90_NOERR) &
               WRITE(*, *) 'NetCDF error in reading observation step information, no.', i
       END DO

       ! Initialize observation interval in file
       delt_obs_file = obs_step1and2(2) - obs_step1and2(1)

       ! observation dimension 
       IF (step <= obsfile_laststep(1)) THEN
          dim_obs_p = dim_state
       ELSE
          dim_obs_p = 0
       END IF

       ! Set dim_obs_p to 0, if we are at the end of a forecast phase without obs.
       IF (.NOT.have_obs) dim_obs_p = 0

       ! Read global observation
       readobs: IF (dim_obs_p > 0) THEN

          ALLOCATE(obs_g(dim_state))

          s = 1
          stat(s) = NF90_INQ_VARID(fileid, 'obs', id_obs)

          WRITE (*,'(8x,a,i6)') &
               '--- Read observation at file position', step / delt_obs_file 

          pos(2) = step/delt_obs_file
          cnt(2) = 1
          pos(1) = 1
          cnt(1) = dim_obs_p
          s = s + 1
          stat(s) = NF90_GET_VAR(fileid, id_obs, obs_g(1:dim_obs_p), start=pos, count=cnt)

          s = s + 1
          stat(s) = NF90_CLOSE(fileid)

          DO i = 1,  s
             IF (stat(i) /= NF90_NOERR) &
                  WRITE(*, *) 'NetCDF error in reading global observation, no.', i
          END DO

       END IF readobs
    ELSE
       ! *** If we generate observations with PDAF or run the twin ***
       ! *** experiment we don't read the observation file         ***

       ALLOCATE(obs_g(dim_state))

       dim_obs_p = dim_state
       obs_g = 0.0
    END IF


    ! For gappy observations initialize index array
    ! and reorder global observation array

    obsgaps: IF (use_obs_mask) THEN

       ! Count observations
       s = 1
       DO i=1, dim_state
          IF (obs_mask(i) == 1) THEN
             s = s + 1
          END IF
       END DO
       dim_obs_p = s - 1

       ALLOCATE(obs_p(dim_obs_p))
       ALLOCATE(obsindx(dim_obs_p))
       obsindx = 0

       ! Initialize index vector and vector of observations
       s = 1
       DO i=1, dim_state
          IF (obs_mask(i) == 1) THEN
             obsindx(s) = i
             obs_p(s) = obs_g(i)
             s = s + 1
          END IF
       END DO

    ELSE

       ALLOCATE(obs_p(dim_obs_p))
       ALLOCATE(obsindx(dim_state))

       DO i=1, dim_state
          obsindx(i) = i
       END DO
       obs_p = obs_g

    END IF obsgaps


! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

    ! *** Count valid observations that lie within the process sub-domain ***

    ! Already set above

    ! *** Initialize vector of observations on the process sub-domain ***

    ! Not required here; it directly is obs_g(1:dim_obs_p)


    ! *** Initialize coordinate array of observations on the process sub-domain ***

    ALLOCATE(ocoord_p(1, dim_obs_p))
    
    ocoord_p(1,:) = REAL(obsindx(1:dim_obs_p))


    ! *** Initialize process local index array                                ***
    ! *** This array holds the information which elements of the state vector ***
    ! *** are used in the observation operation.                              ***

    ALLOCATE(thisobs%id_obs_p(1, dim_obs_p))

    thisobs%id_obs_p(1,:) = obsindx(1:dim_obs_p)


! **********************************************************************
! *** Initialize interpolation coefficients for observation operator ***
! **********************************************************************

    ! Nothing to be done here


! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

    ALLOCATE(ivar_obs_p(dim_obs_p))

    ivar_obs_p = 1 / (rms_obs*rms_obs)


! ****************************************
! *** Gather global observation arrays ***
! ****************************************

    CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
         thisobs%ncoord, cradius, dim_obs)


! *********************************************************
! *** For twin experiment: Read synthetic observations  ***
! *********************************************************

    IF (twin_experiment .AND. filtertype/=100) THEN
       CALL read_syn_obs(file_syntobs, dim_obs, thisobs%obs_f, step_null, 1-mype_filter)
    END IF


! ********************
! *** Finishing up ***
! ********************

    ! Deallocate all local arrays
    DEALLOCATE(obs_g, obs_p, obsindx, ivar_obs_p)

  END SUBROUTINE init_dim_obs_gp



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
  SUBROUTINE obs_op_gp(dim_p, dim_obs, state_p, ostate)

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
       ! Observation operator for observed grid point values
       CALL PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)
    END IF

  END SUBROUTINE obs_op_gp



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
  SUBROUTINE init_dim_obs_l_gp(domain_p, step, dim_obs, dim_obs_l)

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

  END SUBROUTINE init_dim_obs_l_gp



!-------------------------------------------------------------------------------
!> Perform covariance localization for local EnKF
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
  SUBROUTINE localize_covar_gp(dim_p, dim_obs, HP_p, HPH, coords_p)

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

  END SUBROUTINE localize_covar_gp

END MODULE obs_gp_pdafomi
