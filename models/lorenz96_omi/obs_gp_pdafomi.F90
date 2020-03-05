!$Id: obs_gp_pdafomi.F90 346 2020-01-22 08:33:33Z lnerger $
!> PDAF-OMI template observation module 
!!
!! This module handles operations for one data type (called 'module-type' below).
!! 
!! The different subroutines in this module are for the particular handling of
!! a single observation type.
!! The routines are called by the different call-back routines of PDAF or the
!! interface routines in interface_pdafomi. Most of the routines are generic,
!! so that in practice only 2 routines need to be adapted for a particular
!! data type. There are the routines for the initialization of the observation
!! information (init_dim_obs_f) and for the observation operator (obs_op_f).
!!
!! **Using this template:**
!! To be able to distinguish the observation type and the routines in this module,
!! we recommend to rename the module according to the observation module-type.
!! Further,we recommend to replace 'TYPE' in the routine names according to the
!! type of the observation so that they can be identified when 
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
!! 
!! For the particle filter and NETF, and for the localizated NETF these routines
!! are used
!! * likelihood_TYPE \n
!!        Compute global likelihood
!! * likelihood_l_TYPE \n
!!        Compute local likelihood
!!
!! For generating synthetic observations, there are the routines
!! are used
!! * init_obserr_TYPE \n
!!        Initialize vector of observation error standard deviations
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE obs_gp_pdafomi

  USE mod_parallel, &
       ONLY: mype_filter    ! Rank of filter process
  USE PDAFomi_obs_f, &
       ONLY: obs_f          ! Declaration of data type obs_f
  USE PDAFomi_obs_l, &
       ONLY: obs_l          ! Declaration of data type obs_l
 
  IMPLICIT NONE
  SAVE

  ! Variables which are inputs to the module (usually set in init_pdaf)
  LOGICAL :: assim_gp=.true.        !< Whether to assimilate this data type
  REAL    :: rms_obs                !< Observation error standard deviation (for constant errors)

  ! One can declare further variables, e.g. for file names which can
  ! be use-included in init_pdaf() and initialized there.
  CHARACTER(len=110) :: file_obs  ! netcdf file holding observations
  CHARACTER(len=110) :: file_obs_mask  ! ASCII file holding observation mask
  CHARACTER(len=110) :: file_syntobs   ! netcdf file holding synthetic observations
  INTEGER :: obsfile_laststep  ! Last time step in observation file
  INTEGER :: delt_obs_file     ! Observation interval in input file
  LOGICAL :: have_obs          ! Flag whether we consider observations
                               ! at next possible analysis time
  LOGICAL :: use_obs_mask  ! Whether to use a mask for observation gaps
  LOGICAL :: use_maskfile  ! Whether to read mask from a file
  INTEGER :: numobs        ! If not read from file, use this number of obs. (1 to numobs)
  INTEGER :: dx_obs        ! Of not read from file, use this grid point distance of obs.
  INTEGER :: obs_err_type  ! Type of observation error: (0) Gaussian, (1) double-exponential
  INTEGER, ALLOCATABLE :: obs_mask(:) ! Mask array for observation availability


! *** The following two data types are used inside this module. ***
! *** They are declared in PDAFomi_obs_f and PDAFomi_obs_l and  ***
! *** only listed here for reference

! Data type to define the full observations by internally shared variables of the module
!   type obs_f
!           Mandatory variables to be set in init_dim_obs_f
!      INTEGER :: doassim                   ! Whether to assimilate this observation type
!      INTEGER :: dim_obs_p                 ! number of PE-local observations
!      INTEGER :: dim_obs_f                 ! number of full observations
!      INTEGER, ALLOCATABLE :: id_obs_p(:,:) ! indices of observed field in state vector
!      REAL, ALLOCATABLE :: obs_f(:)        ! Full observed field
!      REAL, ALLOCATABLE :: ocoord_f(:,:)   ! Coordinates of full observation vector
!      REAL, ALLOCATABLE :: ivar_obs_f(:)   ! Inverse variance of full observations
!      INTEGER :: disttype                  ! Type of distance computation to use for localization
!      INTEGER :: ncoord                    ! Number of coordinates use for distance computation
!           Variable with predefined values - they can be changed in init_dim_obs_f 
!      LOGICAL :: localfilter=.true.        ! Whether a localized filter is used
!      INTEGER :: obs_err_type=0            ! Type of observation error: (0) Gauss, (1) Laplace
!           Optional variables - they can be set in init_dim_obs_f
!      REAL, ALLOCATABLE :: icoeff_p(:,:)   ! Interpolation coefficients for obs. operator
!      REAL, ALLOCATABLE :: domainsize(:)   ! Size of domain for periodicity (<=0 for no periodicity)
!           Mandatory variable to be set in obs_op_f
!      INTEGER :: off_obs_f                 ! Offset of this observation in overall full obs. vector
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
!! * thisobs\%doassim     - Whether to assimilate this observation type
!! * thisobs\%dim_obs_p   - PE-local number of module-type observations
!! * thisobs\%dim_obs_f   - full number of module-type observations
!! * thisobs\%id_obs_p    - index of module-type observation in PE-local state vector
!! * thisobs\%obs_f       - full vector of module-type observations
!! * thisobs\%ocoord_f    - coordinates of observations in OBS_MOD_F
!! * thisobs\%ivar_obs_f  - full vector of inverse obs. error variances of module-type
!! * thisobs\%disttype    - type of distance computation for localization with this observaton
!! * thisobs\%ncoord      - number of coordinates used for distance computation
!!
!! Optional is the use of
!! * thisobs\%icoeff_p    - Interpolation coefficients for obs. operator (only if interpolation is used)
!! * thisobs\%localfilter - Whether we use a localized filter 
!!                          (default: .true.; only relevant if the model uses domain decomposition)
!! * thisobs\%type_obserr - Type of observation error: (0) Gauss, (1) Laplace
!!                          (default: 0; only relevant for the NETF/LNETF and PF)
!! * thisobs\%domainsize  - Size of domain in each dimension 
!!                          (allocate and use only for Cartesian grid with periodicity)
!!
!! **Adapting the template**
!! In this routine the variables listed above have to be initialized. One
!! can include modules from the model with 'use', e.g. for mesh information.
!! Alternatively one could include these as subroutine arguments
!!
  SUBROUTINE init_dim_obs_f_gp(step, dim_obs_f)

    USE mod_assimilation, &
         ONLY: twin_experiment, filtertype
    USE mod_model, &
         ONLY: dim_state, step_null

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

! *** Arguments ***
    INTEGER, INTENT(in)    :: step       !< Current time step
    INTEGER, INTENT(inout) :: dim_obs_f  !< Dimension of full observation vector

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
    REAL, ALLOCATABLE :: observation_g(:) ! For local filter: global observation vector
    INTEGER, ALLOCATABLE :: obsindx(:)   ! Index array for observations


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    IF (mype_filter==0) &
         WRITE (*,'(8x,a)') 'Assimilate observations - OBS_GP'

    ! Store whether to assimilate this observation type
    IF (assim_gp) thisobs%doassim = 1

    ! Specify type of distance computation
    thisobs%disttype = 1   ! 1=Cartesian with periodicity

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 1

    ! Define size of domain for periodicity
    allocate(thisobs%domainsize(thisobs%ncoord))
    thisobs%domainsize(1) = REAL(dim_state)

    ! Whether we use a local filter (default: .true.) 
    IF (filtertype==6) thisobs%localfilter = .false.  ! (for ESTKF)


! **********************************
! *** Read PE-local observations ***
! **********************************

    IF (.NOT.(filtertype==11 .OR. twin_experiment)) THEN
      ! *** If we don't generate observations with PDAF or run the twin experiment ***
      ! *** we use the observations generated using tools/generate_obs.F90         ***

      ! Retrieve information on observations from file
       s = 1
       stat(s) = NF_OPEN(TRIM(file_obs), NF_NOWRITE, fileid)

      ! Read number of time steps in file
       s = s + 1
       stat(s) = NF_INQ_DIMID(fileid, 'timesteps', id_step)
       s = s + 1
       stat(s) = NF_INQ_DIMLEN(fileid, id_step, nsteps_file)

      ! Read time step information
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, 'step', id_step)

       pos(1) = 1
       cnt(1) = 2
       s = s + 1
       stat(s) = NF_GET_VARA_INT(fileid, id_step, pos, cnt, obs_step1and2)

       pos(1) = nsteps_file
       cnt(1) = 1
       s = s + 1
       stat(s) = NF_GET_VARA_INT(fileid, id_step, pos, cnt, obsfile_laststep)

       DO i = 1,  s
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in reading observation step information, no.', i
       END DO

       ! Initialize observation interval in file
       delt_obs_file = obs_step1and2(2) - obs_step1and2(1)

       ! observation dimension 
       IF (step <= obsfile_laststep) THEN
          dim_obs_p = dim_state
       ELSE
          dim_obs_p = 0
       END IF

       ! Set dim_obs_p to 0, if we are at the end of a forecast phase without obs.
       IF (.NOT.have_obs) dim_obs_p = 0

       ! Read global observation
       readobs: IF (dim_obs_p > 0) THEN

          ALLOCATE(observation_g(dim_state))

          s = 1
          stat(s) = NF_INQ_VARID(fileid, 'obs', id_obs)

          write (*,'(8x,a,i6)') &
               '--- Read observation at file position', step / delt_obs_file

          pos(2) = step/delt_obs_file
          cnt(2) = 1
          pos(1) = 1
          cnt(1) = dim_obs_p
          s = s + 1
          stat(s) = NF_GET_VARA_DOUBLE(fileid, id_obs, pos, cnt, observation_g(1:dim_obs_p))

          s = s + 1
          stat(s) = nf_close(fileid)

          DO i = 1,  s
             IF (stat(i) /= NF_NOERR) &
                  WRITE(*, *) 'NetCDF error in reading global observation, no.', i
          END DO

       END IF readobs
    ELSE
       ! *** If we generate observations with PDAF or run the twin ***
       ! *** experiment we don't read the observation file         ***

       ALLOCATE(observation_g(dim_state))

       dim_obs_p = dim_state
       observation_g = 0.0
    END IF


    ! For gappy observations initialize index array
    ! and reorder global observation array
    ALLOCATE(obsindx(dim_state))

    obsgaps: IF (use_obs_mask) THEN

       obsindx = 0

       s = 1
       DO i=1, dim_state
          IF (obs_mask(i) == 1) THEN
             obsindx(s) = i
             observation_g(s) = observation_g(i)
             s = s + 1
          END IF
       END DO
       dim_obs_p = s - 1

    ELSE

       DO i=1, dim_state
          obsindx(i) = i
       END DO

    END IF obsgaps


! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

    ! *** Count valid observations that lie within the process sub-domain ***

    ! Already set above

    ! *** Initialize vector of observations on the process sub-domain ***

    ! Not required here; it directly is observation_g(1:dim_obs_p)


    ! *** Initialize coordinate array of observations on the process sub-domain ***

    ! Not required here; it directly is REAL(obsindx(1:dim_obs_p))


    ! *** Initialize process local index array                                ***
    ! *** This array holds the information which elements of the state vector ***
    ! *** are used in the observation operation.                              ***

    IF (ALLOCATED(thisobs%id_obs_p)) DEALLOCATE(thisobs%id_obs_p)
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

    dim_obs_f = dim_obs_p

    ! *** Gather full observation vector and corresponding coordinates ***

    ! Allocate full observation arrays
    ! The arrays are deallocated in deallocate_obs in this module
    ALLOCATE(thisobs%obs_f(dim_obs_f))
    ALLOCATE(thisobs%ivar_obs_f(dim_obs_f))
    ALLOCATE(thisobs%ocoord_f(1, dim_obs_f))

    thisobs%obs_f = observation_g(1:dim_obs_f)
    thisobs%ivar_obs_f = ivar_obs_p
    thisobs%ocoord_f(1,:) = REAL(obsindx(1:dim_obs_p))


! *********************************************************
! *** For twin experiment: Read synthetic observations  ***
! *********************************************************

    IF (twin_experiment .AND. filtertype/=11) THEN
       CALL read_syn_obs(file_syntobs, dim_obs_f, thisobs%obs_f, 0, 1-mype_filter)
    END IF


! ********************
! *** Finishing up ***
! ********************

    ! store full and PE-local observation dimension in module variables
    thisobs%dim_obs_p = dim_obs_p
    thisobs%dim_obs_f = dim_obs_f

    ! Clean up arrays
    IF (ALLOCATED(observation_g)) DEALLOCATE(observation_g)
    DEALLOCATE(obsindx, ivar_obs_p)

  END SUBROUTINE init_dim_obs_f_gp



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
  SUBROUTINE obs_op_f_gp(dim_p, dim_obs_f, state_p, obsstate_f, offset_obs)

    USE PDAFomi_obs_op, &
         ONLY: obs_op_f_gridpoint

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
    
       !+++  Choose suitable observation operator from the
       !+++  module PDAFomi_obs_op or implement your own

       ! observation operator for observed grid point values
       CALL obs_op_f_gridpoint(thisobs%localfilter, dim_p, dim_obs_f, thisobs%dim_obs_p, &
            thisobs%dim_obs_f, thisobs%id_obs_p, state_p, obsstate_f, offset_obs)
    END IF

  END SUBROUTINE obs_op_f_gp


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
  SUBROUTINE deallocate_obs_gp()

    USE PDAFomi_obs_f, &
         ONLY: deallocate_obs

    IMPLICIT NONE

    ! Deallocate arrays in full observation type
    CALL deallocate_obs(thisobs)

  END SUBROUTINE deallocate_obs_gp




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
  SUBROUTINE init_obs_f_gp(dim_obs_f, obsstate_f, offset_obs)

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

    IF (thisobs%doassim == 1) THEN
       CALL init_obs_f(thisobs, dim_obs_f, obsstate_f, offset_obs)
    END IF

  END SUBROUTINE init_obs_f_gp



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
!! PE-local  sub-domain. (init_obsvar_l_gp is the 
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
  SUBROUTINE init_obsvar_gp(meanvar, cnt_obs)

    USE PDAFomi_obs_f, &
         ONLY: init_obsvar_f

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: cnt_obs      !< Observation counter
    REAL, INTENT(inout) :: meanvar         !< Mean variance


! *****************************
! *** Compute mean variance ***
! *****************************

    IF (thisobs%doassim == 1) THEN
       CALL init_obsvar_f(thisobs, meanvar, cnt_obs)
    END IF

  END SUBROUTINE init_obsvar_gp



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
  SUBROUTINE init_dim_obs_l_gp(coords_l, lradius, dim_obs_l, offset_obs_l, offset_obs_f)

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

    IF (thisobs%doassim == 1) THEN

       IF (offset_obs_f /= thisobs%off_obs_f) THEN
          WRITE (*,*) 'ERROR: INCONSISTENT ORDER of observation calls in OBS_OP_F and INIT_DIM_OBS_L!'
       END IF


! ********************************************************************
! *** Initialize local observation dimension and local obs. arrays ***
! ********************************************************************

       CALL init_dim_obs_l(thisobs, thisobs_l, coords_l, lradius, dim_obs_l, &
            offset_obs_l, offset_obs_f)
    ELSE
       dim_obs_l = 0
    END IF

  END SUBROUTINE init_dim_obs_l_gp




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
  SUBROUTINE init_obs_l_gp(dim_obs_l, obs_l)

  USE PDAFomi_obs_l, &
       ONLY: init_obs_l

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_obs_l         !< Local dimension of observation vector
    REAL, INTENT(inout) :: obs_l(dim_obs_l)  !< Local observation vector


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

    IF (thisobs%doassim == 1) THEN
       CALL init_obs_l(dim_obs_l, thisobs_l, thisobs, obs_l)
    END IF

  END SUBROUTINE init_obs_l_gp



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
  SUBROUTINE g2l_obs_gp(dim_obs_l, dim_obs_f, obs_f, obs_l)

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

    IF (thisobs%doassim == 1) THEN
       CALL g2l_obs(dim_obs_l, thisobs_l%dim_obs_l, thisobs%dim_obs_f, thisobs_l%id_obs_l, &
            obs_f(thisobs%off_obs_f+1:thisobs%off_obs_f+thisobs%dim_obs_f), &
            thisobs_l%off_obs_l, obs_l)
    END IF

  END SUBROUTINE g2l_obs_gp



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
  SUBROUTINE prodRinvA_l_gp(verbose, dim_obs_l, ncol, locweight, lradius, &
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
    INTEGER :: idummy   ! Dummy variable to prevent compiler warning

    idummy = dim_obs_l


! ***********************
! *** Compute product ***
! ***********************

    IF (thisobs%doassim == 1) THEN
       CALL prodRinvA_l(verbose, thisobs_l%dim_obs_l, ncol, locweight, lradius, sradius, &
            thisobs_l%ivar_obs_l, thisobs_l%distance_l, &
            A_l(thisobs_l%off_obs_l+1 : thisobs_l%off_obs_l+thisobs_l%dim_obs_l, :), &
            C_l(thisobs_l%off_obs_l+1 : thisobs_l%off_obs_l+thisobs_l%dim_obs_l, :))
    END IF

  END SUBROUTINE prodRinvA_l_gp



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
!! domain.  (See init_obsvar_gp for a global variant)
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
  SUBROUTINE init_obsvar_l_gp(meanvar_l, cnt_obs_l)

    USE PDAFomi_obs_l, &
         ONLY: init_obsvar_l

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: cnt_obs_l      !< Observation counter
    REAL, INTENT(inout) :: meanvar_l         !< Mean variance


! ***********************************
! *** Compute local mean variance ***
! ***********************************

    IF (thisobs%doassim == 1) THEN
       CALL init_obsvar_l(thisobs_l, meanvar_l, cnt_obs_l)
    END IF

  END SUBROUTINE init_obsvar_l_gp



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
  SUBROUTINE prodRinvA_gp(ncol, A_p, C_p)

    USE PDAFomi_obs_f, &
         ONLY: prodRinvA

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: ncol        !< Rank of initial covariance matrix
    REAL, INTENT(in   ) :: A_p(:, :)   !< Input matrix
    REAL, INTENT(out)   :: C_p(:, :)   !< Output matrix


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


! ***********************
! *** Compute product ***
! ***********************

       CALL prodRinvA(thisobs%dim_obs_f, ncol, thisobs%ivar_obs_f, &
            A_p(thisobs%off_obs_f+1 : thisobs%off_obs_f+thisobs%dim_obs_f, :), &
            C_p(thisobs%off_obs_f+1 : thisobs%off_obs_f+thisobs%dim_obs_f, :))
    END IF

  END SUBROUTINE prodRinvA_gp



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
  SUBROUTINE add_obs_error_gp(dim_obs, C_p)

    USE PDAFomi_obs_f, &
         ONLY: add_obs_error

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

       CALL add_obs_error(thisobs%dim_obs_f, thisobs%ivar_obs_f, C_p, thisobs%off_obs_f)

    END IF

  END SUBROUTINE add_obs_error_gp



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
  SUBROUTINE init_obscovar_gp(dim_obs, covar, isdiag)

    USE PDAFomi_obs_f, &
         ONLY: init_obscovar

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_obs      !< Number of observations
    REAL, INTENT(out) :: covar(:, :)    !< covariance matrix R
    LOGICAL, INTENT(out) :: isdiag      !< Whether matrix R is diagonal


! ***********************************************
! *** Add observation error covariance matrix ***
! ***********************************************

    IF (thisobs%doassim == 1) THEN
       CALL init_obscovar(thisobs%dim_obs_f, thisobs%ivar_obs_f, thisobs%off_obs_f, &
            covar, isdiag)
    END IF

  END SUBROUTINE init_obscovar_gp



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
  SUBROUTINE likelihood_gp(dim_obs, obs, resid, lhood)

    USE PDAFomi_obs_f, &
         ONLY: likelihood

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
       CALL likelihood(dim_obs, obs, resid, thisobs%ivar_obs_f, lhood, thisobs%obs_err_type)
    END IF

  END SUBROUTINE likelihood_gp



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
  SUBROUTINE likelihood_l_gp(verbose, dim_obs_l, obs_l, resid_l, &
       locweight, lradius, sradius, lhood_l)

    USE PDAFomi_obs_l, &
         ONLY: likelihood_l

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
       CALL likelihood_l(verbose, dim_obs_l, obs_l, resid_l, locweight, lradius,  &
            sradius, thisobs_l%ivar_obs_l, thisobs_l%distance_l, lhood_l, thisobs%obs_err_type)
    END IF

  END SUBROUTINE likelihood_l_gp

END MODULE obs_gp_pdafomi
