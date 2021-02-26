!$Id$
!>  Interface routine to initialize parameters for PDAF and call PDAF_init
!!
!! This routine collects the initialization of variables for PDAF.
!! In addition, the initialization routine PDAF_init is called
!! such that the internal initialization of PDAF is performed.
!! This variant is for the online mode of PDAF.
!!
!! This routine is generic. However, it assumes a constant observation
!! error (rms_obs). Further, with parallelization the local state
!! dimension dim_state_p is used.
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
SUBROUTINE init_pdaf()

  USE mod_parallel_pdaf, &        ! Parallelization variables for assimilation
       ONLY: n_modeltasks, task_id, COMM_filter, COMM_couple, filterpe, &
       mype_world, COMM_model, abort_parallel, MPI_COMM_WORLD, MPIerr, &
       mype_model, mype_filter, npes_filter, writepe, mype_submodel, &
       COMM_filter_echam, mype_filter_echam, npes_filter_echam
  USE mod_assim_pdaf, & ! Variables for assimilation
       ONLY: dim_state, dim_state_p, dim_ens, dim_lag, &
       step_null, off_fields_p, screen, filtertype, subtype, &
       incremental, type_forget, forget, locweight, &
       type_trans, type_sqrt, eff_dim_obs, loctype, &
       DA_couple_type, restart, n_fields, dim_fields_p, dim_fields
  USE mod_assim_atm_pdaf, & ! Variables for assimilation
       ONLY: delt_obs_atm, delt_obs_atm_offset, dp
  USE obs_airt_pdafomi, &
       ONLY: assim_a_airt, rms_obs_airt, lradius_airt, sradius_airt
  USE mo_mpi, ONLY: p_global_comm
  USE timer_pdaf, ONLY: timeit
  USE mo_decomposition, ONLY: dc=>local_decomposition
  USE mo_time_control,  ONLY: get_time_step
  USE output_pdaf, ONLY: write_da, init_output_pdaf

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: i                 ! Counter
  INTEGER :: filter_param_i(7) ! Integer parameter array for filter
  REAL(dp):: filter_param_r(2) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag
  INTEGER :: doexit, steps     ! Not used in this implementation
  REAL(dp):: timenow           ! Not used in this implementation
  INTEGER :: dim_2d_p, dim_3d_p  ! Process-local field dimensions
  INTEGER :: dim_2d_g, dim_3d_g  ! Global field dimensions

! *** External subroutines *** 
  EXTERNAL :: init_ens_pdaf            ! Ensemble initialization
  EXTERNAL :: next_observation_pdaf, & ! Provide time step and model time of next observation
       distribute_state_pdaf, &        ! Routine to distribute a state vector to model fields
       distribute_state_ini_pdaf, &    ! Distribute a state vector to model fields at initial time
       prepoststep_pdaf                ! User supplied pre/poststep routine
  

! ***************************
! ***   Initialize PDAF   ***
! ***************************

  ! Get process-ID in task of model compartment
  CALL MPI_Comm_Rank(p_global_comm, mype_submodel, MPIerr)

  IF (mype_submodel==0) THEN
     WRITE (*,'(1x,a,i5)') 'ECHAM-PDAF: INITIALIZE PDAF, task: ', task_id
  END IF


! ***********************************************************************
! ***   For weakly-coupled assimilation re-define filter communicator ***
! ***********************************************************************

  ! Create a communicator for filter processes for ECHAM only
  CALL MPI_Comm_dup(p_global_comm, COMM_filter_echam, MPIerr)
  CALL MPI_Comm_Size(COMM_filter_echam, npes_filter_echam, MPIerr)
  CALL MPI_Comm_Rank(COMM_filter_echam, mype_filter_echam, MPIerr)

  IF (DA_couple_type == 0) THEN

     ! Set filter communicator to the communicator of ECHAM
     COMM_filter = COMM_filter_echam

     IF (filterpe) THEN
        CALL MPI_Comm_Size(COMM_filter, npes_filter, MPIerr)
        CALL MPI_Comm_Rank(COMM_filter, mype_filter, MPIerr)

        IF (mype_filter_echam==0) THEN
           WRITE (*,'(a)') 'ECHAM-PDAF: Initialize weakly-coupled data assimilation'
        ENDIF
     ENDIF
  ELSE
     IF (filterpe) THEN
        IF (mype_filter_echam==0) THEN
           WRITE (*,'(a)') 'ECHAM-PDAF: Initialize strongly-coupled data assimilation'
        END IF
     END IF
  END IF


! **********************************************************
! ***                  CONTROL OF PDAF                   ***
! ***              used in call to PDAF_init             ***
! **********************************************************

! *** IO options ***
  screen     = 2    ! Write screen output (1) for output, (2) add timings

! *** Filter specific variables
  filtertype = 7    ! Type of filter
                    !   (6) ESTKF
                    !   (7) LESTKF
  dim_ens = n_modeltasks ! Size of ensemble for all ensemble filters
                    ! Number of EOFs to be used for SEEK
  dim_lag = 0       ! Size of lag in smoother
  subtype = 0       ! subtype of filter: 
                    !   ESTKF:
                    !     (0) Standard form of ESTKF
                    !   LESTKF:
                    !     (0) Standard form of LESTKF
  type_trans = 0    ! Type of ensemble transformation
                    !   SEIK/LSEIK and ESTKF/LESTKF:
                    !     (0) use deterministic omega
                    !     (1) use random orthonormal omega orthogonal to (1,...,1)^T
                    !     (2) use product of (0) with random orthonormal matrix with
                    !         eigenvector (1,...,1)^T
                    !   ETKF/LETKF:
                    !     (0) use deterministic symmetric transformation
                    !     (2) use product of (0) with random orthonormal matrix with
                    !         eigenvector (1,...,1)^T
  type_forget = 0   ! Type of forgetting factor in SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
                    !   (0) fixed
                    !   (1) global adaptive
                    !   (2) local adaptive for LSEIK/LETKF/LESTKF
  forget  = 1.0     ! Forgetting factor
  type_sqrt = 0     ! Type of transform matrix square-root
                    !   (0) symmetric square root, (1) Cholesky decomposition
  incremental = 0   ! (1) to perform incremental updating (only in SEIK/LSEIK!)


! **********************************************************
! ***     CONTROL OF USER ROUTINES FOR ASSIMILATION      ***
! **********************************************************

! *** Forecast length (time interval between analysis steps) ***
  delt_obs_atm = 960      ! Number of time steps between analysis/assimilation steps
  delt_obs_atm_offset = 0 ! Offset of time steps until first analysis step

! *** Set assimilation cold start (initialize from covariance matrix)
  restart = .FALSE.   

! *** Set weakly- or strongly-coupled DA
!  This is set in mod_assim_pdaf and can be changed in the namelist file
!  DA_couple_type = 0 ! (0) for weakly- (1) for strongly-coupled DA

! *** Set assimilation variables
  assim_a_airt   = .FALSE.

! *** specifications for observations ***
  rms_obs_airt = 0.5

! *** Localization settings
  locweight = 0     ! Type of localizating weighting
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRADIUS
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
  lradius_airt = 1.0e6         ! Localization radius for air temperature
  sradius_airt = lradius_airt  ! Support radius for localization function

! *** File names
! - are defined in mod_assim_pdaf and modified by the namelist


! *** Read PDAF configuration from namelist ***

  CALL read_config_pdaf()

! Set step_null
  step_null = get_time_step()
  IF (mype_filter_echam==0) WRITE(*,*) 'step_null=',step_null


! ***************************
! *** Define state vector ***
! ***************************

! *** Initialize global and local state vector dimension             ***

  ! *** Define state dimension ***

  ! Process-local
  dim_2d_p = dc%nglat * dc%nglon             ! Size of 2D field
  dim_3d_p = dc%nglat * dc%nglon * dc%nlev   ! Size of 3D field 

  dim_state_p = 6 * dim_3d_p + dim_2d_p  ! Local state dimension

  ! Global
  dim_2d_g = dc%nlat * dc%nlon               ! Size of 2D field
  dim_3d_g = dc%nlat * dc%nlon * dc%nlev     ! Size of 3D field

  dim_state = 6 * dim_3d_g + dim_2d_g    ! Global state dimension


! *** Specify dimension and offset of fields in state vector ***
  
  n_fields = 7  ! Number of model fields in state vector

  ALLOCATE(dim_fields_p(n_fields))
  ALLOCATE(dim_fields(n_fields))
  ALLOCATE(off_fields_p(n_fields))

  ! Process-local field dimensions
  dim_fields_p(1) = dim_3d_p  ! 1 air temperature
  dim_fields_p(2) = dim_2d_p  ! 2 log surface pressure
  dim_fields_p(3) = dim_3d_p  ! 3 vorticity              
  dim_fields_p(4) = dim_3d_p  ! 4 divergence
  dim_fields_p(5) = dim_3d_p  ! 5 specific humidity
  dim_fields_p(6) = dim_3d_p  ! 6 u
  dim_fields_p(7) = dim_3d_p  ! 7 v

  ! Global field dimensions
  dim_fields(1) = dim_3d_g
  dim_fields(2) = dim_2d_g
  dim_fields(3:7) = dim_3d_g

  ! Offsets of fields in process-local state vector
  off_fields_p(1) = 0
  DO i = 2, n_fields
     off_fields_p(i) = off_fields_p(i-1) + dim_fields_p(i-1)
  END DO


! *** Initial Screen output ***

  IF (mype_submodel==0 .AND. task_id==1) CALL init_pdaf_info()

! *** Check ensemble size
  IF (dim_ens /= n_modeltasks) THEN
     WRITE (*,*) 'ERROR: Ensemble size (',dim_ens, &
          ') needs to be identical to the number of model tasks (',n_modeltasks,')'
     CALL abort_parallel()
  END IF


! *****************************************************
! *** Call PDAF initialization routine on all PEs.  ***
! ***                                               ***
! *** For all filters, first the arrays of integer  ***
! *** and real number parameters are initialized.   ***
! *** Subsequently, PDAF_init is called.            ***
! *****************************************************

  filter_param_i(1) = dim_state_p ! State dimension
  filter_param_i(2) = dim_ens     ! Size of ensemble
  filter_param_i(3) = 0           ! Smoother lag (not implemented here)
  filter_param_i(4) = incremental ! Whether to perform incremental analysis
  filter_param_i(5) = type_forget ! Type of forgetting factor
  filter_param_i(6) = type_trans  ! Type of ensemble transformation
  filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
  filter_param_r(1) = forget      ! Forgetting factor

  CALL PDAF_init(filtertype, subtype, 0, &
       filter_param_i, 7,&
       filter_param_r, 2, &
       COMM_model, COMM_filter, COMM_couple, &
       task_id, n_modeltasks, filterpe, init_ens_pdaf, &
       screen, status_pdaf)


! *** Check whether initialization of PDAF was successful ***
  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in initialization of PDAF - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF


! ******************************
! *** Initialize file output ***
! ******************************

 writepe = .FALSE.
 IF (filterpe) THEN
     IF (mype_filter_echam==0) writepe = .TRUE.
 ENDIF

  IF (write_da) THEN
     ! Initialize Netcdf output
     CALL init_output_pdaf(dim_lag, writepe)
  END IF


! **********************************
! *** Prepare ensemble forecasts ***
! **********************************

  IF (mype_submodel==0) THEN
     WRITE (*,'(1x,a,i5)') 'ECHAM-PDAF: INITIALIZE PDAF before barrier, task: ', task_id
  END IF

  CALL timeit(6, 'new')
  CALL MPI_BARRIER(MPI_COMM_WORLD, MPIerr)
  CALL timeit(6, 'old')

  CALL PDAF_get_state(steps, timenow, doexit, next_observation_pdaf, &
       distribute_state_ini_pdaf, prepoststep_pdaf, status_pdaf)


! **********************************************************
! *** Allocate array for effective observation dimension ***
! **********************************************************

  ALLOCATE(eff_dim_obs(dim_2d_g))

END SUBROUTINE init_pdaf
