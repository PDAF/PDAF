!$Id: init_pdaf.F90 2395 2020-10-06 16:46:42Z lnerger $
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
       COMM_filter_fesom, mype_filter_fesom, npes_filter_fesom, MPI_INTEGER
  USE mod_assim_pdaf, &           ! Variables for assimilation
       ONLY: dim_state, dim_state_p, dim_ens, dim_lag, &
       step_null, offset, screen, filtertype, subtype, &
       incremental, type_forget, forget, locweight, &
       type_trans, type_sqrt, eff_dim_obs, loctype, &
       twin_experiment, dim_obs_max, use_global_obs, DA_couple_type, &
       restart
  USE mod_assim_oce_pdaf, &       ! Variables for assimilation - oce-specific
       ONLY: delt_obs_ocn, delt_obs_ocn_offset
  USE obs_SST_CMEMS_pdafomi, &    ! Variables for SST observations
       ONLY: assim_o_sst, path_obs_sst, file_sst_prefix, file_sst_suffix, &
       rms_obs_sst, bias_obs_sst, lradius_sst, sradius_sst, loc_radius_sst, &
       sst_fixed_rmse, sst_exclude_diff, sst_exclude_ice, file_syntobs_sst 
  USE PDAFomi_obs_f, &            ! Routine to initialize domain limits
       ONLY: PDAFomi_get_domain_limits_unstr
  USE output_pdaf, &
       ONLY: write_da, write_ens, init_output_pdaf
  USE g_parfe, &
       ONLY: myDim_nod2D, myDim_nod3D, MPI_COMM_FESOM
  USE o_mesh, &
       ONLY: nod2D, nod3D, coord_nod2D
  USE g_clock, &
       ONLY: timeold
  USE g_rotate_grid, &
       ONLY: r2g
  use timer_pdaf, only: timeit

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: i                 ! Counter
  INTEGER :: filter_param_i(7) ! Integer parameter array for filter
  REAL    :: filter_param_r(2) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag
  INTEGER :: doexit, steps     ! Not used in this implementation
  REAL    :: timenow           ! Not used in this implementation
  REAL, ALLOCATABLE :: gcoords_p(:,:)   ! Array of geographic coordinates

! *** External subroutines *** 
  EXTERNAL :: init_ens_pdaf            ! Ensemble initialization
  EXTERNAL :: next_observation_pdaf, & ! Provide time step and model time of next observation
       distribute_state_pdaf, &        ! Routine to distribute a state vector to model fields
       distribute_state_restart_pdaf, &  ! Routine to distribute a state vector to model fields
       prepoststep_pdaf                ! User supplied pre/poststep routine


! ***************************
! ***   Initialize PDAF   ***
! ***************************

  ! Get process-ID in task of model compartment
  CALL MPI_Comm_Rank(MPI_COMM_FESOM, mype_submodel, MPIerr)

  IF (mype_submodel==0) THEN
!  IF (mype_model==0 .AND. task_id==1) THEN
!     WRITE (*,'(/1x,a)') 'FESOM-PDAF: INITIALIZE PDAF'
     WRITE (*,'(1x,a, i5)') 'FESOM-PDAF: INITIALIZE PDAF, task: ', task_id
  END IF


! ***********************************************************************
! ***   For weakly-coupled assimilation re-define filter communicator ***
! ***********************************************************************

  ! Create a communicator for filter processes for FESOM only
  CALL MPI_Comm_dup(MPI_COMM_FESOM, COMM_filter_fesom, MPIerr)
  CALL MPI_Comm_Size(COMM_filter_fesom, npes_filter_fesom, MPIerr)
  CALL MPI_Comm_Rank(COMM_filter_fesom, mype_filter_fesom, MPIerr)


  IF (DA_couple_type == 0) THEN

     ! Set filter communicator to the communicator of FESOM
     COMM_filter = COMM_filter_fesom !MPI_COMM_FESOM

     IF (filterpe) THEN
        CALL MPI_Comm_Size(COMM_filter, npes_filter, MPIerr)
        CALL MPI_Comm_Rank(COMM_filter, mype_filter, MPIerr)

        IF (mype_filter==0) THEN
           WRITE (*,'(a)') 'FESOM-PDAF: Initialize weakly-coupled data assimilation'
        ENDIF
     ENDIF
  ELSE
     IF (filterpe) THEN
        IF (mype_filter==0) THEN
           WRITE (*,'(a)') 'FESOM-PDAF: Initialize strongly-coupled data assimilation'
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
                    !   (11) GENOBS: Generate synthetic observations
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
  delt_obs_ocn = 480      ! Number of time steps between analysis/assimilation steps
  delt_obs_ocn_offset = 0 ! Offset of time steps until first analysis step

! *** Set assimilation cold start (initialize from covariance matrix)
  restart = .false.   

! *** Set weakly- or strongly-coupled DA
!  This is set in mod_assim_pdaf and can be changed in the namelist file
!  DA_couple_type = 0 ! (0) for weakly- (1) for strongly-coupled DA

! *** Set assimilation variables
  assim_o_sst   = .false.

! *** specifications for observations ***
  ! CMEMS SST
  rms_obs_sst = 0.05 ! error for satellite SST observations
  bias_obs_sst = 0.0    ! observation bias  
  sst_exclude_ice = .false.  ! Exclude SST observations at point with sea ice and T>0
  sst_exclude_diff = 0.0     ! Exclude SST observations if difference from ensemble mean is >sst_exclude_diff
  sst_fixed_rmse = .true.    ! Use the fixed RMS_obs_sst, or varying SST error provided with the data

  ! General
  use_global_obs = 1         ! Use global full obs. or full obs. limited to process domains
  twin_experiment = .false.  ! Whether to run a twin experiment assimilating synthetic observations
  dim_obs_max = 80000        ! Expected maximum number of observations for synthetic obs.

! *** Localization settings
  locweight = 0     ! Type of localizating weighting
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRADIUS
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
  lradius_sst = 1.0  ! Localization radius in meters for SST
  sradius_sst = lradius_sst    ! Support radius for 5th-order polynomial for SST

! *** File names - available as namelist read-in
  path_obs_sst = ''        ! Path to SST observation files
  file_sst_prefix = ''     ! Prefix of file holding SST observations
  file_sst_suffix = '.nc'  ! Suffix of file SST observations


! *** Read PDAF configuration from namelist ***

  CALL read_config_pdaf()


! ***************************
! *** Define state vector ***
! ***************************

! *** Initialize global and local state vector dimension             ***
! *** Fields: ssh, u, v, w, T, S, a_ice, m_ice, m_snow, u_ice, v_ice ***
! *** at the moment the updated variables include: ssh, u, v, w, S, T **

  dim_state   = 5 * nod3d + nod2d + 5 * nod2d
  dim_state_p = 5 * myDim_nod3d + myDim_nod2d + 5 * myDim_nod2d

! *** Specify offset of fields in state vector ***

  ALLOCATE(offset(11))

  offset(1)  = 0                           ! 1 SSH
  offset(2)  = myDim_nod2D                 ! 2 U
  offset(3)  = myDim_nod3D + myDim_nod2D   ! 3 V
  offset(4)  = 2*myDim_nod3D + myDim_nod2D ! 4 W
  offset(5)  = 3*myDim_nod3D + myDim_nod2D ! 5 Temperature
  offset(6)  = 4*myDim_nod3D + myDim_nod2D ! 6 Salinity
  offset(7)  = 5*myDim_nod3D + myDim_nod2D ! 7 aice
  offset(8)  = 5*myDim_nod3D + 2*myDim_nod2D ! 8 mice
  offset(9)  = 5*myDim_nod3D + 3*myDim_nod2D ! 9 msnow
  offset(10) = 5*myDim_nod3D + 4*myDim_nod2D ! 10 uice
  offset(11) = 5*myDim_nod3D + 5*myDim_nod2D ! 11 vice


! *** Initial Screen output ***

  IF (mype_model==0 .AND. task_id==1) CALL init_pdaf_info()

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
     IF (mype_filter==0) writepe = .TRUE.
  ENDIF
  IF (write_da) THEN
     ! Initialize Netcdf output
     CALL  init_output_pdaf(dim_lag, writepe)
  END IF


! ******************************'***
! *** Prepare ensemble forecasts ***
! ******************************'***

  IF (mype_submodel==0) THEN
     WRITE (*,'(1x,a,i5)') 'FESOM-PDAF: INITIALIZE PDAF before barrier, task: ', task_id
  END IF

  call timeit(6, 'new')
  CALL MPI_BARRIER(MPI_COMM_WORLD, MPIerr)
  call timeit(6, 'old')

  IF (restart) THEN
     CALL PDAF_get_state(steps, timenow, doexit, next_observation_pdaf, &
          distribute_state_restart_pdaf, prepoststep_pdaf, status_pdaf)
  ELSE
     CALL PDAF_get_state(steps, timenow, doexit, next_observation_pdaf, &
          distribute_state_pdaf, prepoststep_pdaf, status_pdaf)
  END IF


! ***********************************************************************************
! *** Allocate arrays for effective observation dimension and localization radius ***
! ***********************************************************************************
  
  ALLOCATE(eff_dim_obs(mydim_nod2d))
  ALLOCATE(loc_radius_sst(mydim_nod2d))


! ***************************************
! *** Get domain limiting coordinates ***
! ***************************************

  ! Allocate array for geographic coordinates
  ALLOCATE(gcoords_p(2,myDim_nod2d))

  DO i=1, myDim_nod2d
     ! Get geographic locations of grid points
     CALL r2g(gcoords_p(1,i), gcoords_p(2,i), coord_nod2d(1, i), coord_nod2d(2, i))
  END DO

  CALL PDAFomi_get_domain_limits_unstr(myDim_nod2d, gcoords_p)

  DEALLOCATE(gcoords_p)


END SUBROUTINE init_pdaf
