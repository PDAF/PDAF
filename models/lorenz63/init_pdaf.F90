!$Id$
!BOP
!
! !ROUTINE: init_pdaf - Interface routine to call initialization of PDAF
!
! !INTERFACE:
SUBROUTINE init_pdaf()

! !DESCRIPTION:
! This routine collects the initialization of variables
! for PDAF as well as the call to the initialization
! routine PDAF_init.
!
! This variant is for the Lorenz63 model.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF, &   ! Interface definitions to PDAF core routines
       ONLY: PDAF_init, PDAF_get_state
  USE parser, &
       ONLY: parse
  USE mod_model, &
       ONLY: step_null, dim_state, dt
  USE mod_modeltime, &
       ONLY: total_steps
  USE mod_parallel, &
       ONLY: mype_world, n_modeltasks, task_id, &
       COMM_model, COMM_filter, COMM_couple, filterpe, abort_parallel
  USE mod_assimilation, &
       ONLY: screen, filtertype, subtype, dim_ens, delt_obs, &
       model_error, model_err_amp, incremental, covartype, &
       type_forget, forget, rank_analysis_enkf, &
       file_ini, type_ensinit, seedset, type_trans, &
       type_sqrt, stepnull_means, dim_lag, time, &
       twin_experiment, pf_res_type, init_dt, init_maxtime, pf_noise_type, &
       pf_noise_amp, type_winf, limit_winf
  USE output_netcdf_asml, &
       ONLY: init_netcdf_asml, file_asml, delt_write_asml, write_states, &
       write_stats, write_ens
  USE obs_gp_pdafomi, &
       ONLY: rms_obs, file_obs, use_obs_mask, file_obs_mask, &
       use_maskfile, numobs, dx_obs, obs_err_type, file_syntobs

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: main
! Calls: PDAF_init
! Calls: parse
!EOP

! Local variables
  INTEGER :: filter_param_i(7) ! Integer parameter array for filter
  REAL    :: filter_param_r(3) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag

  ! External subroutines
  EXTERNAL :: init_ens_pdaf    ! Routine for ensemble initialization
  

! ***************************
! ***   Initialize PDAF   ***
! ***************************

  IF (mype_world == 0) THEN
     WRITE (*,'(/1x,a)') 'INITIALIZE PDAF'
  END IF


! **********************************************************
! ***      INITIALIZATION OF PARAMETERS FOR PDAF         ***
! ***             used in call to PDAF_init              ***
! **********************************************************

! *** IO options ***
  screen      = 2   ! Write screen output (1) for output, (2) add timings

! *** Filter specific variables
  filtertype = 6    ! Type of filter
                    !   SEIK (1), EnKF (2), ETKF (4), ESTKF (6), 
                    !   NETF (9), PF (12), GENOBS (100)
  dim_ens = 20      ! Size of ensemble
  dim_lag = 0       ! Size of lag in smoother
  subtype = 0       ! subtype of filter: 
                    !   SEIK:
                    !     (0) mean forecast; new formulation
                    !     (1) mean forecast; old formulation
                    !     (2) fixed error space basis
                    !     (3) fixed state covariance matrix
                    !     (4) SEIK with ensemble transformation
                    !   EnKF:
                    !     (0) analysis for large observation dimension
                    !     (1) analysis for small observation dimension
                    !   ETKF:
                    !     (0) ETKF using T-matrix like SEIK
                    !     (1) ETKF following Hunt et al. (2007)
                    !       There are no fixed basis/covariance cases, as
                    !       these are equivalent to SEIK subtypes 2/3
                    !   ESTKF:
                    !     (0) Standard form of ESTKF
                    !     (2) fixed ensemble perturbations
                    !     (3) fixed state covariance matrix
                    !   NETF:
                    !     (0) Standard form of NETF
                    !   PF:
                    !     (0) Standard form of PF
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
                    !   NETF/LNETF:
                    !     (0) use random orthonormal transformation orthogonal to (1,...,1)^T
                    !     (1) use identity transformation
  incremental = 0   ! (1) to perform incremental updating (only in SEIK/LSEIK!)
  forget  = 1.0     ! Forgetting factor
  type_forget = 0   ! Type of forgetting factor in SEIK/LSEIK
                    ! (0): fixed; (1) global adaptive; (2) local adaptive for LSEIK
  type_sqrt = 0     ! Type of transform matrix square-root
                    !   (0) symmetric square root, (1) Cholesky decomposition
  type_winf = 0     ! NETF/LNETF: Type of weights inflation: (1) use N_eff/N>limit_winf
  limit_winf = 0.0  ! Limit for weights inflation
  rank_analysis_enkf = 0   ! ENKF only: rank to be considered for inversion of HPH
                    ! in analysis step of EnKF; (0) for analysis w/o eigendecomposition
  model_error = .false. ! Whether to apply model error noise
  model_err_amp = 0.1   ! Amplitude of model noise (times dt for error variance)
  pf_res_type = 1   ! Resampling type for particle filter
                    !   (1) probabilistic resampling
                    !   (2) stochastic universal resampling
                    !   (3) residual resampling
  pf_noise_type = 0 ! Type of pertubing noise in PF: (0) no perturbations
                    ! (1) constant stddev, (2) amplitude of stddev relative of ensemble variance
  pf_noise_amp = 0.0 ! Noise amplitude for particle filter


! **********************************************************
! ***         INITIALIZATION OF VARIABLES USED           ***
! ***       IN USER-SUPPLIED (CALL-BACK) ROUTINES        ***
! **********************************************************

! *** Whether to run twin experimnet assimilating synthetic observations ***
  twin_experiment = .false.

! *** IO options ***
  delt_write_asml = 1    ! Output interval for state information in assimilation intervals
  write_states = .TRUE.  ! Whether to write estimates states into the output file
  write_stats  = .TRUE.  ! Whether to write time dependent ensemble statistics (skewness, kurtosis)
  write_ens  = .FALSE.   ! Whether to write full time dependent ensemble stats
  stepnull_means = 1001  ! Step at which the second computation of time mean error is started
                         ! (first computation of mean sis always starting at initial step)

! *** specifications for observations ***
  ! avg. observation error (used for assimilation)
  rms_obs = 2.0      ! This error is the standard deviation 
                     ! for the Gaussian distribution 
  delt_obs = 10      ! Time step interval between analysis/assimilation steps
  use_obs_mask = .FALSE. ! Whether to use observations with gaps
  use_maskfile = .FALSE. ! If a mask is used read it from file
  numobs = 1         ! If not read from file use this number of obs. (1 to numobs)
  dx_obs = 3         ! grid point distance of observations (if not read from file)
  obs_err_type = 0   ! Observation errors: (0) for Gaussian (1) for double-exponential

! *** Filter specific variables ***
  type_ensinit = 'rnd' ! 'eof' for 2nd-order exact sampling from EOFs
                    !    'rnd' for random sampling from true state trajectory
                    !    'ens' for reading an initial ensemble generated from a previous assimilation run
  seedset = 1       ! Index of set of seeds to be used for init (only for 'rnd')
  init_dt = 10      ! Time step interval considered for 'rnd' initialization
  init_maxtime = 5000 ! Maximum time step to pick from for random ensemble initialization
  covartype = 1     ! Definition of factor in covar. matrix used in SEIK
                    ! (0) for (r+1)^-1 (old SEIK); (1): for r^-1 (real ensemble
                    ! covariance matrix) This parameter has also to be set internally
                    ! in PDAF_init

! *** File names ***
  file_asml = 'assim_l63.nc'        ! Output file
  file_ini = 'covar_l63.nc'         ! Initialization file
  file_obs = 'obs_l63.nc'           ! File holding observations
  file_obs_mask = 'obsmask_l63.txt' ! File holding observation mask
  file_syntobs = 'syntobs_l63.nc'   ! File holding synthetic observations generated by PDAF


! ********************************
! ***      INITIALIZATION      ***
! ********************************

! *** Parse command line options   ***
! *** This is optional, but useful ***

  CALL init_pdaf_parse()


! *** Initial Screen output ***
  screen2: IF (mype_world == 0) THEN

     IF (filtertype == 1) THEN
        WRITE (*, '(21x, a)') 'Filter: SEIK'
        IF (subtype == 2) THEN
           WRITE (*, '(6x, a)') '-- fixed error-space basis'
        ELSE IF (subtype == 3) THEN
           WRITE (*, '(6x, a)') '-- fixed state covariance matrix'
        ELSE IF (subtype == 4) THEN
           WRITE (*, '(6x, a)') '-- use ensemble transformation'
        END IF
        WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
        WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
        WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
        IF (model_error) THEN
           WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
        END IF
     ELSE IF (filtertype == 2) THEN
        WRITE (*, '(21x, a)') 'Filter: EnKF'
        WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
        IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
        WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
        IF (model_error) THEN
           WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
        END IF
        IF (rank_analysis_enkf > 0) THEN
           WRITE (*, '(6x, a, i5)') &
                'analysis with pseudo-inverse of HPH, rank:', rank_analysis_enkf
        END IF
     ELSE IF (filtertype == 4) THEN
        WRITE (*, '(21x, a)') 'Filter: ETKF'
        IF (subtype == 0) THEN
           WRITE (*, '(17x, a)') '--> Variant using T-matrix'
        ELSE IF (subtype == 1) THEN
           WRITE (*, '(17x, a)') '--> Variant following Hunt et al. (2007)'
        END IF
        WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
        IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
        WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
        WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
        IF (model_error) THEN
           WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
        END IF
     ELSE IF (filtertype == 6) THEN
        WRITE (*, '(21x, a)') 'Filter: ESTKF'
        WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
        IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
        WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
        WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
        IF (model_error) THEN
           WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
        END IF
     ELSE IF (filtertype == 9) THEN
        WRITE (*, '(21x, a)') 'Filter: NETF'
        WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
        IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
        WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
        WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
        IF (model_error) THEN
           WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
        END IF
        IF (obs_err_type==0) THEN
           WRITE (*, '(6x, a)') 'Gaussian observation errors'
        ELSE IF (obs_err_type==1) THEN
           WRITE (*, '(6x, a)') 'Double-exponential observation errors'
        END IF
     ELSE IF (filtertype == 100) THEN
        WRITE (*, '(6x, a, f5.2)') '-- Generate observations --'
        IF (dim_ens>1) THEN
           WRITE (*, '(14x, a)') 'Use ensemble mean for observations'
           WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
        ELSE
           WRITE (*, '(14x, a)') 'Generate observations from single ensemble state'
        END IF
     ELSE IF (filtertype == 12) THEN
        WRITE (*, '(21x, a)') 'Filter: PF with resampling'
        IF (subtype == 0) THEN
           WRITE (*, '(6x, a)') '-- Standard mode'
        ELSE IF (subtype == 5) THEN
           WRITE (*, '(6x, a)') '-- Offline mode'
        END IF
        WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
        IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
        WRITE (*, '(13x, a, i5)') 'reampling type:', pf_res_type
        WRITE (*, '(17x, a, i5)') 'noise type:', pf_noise_type
        WRITE (*, '(12x, a, f8.3)') 'noise amplitude:', pf_noise_amp
        IF (model_error) THEN
           WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
        END IF
     END IF
     WRITE (*, '(8x, a)') 'File names:'
     WRITE (*, '(11x, a, a)') 'Initialization: ', TRIM(file_ini)
     WRITE (*, '(11x, a, a)') 'Observations:   ', TRIM(file_obs)
     WRITE (*, '(11x, a, a)') 'Output:         ', TRIM(file_asml)
     IF (twin_experiment) &
          WRITE (*, '(/6x, a)') 'Run twin experiment with synthetic observations'
     IF (filtertype==100 .OR. twin_experiment) &
          WRITE (*, '(11x, a, a)') 'File for synthetic observations: ', TRIM(file_syntobs)

  END IF screen2


! *****************************************************
! *** Call filter initialization routine on all PEs ***
! *****************************************************

  whichinit: IF (filtertype == 1) THEN
     ! *** SEIK with init by 2nd order exact sampling ***
     filter_param_i(1) = dim_state   ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
     filter_param_r(1) = forget      ! Forgetting factor

     CALL PDAF_init(filtertype, subtype, step_null, &
          filter_param_i, 7, &
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_pdaf, &
          screen, status_pdaf)
  ELSEIF (filtertype == 2) THEN
     ! *** EnKF with init by 2nd order exact sampling ***
     filter_param_i(1) = dim_state   ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = rank_analysis_enkf ! Maximum rank for matrix inversion
     filter_param_r(1) = forget      ! Forgetting factor
     
     CALL PDAF_init(filtertype, subtype, step_null, &
          filter_param_i, 3, &
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_pdaf, &
          screen, status_pdaf)
  ELSEIF (filtertype == 4) THEN
     ! *** ETKF with init by 2nd order exact sampling ***
     filter_param_i(1) = dim_state   ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = dim_lag     ! Size of lag in smoother
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_r(1) = forget      ! Forgetting factor
     
     CALL PDAF_init(filtertype, subtype, step_null, &
          filter_param_i, 6, &
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_pdaf, &
          screen, status_pdaf)
  ELSEIF (filtertype == 6) THEN
     ! *** ESTKF with init by 2nd order exact sampling ***
     filter_param_i(1) = dim_state   ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = dim_lag     ! Size of lag in smoother
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
     filter_param_r(1) = forget      ! Forgetting factor

     CALL PDAF_init(filtertype, subtype, step_null, &
          filter_param_i, 7, &
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_pdaf, &
          screen, status_pdaf)
  ELSEIF (filtertype == 9) THEN
     ! *** NETF ***
     filter_param_i(1) = dim_state   ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = dim_lag     ! Size of lag in smoother
     filter_param_i(4) = 0           ! Not used for NETF (Whether to perform incremental analysis)
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_i(7) = type_winf   ! Type of weights inflation
     filter_param_r(1) = forget      ! Forgetting factor
     filter_param_r(2) = limit_winf  ! Limit for weights inflation
     
     CALL PDAF_init(filtertype, subtype, step_null, &
          filter_param_i, 7, &
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_pdaf, &
          screen, status_pdaf)
  ELSEIF (filtertype == 100) THEN
     ! *** LETKF with init by 2nd order exact sampling ***
     filter_param_i(1) = dim_state   ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_r(1) = forget      ! Forgetting factor
     
     CALL PDAF_init(filtertype, subtype, step_null, &
          filter_param_i, 2, &
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_pdaf, &
          screen, status_pdaf)
  ELSEIF (filtertype == 12) THEN
     ! *** Particle Filter ***
     filter_param_i(1) = dim_state     ! State dimension
     filter_param_i(2) = dim_ens       ! Size of ensemble
     filter_param_r(1) = pf_noise_amp  ! Noise amplitude
! Optional parameters; you need to re-set the number of parameters if you use them
     filter_param_i(3) = pf_res_type   ! Resampling type
     filter_param_i(4) = pf_noise_type ! Perturbation type
     filter_param_i(5) = type_forget   ! Type of forgetting factor
     filter_param_i(6) = type_winf     ! Type of weights inflation
     filter_param_r(2) = forget        ! Forgetting factor
     filter_param_r(3) = limit_winf    ! Limit for weights inflation

     CALL PDAF_init(filtertype, subtype, step_null, &
          filter_param_i, 6, &
          filter_param_r, 3, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_pdaf, &
          screen, status_pdaf)
  END IF whichinit

! *** Check whether initialization of PDAF was successful ***
  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in initialization of PDAF - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF

  ! Set initial time
  time = time + REAL(step_null) * dt

  ! Initialize netcdf output
  CALL init_netcdf_asml(step_null, dt, dim_state, filtertype, subtype, &
       dim_ens, forget, type_ensinit, rms_obs, delt_obs, total_steps, &
       seedset, stepnull_means, dim_lag)

  ! Initialize mask for observation gaps
  IF (use_obs_mask) THEN
     CALL init_obs_mask(dim_state)
  END IF

  ! Initialize file for synthetic observations
  IF (filtertype==100) THEN
     CALL init_file_syn_obs(dim_state, file_syntobs, 0)
  END IF


END SUBROUTINE init_pdaf
