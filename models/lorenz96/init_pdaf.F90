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
! This variant is for the Lorenz96 model.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE pdaf_interfaces_module, &   ! Interface definitions to PDAF core routines
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
       locweight, cradius, cradius2, sradius, &
       file_ini, type_ensinit, seedset, type_trans, &
       type_sqrt, stepnull_means, dim_lag, time, &
       twin_experiment, pf_res_type, pf_noise_type, pf_noise_amp, &
       type_hyb, hyb_gamma, hyb_kappa, &
       type_winf, limit_winf
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
  REAL    :: filter_param_r(4) ! Real parameter array for filter
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
  filtertype = 1    ! Type of filter
                    !   SEIK (1), EnKF (2), LSEIK (3), ETKF (4), LETKF (5)
                    !   ESTKF (6), LESTKF (7), NETF (9), LNETF (10)
                    !   LKNETF (11), PF (12), GENOBS (100)
  dim_ens = 30      ! Size of ensemble
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
                    !   LSEIK:
                    !     (0) mean forecast;
                    !     (2) fixed error space basis
                    !     (3) fixed state covariance matrix
                    !     (4) LSEIK with ensemble transformation
                    !   ETKF:
                    !     (0) ETKF using T-matrix like SEIK
                    !     (1) ETKF following Hunt et al. (2007)
                    !       There are no fixed basis/covariance cases, as
                    !       these are equivalent to SEIK subtypes 2/3
                    !   LETKF:
                    !     (0) LETKF using T-matrix like SEIK
                    !     (1) LETKF following Hunt et al. (2007)
                    !       There are no fixed basis/covariance cases, as
                    !       these are equivalent to LSEIK subtypes 2/3
                    !   ESTKF:
                    !     (0) Standard form of ESTKF
                    !     (2) fixed ensemble perturbations
                    !     (3) fixed state covariance matrix
                    !   LESTKF:
                    !     (0) Standard form of LESTKF
                    !     (2) fixed ensemble perturbations
                    !     (3) fixed state covariance matrix
                    !   NETF:
                    !     (0) Standard form of NETF
                    !   LNETF:
                    !     (0) Standard form of LNETF
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
  type_forget = 0   ! Type of forgetting factor 
                    ! SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
                    !   (0) fixed
                    !   (1) global adaptive
                    !   (2) local adaptive for LSEIK/LETKF/LESTKF
                    ! NETF/LNETF/PF
                    !   (0) apply inflation on forecast ensemble
                    !   (2) apply inflation on analysis ensemble
  type_sqrt = 0     ! Type of transform matrix square-root
                    !   (0) symmetric square root, (1) Cholesky decomposition
  type_winf = 0     ! NETF/LNETF: Type of weights inflation: (1) use N_eff/N>limit_winf
  limit_winf = 0.0  ! Limit for weights inflation
  rank_analysis_enkf = 0   ! ENKF only: rank to be considered for inversion of HPH
                    ! in analysis step of EnKF; (0) for analysis w/o eigendecomposition
  type_hyb = 0      ! LKNETF: Type of hybrid weight: 
                    !   (0) use fixed hybrid weight hyb_gamma
                    !   (1) use gamma_lin: (1 - N_eff/N_e)*hyb_gamma
                    !   (2) use gamma_alpha: hybrid weight from N_eff/N>=hyb_gamma
                    !   (3) use gamma_ska: 1 - min(s,k)/sqrt(hyb_kappa) with N_eff/N>=hyb_gamma
                    !   (4) use gamma_sklin: 1 - min(s,k)/sqrt(hyb_kappa) >= 1-N_eff/N>=hyb_gamma
  hyb_gamma =  1.0  ! Hybrid filter weight for state (1.0: LETKF, 0.0: LNETF)
  hyb_kappa = 30.0  ! Hybrid norm for using skewness and kurtosis (type_hyb 3 or 4)
  model_error = .false.     ! Whether to apply model error noise
  model_err_amp = 0.1       ! Amplitude of model noise (times dt for error variance)
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

! *** Whether to run twin experiment assimilating synthetic observations ***
  twin_experiment = .false.

! *** IO options ***
  delt_write_asml = 1    ! Output interval for state information in assimilation intervals
  write_states = .TRUE.  ! Whether to write estimates states into the output file
  write_stats  = .TRUE.  ! Whether to write time dependent ensemble statistics (skewness, kurtosis)
  write_ens  = .FALSE.   ! Whether to write full time dependent ensemble stats
  stepnull_means = 3001  ! Step at which the second computation of time mean error is started
                         ! (first computation of mean sis always starting at initial step)

! *** specifications for observations ***
  ! avg. observation error (used for assimilation)
  rms_obs = 1.0      ! This error is the standard deviation 
                     ! for the Gaussian distribution 
  delt_obs = 1       ! Time step interval between analysis/assimilation steps
  use_obs_mask = .FALSE. ! Whether to use observations with gaps
  use_maskfile = .FALSE. ! If a mask is used read it from file
  numobs = dim_state ! If not read from file use this number of obs. (1 to numobs)
  dx_obs = 1         ! grid point distance of observations (if not read from file)
  obs_err_type = 0   ! Observation errors: (0) for Gaussian (1) for double-exponential

! *** Filter specific variables ***
  type_ensinit = 'eof' ! 'eof' for 2nd-order exact sampling from EOFs
                    !    'rnd' for random sampling from true state trajectory
  seedset = 1       ! Index of set of seeds to be used for init (only for 'rnd')
  covartype = 1     ! Definition of factor in covar. matrix used in SEIK
                    ! (0) for (r+1)^-1 (old SEIK); (1): for r^-1 (real ensemble
                    ! covariance matrix) This parameter has also to be set internally
                    ! in PDAF_init
  cradius = 5       ! Cut-off radius in grid points for observation domain 
  locweight = 0     ! Type of localizating weighting
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRADIUS
                    !   (2) use 5th-order polynomial 
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
  sradius = cradius ! Support radius for 5th-order polynomial
                    ! radius for 1/e for exponential weighting

! *** File names ***
  file_asml = 'assimilation.nc' ! Output file
  file_ini = 'covar.nc'         ! Initialization file
  file_obs = 'obs.nc'           ! File holding observations
  file_obs_mask = 'obsmask.txt' ! File holding observation mask
  file_syntobs = 'syntobs.nc'   ! File holding synthetic observations generated by PDAF


! ********************************
! ***      INITIALIZATION      ***
! ********************************

! *** Parse command line options   ***
! *** This is optional, but useful ***

  CALL init_pdaf_parse()

  IF (cradius > dim_state/2) THEN
     cradius = dim_state/2
     WRITE (*,*) 'NOTICE: cradius too large. Reset to dim_state/2'
  END IF
  IF (cradius2 > dim_state/2) THEN
     cradius2 = dim_state/2
     WRITE (*,*) 'NOTICE: cradius2 too large. Reset to dim_state/2'
  END IF
     


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
     ELSE IF (filtertype == 3) THEN
        WRITE (*, '(21x, a)') 'Filter: LSEIK'
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
           WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
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
     ELSE IF (filtertype == 5) THEN
        WRITE (*, '(21x, a)') 'Filter: LETKF'
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
           WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
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
     ELSE IF (filtertype == 7) THEN
        WRITE (*, '(21x, a)') 'Filter: LESTKF'
        WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
        IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
        WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
        WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
        IF (model_error) THEN
           WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
        END IF
     ELSE IF (filtertype == 8) THEN
        WRITE (*, '(21x, a)') 'Filter: LEnKF'
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
     ELSE IF (filtertype == 9) THEN
        WRITE (*, '(21x, a)') 'Filter: NETF'
        WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
        IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
        WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
        WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
        IF (type_winf==1) THEN
           WRITE (*, '(6x, a, f5.2)') 'inflate particle weights so that N_eff/N > ', limit_winf
        END IF
        IF (model_error) THEN
           WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
        END IF
        IF (obs_err_type==0) THEN
           WRITE (*, '(6x, a)') 'Gaussian observation errors'
        ELSE IF (obs_err_type==1) THEN
           WRITE (*, '(6x, a)') 'Double-exponential observation errors'
        END IF
     ELSE IF (filtertype == 10) THEN
        WRITE (*, '(21x, a)') 'Filter: LNETF'
        WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
        IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
        WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
        WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
        IF (type_winf==1) THEN
           WRITE (*, '(6x, a, f5.2)') 'inflate particle weights so that N_eff/N > ', limit_winf
        END IF
        IF (model_error) THEN
           WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
        END IF
        IF (obs_err_type==0) THEN
           WRITE (*, '(6x, a)') 'Gaussian observation errors'
        ELSE IF (obs_err_type==1) THEN
           WRITE (*, '(6x, a)') 'Double-exponential observation errors'
        END IF
     ELSE IF (filtertype == 11) THEN
        WRITE (*, '(21x, a)') 'Filter: LKNETF'
        WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
        IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
        WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
        WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
        WRITE (*, '(12x, a, 2f6.2)') 'hybrid weight gamma: ', hyb_gamma
        WRITE (*, '(11x, a, 2f6.2)') 'hybrid norm kappa: ', hyb_kappa
        IF (model_error) THEN
           WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
        END IF
        IF (obs_err_type==0) THEN
           WRITE (*, '(6x, a)') 'Gaussian observation errors'
        ELSE IF (obs_err_type==1) THEN
           WRITE (*, '(6x, a)') 'Double-exponential observation errors'
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
     ELSE IF (filtertype == 100) THEN
        WRITE (*, '(6x, a, f5.2)') '-- Generate observations --'
        IF (dim_ens>1) THEN
           WRITE (*, '(14x, a)') 'Use ensemble mean for observations'
           WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
        ELSE
           WRITE (*, '(14x, a)') 'Generate observations from single ensemble state'
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
  ELSEIF (filtertype == 3) THEN
     ! *** LSEIK with init by 2nd order exact sampling ***
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
  ELSEIF (filtertype == 5) THEN
     ! *** LETKF with init by 2nd order exact sampling ***
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
  ELSEIF (filtertype == 7) THEN
     ! *** LESTKF with init by 2nd order exact sampling ***
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
  ELSEIF (filtertype == 8) THEN
     ! *** LEnKF with init by 2nd order exact sampling ***
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
  ELSEIF (filtertype == 10) THEN
     ! *** LNETF ***
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
  ELSEIF (filtertype == 11) THEN
     ! *** LKNETF ***
     filter_param_i(1) = dim_state   ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = dim_lag     ! Size of lag in smoother
     filter_param_i(4) = 0           ! Not used for NETF (Whether to perform incremental analysis)
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_i(7) = type_hyb    ! Type of hybrid weight
     filter_param_r(1) = forget      ! Forgetting factor
     filter_param_r(2) = hyb_gamma   ! Hybrid filter weight for state
     filter_param_r(3) = hyb_kappa   ! Normalization factor for hybrid weight 
     
     CALL PDAF_init(filtertype, subtype, step_null, &
          filter_param_i, 7, &
          filter_param_r, 3, &
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
  ELSEIF (filtertype == 100) THEN
     ! *** Observation generation ***
     filter_param_i(1) = dim_state   ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_r(1) = forget      ! Forgetting factor
     
     CALL PDAF_init(filtertype, subtype, step_null, &
          filter_param_i, 2, &
          filter_param_r, 2, &
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
       dim_ens, forget, type_ensinit, cradius, cradius2, &
       locweight, sradius, rms_obs, delt_obs, total_steps, &
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
