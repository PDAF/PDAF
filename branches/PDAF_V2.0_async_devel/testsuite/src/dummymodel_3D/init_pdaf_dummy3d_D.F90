!$Id: init_pdaf_dummy3d_D.F90 783 2009-12-07 10:28:43Z lnerger $
!BOP
!
! !ROUTINE: init_pdaf - Interface routine to call initialization of PDAF
!
! !INTERFACE:
SUBROUTINE init_pdaf(total_steps)

! !DESCRIPTION:
! This routine collects the initialization of variables
! for PDAF as well as the call to the initialization
! routine PDAF_init.
!
! !REVISION HISTORY:
! 2008-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE parser, &
       ONLY: parse
  USE mod_model, &
       ONLY: step_null, dims, dim_l, dt
  USE mod_parallel, &
       ONLY: mype_world, mype_model, n_modeltasks, task_id, &
       COMM_model, COMM_filter, COMM_couple, filterpe
  USE mod_assimilation, &
       ONLY: screen, filtertype, subtype, dim_ens, n_obs, delt_obs, &
       rms_obs, model_error, model_err_amp, incremental, covartype, &
       type_forget, forget, dim_bias, epsilon, rank_analysis_enkf, &
       locweight, local_range, ewidth, int_rediag, filename

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(IN) :: total_steps   ! Number of time steps in experiment

! !CALLING SEQUENCE:
! Called by: main
! Calls: PDAF_init
! Calls: parse
!EOP

! Local variables
  CHARACTER(len=32) :: handle  ! handle for command line parser
  INTEGER :: dim_state_p       ! PE-local state dimension
  INTEGER :: filter_param_i(5) ! Integer parameter array for filter
  REAL    :: filter_param_r(2) ! Real parameter array for filter
  integer :: flag              ! PDAF status flag

  ! External subroutines
  EXTERNAL :: init_seik  ! SEIK ensemble initialization
  EXTERNAL :: init_seek  ! SEEK EOF initialization
  EXTERNAL :: init_enkf  ! EnKF ensemble initialization
  

! ***************************
! ***   Initialize PDAF   ***
! ***************************

  IF (mype_world == 0) THEN
     WRITE (*,'(/1x,a)') 'INITIALIZE PDAF'
  END IF


! **********************************************************
! ***               CONTROL OF PDAF                      ***
! **********************************************************

! *** IO options ***
  screen      = 2   ! Write screen output (1) for output, (2) add timings

! *** specifications for observations ***
  ! avg. observation error (used for assimilation)
  n_obs   = dims(1) * dims(2) ! Dimension of observation vector ("surface")
  rms_obs = 1.0     ! This error is the standard deviation 
                    ! for the Gaussian distribution 
  delt_obs = 2      ! Time step interval between analysis/assimilation steps

! *** Filter specific variables
  filtertype = 0    ! Type of filter
                    !   SEEK (0), SEIK (1), EnKF (2), LSEIK (3), ETKF (4), LETKF (5)
  dim_ens = 101     ! Size of ensemble for SEIK/LSEIK/EnKF/ETKF
                    ! Number of EOFs to be used for SEEK
  subtype = 0       ! subtype of filter: 
                    !   SEEK: 
                    !     (0) evolve normalized modes
                    !     (1) evolve scaled modes with unit U
                    !     (2) fixed basis (V); variable U matrix
                    !     (3) fixed covar matrix (V,U kept static)
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
                    !     (0) ETKF using T-matrix like SEIK
                    !     (1) LETKF following Hunt et al. (2007)
                    !       There are no fixed basis/covariance cases, as
                    !       these are equivalent to LSEIK subtypes 2/3
  int_rediag = 1    ! Interval of analysis steps to perform 
                    !    re-diagonalization in SEEK
  incremental = 0   ! (1) to perform incremental updating (only in SEIK/LSEIK!)
  covartype = 1     ! Definition of factor in covar. matrix used in SEIK
                    ! (0) for (r+1)^-1 (old SEIK); (1): for r^-1 (real ensemble
                    ! covariance matrix) This parameter has also to be set internally
                    ! in PDAF_init
  rank_analysis_enkf = 0   ! rank to be considered for inversion of HPH
                    ! in analysis of EnKF; (0) for analysis w/o eigendecomposition
  type_forget = 0   ! Type of forgetting factor in SEIK/LSEIK
                    ! (0): fixed; (1) global adaptive; (2) local adaptive for LSEIK
  forget  = 1.0     ! Forgetting factor
  epsilon = 1.0E-4  ! epsilon for approx. TLM evolution in SEEK/SEIK
  local_range = 10  ! Range in grid points for observation domain in LSEIK
  locweight = 0     ! Type of localizating weighting
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with EWIDTH
                    !   (2) use 5th-order polynomial
  ewidth = 5.0      ! Inverse relative distance at which the weight
                    ! of the inverse observation variance is 1/e.

! *** File names
  filename = 'output.dat'


! ********************************
! ***      INITIALIZATION      ***
! ********************************

! *** Parse command line options ***

  ! Settings for model and time stepping
  handle = 'dim_bias'                ! Dimension of bias vector
  CALL parse(handle, dim_bias)
  handle = 'model_error'             ! Control application of model error
  CALL parse(handle, model_error)
  handle = 'model_err_amp'           ! Amplitude of model error
  CALL parse(handle, model_err_amp)

  ! Observation settings
  handle = 'n_obs'                   ! Dimension of observation vector
  CALL parse(handle, n_obs)
  handle = 'delt_obs'                ! Time step interval between filter analyses
  CALL parse(handle, delt_obs)
  handle = 'rms_obs'                 ! Assumed uniform RMS error of the observations
  CALL parse(handle, rms_obs)

  ! General settings for PDAF
  handle = 'screen'                  ! set verbosity of PDAF
  CALL parse(handle, screen)
  handle = 'dim_ens'                 ! set ensemble size/rank of covar matrix
  CALL parse(handle, dim_ens)
  handle = 'filtertype'              ! Choose filter algorithm
  CALL parse(handle, filtertype)
  handle = 'subtype'                 ! Set subtype of filter
  CALL parse(handle, subtype)
  handle = 'incremental'             ! Set whether to use incremental updating
  CALL parse(handle, incremental)

  ! Filter-specific settings
  handle = 'epsilon'                 ! Set EPSILON for SEEK
  CALL parse(handle, epsilon)
  handle = 'int_rediag'              ! Time step interval for rediagonalization in SEEK
  CALL parse(handle, int_rediag)
  handle = 'rank_analysis_enkf'      ! Set rank for pseudo inverse in EnKF
  CALL parse(handle, rank_analysis_enkf)
  handle = 'type_forget'             ! Set type of forgetting factor
  CALL parse(handle, type_forget)
  handle = 'forget'                  ! Set forgetting factor
  CALL parse(handle,forget)

  ! Settings for localization in LSEIK/LETKF
  handle = 'local_range'             ! Set range in grid points for observation domain
  CALL parse(handle, local_range)
  handle = 'locweight'               ! Set type of localizating weighting
  CALL parse(handle, locweight)
  handle = 'ewidth'                  ! Set inverse relative distance for weight = 1/e
  CALL parse(handle, ewidth)

  ! Setting for file output
  handle = 'filename'                ! Set name of output file
  CALL parse(handle, filename)


! *** Initial Screen output ***
  screen2: IF (mype_world == 0) THEN

     IF (filtertype == 0) THEN
        WRITE (*, '(/21x, a)') 'Filter: SEEK'
        IF (subtype == 2) THEN
           WRITE (*, '(6x, a)') '-- fixed basis filter with update of matrix U'
           WRITE (*, '(6x, a)') '-- no re-diagonalization of VUV^T'
        END IF
        IF (subtype == 3) THEN
           WRITE (*, '(6x, a)') '-- fixed basis filter & no update of matrix U'
           WRITE (*, '(6x, a)') '-- no re-diagonalization of VUV^T'
        END IF
        WRITE (*, '(13x, a, i5)') 'number of EOFs:', dim_ens
        WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
        WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
        IF ((int_rediag > 0) .AND. ((subtype /= 2) .OR. (subtype /= 3))) &
             WRITE (*, '(10x, a, i4, a)') &
             'Re-diag each ', int_rediag, '-th analysis step'
     ELSE IF (filtertype == 1) THEN
        WRITE (*, '(21x, a)') 'Filter: SEIK'
        IF (subtype == 2) THEN
           WRITE (*, '(6x, a)') '-- fixed error-space basis'
        END IF
        IF (subtype == 3) THEN
           WRITE (*, '(6x, a)') '-- fixed state covariance matrix'
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
        WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
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
        END IF
        IF (subtype == 3) THEN
           WRITE (*,'( 6x, a)') '-- fixed state covariance matrix'
        END IF
        WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
        WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
        WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
        IF (model_error) THEN
           WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
        END IF
     END IF
     
  END IF screen2


! *** Local state dimensions ***

  ! State dimension for my PE-local domain
  dim_state_p = dim_l(1) * dim_l(2) * dim_l(3)


! *****************************************************
! *** Call filter initialization routine on all PEs ***
! *****************************************************

  whichinit: IF (filtertype == 0) THEN
     ! *** SEEK ***
     filter_param_i(1) = dim_state_p ! State dimension
     filter_param_i(2) = dim_ens     ! Number of EOFs
     filter_param_i(3) = int_rediag  ! Interval to perform rediagonalization
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_r(1) = forget      ! Forgetting factor
     filter_param_r(2) = epsilon     ! Epsilon to tangent linear forecast
      
     CALL PDAF_init(filtertype, subtype, step_null, &
          filter_param_i, 5,&
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_seek, &
          screen, flag)
  ELSEIF (filtertype == 1) THEN
     ! *** SEIK with init by 2nd order exact sampling ***
     filter_param_i(1) = dim_state_p ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_r(1) = forget      ! Forgetting factor
     
     CALL PDAF_init(filtertype, subtype, step_null, &
          filter_param_i, 5,&
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_seik, &
          screen, flag)
     ! Initialize in distributed manner
     ! Used for restarting
!         call PDAF_init_dist(filtertype, subtype, step_null, &
!              filter_param_i, 5,&
!              filter_param_r, 2, &
!              COMM_model, COMM_filter, COMM_couple, &
!              task_id, n_modeltasks, filterpe, init_seik, &
!              screen, flag)
  ELSEIF (filtertype == 2) THEN
     ! *** EnKF with Monte Carlo init ***
     filter_param_i(1) = dim_state_p ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = rank_analysis_enkf ! Rank of speudo-inverse in analysis
     filter_param_r(1) = forget      ! Forgetting factor
     
     CALL PDAF_init(filtertype, subtype, step_null, &
          filter_param_i, 5,&
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_enkf, &
          screen, flag)
!         ! Initialize in distributed manner
!         ! Used for restarting
!         call filter_init_dist(filtertype, subtype, step_null, &
!              filter_param_i, 5,&
!              filter_param_r, 2, &
!              COMM_model, COMM_filter, COMM_couple, &
!              task_id, n_modeltasks, filterpe, SEIKEnKF_init_restart, &
!              screen, flag)
  ELSEIF (filtertype == 3) THEN
     ! *** LSEIK with init by 2nd order exact sampling ***
     filter_param_i(1) = dim_state_p ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_r(1) = forget      ! Forgetting factor
     
!      CALL PDAF_init(filtertype, subtype, step_null, &
!           filter_param_i, 5,&
!           filter_param_r, 2, &
!           COMM_model, COMM_filter, COMM_couple, &
!           task_id, n_modeltasks, filterpe, init_seik, &
!           screen, flag)
    ! Initialize in distributed manner
    ! Used for restarting
!         call PDAF_init_dist(filtertype, subtype, step_null, &
!              filter_param_i, 5,&
!              filter_param_r, 2, &
!              COMM_model, COMM_filter, COMM_couple, &
!              task_id, n_modeltasks, filterpe, SEIKEnKF_init_restart, &
!              screen, flag)
  END IF whichinit

END SUBROUTINE init_pdaf
