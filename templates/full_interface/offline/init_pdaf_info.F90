!>  Screen output on assimilation configuration
!!
!! This routine performs a model-sided screen output about
!! the coniguration of the data assimilation system.
!! Using this output is optional. Most of the information
!! is also displayed by PDAF itself when it is initialized
!! in PDAF_init. Not displayed by PDAF is the assimilation
!! interval (delt_obs), which is unknown to PDAF.
!!
!! __Revision history:__
!! * 2011-05 - Lars Nerger - Initial code extracted from init_pdaf
!! * Later revisions - see repository log
!!
SUBROUTINE init_pdaf_info()

  USE PDAF                     ! PDAF interface and parameters
  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: filtertype, subtype, dim_ens, model_error, &
       model_err_amp, forget, rank_ana_enkf, &
       dim_lag, pf_res_type, &
       pf_noise_type, pf_noise_amp, type_hyb, hyb_gamma, hyb_kappa

  IMPLICIT NONE


! *****************************
! *** Initial Screen output ***
! *****************************

  IF (filtertype == PDAF_DA_SEIK) THEN
     WRITE (*, '(21x, a)') 'Filter: SEIK'
     IF (subtype == 2) THEN
        WRITE (*, '(6x, a)') '-- fixed error-space basis'
     ELSE IF (subtype == 10) THEN
        WRITE (*, '(6x, a)') '-- fixed state covariance matrix'
     ELSE IF (subtype == 11) THEN
        WRITE (*, '(6x, a)') '-- use ensemble transformation'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == PDAF_DA_ENKF) THEN
     WRITE (*, '(21x, a)') 'Filter: EnKF'
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
     IF (rank_ana_enkf > 0) THEN
        WRITE (*, '(6x, a, i5)') &
             'analysis with pseudo-inverse of HPH, rank:', rank_ana_enkf
     END IF
  ELSE IF (filtertype == PDAF_DA_LSEIK) THEN
     WRITE (*, '(21x, a)') 'Filter: LSEIK'
     IF (subtype == 2) THEN
        WRITE (*, '(6x, a)') '-- use ensemble transformation'
     ELSEIF (subtype == 10) THEN
        WRITE (*, '(6x, a)') '-- fixed error-space basis'
     ELSE IF (subtype == 11) THEN
        WRITE (*, '(6x, a)') '-- fixed state covariance matrix'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == PDAF_DA_ETKF) THEN
     WRITE (*, '(21x, a)') 'Filter: ETKF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Variant using T-matrix'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(6x, a)') '-- Variant following Hunt et al. (2007)'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == PDAF_DA_LETKF) THEN
     WRITE (*, '(21x, a)') 'Filter: LETKF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Variant using T-matrix'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(6x, a)') '-- Variant following Hunt et al. (2007)'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == PDAF_DA_ESTKF) THEN
     WRITE (*, '(21x, a)') 'Filter: ESTKF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Standard mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == PDAF_DA_LESTKF) THEN
     WRITE (*, '(21x, a)') 'Filter: LESTKF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Standard mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == PDAF_DA_LENKF) THEN
     WRITE (*, '(21x, a)') 'Filter: localized EnKF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Standard mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
     IF (rank_ana_enkf > 0) THEN
        WRITE (*, '(6x, a, i5)') &
             'analysis with pseudo-inverse of HPH, rank:', rank_ana_enkf
     END IF
  ELSE IF (filtertype == PDAF_DA_NETF) THEN
     WRITE (*, '(21x, a)') 'Filter: NETF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Standard mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == PDAF_DA_LNETF) THEN
     WRITE (*, '(21x, a)') 'Filter: LNETF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Standard mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == PDAF_DA_LKNETF) THEN
     WRITE (*, '(21x, a)') 'Filter: LKNETF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- HNK: 2-step LKNETF with NETF before LETKF'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(6x, a)') '-- HKN: 2-step LKNETF with LETKF before NETF'
     ELSE IF (subtype == 2) THEN
        WRITE (*, '(6x, a)') '-- HSync: LKNETF synchronous'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     WRITE (*, '(10x, a, f7.2)') 'forgetting factor:', forget
     IF (type_hyb == 0) THEN
     ELSEIF (type_hyb == 0) THEN
        WRITE (*, '(6x, a)') '-- use fixed hybrid weight hyb_gamma'
     ELSEIF (type_hyb == 1) THEN
        WRITE (*, '(6x, a)') '-- use gamma_lin: (1 - N_eff/N_e)*hyb_gamma'
     ELSEIF (type_hyb == 2) THEN
        WRITE (*, '(6x, a)') '-- use gamma_alpha: hybrid weight from N_eff/N>=hyb_gamma'
     ELSEIF (type_hyb == 3) THEN
        WRITE (*, '(6x, a)') '-- use gamma_ska: 1 - min(s,k)/sqrt(hyb_kappa) with N_eff/N>=hyb_gamma'
     END IF
     WRITE (*, '(8x, a, f7.2)') 'hybrid weight gamma:', hyb_gamma
     WRITE (*, '(10x, a, f7.2)') 'hybrid norm kappa:', hyb_kappa
     IF (model_error) THEN
        WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == PDAF_DA_PF) THEN
     WRITE (*, '(21x, a)') 'Filter: PF with resampling'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Standard mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     WRITE (*, '(13x, a, i5)') 'reampling type:', pf_res_type
     WRITE (*, '(17x, a, i5)') 'noise type:', pf_noise_type
     WRITE (*, '(12x, a, f8.3)') 'noise amplitude:', pf_noise_amp
     IF (model_error) THEN
        WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == PDAF_DA_ENSRF) THEN
     WRITE (*, '(21x, a)') 'Filter: ENSRF/EAKF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- ENSRF with serial observaion processing'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(6x, a)') '-- EAKF with local linear regression'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     WRITE (*, '(10x, a, f7.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == PDAF_DA_GENOBS) THEN
     WRITE (*, '(6x, a, f5.2)') '-- Generate observations --'
     IF (dim_ens>1) THEN
        WRITE (*, '(14x, a)') 'Use ensemble mean for observations'
        WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     ELSE
        WRITE (*, '(14x, a)') 'Generate observations from single ensemble state'
     END IF
  ELSE IF (filtertype == PDAF_DA_3DVAR) THEN
     WRITE (*, '(21x, a)') 'Assimilation using 3D-Var'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Incremental 3D-Var with parameterized covariance matrix'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(6x, a)') '-- 3D ensemble Var using LESTKF for ensemble transformation'
     ELSE IF (subtype == 2) THEN
        WRITE (*, '(6x, a)') '-- 3D ensemble Var using ESTKF for ensemble transformation'
     ELSE IF (subtype == 3) THEN
        WRITE (*, '(6x, a)') '-- Hybrid 3D-Var using LESTKF for ensemble transformation'
     ELSE IF (subtype == 4) THEN
        WRITE (*, '(6x, a)') '-- Hybrid 3D-Var using ESTKF for ensemble transformation'
     END IF
  END IF     
!  IF (twin_experiment) &
!       WRITE (*, '(/6x, a)') 'Run twin experiment with synthetic observations'
!  IF (filtertype==100 .OR. twin_experiment) &
!       WRITE (*, '(11x, a, a)') 'File for synthetic observations: ', TRIM(file_syntobs)

END SUBROUTINE init_pdaf_info
