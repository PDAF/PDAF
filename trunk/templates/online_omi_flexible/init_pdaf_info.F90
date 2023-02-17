!$Id: init_pdaf_info.F90 871 2021-11-22 16:44:34Z lnerger $
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

  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: filtertype, subtype, dim_ens, delt_obs, model_error, &
       model_err_amp, forget, rank_analysis_enkf, &
       dim_lag, twin_experiment, pf_res_type, &
       pf_noise_type, pf_noise_amp, type_hyb, hyb_gamma, hyb_kappa

  IMPLICIT NONE


! *****************************
! *** Initial Screen output ***
! *****************************

  IF (filtertype == 1) THEN
     WRITE (*, '(21x, a)') 'Filter: SEIK'
     IF (subtype == 2) THEN
        WRITE (*, '(6x, a)') '-- fixed error-space basis'
     ELSE IF (subtype == 3) THEN
        WRITE (*, '(6x, a)') '-- fixed state covariance matrix'
     ELSE IF (subtype == 4) THEN
        WRITE (*, '(6x, a)') '-- use ensemble transformation'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 2) THEN
     WRITE (*, '(21x, a)') 'Filter: EnKF'
     IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode'
     END IF
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
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 4) THEN
     WRITE (*, '(21x, a)') 'Filter: ETKF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Variant using T-matrix'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(6x, a)') '-- Variant following Hunt et al. (2007)'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
     IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 5) THEN
     WRITE (*, '(21x, a)') 'Filter: LETKF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Variant using T-matrix'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(6x, a)') '-- Variant following Hunt et al. (2007)'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
     IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 6) THEN
     WRITE (*, '(21x, a)') 'Filter: ESTKF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Standard mode'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
     IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 7) THEN
     WRITE (*, '(21x, a)') 'Filter: LESTKF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Standard mode'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
     IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 8) THEN
     WRITE (*, '(21x, a)') 'Filter: localized EnKF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Standard mode'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode'
     END IF
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
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Standard mode'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
     IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 10) THEN
     WRITE (*, '(21x, a)') 'Filter: LNETF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Standard mode'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
     IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 11) THEN
     WRITE (*, '(21x, a)') 'Filter: LKNETF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- HNK: 2-step LKNETF with NETF before LETKF'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(6x, a)') '-- HKN: 2-step LKNETF with LETKF before NETF'
     ELSE IF (subtype == 4) THEN
        WRITE (*, '(6x, a)') '-- HSync: LKNETF synchronous'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode - HNK: 2-step LKNETF with NETF before LETKF'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
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
     ELSEIF (type_hyb == 4) THEN
        WRITE (*, '(6x, a)') '-- use gamma_sklin: 1 - min(s,k)/sqrt(hyb_kappa) >= 1-N_eff/N>=hyb_gamma'
     END IF
     WRITE (*, '(8x, a, f7.2)') 'hybrid weight gamma:', hyb_gamma
     WRITE (*, '(10x, a, f7.2)') 'hybrid norm kappa:', hyb_kappa
     IF (model_error) THEN
        WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
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
  IF (twin_experiment) &
       WRITE (*, '(/6x, a)') 'Run twin experiment with synthetic observations'
!  IF (filtertype==100 .OR. twin_experiment) &
!       WRITE (*, '(11x, a, a)') 'File for synthetic observations: ', TRIM(file_syntobs)

END SUBROUTINE init_pdaf_info
