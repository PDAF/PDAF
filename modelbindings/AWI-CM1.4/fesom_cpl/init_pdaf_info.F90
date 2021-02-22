!$Id: init_pdaf_info.F90 2249 2020-04-06 06:42:37Z lnerger $
!>   Screen output on assimilation configuration
!!
!! This routine performs a model-sided screen output about
!! the coniguration of the data assimilation system.
!! Using this output is optional. Most of the information
!! is also displayed by PDAF itself when it is initialized
!! in PDAF_init. Not displayed by PDAF is the assimilation
!! interval (delt_obs_ocn), which is unknown to PDAF.
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
SUBROUTINE init_pdaf_info()

  USE mod_assim_pdaf, &           ! Variables for assimilation
       ONLY: filtertype, subtype, dim_ens,  &
       forget, dim_state, twin_experiment
  USE mod_assim_oce_pdaf, &       ! Variables for assimilation - ocean-specific
       ONLY: delt_obs_ocn, file_syntobs

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
  ELSE IF (filtertype == 4) THEN
     WRITE (*, '(21x, a)') 'Filter: ETKF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Variant using T-matrix'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(6x, a)') '-- Variant following Hunt et al. (2007)'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode'
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
  ELSE IF (filtertype == 6) THEN
     WRITE (*, '(21x, a)') 'Filter: ESTKF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Standard mode'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode'
     END IF
  ELSE IF (filtertype == 7) THEN
     WRITE (*, '(a, 21x, a)') 'FESOM-PDAF','Filter: LESTKF'
     IF (subtype == 0) THEN
        WRITE (*, '(a, 6x, a)') 'FESOM-PDAF','-- Standard mode'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 6x, a)') 'FESOM-PDAF','-- Offline mode'
     END IF
  ELSE IF (filtertype == 11) THEN
     WRITE (*, '(a, 6x, a, f5.2)') 'FESOM-PDAF','-- Generate observations --'
     IF (dim_ens>1) THEN
        WRITE (*, '(a, 14x, a)') 'FESOM-PDAF','Use ensemble mean for observations'
     ELSE
        WRITE (*, '(a, 14x, a)') 'FESOM-PDAF','Generate observations from single ensemble state'
     END IF
  END IF     
  WRITE (*, '(a, 14x, a, i5)') 'FESOM-PDAF','ensemble size:', dim_ens
  IF (subtype /= 5) WRITE (*, '(a, 6x, a, i5)') 'FESOM-PDAF','Assimilation interval:', delt_obs_ocn
  WRITE (*, '(a, 10x, a, f5.2)') 'FESOM-PDAF','forgetting factor:', forget
  IF (twin_experiment) &
       WRITE (*, '(/a,6x, a)') 'FESOM-PDAF','Run twin experiment with synthetic observations'
  IF (filtertype==11 .OR. twin_experiment) &
       WRITE (*, '(a,11x, a, a)') 'FESOM-PDAF','File for synthetic observations: ', TRIM(file_syntobs)
  WRITE (*, '(a, 8x, a, 1x, i9)') 'FESOM-PDAF','FESOM state dimension:',dim_state


END SUBROUTINE init_pdaf_info
