!$Id$
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
       forget, dim_state
  USE mod_assim_atm_pdaf, &       ! Variables for assimilation - atmosphere specific
       ONLY: delt_obs_atm

  IMPLICIT NONE


! *****************************
! *** Initial Screen output ***
! *****************************

  IF (filtertype == 1) THEN
     WRITE (*, '(a, 21x, a)') 'ECHAM-PDAF','Filter: SEIK'
     IF (subtype == 2) THEN
        WRITE (*, '(a, 6x, a)') 'ECHAM-PDAF','-- fixed error-space basis'
     ELSE IF (subtype == 3) THEN
        WRITE (*, '(a, 6x, a)') 'ECHAM-PDAF','-- fixed state covariance matrix'
     ELSE IF (subtype == 4) THEN
        WRITE (*, '(a, 6x, a)') 'ECHAM-PDAF','-- use ensemble transformation'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 6x, a)') 'ECHAM-PDAF','-- Offline mode'
     END IF
  ELSE IF (filtertype == 3) THEN
     WRITE (*, '(a, 21x, a)') 'ECHAM-PDAF','Filter: LSEIK'
     IF (subtype == 2) THEN
        WRITE (*, '(a, 6x, a)') 'ECHAM-PDAF','-- fixed error-space basis'
     ELSE IF (subtype == 3) THEN
        WRITE (*, '(a, 6x, a)') 'ECHAM-PDAF','-- fixed state covariance matrix'
     ELSE IF (subtype == 4) THEN
        WRITE (*, '(a, 6x, a)') 'ECHAM-PDAF','-- use ensemble transformation'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 6x, a)') 'ECHAM-PDAF','-- Offline mode'
     END IF
  ELSE IF (filtertype == 4) THEN
     WRITE (*, '(a, 21x, a)') 'ECHAM-PDAF','Filter: ETKF'
     IF (subtype == 0) THEN
        WRITE (*, '(a, 6x, a)') 'ECHAM-PDAF','-- Variant using T-matrix'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(a, 6x, a)') 'ECHAM-PDAF','-- Variant following Hunt et al. (2007)'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 6x, a)') 'ECHAM-PDAF','-- Offline mode'
     END IF
  ELSE IF (filtertype == 5) THEN
     WRITE (*, '(a, 21x, a)') 'ECHAM-PDAF','Filter: LETKF'
     IF (subtype == 0) THEN
        WRITE (*, '(a, 6x, a)') 'ECHAM-PDAF','-- Variant using T-matrix'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(a, 6x, a)') 'ECHAM-PDAF','-- Variant following Hunt et al. (2007)'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 6x, a)') 'ECHAM-PDAF','-- Offline mode'
     END IF
  ELSE IF (filtertype == 6) THEN
     WRITE (*, '(a, 21x, a)') 'ECHAM-PDAF','Filter: ESTKF'
     IF (subtype == 0) THEN
        WRITE (*, '(a, 6x, a)') 'ECHAM-PDAF','-- Standard mode'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 6x, a)') 'ECHAM-PDAF','-- Offline mode'
     END IF
  ELSE IF (filtertype == 7) THEN
     WRITE (*, '(a, 21x, a)') 'ECHAM-PDAF','Filter: LESTKF'
     IF (subtype == 0) THEN
        WRITE (*, '(a, 6x, a)') 'ECHAM-PDAF','-- Standard mode'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 6x, a)') 'ECHAM-PDAF','-- Offline mode'
     END IF
  ELSE IF (filtertype == 100) THEN
     WRITE (*, '(a, 6x, a, f5.2)') 'FESOM-PDAF','-- Generate observations --'
     IF (dim_ens>1) THEN
        WRITE (*, '(a, 14x, a)') 'FESOM-PDAF','Use ensemble mean for observations'
     ELSE
        WRITE (*, '(a, 14x, a)') 'FESOM-PDAF','Generate observations from single ensemble state'
     END IF
  END IF     
  WRITE (*, '(a, 14x, a, i5)') 'ECHAM-PDAF','ensemble size:', dim_ens
  IF (subtype /= 5) WRITE (*, '(a, 6x, a, i5)') 'ECHAM-PDAF','Assimilation interval:', delt_obs_atm
  WRITE (*, '(a, 10x, a, f5.2)') 'ECHAM-PDAF','forgetting factor:', forget
  WRITE (*, '(a, 6x, a, i9)') 'ECHAM-PDAF','ECHAM state dimension:',dim_state

END SUBROUTINE init_pdaf_info
