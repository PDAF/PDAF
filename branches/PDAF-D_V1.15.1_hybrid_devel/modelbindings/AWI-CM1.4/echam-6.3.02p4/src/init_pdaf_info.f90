!$Id: init_pdaf_info.f90 2135 2019-11-22 18:56:29Z lnerger $
!BOP
!
! !ROUTINE: init_pdaf_info - Screen output on assimilation configuration
!
! !INTERFACE:
SUBROUTINE init_pdaf_info()

! !DESCRIPTION:
! This routine performs a model-sided screen output about
! the coniguration of the data assimilation system.
! Using this output is optional. Most of the information
! is also displayed by PDAF itself when it is initialized
! in PDAF_init. Not displayed by PDAF is the assimilation
! interval (delt_obs_atm), which is unknown to PDAF.
!
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mod_assim_pdaf, & ! Variables for assimilation
       ONLY: filtertype, subtype, dim_ens,  &
       forget, delt_obs_ocn, delt_obs_atm, dim_state

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: init_pdaf
!EOP


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
     WRITE (*, '(a, 21x, a)') 'ECHAM-PDAF','Filter: LESTKF'
     IF (subtype == 0) THEN
        WRITE (*, '(a, 6x, a)') 'ECHAM-PDAF','-- Standard mode'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 6x, a)') 'ECHAM-PDAF','-- Offline mode'
     END IF
  ELSE IF (filtertype == 11) THEN
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
  WRITE (*, '(a, 8x, a, i9)') 'ECHAM-PDAF','ECHAM state dimension:',dim_state

END SUBROUTINE init_pdaf_info
