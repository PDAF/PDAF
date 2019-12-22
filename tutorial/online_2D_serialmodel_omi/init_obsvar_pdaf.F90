!$Id$
!BOP
!
! !ROUTINE: init_obsvar_pdaf --- Get mean observation error variance
!
! !INTERFACE:
SUBROUTINE init_obsvar_pdaf(step, dim_obs_p, obs_p, meanvar)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF

! This routine will only be called, if the adaptive
! forgetting factor feature is used. Please note that
! this is an experimental feature.
!
! The routine is called in global filters (like SEIK)
! during the analysis or in local filters (e.g. LSEIK)
! before the loop over local analysis domains 
! by the routine PDAF\_set\_forget that estimates an 
! adaptive forgetting factor.  The routine has to 
! initialize the mean observation error variance.  
! For global filters this should be the global mean,
! while for local filters it should be the mean for the
! PE-local  sub-domain.  (See init\_obsvar\_l_pdaf()
! for a localized variant for local filters.)
!
! Implementation for the 2D online example
! with or without parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_obs_A_pdaf, &
       ONLY: assim_A, init_obsvar_A
  USE mod_obs_B_pdaf, &
       ONLY: assim_B, init_obsvar_B
  USE mod_obs_C_pdaf, &
       ONLY: assim_C, init_obsvar_C

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step          ! Current time step
  INTEGER, INTENT(in) :: dim_obs_p     ! PE-local dimension of observation vector
  REAL, INTENT(in) :: obs_p(dim_obs_p) ! PE-local observation vector
  REAL, INTENT(out)   :: meanvar       ! Mean observation error variance

! !CALLING SEQUENCE:
! Called by: PDAF_set_forget    (as U_init_init_obs_covar)
!EOP

! *** Local variables
  INTEGER :: cnt_obs_f


! *****************************
! *** Compute mean variance ***
! *****************************

  ! Initialize observation counter
  cnt_obs_f = 0

  IF (assim_A) CALL init_obsvar_A(meanvar, cnt_obs_f)
  IF (assim_B) CALL init_obsvar_B(meanvar, cnt_obs_f)
  IF (assim_C) CALL init_obsvar_C(meanvar, cnt_obs_f)

END SUBROUTINE init_obsvar_pdaf
