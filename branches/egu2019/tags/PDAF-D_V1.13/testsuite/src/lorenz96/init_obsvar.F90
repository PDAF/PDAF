!$Id: init_obsvar.F90 1253 2012-01-30 18:56:08Z lnerger $
!BOP
!
! !ROUTINE: init_obsvar --- Get mean observation error variance
!
! !INTERFACE:
SUBROUTINE init_obsvar(step, dim_obs, obs, meanvar)

! !DESCRIPTION:
! User-supplied routine for PDAF (SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF)
! with adaptive forgetting factor. This routine will
! only be called, if the adaptive forgetting factor
! feature is used. Please note that this is an 
! experimental feature.
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
! This variant is for the Lorenz96 mode without
! parallelization.  We assume a diagonal observation
! error covariance matrix with constant variances. 
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: rms_obs

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step       ! Current time step
  INTEGER, INTENT(in) :: dim_obs    ! PE-local dimension of observation vector
  REAL, INTENT(in) :: obs(dim_obs)  ! PE-local observation vector
  REAL, INTENT(out)   :: meanvar    ! Mean observation error variance

! !CALLING SEQUENCE:
! Called by: PDAF_set_forget    (as U_init_init_obs_covar)
!EOP


! *****************************
! *** Compute mean variance ***
! *****************************

  ! Here the mean variance is simply the 
  ! error variance of each single observation.

  meanvar = rms_obs ** 2

END SUBROUTINE init_obsvar
