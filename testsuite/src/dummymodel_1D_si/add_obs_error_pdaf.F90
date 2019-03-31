!$Id: add_obs_error_pdaf.F90 1253 2012-01-30 18:56:08Z lnerger $
!BOP
!
! !ROUTINE: add_obs_error_pdaf --- Add observation error covariance matrix
!
! !INTERFACE:
SUBROUTINE add_obs_error_pdaf(step, dim_obs_p, C_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: EnKF
!
! The routine is called during the analysis step
! by PDAF\_enkf\_analysis_X (X=rlm or rsm).  It 
! has to add the observation error covariance 
! matrix to the provided matrix C_p for the 
! PE-local domain .
! 
! The routine is called by all filter processes.
!
! Version for dummy model with domain decomposition
! with diagonal observation erorr covariance matrix
! and constant variance.
!
! !REVISION HISTORY:
! 2004-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: rms_obs

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step       ! Current time step
  INTEGER, INTENT(in) :: dim_obs_p  ! Dimension of observation vector
  REAL, INTENT(inout) :: C_p(dim_obs_p,dim_obs_p) ! Matrix to that
                                    ! observation covariance R is added

! !CALLING SEQUENCE:
! Called by: PDAF_enkf_analysis_rlm   (as U_add_obs_err)
! Called by: PDAF_enkf_analysis_rsm   (as U_add_obs_err)
!EOP


! *** local variables ***
  INTEGER :: i          ! index of observation component
  REAL :: variance_obs  ! variance of observations


! **********************
! *** INITIALIZATION ***
! **********************

    variance_obs = rms_obs ** 2


! *************************************
! ***   Add observation error       ***
! ***                               ***
! *** Measurements are uncorrelated ***
! *** here, thus R is diagonal      ***
! *************************************

  DO i = 1, dim_obs_p
     C_p(i, i) = C_p(i, i) + variance_obs
  ENDDO

END SUBROUTINE add_obs_error_pdaf
