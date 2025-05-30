!$Id$
!BOP
!
! !ROUTINE: add_obs_error_pdaf --- Add observation error covariance matrix
!
! !INTERFACE:
SUBROUTINE add_obs_error_pdaf(step, dim_obs, C_p)

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
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
!   USE mod_assimilation, &
!        ONLY: rms_obs

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step       ! Current time step
  INTEGER, INTENT(in) :: dim_obs    ! Dimension of observation vector
  REAL, INTENT(inout) :: C_p(dim_obs,dim_obs) ! Matrix to that
                                    ! observation covariance R is added

! !CALLING SEQUENCE:
! Called by: PDAF_enkf_analysis_rlm   (as U_add_obs_err)
! Called by: PDAF_enkf_analysis_rsm   (as U_add_obs_err)
!EOP


! *** local variables ***
!   INTEGER :: i          ! index of observation component
!   REAL :: variance_obs  ! variance of observations


! **********************
! *** INITIALIZATION ***
! **********************


! *************************************
! ***   Add observation error       ***
! *************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE add_obs_error_pdaf.F90: Implement addition of observation error here!'

! C_p = C_p + ?

END SUBROUTINE add_obs_error_pdaf
