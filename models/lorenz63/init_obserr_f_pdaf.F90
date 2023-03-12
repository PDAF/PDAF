!$Id$
!BOP
!
! !ROUTINE: init_obserr_f_pdaf --- Initialize vector of observation errors
!
! !INTERFACE:
SUBROUTINE init_obserr_f_pdaf(step, dim_obs_f, obs_f, obserr_f)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Using for generating observations. The routine 
! provides a vector of observation error standard deviations
! (rms errors) to PDAF to perturb an observed model
! state.
!
! This variant is for the Lorenz63 model without
! parallelization.  We assume a diagonal observation
! error covariance matrix with constant variances. 
!
! !REVISION HISTORY:
! 2019-01 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: rms_obs

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                ! Current time step
  INTEGER, INTENT(in) :: dim_obs_f           ! Full dimension of observation vector
  REAL, INTENT(in)    :: obs_f(dim_obs_f)    ! Full observation vector
  REAL, INTENT(out)   :: obserr_f(dim_obs_f) ! Obervation error stddev

! !CALLING SEQUENCE:
! Called by: PDAF_gen_obs    (as U_init_obserr_f)
!EOP


! ****************************************
! *** Initialize vector of obs. errors ***
! ****************************************

  ! Here we simply use a constant error
  obserr_f = rms_obs
  
END SUBROUTINE init_obserr_f_pdaf
