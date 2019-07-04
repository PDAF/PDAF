!$Id: init_dim_obs_f_pdaf.F90 1249 2012-01-27 13:44:54Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_f_pdaf --- Set full dimension of observation
!
! !INTERFACE:
SUBROUTINE init_dim_obs_f_pdaf(step, dim_obs_f)

! !DESCRIPTION:
! User supplied routine for PDAF (LSEIK):
!
! The routine is called in PDAF\_lseik\_update 
! at the beginning of the analysis step before 
! the loop through all local analysis domains. 
! It has to determine the dimension of the 
! observation vector according to the current 
! time step for all observations required for 
! the analyses in the loop over all local 
! analysis domains on the PE-local state domain.
!
! The routine is called by each filter process.
!
! Implementation for the dummy model with domain 
! decomposition. For simplicity the full 
! observation vector holds a global state.
!
! !REVISION HISTORY:
! 2005-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &
       ONLY: local_dims, dim_state
  USE mod_assimilation, &
       ONLY: dim_obs
  USE mod_parallel, &
       ONLY: mype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step      ! Current time step
  INTEGER, INTENT(out) :: dim_obs_f ! Dimension of full observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_obs)
!EOP


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! dimension for PE-local domain
  dim_obs_f = dim_obs

END SUBROUTINE init_dim_obs_f_pdaf

