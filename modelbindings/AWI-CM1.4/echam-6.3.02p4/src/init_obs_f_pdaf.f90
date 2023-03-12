!$Id: init_obs_f_pdaf.f90 2135 2019-11-22 18:56:29Z lnerger $
!BOP
!
! !ROUTINE: init_obs_f_pdaf --- Initialize observation vector
!
! !INTERFACE:
SUBROUTINE init_obs_f_pdaf(step, dim_obs_f, observation_f)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called in PDAF\_lseik\_update
! before the loop over all local analysis domains
! is entered. It has to provide the full observation 
! vector according to current time step (where 'full' 
! means 'all observations required for the localized 
! analysis on the PE-local domain).  This routine 
! is only used for LSEIK if a globally adaptive 
! forgetting factor is requested, rather than an 
! individual forgetting factor for each analysis 
! domain. This routine has to be implemented 
! consistently with the routines for the full 
! observation dimension and the full observation 
! operator. The forgetting factor will only be 
! globally adaptive, if the full observation vector 
! is the global observation vector.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, ONLY: mype_filter
  USE mo_kind_pdaf, ONLY: dp

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
  INTEGER, INTENT(in) :: dim_obs_f   ! Dimension of full observation vector
  REAL(dp), INTENT(out)   :: observation_f(dim_obs_f) ! Full observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_obs)
! Called by: PDAF_lestkf_update  (as U_init_obs)
! Called by: PDAF_letkf_update   (as U_init_obs)
!EOP


! ******************************************
! *** Initialize full observation vector ***
! ******************************************

  if (mype_filter==0) write (*,*) 'ECHAM-PDAF TEMPLATE init_obs_f_pdaf'

  observation_f = 1.0

END SUBROUTINE init_obs_f_pdaf

