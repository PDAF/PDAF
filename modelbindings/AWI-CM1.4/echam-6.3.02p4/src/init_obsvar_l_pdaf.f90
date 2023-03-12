!$Id: init_obsvar_l_pdaf.f90 2135 2019-11-22 18:56:29Z lnerger $
!BOP
! !ROUTINE: init_obsvar_l_pdaf --- Get local mean observation error variance
!
! !INTERFACE:
SUBROUTINE init_obsvar_l_pdaf(domain, step, dim_obs_l, obs_l, meanvar_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! This routine will only be called, if the 
! local adaptive forgetting factor feature 
! is used. Please note that this is an 
! experimental feature.
!
! The routine is called in the loop over all
! local analysis domains during each analysis
! by the routine PDAF\_set\_forget\_local that 
! estimates a local adaptive forgetting factor.
! The routine has to initialize the mean observation 
! error variance for the current local analysis 
! domain.  (See init_obsvar() for a global variant.)
!
! The routine is executed by all filter processes.
!
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, ONLY: mype_world
  USE mo_kind_pdaf, ONLY: dp

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain        ! Current local analysis domain
  INTEGER, INTENT(in) :: step          ! Current time step
  INTEGER, INTENT(in) :: dim_obs_l     ! Local dimension of observation vector
  REAL(dp), INTENT(in) :: obs_l(dim_obs_l) ! Local observation vector
  REAL(dp), INTENT(out)   :: meanvar_l     ! Mean local observation error variance

! !CALLING SEQUENCE:
! Called by: PDAF_set_forget_local    (as U_init_obsvar_l)
!EOP


! ***********************************
! *** Compute local mean variance ***
! ***********************************

  if (mype_world==768) write (*,*) 'ECHAM-PDAF TEMPLATE init_obsvar_l_pdaf'

  meanvar_l = 1.0

END SUBROUTINE init_obsvar_l_pdaf
