!$Id: init_obsvar_l_pdaf.F90 233 2019-10-01 07:03:44Z lnerger $
!BOP
! !ROUTINE: init_obsvar_l_pdaf --- Get local mean observation error variance
!
! !INTERFACE:
SUBROUTINE init_obsvar_l_pdaf(domain_p, step, dim_obs_l, obs_l, meanvar_l)

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
! Implementation for the 2D online example
! with or without parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_obs_A_pdaf, &
       ONLY: assim_A, init_obsvar_l_A
  USE mod_obs_B_pdaf, &
       ONLY: assim_B, init_obsvar_l_B

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain_p      ! Current local analysis domain
  INTEGER, INTENT(in) :: step          ! Current time step
  INTEGER, INTENT(in) :: dim_obs_l     ! Local dimension of observation vector
  REAL, INTENT(in) :: obs_l(dim_obs_l) ! Local observation vector
  REAL, INTENT(out)   :: meanvar_l     ! Mean local observation error variance

! !CALLING SEQUENCE:
! Called by: PDAF_set_forget_local    (as U_init_obsvar_l)
!EOP

! *** Local variables
  INTEGER :: cnt_obs_l


! ***********************************
! *** Compute local mean variance ***
! ***********************************

  ! Initialize observation counter
  cnt_obs_l = 0

  IF (assim_A) CALL init_obsvar_l_A(meanvar_l, cnt_obs_l)
  IF (assim_B) CALL init_obsvar_l_B(meanvar_l, cnt_obs_l)

END SUBROUTINE init_obsvar_l_pdaf
