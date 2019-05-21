C     $Id$
C!BOP
C     !ROUTINE: init_obsvar_l_pdaf --- Get local mean observation error variance
C
C     !INTERFACE:
      SUBROUTINE init_obsvar_l_pdaf(domain_p, step,
     &      dim_obs_l, obs_l, meanvar_l)

C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C     Used in the filters: LSEIK/LETKF/LESTKF
C
C     This routine will only be called, if the 
C     local adaptive forgetting factor feature 
C     is used. Please note that this is an 
C     experimental feature.
C
C     The routine is called in the loop over all
C     local analysis domains during each analysis
C     by the routine PDAF\_set\_forget\_local that 
C     estimates a local adaptive forgetting factor.
C     The routine has to initialize the mean observation 
C     error variance for the current local analysis 
C     domain.  (See init_obsvar() for a global variant.)
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !USES:
      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: domain_p         ! Current local analysis domain
      INTEGER, INTENT(in) :: step             ! Current time step
      INTEGER, INTENT(in) :: dim_obs_l        ! Local dimension of observation vector
      REAL*8, INTENT(in)  :: obs_l(dim_obs_l) ! Local observation vector
      REAL*8, INTENT(out) :: meanvar_l        ! Mean local observation error variance

C     !CALLING SEQUENCE:
C     Called by: PDAF_set_forget_local    (as U_init_obsvar_l)
C!EOP


C ***********************************
C *** Compute local mean variance ***
C ***********************************

      WRITE (*,*) 'TEMPLATE init_obsvar_l_pdaf.F90:
     &      Set mean observation variance here '

C     this is clearly a dummy value
      meanvar_l = 1.0

      END SUBROUTINE init_obsvar_l_pdaf
