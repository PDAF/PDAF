C     $Id$
C!BOP
C
C     !ROUTINE: init_obsvar_pdaf --- Get mean observation error variance
C
C     !INTERFACE:
      SUBROUTINE init_obsvar_pdaf(step, dim_obs_p,
     &      obs_p, meanvar)

C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C     Used in the filters: SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF

C     This routine will only be called, if the adaptive
C     forgetting factor feature is used. Please note that
C     this is an experimental feature.
C
C     The routine is called in global filters (like SEIK)
C     during the analysis or in local filters (e.g. LSEIK)
C     before the loop over local analysis domains 
C     by the routine PDAF\_set\_forget that estimates an 
C     adaptive forgetting factor.  The routine has to 
C     initialize the mean observation error variance.  
C     For global filters this should be the global mean,
C     while for local filters it should be the mean for the
C     PE-local  sub-domain.  (See init\_obsvar\_l_pdaf()
C     for a localized variant for local filters.)
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !USES:
      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: step             ! Current time step
      INTEGER, INTENT(in) :: dim_obs_p        ! PE-local dimension of observation vector
      REAL*8, INTENT(in)  :: obs_p(dim_obs_p) ! PE-local observation vector
      REAL*8, INTENT(out) :: meanvar          ! Mean observation error variance

C     !CALLING SEQUENCE:
C     Called by: PDAF_set_forget    (as U_init_init_obs_covar)
C!EOP


C     *****************************
C     *** Compute mean variance ***
C     *****************************

      WRITE (*,*) 'TEMPLATE init_obsvar_pdaf.F90:
     &      Set mean observation variance here '

!ML this is clearly a dummy value
      meanvar = 1.0

      END SUBROUTINE init_obsvar_pdaf
