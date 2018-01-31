C     $Id$
C
C     !ROUTINE: g2l_obs_pdaf --- Restrict an obs. vector to local analysis domain
C
C     !INTERFACE:
      SUBROUTINE g2l_obs_pdaf(domain_p, step,
     &      dim_obs_f, dim_obs_l, mstate_f, mstate_l)

C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C     Used in the filters: LSEIK/LETKF/LESTKF
C
C     The routine is called during the analysis step
C     on each of the local analysis domains.
C     It has to restrict the full vector of all 
C     observations required for the loop of localized 
C     analyses on the PE-local domain to the current 
C     local analysis domain.
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !USES:
      USE mod_assimilation,
     &     ONLY: obs, obs_index_l

      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: domain_p   ! Current local analysis domain
      INTEGER, INTENT(in) :: step       ! Current time step
      INTEGER, INTENT(in) :: dim_obs_f  ! Dimension of full PE-local obs. vector
      INTEGER, INTENT(in) :: dim_obs_l  ! Local dimension of observation vector
      REAL*8, INTENT(in)    :: mstate_f(dim_obs_f)   ! Full PE-local obs. vector
      REAL*8, INTENT(out)   :: mstate_l(dim_obs_l)   ! Obs. vector on local domain

C     !CALLING SEQUENCE:
C     Called by: PDAF_lseik_analysis   (as U_g2l_obs)
C     Called by: PDAF_lestkf_analysis  (as U_g2l_obs)
C     Called by: PDAF_letkf_analysis   (as U_g2l_obs)

C     *** local variables ***
      INTEGER :: i              ! Counter


C *******************************************************
C *** Perform localization of some observation vector *** 
C *** to the current local analysis domain.           ***
C *******************************************************

      DO i = 1, dim_obs_l
         mstate_l(i) = mstate_f(obs_index_l(i))
      END DO

      END SUBROUTINE g2l_obs_pdaf
