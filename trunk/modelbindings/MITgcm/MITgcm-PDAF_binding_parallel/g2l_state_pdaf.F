C     $Id$
C
C     !ROUTINE: g2l_state_pdaf --- Restrict a model state to a local analysis domain
C
C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C     Used in the filters: LSEIK/LETKF/LESTKF/LNETF
C
C     The routine is called during the loop over all
C     local analysis domains in PDAF\_lseik\_update
C     before the analysis on a single local analysis
C     domain.  It has to project the full PE-local
C     model state onto the current local analysis
C     domain.
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !INTERFACE:
      SUBROUTINE g2l_state_pdaf(step, domain_p,
     &      dim_p, state_p, dim_l, state_l)

C     !USES:
      USE mod_assimilation,
     &     ONLY: id_lstate_in_pstate

      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: step           ! Current time step
      INTEGER, INTENT(in) :: domain_p       ! Current local analysis domain
      INTEGER, INTENT(in) :: dim_p          ! PE-local full state dimension
      INTEGER, INTENT(in) :: dim_l          ! Local state dimension
      REAL(kind=8), INTENT(in)    :: state_p(dim_p) ! PE-local full state vector
      REAL(kind=8), INTENT(out)   :: state_l(dim_l) ! State vector on local analysis domain

C     !CALLING SEQUENCE:
C     Called by: PDAF_lseik_update    (as U_g2l_state)
C     Called by: PDAF_letkf_update    (as U_g2l_state)
C     Called by: PDAF_lestkf_update   (as U_g2l_state)

C     *** local variables ***
      INTEGER :: i    ! Counter


C     *************************************
C     *** Initialize local state vector ***
C     *************************************

C     Generic initialization using ID_LSTATE_IN_PSTATE set in INIT_DIM_L_PDAF
      DO i = 1, dim_l
         state_l(i) = state_p(id_lstate_in_pstate(i))
      ENDDO

      END SUBROUTINE g2l_state_pdaf
