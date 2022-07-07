C     $Id$
C
C     !ROUTINE: l2g_state_pdaf --- Initialize full state from local analysis
C
C     !DESCRIPTION:
C     User-supplied routine for PDAF.
C     Used in the filters: LSEIK/LETKF/LESTKF/LNETF
C
C     The routine is called during the loop over all
C     local analysis domains in PDAF\_X\_update
C     after the analysis and ensemble transformation
C     on a single local analysis domain. It has to
C     initialize elements of the PE-local full state
C     vector from the provided analysis state vector
C     on the local analysis domain.
C
C     !REVISION HISTORY:
C     2013-02 - Lars Nerger - Initial code
C     Later revisions - see svn log
C
C     !INTERFACE:
      SUBROUTINE l2g_state_pdaf(step, domain_p,
     &      dim_l, state_l, dim_p, state_p)

C     !USES:
      USE mod_assimilation,
     &     ONLY: id_lstate_in_pstate

      IMPLICIT NONE

C     !ARGUMENTS:
      INTEGER, INTENT(in) :: step           ! Current time step
      INTEGER, INTENT(in) :: domain_p       ! Current local analysis domain
      INTEGER, INTENT(in) :: dim_l          ! Local state dimension
      INTEGER, INTENT(in) :: dim_p          ! PE-local full state dimension
      REAL(kind=8), INTENT(in)    :: state_l(dim_l) ! State vector on local analysis domain
      REAL(kind=8), INTENT(inout) :: state_p(dim_p) ! PE-local full state vector 

C     !CALLING SEQUENCE:
C     Called by: PDAF_lestkf_update   (as U_l2g_state)
C     Called by: PDAF_letkf_update    (as U_l2g_state)

C     *** local variables ***
      INTEGER :: i    ! Counter


C     **************************************************
C     *** Initialize elements of global state vector ***
C     **************************************************

C     Generic initialization using ID_LSTATE_IN_PSTATE set in INIT_DIM_L_PDAF
      DO i = 1, dim_l
         state_p(id_lstate_in_pstate(i)) = state_l(i)
      END DO

      END SUBROUTINE l2g_state_pdaf
