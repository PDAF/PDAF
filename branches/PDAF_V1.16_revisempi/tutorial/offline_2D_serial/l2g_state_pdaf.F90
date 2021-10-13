!$Id: l2g_state_pdaf.F90 1369 2013-04-24 16:38:17Z lnerger $
!BOP
!
! !ROUTINE: l2g_state_pdaf --- Initialize full state from local analysis
!
! !INTERFACE:
SUBROUTINE l2g_state_pdaf(step, domain_p, dim_l, state_l, dim_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the loop over all
! local analysis domains in PDAF\_X\_update 
! after the analysis and ensemble transformation 
! on a single local analysis domain. It has to 
! initialize elements of the PE-local full state 
! vector from the provided analysis state vector 
! on the local analysis domain.
!
! Generic implementation using index vector 
! ID_LSTATE_IN_PSTATE.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: id_lstate_in_pstate

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step           ! Current time step
  INTEGER, INTENT(in) :: domain_p       ! Current local analysis domain
  INTEGER, INTENT(in) :: dim_l          ! Local state dimension
  INTEGER, INTENT(in) :: dim_p          ! PE-local full state dimension
  REAL, INTENT(in)    :: state_l(dim_l) ! State vector on local analysis domain
  REAL, INTENT(inout) :: state_p(dim_p) ! PE-local full state vector 

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update    (as U_l2g_state)
! Called by: PDAF_lestkf_update   (as U_l2g_state)
! Called by: PDAF_letkf_update    (as U_l2g_state)
!EOP

! *** local variables ***
  INTEGER :: i                          ! Counter


! **************************************************
! *** Initialize elements of global state vector ***
! **************************************************

  ! Generic initialization using ID_LSTATE_IN_PSTATE set in INIT_DIM_L_PDAF
  DO i = 1, dim_l
     state_p(id_lstate_in_pstate(i)) = state_l(i)
  END DO

END SUBROUTINE l2g_state_pdaf
