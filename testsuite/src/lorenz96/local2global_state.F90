!$Id: local2global_state.F90 836 2010-01-29 17:19:57Z lnerger $
!BOP
!
! !ROUTINE: local2global_state --- Initialize full state from local analysis
!
! !INTERFACE:
SUBROUTINE local2global_state(step, domain, dim_l, state_l, dim, state)

! !DESCRIPTION:
! User-supplied routine for PDAF (LSEIK):
!
! The routine is called during the loop over all
! local analysis domains in PDAF\_lseik\_update 
! after the analysis and ensemble transformation 
! on a single local analysis domain. It has to 
! initialize elements of the PE-local full state 
! vector from the provided analysis state vector 
! on the local analysis domain.
!
! This variant is for the Lorenz96 model without
! parallelization. Here, we simply consider each single 
! grid point as a local analysis domain.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step           ! Current time step
  INTEGER, INTENT(in) :: domain         ! Current local analysis domain
  INTEGER, INTENT(in) :: dim_l          ! Local state dimension
  INTEGER, INTENT(in) :: dim            ! Full state dimension
  REAL, INTENT(in)    :: state_l(dim_l) ! State vector on local analysis domain
  REAL, INTENT(inout) :: state(dim)     ! Full state vector 

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update    (as U_l2g_state)
!EOP


! **************************************************
! *** Initialize elements of global state vector ***
! **************************************************

  ! Here simply the element with index DOMAIN is updated
  state(domain) = state_l(1)

END SUBROUTINE local2global_state
