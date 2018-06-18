!$Id: global2local_state.F90 836 2010-01-29 17:19:57Z lnerger $
!BOP
!
! !ROUTINE: global2local_state --- Restrict a model state to a local analysis domain
!
! !INTERFACE:
SUBROUTINE global2local_state(step, domain, dim, state, dim_l, state_l)

! !DESCRIPTION:
! User-supplied routine for PDAF (LSEIK):
!
! The routine is called during the loop over all
! local analysis domains in PDAF\_lseik\_update
! before the analysis on a single local analysis 
! domain.  It has to project the full PE-local 
! model state onto the current local analysis 
! domain.
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
  INTEGER, INTENT(in) :: dim            ! Full state dimension
  INTEGER, INTENT(in) :: dim_l          ! Local state dimension
  REAL, INTENT(in)    :: state(dim)     ! Full state vector 
  REAL, INTENT(out)   :: state_l(dim_l) ! State vector on local analysis domain

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_g2l_state)
!EOP

! *************************************
! *** Initialize local state vector ***
! *************************************
  
  ! Here simply the element with index DOMAIN
  state_l = state(domain)

END SUBROUTINE global2local_state
