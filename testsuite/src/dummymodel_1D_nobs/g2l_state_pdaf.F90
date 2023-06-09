!$Id: g2l_state_pdaf.F90 1020 2010-07-14 10:06:25Z lnerger $
!BOP
!
! !ROUTINE: g2l_state_pdaf --- Restrict a model state to a local analysis domain
!
! !INTERFACE:
SUBROUTINE g2l_state_pdaf(step, domain, dim_p, state_p, dim_l, state_l)

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
! The routine is called by each filter process.
!
! Implementation for the dummy model with domain 
! decomposition. Here, we simply consider each single 
! grid point as a local analysis domain. In 3D this 
! could be a water column.
!
! !REVISION HISTORY:
! 2005-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step           ! Current time step
  INTEGER, INTENT(in) :: domain         ! Current local analysis domain
  INTEGER, INTENT(in) :: dim_p          ! PE-local full state dimension
  INTEGER, INTENT(in) :: dim_l          ! Local state dimension
  REAL, INTENT(in)    :: state_p(dim_p) ! PE-local full state vector 
  REAL, INTENT(out)   :: state_l(dim_l) ! State vector on local analysis domain

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_g2l_state)
!EOP

! *************************************
! *** Initialize local state vector ***
! *************************************
  
  ! Here simply the element with index DOMAIN
  state_l = state_p(domain)

END SUBROUTINE g2l_state_pdaf
