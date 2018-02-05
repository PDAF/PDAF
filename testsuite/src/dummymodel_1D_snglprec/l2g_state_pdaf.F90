!$Id: l2g_state_pdaf.F90 1711 2016-12-20 10:17:58Z lnerger $
!BOP
!
! !ROUTINE: l2g_state_pdaf --- Initialize full state from local analysis
!
! !INTERFACE:
SUBROUTINE l2g_state_pdaf(step, domain, dim_l, state_l, dim_p, state_p)

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
  INTEGER, INTENT(in) :: dim_l          ! Local state dimension
  INTEGER, INTENT(in) :: dim_p          ! PE-local full state dimension
  REAL, INTENT(in)    :: state_l(dim_l) ! State vector on local analysis domain
  REAL, INTENT(inout) :: state_p(dim_p) ! PE-local full state vector 

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update    (as U_l2g_state)
!EOP


! **************************************************
! *** Initialize elements of global state vector ***
! **************************************************

  ! Here simply the element with index DOMAIN is updated
  state_p(domain) = state_l(1)

END SUBROUTINE l2g_state_pdaf
