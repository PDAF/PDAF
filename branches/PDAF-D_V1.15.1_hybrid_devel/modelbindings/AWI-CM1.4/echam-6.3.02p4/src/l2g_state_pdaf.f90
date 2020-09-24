!$Id: l2g_state_pdaf.f90 2135 2019-11-22 18:56:29Z lnerger $
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
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, ONLY: mype_filter
  USE mo_kind_pdaf, ONLY: dp

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step           ! Current time step
  INTEGER, INTENT(in) :: domain         ! Current local analysis domain
  INTEGER, INTENT(in) :: dim_l          ! Local state dimension
  INTEGER, INTENT(in) :: dim_p          ! PE-local full state dimension
  REAL(dp), INTENT(in)    :: state_l(dim_l) ! State vector on local analysis domain
  REAL(dp), INTENT(inout) :: state_p(dim_p) ! PE-local full state vector 

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update    (as U_l2g_state)
! Called by: PDAF_lestkf_update   (as U_l2g_state)
! Called by: PDAF_letkf_update    (as U_l2g_state)
!EOP


! **************************************************
! *** Initialize elements of global state vector ***
! **************************************************

  if (mype_filter==0) write (*,*) 'ECHAM-PDAF TEMPLATE l2g_state_pdaf'

  state_p = 1.0

END SUBROUTINE l2g_state_pdaf
