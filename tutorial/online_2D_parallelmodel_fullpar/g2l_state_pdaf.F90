!$Id: g2l_state_pdaf.F90 1411 2013-09-25 14:04:41Z lnerger $
!BOP
!
! !ROUTINE: g2l_state_pdaf --- Restrict a model state to a local analysis domain
!
! !INTERFACE:
SUBROUTINE g2l_state_pdaf(step, domain_p, dim_p, state_p, dim_l, state_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the loop over all
! local analysis domains in PDAF\_lseik\_update
! before the analysis on a single local analysis 
! domain.  It has to project the full PE-local 
! model state onto the current local analysis 
! domain.
!
! Implementation for the 2D offline example
! with or without parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step           ! Current time step
  INTEGER, INTENT(in) :: domain_p       ! Current local analysis domain
  INTEGER, INTENT(in) :: dim_p          ! PE-local full state dimension
  INTEGER, INTENT(in) :: dim_l          ! Local state dimension
  REAL, INTENT(in)    :: state_p(dim_p) ! PE-local full state vector 
  REAL, INTENT(out)   :: state_l(dim_l) ! State vector on local analysis domain

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update    (as U_g2l_state)
! Called by: PDAF_letkf_update    (as U_g2l_state)
! Called by: PDAF_lestkf_update   (as U_g2l_state)
!EOP


! *************************************
! *** Initialize local state vector ***
! *************************************
  
  state_l = state_p(domain_p)

END SUBROUTINE g2l_state_pdaf
