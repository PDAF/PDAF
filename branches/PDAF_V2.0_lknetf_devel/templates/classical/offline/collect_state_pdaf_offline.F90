!$Id$
!BOP
!
! !ROUTINE: collect_state_pdaf --- Initialize state vector from model fields
!
! !INTERFACE:
SUBROUTINE collect_state_pdaf(dim_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! For the offline mode of PDAF this routine only needs
! to exist for linking. It is never called.
!
! !REVISION HISTORY:
! 2008-07 - Lars Nerger - Initial code based on online implementation
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! local state vector

! !CALLING SEQUENCE:
! Called by: PDAF_put_state_X   (as U_coll_state)
!EOP
  

! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

  ! Nothing to be done in offline mode.

  
END SUBROUTINE collect_state_pdaf
