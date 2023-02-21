!$Id: collect_state_pdaf_offline.F90 1369 2013-04-24 16:38:17Z lnerger $
!>  Initialize state vector from model fields
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all filters.
!!
!! For the offline mode of PDAF this routine only
!! needs to exist for linking. It is never called.
!!
!! __Revision history:__
!! * 2008-07 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE collect_state_pdaf(dim_p, state_p)

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  !< local state vector
  

! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

  ! Nothing to be done in offline mode.

  
END SUBROUTINE collect_state_pdaf
