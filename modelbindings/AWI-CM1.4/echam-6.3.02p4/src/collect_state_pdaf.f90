!$Id: collect_state_pdaf.f90 2135 2019-11-22 18:56:29Z lnerger $
!BOP
!
! !ROUTINE: collect_state_pdaf --- Initialize state vector from model fields
!
! !INTERFACE:
SUBROUTINE collect_state_pdaf(dim_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/EnKF/SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! This subroutine is called during the forecast 
! phase from PDAF\_put\_state\_X after the 
! propagation of each ensemble member. 
! The supplied state vector has to be initialized
! from the model fields (typically via a module). 
! With parallelization, MPI communication might be 
! required to initialize state vectors for all 
! subdomains on the model PEs. 
!
! The routine is executed by each process that is
! participating in the model integrations.
!
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, ONLY: mype_world
  USE mo_kind_pdaf, ONLY: dp

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL(dp), INTENT(inout) :: state_p(dim_p)  ! local state vector

! !CALLING SEQUENCE:
! Called by: PDAF_put_state_X   (as U_coll_state)
!EOP
  

! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

  if (mype_world==384) write (*,*) 'ECHAM-PDAF TEMPLATE collect_state_pdaf'

!   state_p = ????
  
END SUBROUTINE collect_state_pdaf
