!$Id: collect_state.F90 814 2010-01-27 11:30:16Z lnerger $
!BOP
!
! !ROUTINE: collect_state --- Initialize state vector from model fields
!
! !INTERFACE:
SUBROUTINE collect_state(dim, state)

! !DESCRIPTION:
! User-supplied routine for PDAF (all filters):
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
! This variant is for the Lorenz96 model without
! parallelization. Here, the state vector
! and the model field are identical. Here, the 
! model field is just copied to the state vector.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &
       ONLY: x

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim           ! PE-local state dimension
  REAL, INTENT(inout) :: state(dim)  ! local state vector

! !CALLING SEQUENCE:
! Called by: PDAF_put_state_seek   (as U_coll_state)
! Called by: PDAF_put_state_seik   (as U_coll_state)
! Called by: PDAF_put_state_enkf   (as U_coll_state)
! Called by: PDAF_put_state_lseik   (as U_coll_state)
!EOP
  

! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

  state(:) = x(:)
  
END SUBROUTINE collect_state
