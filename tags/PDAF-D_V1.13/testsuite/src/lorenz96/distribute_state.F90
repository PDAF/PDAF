!$Id: distribute_state.F90 814 2010-01-27 11:30:16Z lnerger $
!BOP
!
! !ROUTINE: distribute_state --- Initialize model fields from state vector
!
! !INTERFACE:
SUBROUTINE distribute_state(dim, state)

! !DESCRIPTION:
! User-supplied routine for PDAF (all filters):
!
! During the forecast phase of the filter this
! subroutine is called from PDAF\_get\_state
! supplying a model state which has to be evolved. 
! The routine has to initialize the fields of the 
! model (typically available through a module) from 
! the state vector of PDAF. With parallelization, 
! MPI communication might be required to 
! initialize all subdomains on the model PEs.
!
! The routine is executed by each process that is
! participating in the model integrations.
!
! This variant is for the Lorenz96 model without
! parallelization. Here, the state vector
! and the model field are identical. Hence, the 
! field array is directly initialized from an 
! ensemble state vector.
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
  INTEGER, INTENT(in) :: dim         ! State dimension
  REAL, INTENT(inout) :: state(dim)  ! State vector

! !CALLING SEQUENCE:
! Called by: PDAF_get_state   (as U_dist_state)
!EOP


! *******************************************
! *** Initialize model fields from state  ***
!********************************************

  x(:) = state(:)


! *******************************
! *** distribute state fields ***
!********************************

  ! Nothing to be done here

END SUBROUTINE distribute_state
