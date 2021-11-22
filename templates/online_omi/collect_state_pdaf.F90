!$Id$
!>  Initialize state vector from model fields
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all filters.
!!
!! This subroutine is called during the forecast 
!! phase from PDAF_put_state_X or PDAF_assimilate_X
!! after the propagation of each ensemble member. 
!! The supplied state vector has to be initialized
!! from the model fields (typically via a module). 
!! With parallelization, MPI communication might be 
!! required to initialize state vectors for all 
!! subdomains on the model PEs. 
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! __Revision history:__
!! * 2004-11 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE collect_state_pdaf(dim_p, state_p)

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  !< local state vector
  
! *** local variables ***


! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE collect_state_pdaf.F90: Implement initialization of state vector here!'

!   state_p = ????
  
END SUBROUTINE collect_state_pdaf
