!$Id: collect_state_dummy3d_D.F90 783 2009-12-07 10:28:43Z lnerger $
!BOP
!
! !ROUTINE: collect_state --- Initialize state vector from model fields
!
! !INTERFACE:
SUBROUTINE collect_state(dim_p, state_p)

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
! For the domain-decomposed dummy model with
! domain-decomposed PDAF, the model fields are just 
! copied to the state vector on all model processes.
!
! !REVISION HISTORY:
! 2004-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &
       ONLY: dim_l, field

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! local state vector

! !CALLING SEQUENCE:
! Called by: PDAF_put_state_seek   (as U_coll_state)
! Called by: PDAF_put_state_seik   (as U_coll_state)
! Called by: PDAF_put_state_enkf   (as U_coll_state)
! Called by: PDAF_put_state_lseik   (as U_coll_state)
!EOP

! *** Local variables ***
  INTEGER :: i, j, k, cnt                ! Counters
  

! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

  cnt = 1
  DO k = 1, dim_l(3)
     DO j = 1, dim_l(2)
        DO i = 1, dim_l(1)
           state_p(cnt) = field(i, j, k)
           cnt = cnt + 1
        END DO
     END DO
  END DO
  
END SUBROUTINE collect_state
