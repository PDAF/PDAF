!$Id: distribute_state_dummy3d_D.F90 783 2009-12-07 10:28:43Z lnerger $
!BOP
!
! !ROUTINE: distribute_state --- Initialize model fields from state vector
!
! !INTERFACE:
SUBROUTINE distribute_state(dim_p, state_p)

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
! For the dummy model and PDAF with domain
! decomposition the state vector and the model
! field are identical. Hence, the field array 
! is directly initialized from an ensemble 
! state vector by each model PE.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &
       ONLY: dim_l, field

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! PE-local state vector

! !CALLING SEQUENCE:
! Called by: PDAF_get_state   (as U_dist_state)
!EOP

! *** Local variables ***
  INTEGER :: i, j, k, cnt                ! Counters


! *******************************************
! *** Initialize model fields from state  ***
! *** Each model PE 0 knows his sub-state ***
!********************************************

  cnt = 1
  DO k = 1, dim_l(3)
     DO j = 1, dim_l(2)
        DO i = 1, dim_l(1)
           field(i, j, k) = state_p(cnt)
           cnt = cnt + 1
        END DO
     END DO
  END DO


! *******************************
! *** distribute state fields ***
!********************************

  ! Nothing to be done here

END SUBROUTINE distribute_state
