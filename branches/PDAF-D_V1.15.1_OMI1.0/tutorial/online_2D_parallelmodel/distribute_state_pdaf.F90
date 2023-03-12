!$Id: distribute_state_pdaf.F90 1589 2015-06-12 11:57:58Z lnerger $
!BOP
!
! !ROUTINE: distribute_state_pdaf --- Initialize model fields from state vector
!
! !INTERFACE:
SUBROUTINE distribute_state_pdaf(dim_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
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
       ONLY: nx_p, ny, field_p

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! PE-local state vector

! !CALLING SEQUENCE:
! Called by: PDAF_get_state      (as U_dist_state)
! Called by: PDAF_assimilate_X   (as U_coll_state)
!EOP

! *** local variables ***
  INTEGER :: j         ! Counters


! *******************************************
! *** Initialize model fields from state  ***
!********************************************

  DO j = 1, nx_p
     field_p(1:ny, j) = state_p(1 + (j-1)*ny : j*ny)
  END DO

END SUBROUTINE distribute_state_pdaf
