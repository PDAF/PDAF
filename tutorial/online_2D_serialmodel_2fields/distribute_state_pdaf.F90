!>  Initialize model fields from state vector
!!
!! User-supplied call-back routine for PDAF.
!!
!! During the forecast phase of the filter this subroutine
!! is called from PDAF_init_forecast or PDAF3_assimilate.
!! supplying a model state, which has to be evolved. 
!! The routine has to initialize the fields of the 
!! model (typically available through a module) from 
!! the state vector of PDAF. With parallelization, 
!! MPI communication might be required to 
!! initialize all subdomains on the model PEs.
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! __Revision history:__
!! * 2004-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE distribute_state_pdaf(dim_p, state_p)

  USE mod_model, &             ! Model variables
       ONLY: nx, ny, field, fieldB
  USE mod_assimilation, &      ! Assimilation variables
       ONLY: fields, id

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  !< PE-local state vector

! *** local variables ***
  INTEGER :: j         ! Counters


! *******************************************
! *** Initialize model fields from state  ***
!********************************************

  ! + For the 2D tutorial model the state vector and
  ! + the model field are identical. Hence, the field
  ! + array is directly initialized from an ensemble 
  ! + state vector by each model PE.

  ! Field
  DO j = 1, nx
     field(1:ny, j) = state_p(fields(id%fieldA)%off + 1 + (j-1)*ny : fields(id%fieldA)%off + j*ny)
  END DO

  ! FieldB
  DO j = 1, nx
     fieldB(1:ny, j) = state_p(fields(id%fieldB)%off + 1 + (j-1)*ny : fields(id%fieldB)%off + j*ny)
  END DO

END SUBROUTINE distribute_state_pdaf
