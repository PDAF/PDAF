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
!! This variant is for the setup that the analysis
!! step is computed on a single process separate from
!! the processes that integrate the ensemble. In this
!! case the domain-decomposed states need to be 
!! distributed according to the domain-decomposition.
!!
!! __Revision history:__
!! * 2004-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE distribute_state_pdaf(dim_p, state_p)

  USE mpi                      ! MPI
  USE mod_model, &             ! Model variables
       ONLY: nx_p, ny, field_p
  USE mod_parallel_model, &    ! Model parallelization variables
       ONLY: COMM_model, mype_model, npes_model
  USE mod_parallel_pdaf, &     ! PDAF parallelization variables
       ONLY: MPIstatus, MPIerr

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  !< PE-local state vector

! *** local variables ***
  INTEGER :: j, rank      ! Counters
  INTEGER :: offset       ! Offset in state vector
  REAL, ALLOCATABLE :: state_part(:) ! local state vector


! *******************************************
! *** Initialize model fields from state  ***
!********************************************

  ! Distribute field on sub-domains to model processes
  IF (mype_model==0) THEN

     ! Directly initialize field on process 0
     DO j = 1, nx_p
        field_p(1:ny, j) = state_p(1 + (j-1)*ny : j*ny)
     END DO

     offset = nx_p*ny
     DO rank = 1, npes_model-1
        CALL MPI_Send(state_p(1+offset : nx_p*ny+offset), nx_p*ny, MPI_DOUBLE_PRECISION, &
             rank, rank, COMM_model, MPIerr)
        offset = offset + nx_p*ny
     END DO
  ELSE
     ALLOCATE(state_part(nx_p*ny))

     CALL MPI_Recv(state_part, nx_p*ny, MPI_DOUBLE_PRECISION, 0, mype_model, &
          COMM_model, MPIstatus, MPIerr)

     ! Initialize field from received state vector on sub-domain
     DO j = 1, nx_p
        field_p(1:ny, j) = state_part(1 + (j-1)*ny : j*ny)
     END DO

     DEALLOCATE(state_part)
  END IF

END SUBROUTINE distribute_state_pdaf
