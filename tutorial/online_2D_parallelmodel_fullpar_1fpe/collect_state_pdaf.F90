!$Id: collect_state_pdaf.F90 1411 2013-09-25 14:04:41Z lnerger $
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
!! For the 2D tutorial model the state vector and
!! the model field are identical. Hence, state vector
!! directly initialized from the model field by
!! each model PE.
!!
!! This variant is for the setup that the analysis
!! step is computed on a single process separate from
!! the processes that integrate the ensemble. In this
!! case the domain-decomposed states need to be 
!! collected into global model fields.
!!
!! __Revision history:__
!! * 2013-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE collect_state_pdaf(dim_p, state_p)

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
  REAL, INTENT(inout) :: state_p(dim_p)  !< local state vector

! *** local variables ***
  INTEGER :: j, rank      ! Counters
  INTEGER :: offset       ! Offset in state vector
  REAL, ALLOCATABLE :: state_part(:) ! local state vector


! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

  ! Collect field on sub-domains on model process 0
  IF (mype_model==0) THEN

     ! Directly initialize local state vector part from field
     DO j = 1, nx_p
        state_p(1 + (j-1)*ny : j*ny) = field_p(1:ny, j)
     END DO

     ! Receive state vector sub-domains from other processes
     offset = nx_p*ny
     DO rank = 1, npes_model-1
        CALL MPI_Recv(state_p(1+offset : nx_p*ny+offset), nx_p*ny, MPI_DOUBLE_PRECISION, &
             rank, rank, COMM_model, MPIstatus, MPIerr)
        offset = offset + nx_p*ny
     END DO
  ELSE
     ALLOCATE(state_part(nx_p*ny))

     ! Initialize local state vector part from field
     DO j = 1, nx_p
        state_part(1 + (j-1)*ny : j*ny) = field_p(1:ny, j)
     END DO

     ! Send state vector on sub-domain
     CALL MPI_Send(state_part, nx_p*ny, MPI_DOUBLE_PRECISION, &
          0, mype_model, COMM_model, MPIerr)

     DEALLOCATE(state_part)
  END IF

END SUBROUTINE collect_state_pdaf
