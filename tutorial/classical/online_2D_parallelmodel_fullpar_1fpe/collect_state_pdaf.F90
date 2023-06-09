!$Id: collect_state_pdaf.F90 1411 2013-09-25 14:04:41Z lnerger $
!BOP
!
! !ROUTINE: collect_state_pdaf --- Initialize state vector from model fields
!
! !INTERFACE:
SUBROUTINE collect_state_pdaf(dim_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! This subroutine is called during the forecast 
! phase from PDAF\_put\_state\_X or PDAF\_assimilate\_X
! after the propagation of each ensemble member. 
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
! 2013-09 - Lars Nerger - Initial code based on online implementation
! Later revisions - see svn log
!
! !USES:
  USE mpi
  USE mod_model, &
       ONLY: nx_p, ny, field_p
  USE mod_parallel_model, &
       ONLY: COMM_model, mype_model, npes_model
  USE mod_parallel_pdaf, &
       ONLY: MPIstatus, MPIerr

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! local state vector
  ! state_p is only allocated if a process participates in the coupling communication

! !CALLING SEQUENCE:
! Called by: PDAF_put_state_X    (as U_coll_state)
! Called by: PDAF_assimilate_X   (as U_coll_state)
!EOP

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
