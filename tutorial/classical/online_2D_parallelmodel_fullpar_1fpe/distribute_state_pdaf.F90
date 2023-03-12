!$Id: distribute_state_pdaf.F90 1411 2013-09-25 14:04:41Z lnerger $
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
  USE mod_parallel_model, &
       ONLY: COMM_model, mype_model, npes_model
  USE mod_parallel_pdaf, &
       ONLY: MPI_DOUBLE_PRECISION, MPIstatus, MPIerr

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! PE-local state vector
  ! state_p is only allocated if a process participates in the coupling communication

! !CALLING SEQUENCE:
! Called by: PDAF_get_state      (as U_dist_state)
! Called by: PDAF_assimilate_X   (as U_coll_state)
!EOP

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
