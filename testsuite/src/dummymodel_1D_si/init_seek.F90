!$Id: init_seek.F90 1099 2011-08-17 14:15:30Z lnerger $
!BOP
!
! !ROUTINE: init_seek --- Initialize state and modes for SEEK
!
! !INTERFACE:
SUBROUTINE init_seek(filtertype, dim_p, rank, state_p, Uinv, &
     eofV_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF (SEEK):
!
! The routine is called when the filter is
! initialized in PDAF\_filter\_init.  It has
! to initialize the state estimate and the 
! approximate covariance matrix in the form
!         $P = V U V^T$
! where P is dim x dim, V is dim x rank, and 
! U is rank x rank. With regard to the 
! parallelization U is global while V is for
! PE-local domain.
!
! The routine is called by all filter processes.
!
! This version is for the dummy model with domain 
! decomposition. This example first initializes 
! the full ensemble on the process with rank 0. 
! Subsequently, the sub-states are distributed 
! according to the domain decomposition.
! (See init\_seik\_dist() for a distributed 
! variant of the initialization.)
!
! !REVISION HISTORY:
! 2004-12 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel, &
       ONLY: mype_filter, npes_filter, COMM_filter, MPI_DOUBLE_PRECISION, &
       MPIerr, MPIstatus
  USE mod_memcount, &
       ONLY: memcount
  USE mod_model, &
       ONLY: local_dims,dim_state

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype          ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p               ! PE-local state dimension
  INTEGER, INTENT(in) :: rank                ! Number of eofs to be used
  REAL, INTENT(out)   :: state_p(dim_p)      ! PE-local model state
  REAL, INTENT(out)   :: Uinv(rank, rank)    ! Inverse of eigenvalue matrix U
  REAL, INTENT(out)   :: eofV_p(dim_p, rank) ! PE-local mode matrix V
  INTEGER, INTENT(inout) :: flag             ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: init_ens_pdaf
! Calls: memcount
! Calls: MPI_send
! Calls: MPI_recv
! Calls: MPI_bcast
!EOP

! *** local variables ***
  INTEGER :: i, col   ! counters
  REAL, ALLOCATABLE :: eofV(:,:) ! global EOF matrix
  REAL, ALLOCATABLE :: state(:)  ! global state vector
  REAL, ALLOCATABLE :: svals(:)  ! singular values
  INTEGER, SAVE :: allocflag = 0 ! Flag for memory counting
  ! Variables and arrays for domain decomposition
  INTEGER :: offset   ! Row-offset according to domain decomposition
  INTEGER :: domain   ! domain counter
  REAL, ALLOCATABLE :: eofV_p_temp(:,:) ! temporary sub-array
  REAL, ALLOCATABLE :: state_p_temp(:)  ! temporary sub-array


! **********************
! *** INITIALIZATION ***
! **********************
  
  ! *** Generate full EOF information on filter-PE 0 ***
  mype0: IF (mype_filter == 0) THEN
     WRITE (*, '(/9x, a)') &
          'Initialize state estimate and EOF decomposed covariance matrix'
     WRITE (*, '(9x, a, i5)') '--- number of EOFs:', rank

     ! allocate memory for temporary fields
     ALLOCATE(svals(rank))
     IF (allocflag == 0) THEN
        ! count allocated memory
        CALL memcount(2, 'r', rank)
     END IF


! *************************************************
! *** Initialize initial state and covar matrix ***
! *************************************************

     ! Very simple initialization for dummy model

     ! Allocate global ensemble and state
     ALLOCATE(eofV(dim_state, rank))
     ALLOCATE(state(dim_state))
     IF (allocflag == 0) THEN
        ! count allocated memory
        CALL memcount(2, 'r', dim_state * rank + dim_state)
     END IF

     ! Just set the entries of the state vector to 2.0
     state(:) = 2.0

     ! Set the initial singular vectors to one
     svals(1 : rank) = 1.0

     ! Set the initial eigenvalue matrix U^-1
     Uinv(:, :) = 0.0
     DO col = 1, rank
        Uinv(col, col) = 1.0 / svals(col)
     END DO

     ! Set the initial ensemble to a part of the identity matrix
     eofV(:, :) = 0.0
     DO col = 1, rank
        eofV(col, col) = 1.0
     END DO

  END IF mype0


! ****************************
! *** Distribute substates ***
! ****************************

  mype0b: IF (mype_filter == 0) THEN
     ! *** Initialize and send sub-state on PE 0 ***

     ! Initialize sub-state and sub_ensemble for PE 0

     ! perform direct initialization of mode matrix for PE 0
     DO col = 1, rank
        DO i = 1, dim_p
           eofV_p(i, col) = eofV(i, col)
        END DO
     END DO
     ! perform direct initialization of state for PE 0
     DO i = 1, dim_p
        state_p(i) = state(i)
     END DO

     ! Define offset in state vectors
     offset = dim_p

     DO domain = 2, npes_filter
        ! Initialize sub-state and sub_ensemble for other PEs
        ! and send sub-arrays

        ! allocate temporary sub-arrays
        ALLOCATE(eofV_p_temp(local_dims(domain), rank))
        ALLOCATE(state_p_temp(local_dims(domain)))
        IF (allocflag == 0) THEN
           ! count allocated memory
           CALL memcount(2, 'r', local_dims(domain) * rank + local_dims(domain))
           allocflag = 1
        END IF

        ! prepare buffer of local mode matrix
        DO col = 1, rank
           DO i = 1, local_dims(domain)
              eofV_p_temp(i, col) = eofV(i + offset, col)
           END DO
        END DO
        ! prepare buffer of local state
        DO i = 1, local_dims(domain)
           state_p_temp(i) = state(i + offset)
        END DO

        ! Send sub-arrays
        CALL MPI_send(eofV_p_temp, rank * local_dims(domain), &
             MPI_DOUBLE_PRECISION, domain - 1, 1, COMM_filter, MPIerr)
        CALL MPI_send(state_p_temp, local_dims(domain), &
             MPI_DOUBLE_PRECISION, domain - 1, 2, COMM_filter, MPIerr)

        DEALLOCATE(eofV_p_temp, state_p_temp)

        ! Increment offset
        offset = offset + local_dims(domain)

     END DO

  ELSE mype0b
     ! *** Receive substate on filter-PEs with rank > 0 ***

     CALL MPI_recv(eofV_p, dim_p * rank, MPI_DOUBLE_PRECISION, &
          0, 1, COMM_filter, MPIstatus, MPIerr)
     CALL MPI_recv(state_p, dim_p, MPI_DOUBLE_PRECISION, &
          0, 2, COMM_filter, MPIstatus, MPIerr)

  END IF mype0b

  ! *** broadcast eigenvalue matrix U to all filter PEs
  CALL MPI_bcast(Uinv, rank * rank, MPI_DOUBLE_PRECISION, &
       0, COMM_filter, MPIerr)


! ****************
! *** clean up ***
! ****************

  IF (mype_filter == 0) THEN
     DEALLOCATE(svals, eofV, state)
  END IF

END SUBROUTINE init_seek
