!$Id: init_seik_dummy3d_D.F90 783 2009-12-07 10:28:43Z lnerger $
!BOP
!
! !ROUTINE: init_seik --- Initialize ensemble for SEIK
!
! !INTERFACE:
SUBROUTINE init_seik(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF (SEIK):
!
! The routine is called when the filter is
! initialized in PDAF\_filter\_init.  It has
! to initialize an ensemble of dim\_ens states
! by exact 2nd order sampling.
! State vectors of the form
!   $x_i = x + sqrt(FAC) eofV (\Omega C^{-1})^T$
! fulfill the condition
!   $P = 1/(FAC)  \sum_{i=1}^{dim\_ens} (x_i - x)(x_i - x)^T$
! The matrix is initialized in the form of
! singular values and singular vectors.
! The routine is called by all filter processes and 
! initializes the ensemble for the local domain.
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
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE timer, &
       ONLY: timeit
  USE mod_memcount, &
       ONLY: memcount
  USE mod_parallel, &
       ONLY: mype_filter, npes_filter, COMM_filter, MPI_DOUBLE_PRECISION, &
       MPIerr, MPIstatus
  USE mod_assimilation, &
       ONLY: covartype
  USE mod_model, &
       ONLY: dims_l_all, dims

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
  REAL, INTENT(out)   :: state_p(dim_p)          ! PE-local model state
  ! It is not necessary to initialize the array 'state_p' for SEIK. 
  ! It is available here only for convenience and can be used freely.
  REAL, INTENT(out)   :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_filter_init    (as U_ens_init)
! Calls: PDAF_seik_omega
! Calls: timeit
! Calls: memcount
! Calls: dgemm (BLAS)
! Calls: MPI_send 
! Calls: MPI_recv
!EOP

! *** local variables ***
  INTEGER :: i, row, col              ! counters
  INTEGER :: dim_state                ! Global state dimension
  INTEGER :: dim_state_l              ! Local state dimension for some domain
  INTEGER, SAVE :: allocflag = 0      ! Flag for memory counting
  REAL, ALLOCATABLE :: ens(:,:)       ! global ensemble
  REAL, ALLOCATABLE :: state(:)       ! global state vector
  REAL, ALLOCATABLE :: eofV(:,:)      ! matrix of eigenvectors V 
  REAL, ALLOCATABLE :: svals(:)       ! singular values
  REAL, ALLOCATABLE :: omega(:,:)     ! Matrix Omega
  INTEGER :: rank     ! Rank of approximated covariance matrix
  REAL :: fac         ! Square-root of dim_ens-1 or dim_ens
  ! variables and arrays for domain decomposition
  INTEGER :: offset   ! Row-offset according to domain decomposition
  INTEGER :: domain   ! domain counter
  REAL,ALLOCATABLE :: ens_p_tmp(:,:) ! Temporary ensemble for some PE-domain
  REAL,ALLOCATABLE :: state_p_tmp(:) ! Temporary state vector for some PE-domain


! **********************
! *** INITIALIZATION ***
! **********************

  ! *** Rank of matrix is ensemble size minus one
  rank = dim_ens - 1
  
  ! *** Generate full ensemble on filter-PE 0 ***
  mype0: IF (mype_filter == 0) THEN
     WRITE (*, '(/9x, a)') 'Generate state ensemble from covariance matrix'
     WRITE (*, '(9x, a)') &
          '--- use rank reduction and 2nd order exact sampling (SEIK type)'
     WRITE (*, '(9x, a, i5)') '--- Ensemble size:  ', dim_ens
     WRITE (*, '(9x, a, i5)') '--- number of EOFs: ', rank

     ! Set global state dimension
     dim_state = dims(1) * dims(2) * dims(3)

     ! allocate memory for temporary fields
     ALLOCATE(eofV(dim_state, rank))
     ALLOCATE(svals(rank))
     ALLOCATE(omega(rank + 1, rank))
     IF (allocflag == 0) THEN
        ! count allocated memory
        CALL memcount(2, 'r', dim_state * rank + rank + rank * (rank + 1))
     END IF


! *************************************************
! *** Initialize initial state and covar matrix ***
! *************************************************

     ! Very simple initialization for dummy model

     ! Allocate global ensemble and state
     ALLOCATE(ens(dim_state, dim_ens))
     ALLOCATE(state(dim_state))
     IF (allocflag == 0) THEN
        ! count allocated memory
        CALL memcount(2, 'r', dim_state * dim_ens + dim_state)
     END IF

     ! Just set the entries of the state vector to 2.0
     state(:) = 2.0

     ! Set the initial singular vectors to one
     svals(1:rank) = 1.0

     ! Set the initial ensemble to a part of the identity matrix
     eofV(:,:) = 0.0
     DO col = 1, rank
        eofV(col, col) = 1.0
     END DO


! *****************************************************
! *** DECOMPOSE COVARIANCE                          ***
! ***                                               ***
! *** P = eofV U eofV^T                             ***
! ***   = eofV C^(-1)^T Omega^T Omega C^(-1) eofV^T ***
! *** where U^(-1) = C C^T                          ***
! ***                                               ***
! *** Since the matrix is already initialized in    ***
! *** decomposed form we directly have the          ***
! *** inverses of C given by the singular values    ***
! *****************************************************


! *************************************************
! *** Generate ensemble of interpolating states ***
! *************************************************

     ! Very simple method here: We generate the full 
     ! ensemble on the filter PE with rank 0. Afterwards
     ! we distribute sub-states to other filter PEs

     CALL timeit(6, 'new')

     WRITE (*, '(9x, a)') '--- generate ensemble of model states'

     ! *** Generate uniform orthogonal matrix OMEGA ***
     CALL PDAF_seik_omega(rank, Omega, 1)

     ! ***      Generate ensemble of states         ***
     ! *** x_i = x + sqrt(FAC) eofV (Omega C^(-1))t ***

     ! A = Omega C^(-1)
     DO col = 1, rank
        DO row = 1, rank + 1
           Omega(row, col) = Omega(row, col) * svals(col)
        END DO
     END DO
      
     ! ens = state+ sqrt(FAC) eofV A^T
     ! FAC depends on the definiton of the factor in the ensemble
     ! covar. matrix which is defined by the variable COVARTYPE
     DO col = 1, rank + 1
        ens(:, col) = state(:)
     END DO

     IF (covartype == 1) THEN
        fac = SQRT(REAL(rank))
     ELSE
        fac = SQRT(REAL(rank + 1))
     END IF
     CALL dgemm('n', 't', dim_state, rank + 1, rank, &
          fac, eofV, dim_state, Omega, rank + 1, &
          1.0, ens, dim_state)

     CALL timeit(6, 'old')
     
  END IF mype0
  

! ****************************
! *** Distribute substates ***
! ****************************

  mype0b: IF (mype_filter == 0) THEN
     ! *** Initialize and send sub-state on PE 0 ***

     ! Initialize sub_ensemble for PE 0
     DO col = 1, dim_ens
        DO i=1, dim_p
           ens_p(i, col) = ens(i, col)
        END DO
     END DO
           ! perform reordering of state for PE 0
!            DO i = 1, dim_p
!               state_p(i) = state(i)
!            END DO

     ! Define offset in state vectors
     offset = dim_p

     DO domain = 2, npes_filter

        ! Initialize sub_ensemble for other PEs and send sub-arrays

        ! Allocate temporary sub-arrays
        dim_state_l = dims_l_all(1, domain) * dims_l_all(2, domain) &
             * dims_l_all(3, domain)
        ALLOCATE(ens_p_tmp(dim_state_l, dim_ens))
!            ALLOCATE(state_p_tmp(dim_state_l))
        IF (allocflag ==0) THEN
           ! count allocated memory
           CALL memcount(2, 'r', dim_state_l * dim_ens)
           allocflag = 1
        END IF

        ! perform reordering of mode matrix
        DO col = 1, dim_ens
           DO i = 1, dim_state_l
              ens_p_tmp(i, col) = ens(i + offset, col)
           END DO
        END DO
           ! perform reordering of state
!            DO i = 1, dim_state_l
!               state_p_tmp(i) = state(i + offset)
!            END DO

        ! Send sub-arrays
        CALL MPI_send(ens_p_tmp, dim_ens * dim_state_l, &
             MPI_DOUBLE_PRECISION, domain - 1, 1, COMM_filter, MPIerr)
!            CALL MPI_send(state_p_tmp, dim_state_l, &
!                 MPI_DOUBLE_PRECISION, domain - 1, 2, COMM_filter, MPIerr)

!            DEALLOCATE(ens_p_tmp, state_p_tmp)
        DEALLOCATE(ens_p_tmp)

        ! Increment offset
        offset = offset + dim_state_l

     END DO

  ELSE mype0b
     ! *** Receive substate on filter-PEs with rank > 0 ***

     CALL MPI_recv(ens_p, dim_p * dim_ens, MPI_DOUBLE_PRECISION, &
          0, 1, COMM_filter, MPIstatus, MPIerr)
!      CALL MPI_recv(state_p, dim_p, MPI_DOUBLE_PRECISION, &
!           0, 2, COMM_filter, MPIstatus, MPIerr)
     
  END IF mype0b


! ****************
! *** clean up ***
! ****************

  IF (mype_filter == 0) THEN
     DEALLOCATE(svals, eofV, omega)
     DEALLOCATE(ens, state)
  END IF

END SUBROUTINE init_seik
