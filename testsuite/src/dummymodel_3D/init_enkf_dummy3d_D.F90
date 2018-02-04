!$Id: init_enkf_dummy3d_D.F90 783 2009-12-07 10:28:43Z lnerger $
!BOP
!
! !ROUTINE: init_enkf --- Initialize ensemble for EnKF
!
! !INTERFACE:
SUBROUTINE init_enkf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF (EnKF):
!
! The routine is called when the filter is
! initialized in PDAF\_filter\_init.  It has
! to initialize an ensemble of states from a 
! mean state with prescribed covariance by 
! random sampling as the transformation of a 
! set of independent random numbers. The matrix 
! is initialized in the form of singular values
! and singular vectors. Based on this 
! information, the random  ensemble is generated.
!
! The routine is called by all filter processes and
! initializes the ensemble for the PE-local domain.
!
! This version is for the dummy model with domain 
! decomposition. 
!
! !REVISION HISTORY:
! 2004-12 - Lars Nerger - Initial code
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
  USE mod_model, &
       ONLY: dims_l_all, dims

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype             ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                  ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                ! Size of ensemble
  REAL, INTENT(out)   :: state_p(dim_p)         ! PE-local model state
  ! It is not necessary to initialize the array 'state' for EnKF. 
  ! It is available here only for convenience and can be used freely.
  REAL, INTENT(in)    :: Uinv(1, 1)             ! Array not referenced for EnKF
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)  ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_filter_init   (as U_ens_init)
! Calls: PDAF_seik_omega
! Calls: timeit
! Calls: memcount
! Calls: dgemm (BLAS)
! Calls: MPI_send 
! Calls: MPI_recv
!EOP

! *** local variables ***
  INTEGER :: i, j, member, col    ! Counters
  INTEGER :: dim_state            ! Global state dimension
  INTEGER :: dim_state_l          ! Local state dimension for some domain
  INTEGER, SAVE :: allocflag = 0  ! Flag for memory counting
  REAL, ALLOCATABLE :: ens(:,:)   ! Global ensemble
  REAL, ALLOCATABLE :: state(:)   ! Global state vector
  REAL, ALLOCATABLE :: eofV(:,:)  ! Matrix of eigenvectors V 
  REAL, ALLOCATABLE :: svals(:)   ! Singular values
  REAL, ALLOCATABLE :: Omega(:,:) ! Random matrix
  INTEGER :: omegatype            ! Type of matrix OMEGA
  INTEGER :: iseed(4)             ! Seed array for random number routine
  INTEGER :: rank                 ! Rank of approximated covariance matrix
  REAL :: randval                 ! Value of random number
  REAL :: norm                    ! Norm for ensemble transformation
  ! variables and arrays for domain decomposition
  INTEGER :: offset   ! Row-offset according to domain decomposition
  INTEGER :: domain   ! domain counter
  REAL, ALLOCATABLE :: ens_p_tmp(:,:)  ! temporary sub-array
  REAL, ALLOCATABLE :: state_p_tmp(:)  ! temporary sub-array



! **********************
! *** INITIALIZATION ***
! **********************

  ! *** Rank of matrix is ensemble size minus one
  rank = dim_ens - 1
  
  ! *** Generate full ensemble on filter-PE 0 ***
  mype0: IF (mype_filter == 0) THEN
     WRITE (*, '(/9x, a)') 'Generate state ensemble from covariance matrix'
     WRITE (*, '(9x, a)') '--- use random stochastic ensemble (EnKF type)'
     WRITE (*, '(9x, a, i5)') '--- Ensemble size:', dim_ens

     ! Set global state dimension
     dim_state = dims(1) * dims(2) * dims(3)

     ! allocate memory for temporary fields
     ALLOCATE(eofV(dim_state, rank))
     ALLOCATE(svals(rank))
     IF (allocflag == 0) THEN
        ! count allocated memory
        CALL memcount(2, 'r', dim_state * rank + rank)
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
     svals(1 : rank) = 1.0

     ! Set the initial ensemble to a part of the identity matrix
     eofV(:, :) = 0.0
     DO col = 1, rank
        eofV(col, col) = 1.0
     END DO


! ********************************************
! *** GENERATE ENSEMBLE                    ***
! *** as  _                              T ***
! *** X = X + norm eofV diag(svals) Omega  ***
! ********************************************

     WRITE (*, '(9x, a)') '--- Generate ensemble'

     CALL timeit(6, 'new')

     WRITE (*, '(13x, a, i5)') 'Rank of covar matrix = ', rank

     ! rescale eigenvectors
     DO j = 1, rank
        DO i = 1, dim_state
           eofV(i, j) = eofV(i, j) * svals(j)
        END DO
     END DO

     ! Initialize ISEED
     iseed(1) = 1
     iseed(2) = 2
     iseed(3) = 3
     iseed(4) = 55

     ! Generate random transformation matrix OMEGA
     ALLOCATE(Omega(dim_ens, rank))
     IF (allocflag == 0) THEN
        ! count allocated memory
        CALL memcount(2, 'r', dim_ens * rank)
     END IF

     omegatype = 6 ! Use unit columns orthogonal to (1,..,1)^T in OMEGA
     CALL PDAF_enkf_omega(iseed, rank, dim_ens, Omega, norm, omegatype, 1)

     ! generate random states
     DO member = 1, dim_ens
        ens(:, member) = state(:)
     END DO

     CALL dgemm('n', 't', dim_state, dim_ens, rank, &
          norm, eofV, dim_state, Omega, dim_ens, &
          1.0, ens, dim_state)
     
     CALL timeit(6, 'old')
     
  END IF mype0


! ****************************
! *** Distribute substates ***
! ****************************

  mype0b: IF (mype_filter == 0) THEN
     ! *** Initialize and send sub-state on PE 0 ***

     ! Initialize  sub_ensemble for PE 0
     DO col = 1, dim_ens
        DO i = 1, dim_p
           ens_p(i, col) = ens(i, col)
        END DO
     END DO
!            ! perform reordering of state for PE 0
!            DO i = 1, dim_p
!               state_p(i) = state(i)
!            END DO

     ! Define offset in state vectors
     offset = dim_p

     DO domain = 2, npes_filter

!       write (*,*) 'Domain: ',domain,' Offset: ',offset

        ! Initialize and sub_ensemble for other PEs and send sub-arrays

        ! allocate temporary sub-arrays
        dim_state_l = dims_l_all(1, domain) * dims_l_all(2, domain) &
             * dims_l_all(3, domain)
        ALLOCATE(ens_p_tmp(dim_state_l, dim_ens))
!         ALLOCATE(state_p_tmp(dim_state_l))
        IF (allocflag == 0) THEN
           ! count allocated memory
           CALL memcount(2, 'r', dim_state_l * dim_ens)
           allocflag = 1
        END IF

        ! perform reordering of mode matrix
        DO col = 1,dim_ens
           DO i = 1, dim_state_l
              ens_p_tmp(i, col) = ens(i + offset, col)
           END DO
        END DO
!            ! perform reordering of state
!            DO i = 1, dim_state_l
!               state_p_tmp(i) = state(i + offset)
!            END DO

        ! Send sub-arrays
        CALL MPI_send(ens_p_tmp, dim_ens * dim_state_l, &
             MPI_DOUBLE_PRECISION, domain - 1, 1, COMM_filter, MPIerr)
!            CALL MPI_send(state_p_tmp, dim_state_l, &
!                 MPI_DOUBLE_PRECISION, domain - 1, 2, COMM_filter, MPIerr)

!         DEALLOCATE(ens_p_tmp, state_p_tmp)
        DEALLOCATE(ens_p_tmp)

        ! Increment offset
        offset = offset + dim_state_l
        
     END DO

  ELSE mype0b
     ! *** Receive substate on filter-PEs with rank > 0 ***

!     write (*,*) 'RECV: domain ',mype_filter+1,'dim_p',dim_p,'dim_ens',dim_ens
     CALL MPI_recv(ens_p, dim_p * dim_ens, MPI_DOUBLE_PRECISION, &
          0, 1, COMM_filter, MPIstatus, MPIerr)
!      CALL MPI_recv(state_p, dim_p, MPI_DOUBLE_PRECISION, &
!           0, 2, COMM_filter, MPIstatus, MPIerr)

  END IF mype0b


! ****************
! *** clean up ***
! ****************

  IF (mype_filter == 0) THEN
     DEALLOCATE(Omega)
     DEALLOCATE(svals, eofV)
     DEALLOCATE(ens, state)
  END IF

END SUBROUTINE init_enkf
