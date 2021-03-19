!$Id$
!BOP
!
! !ROUTINE: init_3dvar_pdaf --- Initialize ensemble
!
! !INTERFACE:
SUBROUTINE init_3dvar_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in: 3D-Var
!
! The routine is called when the filter is
! initialized in PDAF\_filter\_init. It has
! to initialize the initial state estimate
! for the assimilation. In addition it 
! initializes the square-root of P that is 
! used for 3D-Var (in this example we use
! explicitly a matrix holding the square-root
! of P, which is given by the scaled ensemble
! perturbations.
!
! The initialization is done here analogously
! to the ensemble initialization using second-
! order exact sampling as done in init_seik_pdaf.
! The internal ensemble-initialization is done 
! here for a sample size dim_cvec (instead dim_ens)
!
! !REVISION HISTORY:
! 2021-03 - Lars Nerger - Initial code cased on init_seik
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
       ONLY: local_dims, dim_state
  USE mod_assimilation, &
       ONLY: Vmat_p, dim_cvec, subtype
  USE PDAF_interfaces_module, &
       ONLY: PDAF_sampleens

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state
  ! It is not necessary to initialize the array 'state_p' for 3D-Var
  ! It is available here only for convenience and can be used freely.
  REAL, INTENT(inout) :: Uinv(1,1)               ! Array not referenced for 3D-Var
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_init    (as U_init_ens)
! Calls: PDAF_sampleens
! Calls: timeit
! Calls: memcount
! Calls: dgemm (BLAS)
! Calls: MPI_send 
! Calls: MPI_recv
!EOP

! *** local variables ***
  INTEGER :: i, col                   ! counters
  INTEGER, SAVE :: allocflag = 0      ! Flag for memory counting
  REAL, ALLOCATABLE :: ens(:,:)       ! global ensemble
  REAL, ALLOCATABLE :: state(:)       ! global state vector
  REAL, ALLOCATABLE :: eofV(:,:)      ! matrix of eigenvectors V 
  REAL, ALLOCATABLE :: svals(:)       ! singular values
  INTEGER :: rank                     ! Rank of approximated covariance matrix
  ! variables and arrays for domain decomposition
  INTEGER :: offset                   ! Row-offset according to domain decomposition
  INTEGER :: domain                   ! domain counter
  REAL,ALLOCATABLE :: ens_p_tmp(:,:)  ! Temporary ensemble for some PE-domain
  REAL,ALLOCATABLE :: Vmat_p_tmp(:,:) ! Temporary Vmat for some PE-domain
  REAL :: fact                        ! Scaling factor


! **********************
! *** INITIALIZATION ***
! **********************

  ! *** Rank of matrix is ensemble size minus one
  rank = dim_cvec - 1

  ! Allocate matrix holding B^1/2 (from mod_assimilation)
  ALLOCATE(Vmat_p(dim_p, dim_cvec))
  
  ! *** Generate full ensemble on filter-PE 0 ***
  mype0: IF (mype_filter == 0) THEN
     WRITE (*, '(/9x, a)') 'Initialize state and B^1/2 for 3D-Var'
     WRITE (*, '(9x, a)') '--- generate B^1/2 from covariance matrix'
     WRITE (*, '(9x, a)') &
          '--- use rank reduction and 2nd order exact sampling (SEIK type)'
     IF (filtertype==11) &
          WRITE (*, '(9x, a)') '+++ Use true initial state for generating synthetic observations'
     WRITE (*, '(9x, a, i5)') '--- number of EOFs:    ', rank
     WRITE (*, '(9x, a, i5)') '--- columns in B^1/2:  ', dim_cvec

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

! We generate the ensemble by 2nd order exact sampling
! from a covariance matrix. The covariance matrix is
! decomposed as 
!      P = eofV svals eofV^T                            
! which can be obtained e.g. using the tool PDAF_eofcovar.
! The arrays eofV and svals are then used in PDAF_sampleens
! to generate the ensemble array 'ens'. In addtion
! PDAF_sampleens needs to ensemble mean state aroudn which
! the ensmeble is sampled.
!
! For the dummy model, we use a very simple initialization:
! We just specify the ensemble mean state ('state'), the
! singular values ('svals'), and the eigentvectors ('eof')
! before we call PDAF_sampleens.

     ! Allocate global ensemble and state
     ALLOCATE(ens(dim_state, dim_cvec))
     ALLOCATE(state(dim_state))
     IF (allocflag == 0) THEN
        ! count allocated memory
        CALL memcount(2, 'r', dim_state * dim_cvec + dim_state)
     END IF

     ! Just set the entries of the state vector to 2.0
     IF (filtertype/=11) THEN
        ! For applying data assimilation
        state(:) = 2.0
     ELSE
        ! For generating synthetic observations
        state(:) = 1.0
        IF (mype_filter == 0) &
             WRITE (*, '(9x, a)') '--- initialize ensemble mean for observation generation'
     END IF

     ! Set the initial singular vectors to one
     svals(1:rank) = 1.0

     ! Set the initial ensemble to a part of the identity matrix
     eofV(:,:) = 0.0
     DO col = 1, rank
        eofV(col, col) = 1.0
     END DO


! *************************************************
! *** Generate ensemble of interpolating states ***
! *************************************************

     ! Very simple method here: We generate the full 
     ! ensemble on the filter PE with rank 0. Afterwards
     ! we distribute sub-states to other filter PEs

     CALL timeit(6, 'new')

     WRITE (*, '(9x, a)') '--- generate ensemble of states'

     ! Generate ensemble using PDAF sampling routine
     CALL PDAF_SampleEns(dim_state, dim_cvec, eofV, svals, state, &
          ens, 1, flag)

     CALL timeit(6, 'old')

  END IF mype0
  

! ****************************
! *** Distribute substates ***
! ****************************

  ! *********************************************
  ! *** Distribute ensmeble to prepare B^1/2  ***
  ! *********************************************

  mype0b: IF (mype_filter == 0) THEN
     ! *** Initialize and send sub-state on PE 0 ***

     ! Initialize sub-ensemble for PE 0
     DO col = 1, dim_cvec
        DO i=1, dim_p
           Vmat_p(i, col) = ens(i, col)
        END DO
     END DO

     ! Define offset in state vectors
     offset = dim_p

     WRITE (*, '(9x, a)') '--- Distribute substates'

     DO domain = 2, npes_filter
        ! Initialize sub-ensemble for other PEs and send sub-arrays

        ! Allocate temporary buffer array
        ALLOCATE(Vmat_p_tmp(local_dims(domain), dim_cvec))
        IF (allocflag ==0) THEN
           ! count allocated memory
           CALL memcount(2, 'r', local_dims(domain) * dim_cvec)
           allocflag = 1
        END IF

        ! Initialize MPI buffer for local ensemble
        DO col = 1, dim_cvec
           DO i = 1, local_dims(domain)
              Vmat_p_tmp(i, col) = ens(i + offset, col)
           END DO
        END DO

        ! Send sub-arrays
        CALL MPI_send(Vmat_p_tmp, dim_cvec * local_dims(domain), &
             MPI_DOUBLE_PRECISION, domain - 1, 1, COMM_filter, MPIerr)

        DEALLOCATE(Vmat_p_tmp)

        ! Increment offset
        offset = offset + local_dims(domain)

     END DO

  ELSE mype0b
     ! *** Receive ensemble substates on filter-PEs with rank > 0 ***

     CALL MPI_recv(Vmat_p, dim_p * dim_cvec, MPI_DOUBLE_PRECISION, &
          0, 1, COMM_filter, MPIstatus, MPIerr)
     
  END IF mype0b


  ! *****************************************
  ! *** Distribute the state estimate     ***
  ! *** (stored in ens_p returned to PDAF ***
  ! *****************************************

  mype0c: IF (mype_filter == 0) THEN
     ! *** Initialize and send sub-state on PE 0 ***

     ! Initialize sub-ensemble for PE 0
     DO col = 1, dim_ens
        DO i=1, dim_p
           ens_p(i, col) = state(i)
        END DO
     END DO

     ! Define offset in state vectors
     offset = dim_p

     DO domain = 2, npes_filter
        ! Initialize sub-ensemble for other PEs and send sub-arrays

        ! Allocate temporary buffer array
        ALLOCATE(ens_p_tmp(local_dims(domain), dim_ens))
        IF (allocflag ==0) THEN
           ! count allocated memory
           CALL memcount(2, 'r', local_dims(domain) * dim_ens)
           allocflag = 1
        END IF

        ! Initialize MPI buffer for local ensemble
        DO col = 1, dim_ens
           DO i = 1, local_dims(domain)
              ens_p_tmp(i, col) = state(i + offset)
           END DO
        END DO

        ! Send sub-arrays
        CALL MPI_send(ens_p_tmp, dim_ens * local_dims(domain), &
             MPI_DOUBLE_PRECISION, domain - 1, 1, COMM_filter, MPIerr)

        DEALLOCATE(ens_p_tmp)

        ! Increment offset
        offset = offset + local_dims(domain)

     END DO

  ELSE mype0c
     ! *** Receive ensemble substates on filter-PEs with rank > 0 ***

     CALL MPI_recv(ens_p, dim_p * dim_ens, MPI_DOUBLE_PRECISION, &
          0, 1, COMM_filter, MPIstatus, MPIerr)
     
  END IF mype0c


! ****************************************
! *** Initialize B^1/2 for 3D-Var      ***
! *** stored as Vmat_p outside of PDAF ***
! ****************************************

  WRITE (*, '(9x, a)') '--- initialize B^1/2'

  ! Here, we simply use the scaled ensemble perturbations
  DO col = 1, dim_cvec
     Vmat_p(:,col) = Vmat_p(:,col) - ens_p(:,1)
  END DO

  fact = 1.0/SQRT(REAL(dim_cvec-1))

  Vmat_p = Vmat_p * fact


! ****************
! *** clean up ***
! ****************

  IF (mype_filter == 0) THEN
     DEALLOCATE(svals, eofV)
     DEALLOCATE(ens, state)
  END IF

END SUBROUTINE init_3dvar_pdaf
