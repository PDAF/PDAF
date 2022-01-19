!$Id: init_ens_offline.F90 1366 2013-04-24 16:25:05Z lnerger $
!>  Initialize ensemble
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all ensemble filters.
!!
!! The routine is called when the filter is
!! initialized in PDAF_filter_init.  It has
!! to initialize an ensemble of dim_ens states.
!!
!! The routine is called by all filter processes and 
!! initializes the ensemble for the PE-local domain.
!!
!! Implementation for the 2D offline example
!! with parallelization.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_ens_offline(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

  USE mpi                    ! MPI
  USE mod_parallel, &        ! Parallelization
       ONLY: mype_filter, npes_filter, COMM_filter, MPIerr, MPIstatus
  USE mod_assimilation, &    ! Assimilation variables
       ONLY: nx, ny, dim_state, local_dims

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: filtertype                !< Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                     !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                   !< Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)            !< PE-local model state
  !< (It is not necessary to initialize the array 'state_p' for ensemble filters.
  !< It is available here only for convenience and can be used freely.)
  REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) !< Array not referenced for ensemble filters
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)     !< PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                   !< PDAF status flag

! *** local variables ***
  INTEGER :: i, j, col, member        ! Counters
  REAL, ALLOCATABLE :: field(:,:)     ! global model field
  REAL, ALLOCATABLE :: ens(:,:)       ! global ensemble array
  CHARACTER(len=2) :: ensstr          ! String for ensemble member
  ! variables and arrays for domain decomposition
  INTEGER :: offset                   ! Row-offset according to domain decomposition
  INTEGER :: domain                   ! domain counter
  REAL,ALLOCATABLE :: ens_p_tmp(:,:)  ! Temporary ensemble for some PE-domain


! **********************
! *** INITIALIZATION ***
! **********************

  mype0a: IF (mype_filter==0) THEN

     WRITE (*, '(/9x, a)') 'Initialize state ensemble'
     WRITE (*, '(9x, a)') '--- read ensemble from files'
     WRITE (*, '(9x, a, i5)') '--- Ensemble size:  ', dim_ens

  END IF mype0a


! ********************************
! *** Read ensemble from files ***
! ********************************

  mype0b: IF (mype_filter==0) THEN
  
     ! *** Generate full ensemble on filter-PE 0 ***

     ! allocate memory for temporary fields
     ALLOCATE(field(ny, nx))
     ALLOCATE(ens(dim_state, dim_ens))

     DO member = 1, dim_ens
        WRITE (ensstr, '(i1)') member
        OPEN(11, file = '../inputs_offline/ens_'//TRIM(ensstr)//'.txt', status='old')

        DO i = 1, ny
           READ (11, *) field(i, :)
        END DO
        DO j = 1, nx
           ens(1 + (j-1)*ny : j*ny, member) = field(1:ny, j)
        END DO

        CLOSE(11)
     END DO
  END IF mype0b


! ****************************
! *** Distribute substates ***
! ****************************

  mype0c: IF (mype_filter == 0) THEN
     ! *** Initialize and send sub-state on PE 0 ***

     ! Initialize sub-ensemble for PE 0
     DO col = 1, dim_ens
        DO i=1, dim_p
           ens_p(i, col) = ens(i, col)
        END DO
     END DO

     ! Define offset in state vectors
     offset = local_dims(1)

     DO domain = 2, npes_filter
        ! Initialize sub-ensemble for other PEs and send sub-arrays

        ! Allocate temporary buffer array
        ALLOCATE(ens_p_tmp(local_dims(domain), dim_ens))

        ! Initialize MPI buffer for local ensemble
        DO col = 1, dim_ens
           DO i = 1, local_dims(domain)
              ens_p_tmp(i, col) = ens(i + offset, col)
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


! ****************
! *** clean up ***
! ****************

  IF (mype_filter==0) DEALLOCATE(field, ens)

END SUBROUTINE init_ens_offline
