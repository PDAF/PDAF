!$Id$
!BOP
!
! !ROUTINE: init_3dvar_eof --- Initialize ensemble from EOF decomposition
!
! !INTERFACE:
SUBROUTINE init_3dvar_eof(dim, dim_ens, state, ens, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF (all filters):
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
! This version is for the Lorenz96 model
! without parallelization.
!
! !REVISION HISTORY:
! 2021-03 - Lars Nerger - Initial code based on init_ens_eof
! Later revisions - see svn log
!
! !USES:
  USE mod_memcount, &
       ONLY: memcount
  USE mod_assimilation, &
       ONLY: covartype, file_ini, Vmat, dim_cvec

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim                 ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens             ! Size of ensemble
  REAL, INTENT(inout) :: state(dim)          ! PE-local model state
  ! It is not necessary to initialize the array 'state' for SEIK. 
  ! It is available here only for convenience and can be used freely.
  REAL, INTENT(out)   :: ens(dim, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag             ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: init_3dvar
! Calls: timeit
! Calls: memcount
! Calls: PDAF_SampleEns
!EOP

! *** local variables ***
  INTEGER :: i, s, row, col       ! counters
  INTEGER, SAVE :: allocflag = 0  ! Flag for memory counting
  REAL, ALLOCATABLE :: eofV(:,:)  ! matrix of eigenvectors V 
  REAL, ALLOCATABLE :: svals(:)   ! singular values
  REAL :: invdim_ens              ! Inverse ensemble size
  REAL :: fact                    ! Scaling factor
  INTEGER :: dim_file             ! State dimension in file
  INTEGER :: rank                 ! Rank of approximated covariance matrix
  INTEGER :: rank_file            ! Rank of covariance matrix stored in file
  INTEGER :: stat(50)             ! Array for status flag
  INTEGER :: fileid               ! ID for NetCDF file
  INTEGER :: id_svals, id_eofV    ! IDs for fields
  INTEGER :: id_state             ! ID for field
  INTEGER :: id_dim               ! ID for dimension
  INTEGER :: pos(2)               ! Position index for writing
  INTEGER :: cnt(2)               ! Count index for writing


! **********************
! *** INITIALIZATION ***
! **********************

  ! *** Rank of matrix is ensemble size minus one
  rank = dim_cvec - 1
  
  ! *** Generate full ensemble on filter-PE 0 ***
  WRITE (*, '(9x, a)') 'Initialize state and B^1/2 for 3D-Var'
  WRITE (*, '(9x, a)') &
       '--- use rank reduction and 2nd order exact sampling (SEIK type)'
  WRITE (*, '(9x, a, i5)') '--- number of EOFs:    ', rank
  WRITE (*, '(9x, a, i5)') '--- columns in B^1/2:  ', dim_cvec

  ! Initialize numbers 
  invdim_ens = 1.0 / REAL(dim_cvec)

  ! allocate memory for temporary fields
  ALLOCATE(eofV(dim, rank))
  ALLOCATE(svals(rank))
  IF (allocflag == 0) THEN
     ! count allocated memory
     CALL memcount(2, 'r', dim * rank + rank)
  END IF

  ! Allocate matrix holding B^1/2 (from mod_assimilation)
  ALLOCATE(Vmat(dim, dim_cvec))


! *************************************************
! *** Initialize initial state and covar matrix ***
! *************************************************
  
  WRITE(*,'(9x,a,a)') '--- Reading covariance information from ', TRIM(file_ini)

  s = 1
  stat(s) = NF_OPEN(file_ini, NF_NOWRITE, fileid)

  ! Read size of state vector
  s = s + 1
  stat(s) = NF_INQ_DIMID(fileid, 'dim_state', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(fileid, id_dim, dim_file)

  ! Read rank stored in file
  s = s + 1
  stat(s) = NF_INQ_DIMID(fileid, 'rank', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(fileid, id_dim, rank_file)

  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in reading dimensions from init file, no.', i
  END DO

  ! Check consistency of dimensions
  checkdim: IF (dim == dim_file .AND. rank_file >= rank) THEN

     ! Inquire IDs for mean state, singular vectors and values
     s = 1
     stat(s) = NF_INQ_VARID(fileid, 'meanstate', id_state)
     s = s + 1
     stat(s) = NF_INQ_VARID(fileid, 'u_svd', id_eofV)
     s = s + 1
     stat(s) = NF_INQ_VARID(fileid, 'sigma', id_svals)

     ! Read initialization information
     s = s + 1
     stat(s) = NF_GET_VAR_DOUBLE(fileid, id_state, state)

     pos(2) = 1
     cnt(2) = rank
     pos(1) = 1
     cnt(1) = dim
     s = s + 1
     stat(s) = NF_GET_VARA_DOUBLE(fileid, id_eofV, pos, cnt, eofV)

     s = s + 1
     stat(s) = NF_GET_VARA_DOUBLE(fileid, id_svals, 1, rank, svals)

     s = s + 1
     stat(s) = nf_close(fileid)

     DO i = 1,  s
        IF (stat(i) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in reading initialization file, no.', i
     END DO


! *****************************************
! *** Generate ensemble of model states ***
! *****************************************

     WRITE (*, '(9x, a)') '--- generate ensemble of states'
     
     ! Use PDAF routine to generate ensemble from covariance matrix
     CALL PDAF_SampleEns(dim, dim_cvec, eofV, svals, state, Vmat, 1, flag)
     
  ELSE

      ! *** Rank stored in file is smaller than requested EOF rank ***
     WRITE(*,*) 'Rank stored in file is smaller than requested EOF rank'

     stat(s) = nf_close(fileid)
     STOP

  END IF checkdim


! **********************************************
! *** Initialize square-root of P for 3D-Var ***
! **********************************************

  WRITE (*, '(9x, a)') 'Initialize B^1/2'

  ! Here, we simply use the scaled ensemble perturbations

  ! Compute ensemble mean
  state = 0.0
  DO col = 1, dim_cvec
     DO i = 1, dim
        state(i) = state(i) + Vmat(i, col)
     END DO
  END DO
  state(:) = invdim_ens * state(:)

  ! Initialize ensemble perturbations
  DO col = 1, dim_cvec
     Vmat(:,col) = Vmat(:,col) - state(:)
  END DO

  fact = 1.0/SQRT(REAL(dim_cvec-1))

  Vmat = Vmat * fact


! ******************************************
! *** Initialize ensemble array for PDAF ***
! ******************************************

  ens(:,1) = state(:)


! ****************
! *** clean up ***
! ****************

  DEALLOCATE(svals, eofV)

END SUBROUTINE init_3dvar_eof
