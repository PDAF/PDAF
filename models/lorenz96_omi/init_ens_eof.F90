!$Id$
!BOP
!
! !ROUTINE: init_ens_eof --- Initialize ensemble from EOF decomposition
!
! !INTERFACE:
SUBROUTINE init_ens_eof(dim, dim_ens, state, ens, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF (all filters):
!
! The routine is called by init_seik. It 
! initializes an ensemble of dim\_ens states
! by exact 2nd order sampling.
! State vectors of the form
!   $x_i = x + sqrt(dim_ens-1) eofV (\Omega C^{-1})^T$
! fulfill the condition
!   $P = 1/(dim_ens-1)  \sum_{i=1}^{dim\_ens} (x_i - x)(x_i - x)^T$
! The matrix is initialized in the form of
! singular values and singular vectors. Here, the
! subroutine PDAF_SampleEns is used to generate
! the ensemble states.
!
! This version is for the Lorenz96 model
! without parallelization.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE netcdf
  USE mod_memcount, &
       ONLY: memcount
  USE mod_assimilation, &
       ONLY: file_ini

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim                 ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens             ! Size of ensemble
  REAL, INTENT(inout) :: state(dim)          ! PE-local model state
  ! It is not necessary to initialize the array 'state' for SEIK. 
  ! It is available here only for convenience and can be used freely.
  REAL, INTENT(out)   :: ens(dim, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag             ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: init_ens
! Calls: timeit
! Calls: memcount
! Calls: PDAF_SampleEns
!EOP

! *** local variables ***
  INTEGER :: i, s                 ! counters
  INTEGER, SAVE :: allocflag = 0  ! Flag for memory counting
  REAL, ALLOCATABLE :: eofV(:,:)  ! matrix of eigenvectors V 
  REAL, ALLOCATABLE :: svals(:)   ! singular values
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
  rank = dim_ens - 1
  
  ! *** Generate full ensemble on filter-PE 0 ***
  WRITE (*, '(9x, a)') '--- generate ensemble from covariance matrix'
  WRITE (*, '(9x, a)') &
       '--- use rank reduction and 2nd order exact sampling (SEIK type)'
  WRITE (*, '(9x, a, i5)') '--- Ensemble size:  ', dim_ens
  WRITE (*, '(9x, a, i5)') '--- number of EOFs: ', rank

  ! allocate memory for temporary fields
  ALLOCATE(eofV(dim, rank))
  ALLOCATE(svals(rank))
  IF (allocflag == 0) THEN
     ! count allocated memory
     CALL memcount(2, 'r', dim * rank + rank)
  END IF


! *************************************************
! *** Initialize initial state and covar matrix ***
! *************************************************
  
  WRITE(*,'(9x,a,a)') '--- Reading covariance information from ', TRIM(file_ini)

  s = 1
  stat(s) = NF90_OPEN(file_ini, NF90_NOWRITE, fileid)

  ! Read size of state vector
  s = s + 1
  stat(s) = NF90_INQ_DIMID(fileid, 'dim_state', id_dim)
  s = s + 1
  stat(s) = NF90_Inquire_dimension(fileid, id_dim, len=dim_file)

  ! Read rank stored in file
  s = s + 1
  stat(s) = NF90_INQ_DIMID(fileid, 'rank', id_dim)
  s = s + 1
  stat(s) = NF90_Inquire_dimension(fileid, id_dim, len=rank_file)

  DO i = 1,  s
     IF (stat(i) /= NF90_NOERR) &
          WRITE(*, *) 'NetCDF error in reading dimensions from init file, no.', i
  END DO

  ! Check consistency of dimensions
  checkdim: IF (dim == dim_file .AND. rank_file >= rank) THEN

     ! Inquire IDs for mean state, singular vectors and values
     s = 1
     stat(s) = NF90_INQ_VARID(fileid, 'meanstate', id_state)
     s = s + 1
     stat(s) = NF90_INQ_VARID(fileid, 'u_svd', id_eofV)
     s = s + 1
     stat(s) = NF90_INQ_VARID(fileid, 'sigma', id_svals)

     ! Read initialization information
     s = s + 1
     stat(s) = NF90_GET_VAR(fileid, id_state, state)

     pos(2) = 1
     cnt(2) = rank
     pos(1) = 1
     cnt(1) = dim
     s = s + 1
     stat(s) = NF90_GET_VAR(fileid, id_eofV, eofV, start=pos, count=cnt)

     pos(1) = 1
     cnt(1) = rank
     s = s + 1
     stat(s) = NF90_GET_VAR(fileid, id_svals, svals, start=pos(1:1), count=cnt(1:1))

     s = s + 1
     stat(s) = NF90_CLOSE(fileid)

     DO i = 1,  s
        IF (stat(i) /= NF90_NOERR) &
             WRITE(*, *) 'NetCDF error in reading initialization file, no.', i
     END DO


! *****************************************
! *** Generate ensemble of model states ***
! *****************************************

     WRITE (*,'(9x, a)') '--- generate state ensemble'
     
     ! Use PDAF routine to generate ensemble from covariance matrix
     CALL PDAF_SampleEns(dim, dim_ens, eofV, svals, state, ens, 1, flag)
     
  ELSE

      ! *** Rank stored in file is smaller than requested EOF rank ***
     WRITE(*,*) 'Rank stored in file is smaller than requested EOF rank'

     stat(s) = NF90_CLOSE(fileid)
     STOP

  END IF checkdim


! ****************
! *** clean up ***
! ****************

  DEALLOCATE(svals, eofV)

END SUBROUTINE init_ens_eof
