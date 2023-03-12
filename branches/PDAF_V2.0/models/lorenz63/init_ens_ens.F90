!$Id: init_ens_eof.F90 178 2019-07-03 15:07:04Z lnerger $
!BOP
!
! !ROUTINE: init_ens_ens --- Read initial ensemble from previous DA run
!
! !INTERFACE:
SUBROUTINE init_ens_ens(dim, dim_ens, state, ens, flag)

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
! This version is for the Lorenz63 model
! without parallelization.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE netcdf
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
  INTEGER :: dim_file             ! State dimension in file
  INTEGER :: dim_ens_file            ! Rank of covariance matrix stored in file
  INTEGER :: stat(50)             ! Array for status flag
  INTEGER :: fileid               ! ID for NetCDF file
  INTEGER :: id_ens               ! ID for field
  INTEGER :: id_dim               ! ID for dimension
  INTEGER :: pos(3)               ! Position index for reading
  INTEGER :: cnt(3)               ! Count index for reading


! **********************
! *** INITIALIZATION ***
! **********************

  ! *** Generate full ensemble on filter-PE 0 ***
  WRITE (*, '(9x, a)') '--- read ensemble from files'
  WRITE (*, '(9x, a, i5)') '--- Ensemble size:  ', dim_ens


! *****************************************************************
! *** Initialize initial ensemble from assimilation output file ***
! *****************************************************************
  
  WRITE(*,'(9x,a,a)') '--- Reading ensemble from ', TRIM(file_ini)

  s = 1
  stat(s) = NF90_OPEN(file_ini, NF90_NOWRITE, fileid)

! Read size of state vector
  s = s + 1
  stat(s) = NF90_INQ_DIMID(fileid, 'dim_state', id_dim)
  s = s + 1
  stat(s) = NF90_Inquire_dimension(fileid, id_dim, len=dim_file)

  ! Read rank stored in file
  s = s + 1
  stat(s) = NF90_INQ_DIMID(fileid, 'dim_ens', id_dim)
  s = s + 1
  stat(s) = NF90_Inquire_dimension(fileid, id_dim, len=dim_ens_file)

  DO i = 1,  s
     IF (stat(i) /= NF90_NOERR) &
          WRITE(*, *) 'NetCDF error in reading dimensions from init file, no.', i
  END DO

  ! Check consistency of dimensions
  checkdim: IF (dim == dim_file .AND. dim_ens_file >= dim_ens) THEN

     ! Inquire IDs for mean state, singular vectors and values
     s = 1
     stat(s) = NF90_INQ_VARID(fileid, 'ens_ana', id_ens)

     ! Read initialization information
     pos(3) = 1
     cnt(3) = 1
     pos(2) = 1
     cnt(2) = dim_ens
     pos(1) = 1
     cnt(1) = dim
     s = s + 1
     stat(s) = NF90_GET_VAR(fileid, id_ens, ens, start=pos(1:3), count=cnt(1:3))

     s = s + 1
     stat(s) = NF90_CLOSE(fileid)

     DO i = 1,  s
        IF (stat(i) /= NF90_NOERR) &
             WRITE(*, *) 'NetCDF error in reading initialization file, no.', i
     END DO
     
  ELSE

      ! *** Rank stored in file is smaller than requested EOF rank ***
     WRITE(*,*) 'Ensemble stored in file is smaller than requested ensemble size'

     stat(s) = NF90_CLOSE(fileid)
     STOP

  END IF checkdim

END SUBROUTINE init_ens_ens
