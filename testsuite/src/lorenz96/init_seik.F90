!$Id: init_seik.F90 841 2010-02-01 12:46:02Z lnerger $
!BOP
!
! !ROUTINE: init_seik --- Initialize ensemble for SEIK
!
! !INTERFACE:
SUBROUTINE init_seik(filtertype, dim, dim_ens, state, Uinv, &
     ens, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF (SEIK):
!
! The routine is called when the filter is
! initialized in PDAF\_filter\_init.  It has
! to initialize an ensemble of dim\_ens model 
! states. Here, we supply two methods to 
! initialize the ensemble:
! 1. By second-order exact sampling. 
!    This follows Pham (2001) and was used
!    in most of our papers.
! 2. By random sampling form a long state
!    trajectory. This method is often described
!    In papers on the EnKF and the ensemble
!    square-root filters.
!
! This version is for the Lorenz96 model
! without parallelization.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE timer, &
       ONLY: timeit
  USE mod_memcount, &
       ONLY: memcount
  USE mod_assimilation, &
       ONLY: covartype, file_ini, type_ensinit

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype          ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim                 ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens             ! Size of ensemble
  REAL, INTENT(inout) :: state(dim)          ! PE-local model state
  ! It is not necessary to initialize the array 'state' for SEIK. 
  ! It is available here only for convenience and can be used freely.
  REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
  REAL, INTENT(out)   :: ens(dim, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag             ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_filter_init    (as U_ens_init)
! Calls: PDAF_seik_omega
! Calls: timeit
! Calls: memcount
! Calls: dgemm (BLAS)
!EOP

! *** local variables ***
!   INTEGER :: i, s, row, col       ! counters
!   INTEGER, SAVE :: allocflag = 0  ! Flag for memory counting
!   REAL, ALLOCATABLE :: eofV(:,:)  ! matrix of eigenvectors V 
!   REAL, ALLOCATABLE :: svals(:)   ! singular values
!   REAL, ALLOCATABLE :: omega(:,:) ! Matrix Omega
!   REAL :: fac                     ! Square-root of dim_eof+1 or dim_eof
!   INTEGER :: dim_file             ! State dimension in file
!   INTEGER :: rank                 ! Rank of approximated covariance matrix
!   INTEGER :: rank_file            ! Rank of covariance matrix stored in file
!   INTEGER :: stat(50)             ! Array for status flag
!   INTEGER :: fileid               ! ID for NetCDF file
!   INTEGER :: id_svals, id_eofV    ! IDs for fields
!   INTEGER :: id_state             ! ID for field
!   INTEGER :: id_dim               ! ID for dimension


! **********************
! *** INITIALIZATION ***
! **********************
  
  ! *** Generate full ensemble ***
  WRITE (*, '(/9x, a)') 'Initialize state ensemble'

  CALL timeit(6, 'new')

  IF (TRIM(type_ensinit) == 'eof') THEN
     ! Initialize by 2nd-order exact sampling from EOFs
     CALL init_ens_eof(dim, dim_ens, state, ens, flag)
  ELSE IF (TRIM(type_ensinit) == 'rnd') THEN
     ! Initialize by random sampling from state trajectory
     CALL init_ens_rnd(dim, dim_ens, state, ens, flag)
  END IF

  CALL timeit(6, 'old')


END SUBROUTINE init_seik
