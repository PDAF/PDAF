!$Id: init_seik_dist.F90 1018 2010-07-14 09:25:43Z lnerger $
!BOP
!
! !ROUTINE: init_seik_dist --- Initialize distributed ensemble for SEIK
!
! !INTERFACE:
SUBROUTINE init_seik_dist(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF (SEIK):
!
! The routine initializes an ensemble of dim\_ens states
! by exact 2nd order sampling.
! States of the form
!         $x_i = x + sqrt(FAC) eofV (\Omega C^{-1})^T$
! fulfill the condition
!     $P = 1/(FAC)  \sum_{i=1}^{dim\_ens} (x_i - x)(x_i - x)^T$
! The matrix is initialized in the form of
! singular values and singular vectors.
!
! The routine is called by all filter processes and 
! initializes the ensemble for the local domain.
!
! This version is for the dummy model with domain 
! decomposition. This example performs the initialization 
! in a distributed form. Each ensemble of local sub-states 
! is directly initialized. (See init\_seik() for a 
! different variant of initialization.)
!
! !REVISION HISTORY:
! 2006-05 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE timer, &
       ONLY: timeit
  USE mod_memcount, &
       ONLY: memcount
  USE mod_parallel, &
       ONLY: mype_filter
  USE mod_assimilation, &
       ONLY: covartype
  USE mod_model, &
       ONLY: local_dims, dim_state

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize  
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
  REAL, INTENT(out) :: state_p(dim_p)            ! PE-local model state
  REAL, INTENT(out) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
  REAL, INTENT(out) :: ens_p(dim_p,dim_ens)      ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_filter_init (as U_ens_init)
! Calls: PDAF_seik_omega
! Calls: timeit
! Calls: memcount
! Calls: dgemm (BLAS)
!EOP

! *** local variables ***
  INTEGER :: i, row, col           ! counters
  INTEGER, SAVE :: allocflag = 0   ! Flag for memory counting
  REAL, ALLOCATABLE :: eofV_p(:,:) ! Matrix of PE-local eigenvectors V 
  REAL, ALLOCATABLE :: svals(:)    ! singular values
  REAL, ALLOCATABLE :: omega(:,:)  ! Matrix Omega
  INTEGER :: rank     ! Rank of approximated covariance matrix
  REAL :: fac         ! Square-root of dim_ens-1 or dim_ens
  ! variables and arrays for domain decomposition
  INTEGER :: offset   ! Row-offset according to domain decomposition


! **********************
! *** INITIALIZATION ***
! **********************

  ! *** Rank of matrix is ensemble size minus one
  rank = dim_ens - 1
  
  ! *** Generate local ensemble on current PE ***
  mype0: IF (mype_filter == 0) THEN
     WRITE (*, '(/9x, a)') 'Generate state ensemble from covariance matrix'
     WRITE (*, '(9x, a)') &
          '--- use rank reduction and 2nd order exact sampling (SEIK type)'
     WRITE (*, '(9x, a, i5)') '--- number of EOFs:', rank
  ENDIF mype0

  ! allocate memory for temporary fields
  ALLOCATE(eofV_p(dim_p, rank))
  ALLOCATE(svals(rank))
  ALLOCATE(omega(rank + 1, rank))
  IF (allocflag == 0) THEN
     ! count allocated memory
     CALL memcount(2, 'r', dim_p * rank + rank + rank * (rank + 1))
  END IF
  

! *************************************************
! *** Initialize initial state and covar matrix ***
! *************************************************

  ! Very simple initialization for dummy model
   
  ! Just set the entries of the state vector to 2.0
  state_p(1 : dim_p) = 2.0

  ! Set the initial singular vectors to one
  svals(1 : rank) = 1.0

  ! Set the initial ensemble to a part of the identity matrix
  ! Which columns need to be set to one depends on the process rank
  eofV_p(:, :) = 0.0
  offset = 0   ! Compute column offset
  DO i = 1, mype_filter
     offset = offset + local_dims(mype_filter)
  END DO
  DO col = 1 + offset, MIN(local_dims(mype_filter + 1) + offset, rank)
     eofV_p(col - offset, col) = 1.0
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
  
  WRITE (*, '(9x, a)') '--- generate ensemble of interpolating states'

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

  ! state_ens = state+ sqrt(FAC) eofV A^T
  ! FAC depends on the definiton of the factor in the ensemble
  ! covar. matrix which is defined by the variable COVARTYPE
  DO col = 1, rank + 1
     ens_p(1 : dim_p, col) = state_p(1 : dim_p)
  END DO

  IF (covartype == 1) THEN
     fac = SQRT(REAL(rank))
  ELSE
     fac = SQRT(REAL(rank + 1))
  END IF
  CALL dgemm('n', 't', dim_p, rank + 1, rank, &
       fac, eofV_p, dim_p, Omega, rank + 1, &
       1.0, ens_p, dim_p)

  CALL timeit(6, 'old')


! ****************
! *** clean up ***
! ****************

  DEALLOCATE(svals, eofV_p, omega)

END SUBROUTINE init_seik_dist
