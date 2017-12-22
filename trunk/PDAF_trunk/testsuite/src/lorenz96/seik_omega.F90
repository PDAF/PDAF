!$Id: seik_omega.F90 1160 2011-09-14 09:32:08Z lnerger $
!BOP
!
! !ROUTINE: seik_omega - Generate random matrix with special properties
!
! !INTERFACE:
SUBROUTINE seik_omega(rank, omega, omegatype)

! !DESCRIPTION:
! This routine is a copy of the routine PDAF_seik_omega
! from PDAF. It is extended to use different seed set for 
! the generation of random numbers, specified by 'seedset'
! for mod_assimilation. In addition, the routine is 
! modified to run outside of PDAF. 
!
! Generate a transformation matrix OMEGA for
! the generation and transformation of the 
! ensemble in the SEIK and LSEIK filter.
! Generated is a uniform orthogonal matrix OMEGA
! with R columns orthonormal in $R^{r+1}$
! and orthogonal to (1,...,1)' by iteratively 
! applying the Householder matrix onto random 
! vectors distributed uniformly on the unit sphere.
!
! This version initializes at each iteration step
! the whole Householder matrix and subsequently
! computes Omega using DGEMM from BLAS. All fields are 
! allocated once at their maximum required size.
! (On SGI O2K this is about a factor of 2.5 faster
! than the version applying BLAS DDOT, but requires
! more memory.)
!
! For omegatype/=1 a deterministic omega is computed
! where the Housholder matrix of (1,...,1)' is operated
! on an identity matrix.
!
! !REVISION HISTORY:
! 2002-01 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_memcount, &
       ONLY: memcount
  USE mod_assimilation, &
       ONLY: seedset

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: rank      ! Approximated rank of covar matrix
  REAL, INTENT(inout) :: omega(rank+1, rank) ! Matrix Omega
  INTEGER, INTENT(in) :: omegatype ! Select type of omega:
                                   !     (1) generated from random vectors
                                   ! (other) generated from deterministic vectors

! !CALLING SEQUENCE:
! Called by: init_ens_eof
! Calls: dlarnv (LAPACK)
! Calls: dgemm  (BLAS)
!EOP

!  *** local variables ***
  INTEGER :: iter, col, row        ! counters
  INTEGER :: iseed(4)              ! seed array for random number routine
  REAL :: norm                     ! norm of random vector
  INTEGER :: pflag                 ! pointer flag
  INTEGER, SAVE :: first = 1       ! flag for init of random number seed
  REAL :: rndval                   ! temporary value for init of Householder matrix
  REAL :: rndnum                   ! Value of randum entry
  INTEGER, SAVE :: allocflag = 0   ! Flag for dynamic allocation
  REAL, ALLOCATABLE :: rndvec(:)   ! vector of random numbers
  REAL, ALLOCATABLE :: house(:,:)  ! Row of the Householder matrix
  REAL, POINTER :: omega_iter(:,:)         ! Pointer to temporary Omega field
  REAL, POINTER :: omega_itermin1(:,:)     ! Pointer to temporary Omega field
  REAL, ALLOCATABLE, TARGET :: temp1(:,:)  ! fields holding temporary Omega
  REAL, ALLOCATABLE, TARGET :: temp2(:,:)  ! fields holding temporary Omega


! **********************
! *** INITIALIZATION ***
! **********************

  randomega: IF (omegatype == 1) THEN
     ! *** Generate omega by random vectors ***

     WRITE (*, '(9x,a,i3)') &
          '--- Compute random Omega using seed set no. ', seedset

     ! allocate fields
     ALLOCATE(rndvec(rank))
     ALLOCATE(house(rank + 1, rank))
     ALLOCATE(temp1(rank, rank), temp2(rank, rank))
     IF (allocflag == 0) THEN
        ! count allocated memory
        CALL memcount(2, 'r', rank + (rank + 1) * rank + 2 * rank**2)
        allocflag = 1
     END IF

     ! set pointers
     omega_itermin1 => temp1
     omega_iter     => temp2
     pflag = 0

     ! Initialized seed for random number routine
     IF (seedset == 2) THEN
        iseed(1)=1
        iseed(2)=5
        iseed(3)=7
        iseed(4)=9
     ELSE IF (seedset == 3) THEN
        iseed(1)=2
        iseed(2)=5
        iseed(3)=7
        iseed(4)=9
     ELSE IF (seedset == 4) THEN
        iseed(1)=1
        iseed(2)=6
        iseed(3)=7
        iseed(4)=9
     ELSE IF (seedset == 5) THEN
        iseed(1)=1
        iseed(2)=5
        iseed(3)=8
        iseed(4)=9
     ELSE IF (seedset == 6) THEN
        iseed(1)=2
        iseed(2)=5
        iseed(3)=8
        iseed(4)=9
     ELSE IF (seedset == 7) THEN
        iseed(1)=2
        iseed(2)=6
        iseed(3)=8
        iseed(4)=9
     ELSE IF (seedset == 8) THEN
        iseed(1)=2
        iseed(2)=6
        iseed(3)=8
        iseed(4)=11
     ELSE IF (seedset == 9) THEN
        iseed(1)=3
        iseed(2)=6
        iseed(3)=8
        iseed(4)=11
     ELSE IF (seedset == 10) THEN
        iseed(1)=3
        iseed(2)=7
        iseed(3)=8
        iseed(4)=11
     ELSE IF (seedset == 11) THEN
        iseed(1)=13
        iseed(2)=7
        iseed(3)=8
        iseed(4)=11
     ELSE IF (seedset == 12) THEN
        iseed(1)=13
        iseed(2)=11
        iseed(3)=8
        iseed(4)=11
     ELSE IF (seedset == 13) THEN
        iseed(1)=13
        iseed(2)=13
        iseed(3)=8
        iseed(4)=11
     ELSE IF (seedset == 14) THEN
        iseed(1)=13
        iseed(2)=13
        iseed(3)=17
        iseed(4)=11
     ELSE IF (seedset == 15) THEN
        iseed(1)=13
        iseed(2)=13
        iseed(3)=19
        iseed(4)=11
     ELSE IF (seedset == 16) THEN
        iseed(1)=15
        iseed(2)=13
        iseed(3)=19
        iseed(4)=11
     ELSE IF (seedset == 17) THEN
        iseed(1)=15
        iseed(2)=135
        iseed(3)=19
        iseed(4)=11
     ELSE IF (seedset == 18) THEN
        iseed(1)=19
        iseed(2)=135
        iseed(3)=19
        iseed(4)=11
     ELSE IF (seedset == 19) THEN
        iseed(1)=19
        iseed(2)=135
        iseed(3)=19
        iseed(4)=17
     ELSE IF (seedset == 20) THEN
        iseed(1)=15
        iseed(2)=15
        iseed(3)=47
        iseed(4)=17
     ELSE
        ! This is the seed used in PDAF_seik_omega
        iseed(1) = 1000
        iseed(2) = 2034
        iseed(3) = 0
        iseed(4) = 3
     END IF


! ***************************************
! *** First step of iteration         ***  
! *** Determine omega_iter for iter=1 ***
! ***************************************

     ! Get random number [-1,1]
     CALL dlarnv(2, iseed, 1, rndvec(1))

     IF (rndvec(1) >= 0.0) THEN
        omega_itermin1(1, 1) = +1.0
     ELSE
        omega_itermin1(1, 1) = -1.0
     END IF

! *****************
! *** Iteration ***
! *****************

     iteration: DO iter = 2, rank

! *** Initialize new random vector ***
      
        ! Get random vector of dimension DIM (elements in [-1,1])
        CALL dlarnv(2, iseed, iter, rndvec(1:iter))

        ! Normalize random vector
        norm = 0.0
        DO col = 1, iter
           norm = norm + rndvec(col)**2
        END DO
        norm = SQRT(norm)
        
        DO col = 1, iter
           rndvec(col) = rndvec(col) / norm
        END DO

! *** Compute Householder matrix ***

        ! First ITER-1 rows
        rndval = 1.0 / (ABS(rndvec(iter)) + 1.0)
        housecol: DO col = 1, iter - 1
           houserow: DO row = 1,iter - 1
              house(row, col) = - rndvec(row) * rndvec(col) * rndval
           END DO houserow
        END DO housecol
        
        DO col = 1, iter - 1
           house(col, col) = house(col, col) + 1.0
        END DO

        ! Last row
        housecol2: DO col = 1, iter - 1
           house(iter, col) = - (rndvec(iter) + SIGN(1.0, rndvec(iter))) &
                * rndvec(col) * rndval
        END DO housecol2

! *** Compute omega on this iteration stage ***

        ! First iter-1 columns
        CALL dgemm ('n', 'n', iter, iter - 1, iter - 1, &
             1.0, house, rank + 1, omega_itermin1, rank, &
             0.0, omega_iter, rank)

        ! Final column
        DO row = 1, iter
           omega_iter(row, iter) = rndvec(row)
        END DO

! *** Adjust pointers to temporal OMEGA fields ***

        IF (pflag == 0) THEN
           omega_itermin1 => temp2
           omega_iter     => temp1
           pflag = 1
        ELSE IF (pflag == 1) THEN
           omega_itermin1 => temp1
           omega_iter     => temp2
           pflag = 0
        END IF

     END DO iteration


! ********************************************
! ***            Final step                ***
! *** Projecting orthogonal to (1,...,1)^T ***
! ********************************************
    
     rndvec(1) = 1.0 / SQRT(REAL(rank + 1))

! *** Compute Householder matrix ***

     ! First r rows
     rndval = - rndvec(1) * rndvec(1) / (rndvec(1) + 1.0)
     housecolb: DO col = 1, rank
        houserowb: DO row = 1, rank
           house(row, col) = rndval
        END DO houserowb
     END DO housecolb
     
     DO col = 1, rank
        house(col, col) = house(col, col) + 1.0
     END DO
     
     ! Last row
     rndval = - (rndvec(1) + 1.0) * rndvec(1) / (rndvec(1) + 1.0)
     housecolc: DO col = 1, rank
        house(rank + 1, col) = rndval
     END DO housecolc

! *** Compute omega ***

     ! First iter-1 columns
     CALL dgemm ('n', 'n', rank + 1, rank, rank, &
          1.0, house, rank + 1, omega_itermin1, rank, &
          0.0, omega, rank + 1)

! *** CLEAN UP ***

     NULLIFY(omega_itermin1, omega_iter)
     DEALLOCATE(temp1, temp2)
     DEALLOCATE(rndvec, house)

  ELSE randomega
     ! *** Generate Omega by deterministic vectors ***
     ! *** given by last step of upper recursion   ***
     ! *** with omega_itermin1 being the identity  ***

!     if (mype==0) write (*,'(8x,a)') '--- Compute deterministic Omega'

     rndnum = 1.0 / SQRT(REAL(rank + 1))

! *** Compute Householder matrix ***

     ! First r rows
     rndval = - rndnum * rndnum / (rndnum + 1.0)
     omegacolb: DO col = 1, rank
        omegarowb: DO row = 1, rank
           omega(row, col) = rndval
        END DO omegarowb
     END DO omegacolb
     
     DO col = 1, rank
        omega(col, col) = omega(col, col) + 1.0
     END DO

     ! Last row
     rndval = - (rndnum + 1.0) * rndnum / (rndnum + 1.0)
     omegacolc: DO col = 1, rank
        omega(rank + 1, col) = rndval
     END DO omegacolc

  END IF randomega

END SUBROUTINE seik_omega
