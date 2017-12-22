! Copyright (c) 2004-2016 Lars Nerger
!
! This file is part of PDAF.
!
! PDAF is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version
! 3 of the License, or (at your option) any later version.
!
! PDAF is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with PDAF.  If not, see <http://www.gnu.org/licenses/>.
!
!$Id: PDAF_generate_rndmat.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_generate_rndmat - Generate random matrix with special properties
!
! !INTERFACE:
SUBROUTINE PDAF_generate_rndmat(dim, rndmat, mattype)

! !DESCRIPTION:
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
! computes Omega using GEMM from BLAS. All fields are 
! allocated once at their maximum required size.
! (On SGI O2K this is about a factor of 2.5 faster
! than the version applying BLAS DDOT, but requires
! more memory.)
!
! For omegatype=0 a deterministic omega is computed
! where the Housholder matrix of (1,...,1)' is operated
! on an identity matrix.
!
! !  This is a core routine of PDAF and 
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2002-01 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_memcounting, &
       ONLY: PDAF_memcount

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim       ! Size of matrix mat
  REAL, INTENT(out)   :: rndmat(dim, dim) ! Matrix
  INTEGER, INTENT(in) :: mattype   ! Select type of random matrix:
                                   !   (1) orthonormal random matrix
                                   !   (2) orthonormal with eigenvector (1,...,1)^T

! !CALLING SEQUENCE:
! Called by: PDAF_seik_omega
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: gemvTYPE (BLAS; dgemv or sgemv dependent on precision)
! Calls: larnvTYPE (BLAS; dlarnv or slarnv dependent on precision)
!EOP

!  *** local variables ***
  INTEGER :: iter, col, row          ! counters
  INTEGER :: i, j, k                 ! counters
  INTEGER :: seedset = 1             ! Choice of seed set for random numbers
  INTEGER :: dimrnd                  ! Size of random matrix to be generation at first part
  INTEGER, SAVE :: iseed(4)          ! seed array for random number routine
  REAL :: norm                       ! norm of random vector
  INTEGER :: pflag                   ! pointer flag
  INTEGER, SAVE :: first = 1         ! flag for init of random number seed
  REAL :: rndval                     ! temporary value for init of Householder matrix
  INTEGER, SAVE :: allocflag = 0     ! Flag for dynamic allocation
  REAL, ALLOCATABLE :: rndvec(:)     ! vector of random numbers
  REAL, ALLOCATABLE :: h_rndvec(:)   ! vector of random numbers
  REAL, ALLOCATABLE :: house(:,:)    ! Householder matrix
  REAL, ALLOCATABLE :: matUBB(:,:)   ! Temporary matrix
  REAL, POINTER :: mat_iter(:,:)     ! Pointer to temporary random array
  REAL, POINTER :: mat_itermin1(:,:) ! Pointer to temporary random array
  REAL, POINTER :: matU(:,:)         ! Pointer to temporary array
  REAL, POINTER :: matUB(:,:)        ! Pointer to temporary array
  REAL, POINTER :: matB(:,:)         ! Pointer to temporary array
  REAL, ALLOCATABLE, TARGET :: temp1(:,:)  ! Target array
  REAL, ALLOCATABLE, TARGET :: temp2(:,:)  ! Target array


! **********************
! *** INITIALIZATION ***
! **********************

  ! Determine size of matrix build through householder reflections
  randomega: IF (mattype == 1) THEN
     ! Random orthonormal matrix
     dimrnd = dim
  ELSE
     ! Random orthonormal matrix with eigenvector (1,...,1)^T
     dimrnd = dim - 1
  END IF randomega


! ******************************************
! *** Generate orthonormal random matrix ***
! ******************************************

  ! allocate fields
  ALLOCATE(rndvec(dim))
  ALLOCATE(house(dim + 1, dim))
  ALLOCATE(temp1(dim, dim), temp2(dim, dim))
  IF (allocflag == 0) THEN
     ! count allocated memory
     CALL PDAF_memcount(3, 'r', dim + (dim + 1) * dim + 2 * dim**2)
     allocflag = 1
  END IF

  ! set pointers
  mat_itermin1 => temp1
  mat_iter     => temp2
  pflag = 0

  ! Initialized seed for random number routine
  IF (first == 1) THEN
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
        ! Standard seed
        iseed(1) = 1000
        iseed(2) = 2034
        iseed(3) = 0
        iseed(4) = 3
     END IF
     first = 2
  END IF


! *** First step of iteration       ***  
! *** Determine mat_iter for iter=1 ***

  ! Get random number [-1,1]
  CALL larnvTYPE(2, iseed, 1, rndvec(1))
  
  IF (rndvec(1) >= 0.0) THEN
     mat_itermin1(1, 1) = +1.0
  ELSE
     mat_itermin1(1, 1) = -1.0
  END IF

! *** Iteration ***

  iteration: DO iter = 2, dimrnd

! Initialize new random vector
      
     ! Get random vector of dimension DIM (elements in [-1,1])
     CALL larnvTYPE(2, iseed, iter, rndvec(1:iter))

     ! Normalize random vector
     norm = 0.0
     DO col = 1, iter
        norm = norm + rndvec(col)**2
     END DO
     norm = SQRT(norm)
        
     DO col = 1, iter
        rndvec(col) = rndvec(col) / norm
     END DO

! Compute Householder matrix

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

! Compute matrix on this iteration stage

     ! First iter-1 columns
     CALL gemmTYPE ('n', 'n', iter, iter - 1, iter - 1, &
          1.0, house, dim + 1, mat_itermin1, dim, &
          0.0, mat_iter, dim)

     ! Final column
     DO row = 1, iter
        mat_iter(row, iter) = rndvec(row)
     END DO

! Adjust pointers to temporal OMEGA fields

     IF (pflag == 0) THEN
        mat_itermin1 => temp2
        mat_iter     => temp1
        pflag = 1
     ELSE IF (pflag == 1) THEN
        mat_itermin1 => temp1
        mat_iter     => temp2
        pflag = 0
     END IF

  END DO iteration


! ****************************************************
! *** Ensure eigenvector (1,...1,)^T for mattype=2 ***
! ****************************************************

  mattype2: IF (mattype == 1) THEN

     ! *** Generation of random matrix completed for mattype=1
     rndmat = mat_itermin1

  ELSE mattype2

     ! *** Complete generation of random matrix with eigenvector
     ! *** (1,...,1)^T by transformation with a basis that
     ! *** includes (1,...,1)^T. (We follow the description 
     ! *** Sakov and Oke, MWR 136, 1042 (2008)).

     NULLIFY(mat_iter, mat_itermin1)

     ALLOCATE(h_rndvec(dim))

! *** Complete initialization of random matrix with eigenvector ***
! *** (1,...,1)^T in the basis that includes (1,...,1)^T        ***

     IF (pflag == 0) THEN
        matU   => temp1
        matUB => temp2
     ELSE
        matU   => temp2
        matUB => temp1
     END if

     matUB(:,:) = 0.0
     matUB(1,1) = 1.0
     DO col = 2, dim
        DO row = 2, dim
           matUB(row, col) = matU(row - 1, col - 1)
        END DO
     END DO
     NULLIFY(matU)

! *** Generate orthonormal basis including (1,...,1)^T as leading vector ***
! *** We again use houesholder reflections.                              ***

     IF (pflag == 0) THEN
        matB => temp1
     ELSE
        matB => temp2
     END IF

     ! First column
     DO row = 1, dim
        matB(row, 1) = 1.0 / SQRT(REAL(dim))
     END DO

     ! columns 2 to dim
     buildB: DO col = 2, dim

        ! Get random vector of dimension DIM (elements in [0,1])
        CALL larnvTYPE(1, iseed, dim, rndvec)

        loopcols: DO i = 1, col - 1
           DO j = 1, dim
              DO k = 1, dim
                 house(k, j) = - matB(k,i) * matB(j,i)
              END DO
           END DO
           DO j = 1, dim
              house(j, j) = house(j, j) + 1.0
           END DO

           ! Apply house to random vector
           CALL gemvTYPE ('n', dim, dim, &
                1.0, house, dim+1, rndvec, 1, &
                0.0, h_rndvec, 1)
           rndvec = h_rndvec

        END DO loopcols

        ! Normalize vector
        norm = 0.0
        DO i = 1, iter
           norm = norm + h_rndvec(i)**2
        END DO
        norm = SQRT(norm)
        
        DO i = 1, iter
           h_rndvec(i) = h_rndvec(i) / norm
        END DO

        ! Inialize column of matB
        matB(:, col) = h_rndvec

     END DO buildB


! *** Final step: Transform random matrix  ***
! *** rndmat = matB matUB matB^T  ***

     ALLOCATE(matUBB(dim, dim))

     ! matUB * matB^T
     CALL gemmTYPE ('n', 't', dim, dim, dim, &
          1.0, matUB, dim, matB, dim, &
          0.0, matUBB, dim)

     ! matB * matUB * matB^T
     CALL gemmTYPE ('n', 'n', dim, dim, dim, &
          1.0, matB, dim, matUBB, dim, &
          0.0, rndmat, dim)

! *** CLEAN UP ***

     NULLIFY(matUB, matB)
     DEALLOCATE(matUBB)
     DEALLOCATE(h_rndvec)

  END IF mattype2


! ************************
! *** General clean up ***
! ************************

  DEALLOCATE(temp1, temp2)
  DEALLOCATE(rndvec, house)

END SUBROUTINE PDAF_generate_rndmat
