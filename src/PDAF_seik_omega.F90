! Copyright (c) 2004-2025 Lars Nerger
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
!
!> Generate random matrix with special properties
!!
!! Generate a transformation matrix OMEGA for
!! the generation and transformation of the 
!! ensemble in the SEIK and LSEIK filter.
!! Generated is a uniform orthogonal matrix OMEGA
!! with R columns orthonormal in $R^{r+1}$
!! and orthogonal to (1,...,1)' by iteratively 
!! applying the Householder matrix onto random 
!! vectors distributed uniformly on the unit sphere.
!!
!! This version initializes at each iteration step
!! the whole Householder matrix and subsequently
!! computes Omega using GEMM from BLAS. All fields are 
!! allocated once at their maximum required size.
!! (On SGI O2K this is about a factor of 2.5 faster
!! than the version applying BLAS DDOT, but requires
!! more memory.)
!!
!! For Omegatype=0 a deterministic Omega is computed
!! where the Housholder matrix of (1,...,1)' is operated
!! on an identity matrix.
!!
!! !  This is a core routine of PDAF and 
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2002-01 - Lars Nerger - Initial code
!! * Later revisions - see svn log
!!
SUBROUTINE PDAF_seik_Omega(rank, Omega, Omegatype, screen)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_mod_filtermpi, &
       ONLY: mype

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: rank                !< Approximated rank of covar matrix
  REAL, INTENT(inout) :: Omega(rank+1, rank) !< Matrix Omega
  INTEGER, INTENT(in) :: Omegatype           !< Select type of Omega:
                                             !<   (1) generated from random vectors
                                             !<   (0) generated from deterministic vectors
                                             !< (other) product of matrix from (2) with
                                             !<      orthonormal random matrix orthogonal (1....1)T
  INTEGER, INTENT(in) :: screen              !< Verbosity flag

!  *** local variables ***
  INTEGER :: col, row              ! counters
  REAL :: rndval                   ! temporary value for init of Householder matrix
  REAL :: rndnum                   ! Value of randum entry
  REAL, ALLOCATABLE :: house(:,:)  ! Householder matrix
  REAL, POINTER :: rndmat(:,:)         ! Pointer to temporary Omega field


  randOmega: IF (Omegatype == 0) THEN 
! *************************************************
! *** Generate deterministic Omega as           ***
! *** Householder matrix associated with the    ***
! *** vector  1/sqrt(rank) (1,...,1)^T          ***
! *************************************************

     IF (mype == 0 .AND. screen > 0) &
          WRITE (*,'(a, 5x, a)') 'PDAF','--- Compute deterministic Omega'

     rndnum = 1.0 / SQRT(REAL(rank + 1))

     ! First r rows
     rndval = - rndnum * rndnum / (rndnum + 1.0)
     Omegacolb: DO col = 1, rank
        Omegarowb: DO row = 1, rank
           Omega(row, col) = rndval
        END DO Omegarowb
     END DO Omegacolb
     
     DO col = 1, rank
        Omega(col, col) = Omega(col, col) + 1.0
     END DO

     ! Last row
     rndval = - (rndnum + 1.0) * rndnum / (rndnum + 1.0)
     Omegacolc: DO col = 1, rank
        Omega(rank + 1, col) = rndval
     END DO Omegacolc

  ELSEIF (Omegatype == 1) THEN randOmega
! ****************************************
! *** Generate Omega by random vectors ***
! ****************************************

     IF (mype == 0 .AND. screen > 0) &
          WRITE (*,'(a, 5x, a)') 'PDAF','--- Compute random Omega'

! *** Initialization ***

     ! Allocate fields
     ALLOCATE(house(rank + 1, rank))
     ALLOCATE(rndmat(rank, rank))

! *** Initialize orthonormal random matrix of size rank*rank ***

     CALL PDAF_generate_rndmat(rank, rndmat, 1)

! *** Project rndmat orthogonal to (1,...,1)^T ***

     ! *** Compute Householder matrix ***

     rndnum = 1.0 / SQRT(REAL(rank + 1))

     ! First r rows
     rndval = - rndnum * rndnum / (rndnum + 1.0)
     housecol: DO col = 1, rank
        houserow: DO row = 1, rank
           house(row, col) = rndval
        END DO houserow
     END DO housecol
     
     DO col = 1, rank
        house(col, col) = house(col, col) + 1.0
     END DO
     
     ! Last row
     rndval = - (rndnum + 1.0) * rndnum / (rndnum + 1.0)
     housecolb: DO col = 1, rank
        house(rank + 1, col) = rndval
     END DO housecolb

     ! *** Complete Omega: house * rndmat ***

     CALL gemmTYPE ('n', 'n', rank + 1, rank, rank, &
          1.0, house, rank + 1, rndmat, rank, &
          0.0, Omega, rank + 1)

! *** CLEAN UP ***

     DEALLOCATE(house)
     DEALLOCATE(rndmat)

  ELSE randOmega
! *** Generate Omega as a product of a deterministic  ***
! *** transformation with an orthonormal random       ***
! *** matrix that preserves the mean.                 ***
! *** 1. The deterministic matrix matrix given by the ***
! *** householder matrix from Omegatype=0.            ***
! *** 2. The random matrix is generated analogously   ***
! *** to Omegatype=1 followed by a transformation to  ***
! *** ensure the (1,....,1)^T is an eigenvector of    ***
! *** the matrix.                                     ***

     IF (mype == 0 .AND. screen > 0) &
          WRITE (*,'(a, 5x, a)') 'PDAF','--- Compute product Omega'

! *** Initialization ***

    ! Allocate fields
     ALLOCATE(house(rank + 1, rank))
     ALLOCATE(rndmat(rank, rank))

! *** 1. Deterministic part:                            ***
! *** Compute Householder matrix associated with the    ***
! *** vector  1/sqrt(rank) (1,...,1)^T                  ***
! *** (this is the transformation used for Omegatype=0) ***

     rndnum = 1.0 / SQRT(REAL(rank + 1))

     ! First r rows
     rndval = - rndnum * rndnum / (rndnum + 1.0)
     housecolc: DO col = 1, rank
        houserowc: DO row = 1, rank
           house(row, col) = rndval
        END DO houserowc
     END DO housecolc
     
     DO col = 1, rank
        house(col, col) = house(col, col) + 1.0
     END DO

     ! Last row
     rndval = - (rndnum + 1.0) * rndnum / (rndnum + 1.0)
     housecalc: DO col = 1, rank
        house(rank + 1, col) = rndval
     END DO housecalc

! *** 2. Random part: 
! *** Initialize orthonormal random matrix of size rank*rank 
! *** with eigenvector (1,...,1)^T

     CALL PDAF_generate_rndmat(rank, rndmat, 2)

! *** 3. Multiply deterministic and random parts: 

     CALL gemmTYPE ('n', 'n', rank+1, rank, rank, &
          1.0, house, rank+1, rndmat, rank, &
          0.0, Omega, rank+1)

! *** CLEAN UP ***

     DEALLOCATE(house, rndmat)

  END IF randOmega

END SUBROUTINE PDAF_seik_Omega
