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
!$Id: PDAF-D_estkf_AOmega.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_estkf_AOmega --- Operate matrix Omega on A as AOmega
!
! !INTERFACE:
SUBROUTINE PDAF_estkf_AOmega(dim, dim_ens, A)

! !DESCRIPTION:
! Operate matrix Omega on another matrix as
!         $A = A Omega$.
!
! Omega is a dim_ens x (dim_ens-1) matrix with matximum
! rank and zero column sums. It is computed by the 
! Householder reflection associate with the vector
! (N-1)^(-1) (1,...,1)^T.
!
! The values of Omega are
!    1 - 1 / (N (1/sqrt(N) + 1)) for i=j
!    - 1 / (N (1/sqrt(N) + 1)) for i/=j, i<N
!    - 1 / sqrt(N) for i=N
!
! In this routine the product A Omega is implemented as
! operations:
! 1. Compute the row sums of A
! 2. Normalize row sums by 1/(sqrt(N) + N)
! 3. Subtract value of last column multiplied by 1/(1+sqrt(N))
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2011-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim               ! dimension of states
  INTEGER, INTENT(in) :: dim_ens           ! Size of ensemble
  REAL, INTENT(inout) :: A(dim, dim_ens)   ! Input/output matrix

! !CALLING SEQUENCE:
! Called by: PDAF_estkf_analysis
! Calls PDAF_memcount
!EOP
  
! *** local variables ***
  INTEGER :: row, col  ! Counters
  REAL :: normsum      ! Normalization for row sum
  REAL :: normlast     ! Normalization for last column
  REAL :: val          ! Temporary variable
  INTEGER, SAVE :: allocflag = 0  ! Flag for dynamic allocation
  REAL, ALLOCATABLE :: rowsums(:) ! Row sums of A

!$OMP threadprivate(allocflag)


! **********************
! *** INITIALIZATION ***
! **********************

  ALLOCATE(rowsums(dim))
  IF (allocflag == 0) THEN
     ! count allocated memory
     CALL PDAF_memcount(3, 'r', dim)
     allocflag = 1
  END IF
  rowsums   = 0.0

  ! Initialize normalization values
  normsum = 1.0 / REAL(dim_ens) / (1.0/SQRT(REAL(dim_ens))+1.0)
  normlast = 1.0 / (1.0 + SQRT(REAL(dim_ens)))


  ! *** Compute row sums of A ***
  DO col = 1, dim_ens
     DO row = 1, dim
        rowsums(row) = rowsums(row) + A(row, col)
     END DO
  END DO

  ! Scale by NORMSUM
  rowsums = normsum * rowsums

  ! Substract scale value for last column
  DO row = 1, dim
     val = A(row, dim_ens) * normlast
     rowsums(row) = rowsums(row) + val
  END DO


! **********************************************
! ***  Operate Omega on A                    ***
! **********************************************

  DO col = 1, dim_ens - 1
     DO row = 1, dim
        A(row, col) = A(row, col) - rowsums(row)
     END DO
  END DO
  
  DO row = 1, dim
     A(row, dim_ens) = 0.0
  END DO


! ********************
! *** FINISHING UP ***
! ********************

  DEALLOCATE(rowsums)

END SUBROUTINE PDAF_estkf_AOmega
