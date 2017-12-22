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
!$Id: PDAF-D_seik_uinv.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_seik_Uinv - Initialize matrix Uinv from matrix T
!
! !INTERFACE:
SUBROUTINE PDAF_seik_Uinv(rank, Uinv)

! !DESCRIPTION:
! Initialize matrix Uinv by
! $U^{-1} = FAC\ T^T T$
! where $FAC$ = rank+1 for a covariance matrix with factor
! (rank+1)$^{-1}$ and $FAC$ = rank for a covariance matrix
! with factor rank$^{-1}$.
!
! There are two proposed forms of T (ensemble size N):\\
! typeT=0: diag(T)=1-1/N; nondiag(T)=-1/N; 
!          last row= -1/N\\
! typeT=1: diag(T)=1; nondiag(T)=0; last row = -1\\
! We typically use TypeT=0, but both variants are implemented.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2002-01 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filter, &
       ONLY: Nm1vsN

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: rank             ! Rank of initial covariance matrix
  REAL, INTENT(inout) :: Uinv(rank, rank) ! Inverse of matrix U
!EOP
  
! *** local variables ***
  INTEGER :: row, col       ! counters
  INTEGER :: typeT = 0      ! Choose type of T
  REAL :: rdivrp1, r2divrp1 ! scaling factors for Uinv


! ***********************
! *** Initialize Uinv ***
! ***********************

  ttype: IF (typeT == 0) THEN

     ! Scaling factors
     IF (Nm1vsN == 1) THEN
        ! For ensemble covariance matrix with factor (N-1)^-1
        rdivrp1 = REAL(rank) / REAL(rank + 1)
        r2divrp1 = REAL(rank)**2 / REAL(rank + 1)
     ELSE
        ! For ensemble covariance matrix with factor N^-1
        rdivrp1 = 1
        r2divrp1 = REAL(rank)
     END IF

     DO col = 1, rank
        ! non-diagonal elements - upper triangle
        DO row = 1, col - 1
           Uinv(row, col) = - rdivrp1
        END DO
        ! diagonal
        Uinv(col, col) = r2divrp1
        ! non-diagonal elements - lower triangle
        DO row = col + 1, rank
           Uinv(row, col) = -rdivrp1
        END DO
     END DO

  ELSE ttype

     ! Scaling factors
     IF (Nm1vsN == 1) THEN
        ! For ensemble covariance matrix with factor (N-1)^-1
        rdivrp1 = REAL(rank) / REAL(rank + 1)
     ELSE
        ! For ensemble covariance matrix with factor N^-1
        rdivrp1 = 1
     END IF

     DO col = 1, rank
        ! non-diagonal elements - upper triangle
        DO row = 1, col - 1
           Uinv(row, col) = rdivrp1
        END DO
        ! diagonal
        Uinv(col, col) = 2.0 * rdivrp1
        ! non-diagonal elements - lower triangle
        DO row = col + 1, rank
           Uinv(row, col) = rdivrp1
        END DO
     END DO
     
  END IF ttype

! ********************
! *** FINISHING UP ***
! ********************

END SUBROUTINE PDAF_seik_Uinv
