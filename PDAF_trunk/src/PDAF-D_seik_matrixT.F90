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
!$Id: PDAF-D_seik_matrixT.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_seik_matrixT --- Operate matrix T on A as AT
!
! !INTERFACE:
SUBROUTINE PDAF_seik_matrixT(dim, dim_ens, A)

! !DESCRIPTION:
! Operate matrix T on another matrix as
!         $A = A T$.
!
! T is a dim_ens x (dim_ens-1) matrix with zero column sums.
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
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim               ! dimension of states
  INTEGER, INTENT(in) :: dim_ens           ! Size of ensemble
  REAL, INTENT(inout) :: A(dim, dim_ens)   ! Input/output matrix

! !CALLING SEQUENCE:
! Called by: PDAF_seik_analysis
! Called by: PDAF_seik_analysis_newT
! Calls PDAF_memcount
!EOP
  
! *** local variables ***
  INTEGER :: row, col  ! counters
  INTEGER :: typeT = 0 ! Which type of T
  REAL :: invdimens    ! Inverse of ensemble size
  INTEGER, SAVE :: allocflag = 0  ! Flag for dynamic allocation
  REAL, ALLOCATABLE :: rowmean(:) ! Mean values of rows of A

!$OMP THREADPRIVATE(allocflag)


! **********************
! *** INITIALIZATION ***
! **********************

  ALLOCATE(rowmean(dim))
  IF (allocflag == 0) THEN
     ! count allocated memory
     CALL PDAF_memcount(3, 'r', dim)
     allocflag = 1
  END IF
  rowmean   = 0.0
  invdimens = 1.0 / REAL(dim_ens)

  IF (typeT == 0) THEN

     ! *** Compute row means of A ***
     DO col = 1, dim_ens
        DO row = 1, dim
           rowmean(row) = rowmean(row) + A(row, col)
        END DO
     END DO
     rowmean = invdimens * rowmean

  ELSE

     ! *** Get last column of A ***
     DO row = 1, dim
        rowmean(row) = A(row, dim_ens)
     END DO
     
  END IF


! **********************************************
! ***  Operate T on A                        ***
! ***                                        ***
! *** v^TT = (v_1-mean(v), ... ,v_r-mean(v)) ***
! *** with v = (v_1,v_2, ... ,r_(r+1))       ***
! **********************************************

  DO col = 1, dim_ens - 1
     DO row = 1, dim
        A(row, col) = A(row, col) - rowmean(row)
     END DO
  END DO
  
  DO row = 1, dim
     A(row, dim_ens) = 0.0
  END DO


! ********************
! *** FINISHING UP ***
! ********************

  DEALLOCATE(rowmean)

END SUBROUTINE PDAF_seik_matrixT
