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
!$Id: PDAF-D_etkf_Tright.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_etkf_Tright --- Operate matrix T from right to some matrix
!
! !INTERFACE:
SUBROUTINE PDAF_etkf_Tright(dim, dim_ens, A)

! !DESCRIPTION:
! Operate matrix T on another matrix as
!         $A <- A T$.
!
! T is a symmetric dim_ens x dim_ens matrix with zero column 
! sums defined as:  
!            diag(T)=1-1/dim_ens; nondiag(T)=-1/dim_ens
! Applied from the right to a matrix A, T subtracts the 
! mean over the columns of A from this matrix. In this routine,
! T is not explicitly constructed, but its operation on the 
! matrix A is implemented.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2009-07 - Lars Nerger - Initial code
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
! Called by: PDAF_etkf_analysis
! Called by: PDAF_etkf_analysis_T
! Calls PDAF_memcount
!EOP
  
! *** local variables ***
  INTEGER :: row, col  ! counters
  REAL :: invdimens    ! Inverse of ensemble size
  INTEGER, SAVE :: allocflag = 0  ! Flag for dynamic allocation
  REAL, ALLOCATABLE :: rowmean(:) ! Mean values of rows of A


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

  ! *** Compute row means of A ***
  DO col = 1, dim_ens
     DO row = 1, dim
        rowmean(row) = rowmean(row) + A(row, col)
     END DO
  END DO
  rowmean = invdimens * rowmean


! **********************************************
! ***  Operate T on A                        ***
! ***                                        ***
! *** v^TT = (v_1-mean(v), ... ,v_r-mean(v)) ***
! *** with v = (v_1,v_2, ... ,r_N)           ***
! **********************************************

  DO col = 1, dim_ens
     DO row = 1, dim
        A(row, col) = A(row, col) - rowmean(row)
     END DO
  END DO
  

! ********************
! *** FINISHING UP ***
! ********************

  DEALLOCATE(rowmean)

END SUBROUTINE PDAF_etkf_Tright
