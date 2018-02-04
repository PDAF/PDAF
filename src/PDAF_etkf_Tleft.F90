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
!$Id: PDAF_etkf_Tleft.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_etkf_Tleft --- Operate matrix T from left to some matrix
!
! !INTERFACE:
SUBROUTINE PDAF_etkf_Tleft(dim_ens, dim, A)

! !DESCRIPTION:
! Operate matrix T on another matrix as
!                 B = T A\\
! \\
! T is a symmetric dim_ens x dim_ens matrix with zero column 
! sums defined as:  
!            diag(T)=1-1/dim_ens; nondiag(T)=-1/dim_ens
! Applied from the left to a matrix A, T subtracts the 
! mean over the rows of A from this matrix. In this routine, T
! is not explicitly constructed, but its operation on the 
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
  INTEGER, INTENT(in) :: dim_ens           ! Rank of initial covariance matrix
  INTEGER, INTENT(in) :: dim               ! Number of columns in A and B
  REAL, INTENT(inout) :: A(dim_ens, dim)   ! Input/output matrix

! !CALLING SEQUENCE:
! Called by: User-provided prepoststep routines for ETKF
! Called by: PDAF_etkf_analysis_T
! Called by: PDAF_letkf_analysis_T
! Calls: PDAF_memcount
!EOP
  
! *** local variables ***
  INTEGER :: row, col  ! counters
  REAL :: invdimens    ! Inversize of ensemble size
  INTEGER, SAVE :: allocflag = 0  ! Flag for dynamic allocation
  REAL, ALLOCATABLE :: colmean(:) ! Mean values of columns of A


! **********************
! *** INITIALIZATION ***
! **********************

  ALLOCATE(colmean(dim))
  IF (allocflag == 0) THEN
     ! count allocated memory
     CALL PDAF_memcount(3, 'r', dim)
     allocflag = 1
  END IF
  colmean = 0.0
  invdimens = 1.0 / REAL(dim_ens)

! *** Compute column means of A ***
  DO col = 1, dim
     DO row = 1, dim_ens
        colmean(col) = colmean(col) + A(row, col)
     END DO
  END DO
  colmean = invdimens * colmean


! ****************************************************
! ***  Operate T on A                              ***
! ***                                              ***
! *** Tv_1 = (v_11-mean(v_1), ... ,v_r1-mean(v_1)) ***
! *** with v_1 = (v_11,v_21, ... ,v_N)             ***
! ****************************************************

  ! first DIM rows
  DO col = 1, dim
     DO row = 1, dim_ens
        A(row, col) = A(row, col) - colmean(col)
     END DO
  END DO


! ********************
! *** FINISHING UP ***
! ********************

  DEALLOCATE(colmean)

END SUBROUTINE PDAF_etkf_Tleft
