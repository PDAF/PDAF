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
!> Subtract column-wise means from an array
!!
!! This routine subtracts the mean value of each row from
!! a matrix. Usually the array is the transpose of an ensemble
!! matrix. This is used e.g. in the ETKF using a formulation
!! with T-matrix. This is the multiplication from the left: B = T A
!! where T is a symmetric dim_ens x dim_ens matrix with zero column 
!! sums defined as:  
!!            diag(T)=1-1/dim_ens; nondiag(T)=-1/dim_ens
!!
!! This is a copy of PDAF_etkf_Tleft with a clearer name.
!!
!! ! This is a core routine of PDAF and
!!   should not be changed by the user   !
!!
!! __Revision history:__
!! * 2024-11 - Lars Nerger - Initial code as copy of PDAF_etkf_Tleft
!! * Later revisions - see repository log
!!
SUBROUTINE PDAF_subtract_colmean(dim_ens, dim, A)

  USE PDAF_memcounting, &
       ONLY: PDAF_memcount

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim_ens           ! Rank of initial covariance matrix
  INTEGER, INTENT(in) :: dim               ! Number of columns in A and B
  REAL, INTENT(inout) :: A(dim_ens, dim)   ! Input/output matrix
  
! *** Local variables ***
  INTEGER :: row, col             ! counters
  REAL :: invdimens               ! Inverse of ensemble size
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

END SUBROUTINE PDAF_subtract_colmean
