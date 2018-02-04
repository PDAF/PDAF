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
!$Id: PDAF_seik_TtimesA.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_seik_TtimesA() --- Operate matrix T on some matrix
!
! !INTERFACE:
SUBROUTINE PDAF_seik_TtimesA(rank, dim_col, A, B)

! !DESCRIPTION:
! Operate matrix T on another matrix as
!                 B = T A\\
! \\
! T is a dim_ens x (dim_ens-1) matrix with zero column sums.
! There are two proposed forms of T (ensemble size N):\\
! typeT=0: diag(T)=1-1/N; nondiag(T)=-1/N; 
!          last row= -1/N\\
! typeT=1: diag(T)=1; nondiag(T)=0; last row = -1\\
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
  INTEGER, INTENT(in) :: rank               ! Rank of initial covariance matrix
  INTEGER, INTENT(in) :: dim_col            ! Number of columns in A and B
  REAL, INTENT(in)    :: A(rank, dim_col)   ! Input matrix
  REAL, INTENT(out)   :: B(rank+1, dim_col) ! Output matrix (TA)

! !CALLING SEQUENCE:
! Called by: User-provided prepoststep routines for SEIK and LSEIK
! Called by: PDAF_seik_analysis_newT
! Called by: PDAF_seik_resample_newT
! Called by: PDAF_lseik_analysis
! Called by: PDAF_lseik_resample
! Calls: PDAF_memcount
!EOP
  
! *** local variables ***
  INTEGER :: row, col  ! counters
  INTEGER :: typeT = 0 ! Which type of T
  REAL :: invdimens    ! Inversize of ensemble size
  INTEGER, SAVE :: allocflag = 0  ! Flag for dynamic allocation
  REAL, ALLOCATABLE :: colmean(:) ! Mean values of columns of A

!$OMP THREADPRIVATE(allocflag)


! **********************
! *** INITIALIZATION ***
! **********************

  ALLOCATE(colmean(dim_col))
  IF (allocflag == 0) THEN
     ! count allocated memory
     CALL PDAF_memcount(3, 'r', dim_col)
     allocflag = 1
  END IF
  colmean = 0.0
  invdimens = -1.0 / REAL(rank + 1)

  whichT: IF (typeT == 0) THEN

! *** Compute column means of A ***
    DO col = 1, dim_col
       DO row = 1, rank
          colmean(col) = colmean(col) + invdimens * A(row, col)
       END DO
    END DO

  END IF whichT


! ****************************************************
! ***  Operate T on A                              ***
! ***                                              ***
! *** Tv_1 = (v_11-mean(v_1), ... ,v_r1-mean(v_1)) ***
! *** with v_1 = (v_11,v_21, ... ,v_N1  )          ***
! ****************************************************

  ! first DIM rows
  DO col = 1, dim_col
     DO row = 1, rank
        B(row, col) = A(row, col) + colmean(col)
     END DO
  END DO

  ! row RANK+1
  DO col = 1, dim_col
     B(rank + 1, col) = colmean(col)
  END DO


! ********************
! *** FINISHING UP ***
! ********************

  DEALLOCATE(colmean)

END SUBROUTINE PDAF_seik_TtimesA
