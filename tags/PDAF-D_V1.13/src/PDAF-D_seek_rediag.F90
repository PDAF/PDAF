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
!$Id: PDAF-D_seek_rediag.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_seek_rediag --- Perform rediagonalization of P in SEEK
!
! !INTERFACE:
SUBROUTINE PDAF_seek_rediag(dim_p, dim_eof, Uinv, eofV_p, subtype, &
     screen, flag)

! !DESCRIPTION:
! Re-orthogonalization of the modes V of the
! low-rank approximated covariance matrix in
! its decomposed form P = V U V$^T$.
!
! Compute eigenmodes of the matrix B = L$^T$ L = C D C$^T$
! where L = V U$^{1/2}$ (from Cholesky decomposition)
! and get new modes V as V = L C D$^{-1/2}$,
! $D = diag(\lambda_1,...,\lambda_r),\ \lambda_i > \lambda_i+1$
! and C matrix of corresponding eigenvectors.
! The new U is given by the matrix D.
!
! Variant for domain decomposed states.
!
! New version to compute matrix B. More efficient for
! dim $>>$ dim\_eof
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2003-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: mype, MPIerr, COMM_filter, MPI_SUM, MPI_REALTYPE

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p    ! PE-Local state dimension
  INTEGER, INTENT(in) :: dim_eof  ! Number of EOFs
  REAL, INTENT(inout) :: Uinv(dim_eof,dim_eof) ! Inverse of matrix U
  REAL, INTENT(inout) :: eofV_p(dim_p,dim_eof) ! PE-local matrix V
  INTEGER, INTENT(in) :: subtype  ! Filter subtype
  INTEGER, INTENT(in) :: screen   ! Verbosity flag
  INTEGER, INTENT(inout) :: flag  ! Status Flag

! !CALLING SEQUENCE:
! Called by: PDAF_seek_update
! Calls: PDAF_timeit
! Calls: PDAF_memcount
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: potrfTYPE (LAPACK; dpotrf ot spotrf dependent on precision)
! Calls: gesvTYPE (LAPACK; dgesv or sgesv dependent on precision)
! Calls: MPI_allreduce (MPI)
!EOP
       
! *** local variables ***
  INTEGER :: row, col              ! counters
  INTEGER, SAVE :: allocflag = 0   ! Flag whether first time allocation is done
  INTEGER, ALLOCATABLE :: ipiv(:)  ! vector of pivot indices for SGESV
  REAL, ALLOCATABLE :: ev(:)       ! vector of eigenvalues
  REAL, ALLOCATABLE :: rwork(:)    ! workarray for eigenproblem
  REAL, ALLOCATABLE :: U(:,:)      ! temporary for covar matrix
  REAL, ALLOCATABLE :: L(:,:)      ! covariance matrix
  REAL, ALLOCATABLE :: Temp1(:,:)  ! temporary for covar matrix
  REAL, ALLOCATABLE :: B(:,:)      ! temporary for covar matrix
  INTEGER :: syev_info, gesv_info, potrf_info  ! Info flags for LAPACK calls


! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a)') 'PDAF', 'Re-orthogonalize covariance matrix modes'
  END IF


! **************************************
! *** Compute matrix B = A^T V^T V A ***
! **************************************

  CALL PDAF_timeit(20, 'new')
  ! *** Get U by inversion of Uinv ***
  ALLOCATE(U(dim_eof, dim_eof))
  ALLOCATE(ipiv(dim_eof))
  IF (allocflag == 0) CALL PDAF_memcount(4, 'r', dim_eof * dim_eof)
  IF (allocflag == 0) CALL PDAF_memcount(4, 'i', dim_eof)

  ! Initialize U
  U = 0.0
  DO row = 1, dim_eof
     U(row, row) = 1.0
  END DO

  ! call solver
  CALL PDAF_timeit(32, 'new')
  CALL gesvTYPE(dim_eof, dim_eof, Uinv, dim_eof, ipiv, &
       U, dim_eof, gesv_info)
  CALL PDAF_timeit(32, 'old')

  DEALLOCATE(ipiv)

  ! Check if inversion was successful
  INVok: IF (gesv_info == 0) THEN

     ! *** Cholesky decomposition of U: U = A A^T ***
     CALL PDAF_timeit(33, 'new')
     CALL potrfTYPE('l', dim_eof, U, dim_eof, potrf_info)
     CALL PDAF_timeit(33, 'old')

     ! *** set upper elements to zero ***
     DO col = 2, dim_eof
        DO row = 1, col - 1
           U(row, col) = 0.0
        END DO
     END DO

     !*** Compute B = A^T V^T V A ***
     ALLOCATE(Temp1(dim_eof, dim_eof))
     ALLOCATE(B(dim_eof, dim_eof))
     IF (allocflag == 0) CALL PDAF_memcount(4, 'r', 2 * dim_eof * dim_eof)

     ! local V^T V
     CALL gemmTYPE('t', 'n', dim_eof, dim_eof, dim_p, &
          1.0, eofV_p, dim_p, eofV_p, dim_p, &
          0.0, Temp1, dim_eof)
     CALL PDAF_timeit(20, 'old')

     CALL PDAF_timeit(21, 'new')
     CALL MPI_allreduce(Temp1, B, dim_eof * dim_eof, MPI_REALTYPE, &
          MPI_SUM, COMM_filter, MPIerr)
     CALL PDAF_timeit(21, 'old')

     CALL PDAF_timeit(20, 'new')
     ! (V^T V) A (A stored in U)
     CALL gemmTYPE('n', 'n', dim_eof, dim_eof, dim_eof, &
          1.0, B, dim_eof, U, dim_eof, &
          0.0, Temp1, dim_eof)

     ! B = A^T (V^T V A) (A stored in U)
     CALL gemmTYPE('t', 'n', dim_eof, dim_eof, dim_eof, &
          1.0, U, dim_eof, Temp1, dim_eof, &
          0.0, B, dim_eof)
     CALL PDAF_timeit(31, 'old')


! *******************************
! *** Eigendecomposition of B ***
! ***       B = C D C^T       ***
! *******************************

     ALLOCATE(ev(dim_eof))
     ALLOCATE(rwork(3 * dim_eof))
     IF (allocflag == 0) CALL PDAF_memcount(4, 'r', 4 * dim_eof)
  
     CALL syevTYPE('v', 'u', dim_eof, B, dim_eof, &
          ev, rwork, 3 * dim_eof, syev_info)
     CALL PDAF_timeit(20, 'old')
     DEALLOCATE(rwork)

     ! check if eigendecomposition was successful
     EVPok: IF (syev_info == 0) THEN
        ! Eigendecomposition OK, continue
      
      ! *** Reorder matrix of eigenvectors ***
!       eof_mid = floor(real(dim_eof)/2.0)

!       do col=1,eof_mid
!         rwork(1:dim_eof) = U(:,col)
!         U(:,col) = U(:,dim_eof-col+1)
!         U(:,dim_eof-col+1) = rwork(1:dim_eof)
!       end do


! *****************************************************
! *** Initialize covar with re-orthogonalized modes ***
! *****************************************************

        CALL PDAF_timeit(22, 'new')
        ! AC = A C (AC stored in Temp1, A stored in U, C stored in B) 
        CALL gemmTYPE('n', 'n', dim_eof, dim_eof, dim_eof, &
             1.0, U, dim_eof, B, dim_eof, &
             0.0, Temp1, dim_eof)

        ! initialize L from V
        ALLOCATE(L(dim_p, dim_eof))
        IF (allocflag == 0) CALL PDAF_memcount(4, 'r', dim_eof * dim_p)

        L = eofV_p

        ! V = L AC (AC stored in Temp1)
        CALL gemmTYPE('n', 'n', dim_p, dim_eof, dim_eof, &
             1.0, L, dim_p, Temp1, dim_eof, &
             0.0, eofV_p, dim_p)
        DEALLOCATE(L, U, Temp1, B)

        unitmodes: IF (subtype /= 1) THEN
           ! use eigenvalues in U with unit modes in V
           !    U = diag(lambda_1,...,lambda_r)

           IF (mype == 0 .AND. screen > 0) THEN
              WRITE (*, '(a, 5x, a)') 'PDAF', '--- Use normalized modes'
           END IF

           Uinv = 0.0
           DO row = 1, dim_eof
              Uinv(row, row) = 1.0 / ev(row)
              ev(row) = SQRT(1.0 / ev(row))
           END DO

           ! Rescale modes V
           DO col = 1, dim_eof
              DO row = 1, dim_p
                 eofV_p(row, col) = eofV_p(row, col) * ev(col)
              END DO
           END DO
        ELSE unitmodes
           ! use unit U with modes in V scaled by eigenvalues

           IF (mype==0 .AND. screen > 0) THEN
              WRITE (*, '(a, 5x, a)') 'PDAF', '--- Use unit U'
           END IF

           Uinv = 0.0
           DO row = 1, dim_eof
              Uinv(row, row) = 1.0
           END DO
        END IF unitmodes
        CALL PDAF_timeit(22, 'old')
      
     ELSE
        ! eigendecomposition failed
        IF (mype == 0) WRITE (*, '(/5x, a/)') &
             'PDAF-ERROR(1): Problem with EOF decomposition of matrix Lt L !!!'
        flag = 1
        
     ENDIF EVPok

     DEALLOCATE(ev)

  ELSE INVok
     ! Inversion failed
     IF (mype == 0) WRITE (*, '(/5x, a/)') &
          'PDAF-ERROR(2): Problem with inversion of Uinv !!!'
     flag = 2

  ENDIF INVok


! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1
  
END SUBROUTINE PDAF_seek_rediag
