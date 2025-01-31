! Copyright (c) 2004-2024 Lars Nerger
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
!$Id$
!BOP
!
! !ROUTINE: PDAF_seik_resample_newT --- Perform ensemble transformation in SEIK
!
! !INTERFACE:
SUBROUTINE PDAF_seik_resample_newT(subtype, dim_p, dim_ens, rank, &
     Uinv, state_p, ens_p, type_sqrt, screen, flag)

! !DESCRIPTION:
! Routine for ensemble transformation in the SEIK filter.
! The routine transforms a forecast ensemble to represent
! the analysis state und the analysis covariance
! matrix given in factored form P = L U L$^T$.
!
! This variant is for the more efficient implementation 
! of XT. Thus ens\_p contains the real state ensemble not
! the error modes L=XT.\\
! In addition, this variant uses a block formulation for 
! the resampling, which reduces the memory requirements.
!
! Variant for PDAF with domain-decomposition.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2003-11 - Lars Nerger - Initial code
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
  USE PDAF_mod_filter, &
       ONLY: Nm1vsN, type_trans
  USE PDAF_mod_filtermpi, &
       ONLY: mype
  USE PDAF_mod_filter, &
       ONLY: debug

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: subtype      ! Filter subtype
  INTEGER, INTENT(in) :: dim_p        ! PE-local dimension of model state
  INTEGER, INTENT(in) :: dim_ens      ! Size of ensemble
  INTEGER, INTENT(in) :: rank         ! Rank of initial covariance matrix
  REAL, INTENT(inout) :: Uinv(rank, rank)      ! Inverse of matrix U
  REAL, INTENT(inout) :: state_p(dim_p)        ! PE-local model state
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) ! PE-local state ensemble
  INTEGER, INTENT(in) :: type_sqrt    ! Type of square-root of A
                                      ! (0): symmetric sqrt; (1): Cholesky decomposition
  INTEGER, INTENT(in) :: screen       ! Verbosity flag
  INTEGER, INTENT(inout) :: flag      ! Status flag

! !CALLING SEQUENCE:
! Called by: PDAF_seik_update
! Calls: PDAF_seik_Omega
! Calls: PDAF_seik_TtimesA
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: syevTYPE (LAPACK; dsyev or ssyev dependent on precision)
! Calls: potrfTYPE (LAPACK; dpotrf or spotrf dependent on precision)
! Calls: trtrsTYPE (LAPACK; dtrtrs or strtrs dependent on precision)
!EOP

! *** local variables ***
  INTEGER :: i, j, row, col           ! Counters
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  INTEGER :: lib_info                 ! Status flags for library calls
  INTEGER :: ldwork                   ! Size of work array for SYEV
  INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
  REAL :: fac                         ! Temporary variable: sqrt(dim_ens) or sqrt(rank)
  REAL :: rdim_ens                    ! Size of ensemble in real format
  LOGICAL :: storeOmega = .FALSE.     ! Store matrix Omega instead of recomputing it
  LOGICAL, SAVE :: firsttime = .TRUE. ! indicates first call to resampling
  REAL, ALLOCATABLE :: Omega(:,:)            ! Orthogonal matrix Omega
  REAL, ALLOCATABLE :: OmegaT(:,:)           ! Transpose of Omega
  REAL, SAVE, ALLOCATABLE :: OmegaTsave(:,:) ! Saved transpose of Omega
  REAL, ALLOCATABLE :: TA(:,:)               ! Temporary matrix
  REAL, ALLOCATABLE :: ens_blk(:,:)          ! Temporary blocked state ensemble
  REAL, ALLOCATABLE :: tempUinv(:,:)         ! Temporary matrix Uinv
  REAL, ALLOCATABLE :: Ttrans(:,:)           ! Temporary matrix T^T
  REAL, ALLOCATABLE :: svals(:)   ! Singular values of Uinv
  REAL, ALLOCATABLE :: work(:)    ! Work array for SYEV
  REAL, ALLOCATABLE :: Usqrt(:,:) ! Temporary for square-root of U


! **********************
! *** INITIALIZATION ***
! **********************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_seik_resample -- START'

  IF (mype == 0 .AND. screen > 0) THEN
     IF (subtype /= 3) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', 'Transform state ensemble'
     ELSE
        WRITE (*, '(a, 5x, a)') 'PDAF', 'Transform state ensemble for fixed ensemble case'
     END IF
     IF (type_sqrt == 1) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- use Cholesky square-root of U'
     ELSE
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- use symmetric square-root of U'
     END IF
  END IF

  CALL PDAF_timeit(20, 'new')
  CALL PDAF_timeit(32, 'new')


! ************************************
! *** Compute square-root of U     ***
! ************************************

  ! initialize Uinv for internal use
  ALLOCATE(tempUinv(rank, rank))
  IF (allocflag == 0) CALL PDAF_memcount(4, 'r', rank ** 2)
  IF (subtype /= 3) THEN
     tempUinv(:,:) = Uinv(:,:)
  ELSE
     rdim_ens = REAL(dim_ens)

     ! Initialize matrix T^T
     ALLOCATE(Ttrans(rank, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(4, 'r', dim_ens * rank)
     DO i = 1, rank
        DO j = 1, dim_ens
           Ttrans(i, j) = -1.0 / rdim_ens
        END DO
     END DO
     DO i = 1, rank
        Ttrans(i, i) = Ttrans(i, i) + 1.0
     END DO

     IF (Nm1vsN == 1) THEN
        ! Use factor (N-1)
        fac = dim_ens - 1
     ELSE
        ! Use factor N
        fac = dim_ens
     END IF

     CALL gemmTYPE('n', 't', rank, rank, dim_ens, &
          fac, Ttrans, rank, Ttrans, rank, &
          0.0, tempUinv, rank)
     DEALLOCATE(Ttrans)
  END IF

  ! Usqrt is allocated with dim_ens cols, because this is 
  ! required further below. Now only rank columns are used
  ALLOCATE(Usqrt(rank, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens * rank)

  typesqrtU: IF (type_sqrt == 1) THEN
     ! Compute square-root by Cholesky-decomposition

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_seik_resample:', debug, &
          '  Compute Cholesky decomposition of U^-1'

     CALL potrfTYPE('l', rank, tempUinv, rank, lib_info)

  ELSE
     ! Compute symmetric square-root by SVD of Uinv

     ALLOCATE(svals(rank))
     ALLOCATE(work(3 * rank))
     ldwork = 3 * rank
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3 * rank)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_seik_resample:', debug, &
          '  Compute eigenvalue decomposition of U^-1'

     ! Compute SVD of Uinv
     CALL syevTYPE('v', 'l', rank, Uinv, rank, svals, work, ldwork, lib_info)

     DEALLOCATE(work)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_seik_resample:', debug, '  eigenvalues', svals

     DO col = 1, rank
        DO row = 1, rank
           Usqrt(row, col) = Uinv(row, col) / SQRT(svals(col))
        END DO
     END DO

     CALL gemmTYPE('n', 't', rank, rank, rank, &
          1.0, Usqrt, rank, Uinv, rank, &
          0.0, tempUinv, rank)
     DEALLOCATE(svals)

  END IF typesqrtU

  CALL PDAF_timeit(32, 'old')
     

  ! check if computation of square-root was successful
  CholeskyOK: IF (lib_info == 0) THEN
     ! Decomposition OK, continue

  
! *************************************************
! *** Generate ensemble of interpolating states ***
! *************************************************

     ! allocate fields
     ALLOCATE(OmegaT(rank, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(4, 'r', dim_ens * rank)
     
     IF (storeOmega .AND. allocflag == 0) THEN
        ALLOCATE(OmegaTsave(rank, dim_ens))
        CALL PDAF_memcount(4, 'r', dim_ens * rank)
     END IF

     CALL PDAF_timeit(33, 'new')
     Omega_store: IF (storeOmega) THEN
        first: IF (firsttime) THEN
           ! *** At first call to SEIK_RESAMPLE initialize   ***
           ! *** the matrix Omega in SEIK_Omega and store it ***

           ALLOCATE(Omega(dim_ens, rank))
           IF (allocflag == 0) CALL PDAF_memcount(4, 'r', dim_ens * rank)

           ! *** Generate uniform orthogonal matrix OMEGA ***
           CALL PDAF_seik_Omega(rank, Omega, type_trans, screen)

           ! transpose Omega
           IF (type_sqrt == 1) THEN
              OmegaT = TRANSPOSE(Omega)
              ! store transposed Omega
              OmegaTsave = OmegaT
           ELSE
              Usqrt = TRANSPOSE(Omega)
              ! store transposed Omega
              OmegaTsave = Usqrt
           END IF

           firsttime  = .FALSE.
      
           DEALLOCATE(Omega)

        ELSE first
           IF (mype == 0 .AND. screen > 0) &
                WRITE (*, '(a, 5x, a)') 'PDAF', '--- use stored Omega'
           IF (type_sqrt == 1) THEN
              OmegaT = OmegaTsave
           ELSE
              Usqrt = OmegaTsave
           END IF
        END IF first

     ELSE Omega_store

        ! *** Initialize the matrix Omega in SEIK_Omega ***
        ! *** each time SEIK_RESAMPLE is called         ***

        ALLOCATE(Omega(dim_ens, rank))
        IF (allocflag == 0) CALL PDAF_memcount(4, 'r', dim_ens * rank)
      
        ! *** Generate uniform orthogonal matrix OMEGA ***
        CALL PDAF_seik_Omega(rank, Omega, type_trans, screen)
        
        ! transpose Omega
        IF (type_sqrt == 1) THEN
           OmegaT = TRANSPOSE(Omega)
        ELSE
           Usqrt = TRANSPOSE(Omega)
        END IF
        
        DEALLOCATE(Omega)

     END IF Omega_store
     IF (debug>0) THEN
        IF (type_sqrt == 1) THEN
           WRITE (*,*) '++ PDAF-debug PDAF_seik_update:', debug, '  Omega^T', OmegaT
        ELSE
           WRITE (*,*) '++ PDAF-debug PDAF_seik_update:', debug, '  Omega^T', Usqrt
        END IF
     END IF

     CALL PDAF_timeit(33, 'old')


     ! ***     Generate ensemble of states                             ***
     ! *** x_i = x + sqrt(FAC) L (Omega C^(-1))t                       ***
     ! *** Here FAC depends on the use definition of the covariance    ***
     ! *** matrix P using a factor (r+1)^-1 or r^-1.                   ***
    
     CALL PDAF_timeit(34, 'new')
     IF (type_sqrt == 1) THEN
        ! A = (Omega C^(-1)) by solving Ct A = OmegaT for A
        CALL trtrsTYPE('L', 'T', 'N', rank, dim_ens, &
          tempUinv, rank, OmegaT, rank, lib_info)
     ELSE
        ! TMP_UINV already contains matrix C (no more inversion)

        CALL gemmTYPE('n', 'n', rank, dim_ens, rank, &
             1.0, tempUinv, rank, Usqrt, rank, &
             0.0, OmegaT, rank)

        lib_info = 0

     END IF
     CALL PDAF_timeit(34, 'old')

     ! check if solve was successful
     solveOK: IF (lib_info == 0) THEN
        ! Solve for A OK, continue
      
        CALL PDAF_timeit(35, 'new')

        ! *** T A' (A' stored in OmegaT) ***
        ALLOCATE(TA(dim_ens, dim_ens))
        IF (allocflag == 0) CALL PDAF_memcount(4, 'r', dim_ens**2)

        CALL PDAF_seik_TtimesA(rank, dim_ens, OmegaT, TA)

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_seik_resample:', debug, '  transform', TA

        CALL PDAF_timeit(35, 'old')
        CALL PDAF_timeit(20, 'old')

        ! *** Block formulation for resampling
        maxblksize = 200
        IF (mype == 0 .AND. screen > 0) &
             WRITE (*, '(a, 5x, a, i5)') 'PDAF', '--- use blocking with size ', maxblksize
        
        ALLOCATE(ens_blk(maxblksize, dim_ens))
        IF (allocflag == 0) CALL PDAF_memcount(4, 'r', maxblksize * dim_ens)

        blocking: DO blklower = 1, dim_p, maxblksize
           
           blkupper = MIN(blklower + maxblksize - 1, dim_p)

           ! Store old state ensemble
           CALL PDAF_timeit(21, 'new')
           DO col = 1, dim_ens
              ens_blk(1 : blkupper - blklower + 1, col) &
                   = ens_p(blklower : blkupper, col)
           END DO
           
           DO col = 1,dim_ens
              ens_p(blklower : blkupper, col) = state_p(blklower : blkupper)
           END DO
           CALL PDAF_timeit(21, 'old')

           ! *** X = state+ sqrt(FAC) state_ens T A^T (A^T stored in OmegaT) ***
           CALL PDAF_timeit(22, 'new')

           IF (Nm1vsN == 1) THEN
              ! Use factor (N-1)^-1
              fac = SQRT(REAL(dim_ens - 1))
           ELSE
              ! Use factor N^-1
              fac = SQRT(REAL(dim_ens))
           END IF

           CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
                fac, ens_blk(1, 1), maxblksize, TA(1, 1), dim_ens, &
                1.0, ens_p(blklower, 1), dim_p)
           CALL PDAF_timeit(22, 'old')
        END DO blocking

        DEALLOCATE(ens_blk, TA)

     ELSE SolveOK

        ! Solve for A failed
        WRITE (*, '(/5x, a/)') 'PDAF-ERROR(2): Problem with solve for A in SEIK_RESAMPLE !!!'
        flag = 2

        CALL PDAF_timeit(20, 'old')
     ENDIF SolveOK

     DEALLOCATE(OmegaT)

  ELSE CholeskyOK

     ! eigendecomposition failed
     WRITE (*, '(/5x, a/)') 'PDAF-ERROR(3): Problem with Cholesky decomposition of Uinv !!!'
     flag = 3

     CALL PDAF_timeit(20, 'old')

  ENDIF CholeskyOK


! ****************
! *** clean up ***
! ****************

  DEALLOCATE(Usqrt)
  DEALLOCATE(tempUinv)

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_seik_resample -- END'

END SUBROUTINE PDAF_seik_resample_newT

