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
! !ROUTINE: PDAF_lseik_resample --- Perform LSEIK ensemble transformation
!
! !INTERFACE:
SUBROUTINE PDAF_lseik_resample(domain_p, subtype, dim_l, dim_ens, &
     rank, Uinv_l, state_l, ens_l, OmegaT_in, type_sqrt, screen, flag)

! !DESCRIPTION:
! Routine for ensemble transformation in the 
! LSEIK filter. The routine generates a local
! ensemble of states that represents the local
! analysis state und the local analysis 
! covariance matrix given in factored 
! form P = L U L$^T$.
!
! Variant for domain decomposition. This variant is 
! also using the more efficient implementation of XT. 
! Thus ens\_l contains the real state ensemble not
! the error modes L=XT.
! In addition this variant uses a block formulation for 
! the resampling, which reduces the memory requirements.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
 !REVISION HISTORY:
! 2005-09 - Lars Nerger - Initial code
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
       ONLY: Nm1vsN, debug
  USE PDAF_mod_filtermpi, &
       ONLY: mype
#if defined (_OPENMP)
  USE omp_lib, &
       ONLY: omp_get_num_threads, omp_get_thread_num
#endif

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain_p  ! Current local analysis domain
  INTEGER, INTENT(in) :: subtype   ! Specification of filter subtype
  INTEGER, INTENT(in) :: dim_l     ! State dimension on local analysis domain
  INTEGER, INTENT(in) :: dim_ens   ! Size of ensemble
  INTEGER, INTENT(in) :: rank      ! Rank of initial covariance matrix
  REAL, INTENT(inout) :: Uinv_l(rank, rank)       ! Inverse of matrix U
  REAL, INTENT(inout) :: state_l(dim_l)           ! Local model state
  REAL, INTENT(inout) :: ens_l(dim_l, dim_ens)    ! Local state ensemble
  REAL, INTENT(inout) :: OmegaT_in(rank, dim_ens) ! Matrix Omega
  INTEGER, INTENT(in) :: type_sqrt ! Type of square-root of A
                                   ! (0): symmetric sqrt; (1): Cholesky decomposition
  INTEGER, INTENT(in) :: screen    ! Verbosity flag
  INTEGER, INTENT(inout) :: flag   ! Status flag

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update
! Calls: PDAF_seik_Omega
! Calls: PDAF_seik_TtimesA
! Calls: PDAF_timeit
! Calls: PDAF_memcount
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
  REAL :: fac                         ! Temporary variable sqrt(dim_ens) or sqrt(rank)
  REAL    :: rdim_ens                 ! Inverse ensemble size as real
  INTEGER, SAVE :: lastdomain = -1    ! store domain index
  LOGICAL, SAVE :: screenout = .true. ! Whether to print information to stdout
  REAL, ALLOCATABLE :: OmegaT(:,:)    ! Transpose of Omega
  REAL, ALLOCATABLE :: TA(:,:)        ! Temporary matrix
  REAL, ALLOCATABLE :: ens_block(:,:) ! Temporary blocked state ensemble
  REAL, ALLOCATABLE :: tmpUinv_l(:,:) ! Temporary matrix Uinv
  REAL, ALLOCATABLE :: Ttrans(:,:)    ! Temporary matrix T^T
  REAL, ALLOCATABLE :: svals(:)       ! Singular values of Uinv
  REAL, ALLOCATABLE :: work(:)        ! Work array for SYEV
  INTEGER, SAVE :: mythread, nthreads ! Thread variables for OpenMP

!$OMP THREADPRIVATE(mythread, nthreads, lastdomain, allocflag, screenout)


! *******************
! *** Preparation ***
! *******************

  CALL PDAF_timeit(51, 'new')

#if defined (_OPENMP)
  nthreads = omp_get_num_threads()
  mythread = omp_get_thread_num()
#else
  nthreads = 1
  mythread = 0
#endif

  ! Control screen output
  IF (lastdomain<domain_p .AND. lastdomain>-1) THEN
     screenout = .false.
  ELSE
     screenout = .true.

     ! In case of OpenMP, let only thread 0 write output to the screen
     IF (mythread>0) screenout = .false.
  END IF

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lseik_resample -- START'


! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype == 0 .AND. screen > 0 .AND. screenout) THEN
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

  ! allocate fields
  ALLOCATE(OmegaT(rank, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(4, 'r', rank * dim_ens)


! ************************************
! *** Compute square-root of U     ***
! ************************************

  ! initialize Uinv for internal use
  ALLOCATE(tmpUinv_l(rank, rank))
  IF (allocflag == 0) CALL PDAF_memcount(4, 'r', rank**2)
  IF (subtype /= 3) THEN
     tmpUinv_l(:, :) = Uinv_l(:, :)
  ELSE
     rdim_ens = REAL(dim_ens)

     ! Initialize matrix T^T
     ALLOCATE(Ttrans(rank, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(4, 'r', rank * dim_ens)
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
          0.0, tmpUinv_l, rank)
     DEALLOCATE(Ttrans)
  END IF

  typesqrtU: IF (type_sqrt == 1) THEN
     ! Compute square-root by Cholesky-decomposition

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lseik_resample:', debug, &
          '  Compute Cholesky decomposition of U^-1_l'

     CALL potrfTYPE('l', rank, tmpUinv_l, rank, lib_info)

  ELSE
     ! Compute symmetric square-root by SVD of Uinv

     ALLOCATE(svals(rank))
     ALLOCATE(work(3 * rank))
     ldwork = 3 * rank
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 4 * rank)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lseik_resample:', debug, &
          '  Compute eigenvalue decomposition of U^-1_l'

     ! Compute SVD of Uinv
     CALL syevTYPE('v', 'l', rank, Uinv_l, rank, svals, work, ldwork, lib_info)

     DEALLOCATE(work)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lseik_resample:', debug, '  eigenvalues', svals

     ! Use OmegaT as temporary array
     DO col = 1, rank
        DO row = 1, rank
           OmegaT(row, col) = Uinv_l(row, col) / SQRT(svals(col))
        END DO
     END DO

     CALL gemmTYPE('n', 't', rank, rank, rank, &
          1.0, OmegaT, rank, Uinv_l, rank, &
          0.0, tmpUinv_l, rank)

     DEALLOCATE(svals)

  END IF typesqrtU

  CALL PDAF_timeit(32, 'old')


  ! check if Cholesky decomposition was successful
  CholeskyOK: IF (lib_info == 0) THEN
     ! Decomposition OK, continue

  
! *************************************************
! *** Generate ensemble of interpolating states ***
! *************************************************

     ! ***     Generate ensemble of states                             ***
     ! *** x_i = x + sqrt(FAC) X T (Omega C^(-1))t                     ***
     ! *** Here FAC depends on the use definition of the covariance    ***
     ! *** matrix P using a factor (r+1)^-1 or r^-1.                   ***
    
     CALL PDAF_timeit(34, 'new')
     IF (type_sqrt == 1) THEN
        ! Initialize the matrix Omega from argument OmegaT_in
        OmegaT = OmegaT_in

        ! A = (Omega C^(-1)) by solving Ct A = OmegaT for A
        CALL trtrsTYPE('L', 'T', 'N', rank, dim_ens, &
             tmpUinv_l, rank, OmegaT, rank, lib_info)
     ELSE
        ! TMP_UINV already contains matrix C (no more inversion)

        CALL gemmTYPE('n', 'n', rank, dim_ens, rank, &
             1.0, tmpUinv_l, rank, OmegaT_in, rank, &
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
             WRITE (*,*) '++ PDAF-debug PDAF_lseik_resample:', debug, '  transform', TA

        CALL PDAF_timeit(35, 'old')
        CALL PDAF_timeit(20, 'old')

        ! *** Block formulation for resampling
        maxblksize = 200
        IF (mype == 0 .AND. screen > 0 .AND. screenout) &
             WRITE (*, '(a, 5x, a, i5)') 'PDAF', '--- use blocking with size ', maxblksize

        ALLOCATE(ens_block(maxblksize, dim_ens))
        IF (allocflag == 0) CALL PDAF_memcount(4, 'r', maxblksize * dim_ens)

        blocking: DO blklower = 1, dim_l, maxblksize

           blkupper = MIN(blklower + maxblksize - 1, dim_l)

           ! Store old state ensemble
           CALL PDAF_timeit(21, 'new')
           DO col = 1, dim_ens
              ens_block(1 : blkupper - blklower + 1, col) &
                   = ens_l(blklower : blkupper, col)
           END DO

           DO col = 1, dim_ens
              ens_l(blklower : blkupper, col) = state_l(blklower : blkupper)
           END DO
           CALL PDAF_timeit(21, 'old')

           ! *** X = state+ sqrt(FAC) state_ens T A^T (A^T stored in OmegaT) ***
           ! *** Here FAC depends on the use definition of the covariance    ***
           ! *** matrix P using a factor (r+1)^-1 or r^-1.                   ***
           CALL PDAF_timeit(22, 'new')

           IF (Nm1vsN == 1) THEN
              ! Use factor (N-1)^-1
              fac = SQRT(REAL(rank))
           ELSE
              ! Use factor N^-1
              fac = SQRT(REAL(dim_ens))
           END IF

           CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
                fac, ens_block(1, 1), maxblksize, TA(1, 1), dim_ens, &
                1.0, ens_l(blklower, 1), dim_l)
           CALL PDAF_timeit(22, 'old')
        END DO blocking
        
        DEALLOCATE(ens_block, TA)

     ELSE SolveOK

        ! Solve for A failed
        WRITE (*, '(/5x, a, i10, a/)') &
             'PDAF-ERROR(2): Problem with solve for A in SEIK_RESAMPLE - domain ', &
             domain_p, ' !!!'
        flag = 2

        CALL PDAF_timeit(20, 'old')

     ENDIF SolveOK

  ELSE CholeskyOK

     ! eigendecomposition failed
     IF (type_sqrt == 1) THEN
        WRITE (*, '(/5x, a, i10, a/)') &
             'PDAF-ERROR(1): Problem with Cholesky decomposition of Uinv - domain ', &
             domain_p, ' !!!'
     ELSE
        WRITE (*, '(/5x, a, i10, a/)') &
             'PDAF-ERROR(1): Problem with eigenvalue decomposition of Uinv - domain ', &
             domain_p, ' !!!'
     END IF
     flag = 1

  ENDIF CholeskyOK
  CALL PDAF_timeit(51, 'old')


! ****************
! *** clean up ***
! ****************

  DEALLOCATE(tmpUinv_l)
  DEALLOCATE(OmegaT)

  IF (allocflag == 0) allocflag = 1

  ! Store domain index
  lastdomain = domain_p

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lseik_resample -- END'

END SUBROUTINE PDAF_lseik_resample

