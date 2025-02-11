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
!
!> Perform LSEIK analysis step
!!
!! Analysis step of the LSEIK filter
!! with adaptive forgetting factor.
!!
!! Variant for domain decomposed states.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2005-09 - Lars Nerger - Initial code
!! * Later revisions - see svn log
!!
MODULE PDAF_lseik_analysis

CONTAINS
SUBROUTINE PDAFlseik_analysis(domain_p, step, dim_l, dim_obs_l, dim_ens, &
     rank, state_l, Uinv_l, ens_l, HL_l, HXbar_l, &
     obs_l, state_inc_l, forget, &
     U_prodRinvA_l, incremental, screen, debug, flag)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: mype
  USE PDAF_analysis_utils, &
       ONLY: PDAF_seik_matrixT, PDAF_seik_TtimesA, PDAF_seik_Uinv
#if defined (_OPENMP)
  USE omp_lib, &
       ONLY: omp_get_num_threads, omp_get_thread_num
#endif

  IMPLICIT NONE

! *** Arguments ***
!  Variable naming scheme:
!    suffix _p: Denotes a full variable on the PE-local domain
!    suffix _l: Denotes a local variable on the current analysis domain
  INTEGER, INTENT(in) :: domain_p    !< Current local analysis domain
  INTEGER, INTENT(in) :: step        !< Current time step
  INTEGER, INTENT(in) :: dim_l       !< State dimension on local analysis domain
  INTEGER, INTENT(in) :: dim_obs_l   !< Size of obs. vector on local ana. domain
  INTEGER, INTENT(in) :: dim_ens     !< Size of ensemble 
  INTEGER, INTENT(in) :: rank        !< Rank of initial covariance matrix
  REAL, INTENT(inout) :: state_l(dim_l)           !< State on local analysis domain
  REAL, INTENT(inout) :: Uinv_l(rank, rank)       !< Inverse of matrix U
  REAL, INTENT(in) :: ens_l(dim_l, dim_ens)       !< Local state ensemble
  REAL, INTENT(inout) :: HL_l(dim_obs_l, dim_ens) !< Local observed state ensemble (perturbation)
  REAL, INTENT(in) :: HXbar_l(dim_obs_l)          !< Local observed ensemble mean
  REAL, INTENT(in) :: obs_l(dim_obs_l)            !< Local observation vector
  REAL, INTENT(in) :: state_inc_l(dim_l)          !< Local state increment
  REAL, INTENT(inout) :: forget      !< Forgetting factor
  INTEGER, INTENT(in) :: incremental !< Control incremental updating
  INTEGER, INTENT(in) :: screen      !< Verbosity flag
  INTEGER, INTENT(in) :: debug       !< Flag for writing debug output
  INTEGER, INTENT(inout) :: flag     !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA_l          !< Provide product R^-1 A for local analysis domain

       
! *** local variables ***
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  INTEGER, SAVE :: lastdomain = -1     ! store domain index
  LOGICAL, SAVE :: screenout = .true.  ! Whether to print information to stdout
  REAL, ALLOCATABLE :: RiHL_l(:,:)     ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: Uinv_inc(:,:)   ! local Uinv
  REAL, ALLOCATABLE :: innov_l(:)      ! observation innovation
  REAL, ALLOCATABLE :: RiHLd_l(:)      ! local RiHLd
  REAL, ALLOCATABLE :: TRiHLd_l(:,:)   ! Temporary vector for analysis 
  REAL, ALLOCATABLE :: Uinv_l_tmp(:,:) ! Temporary storage of Uinv
  INTEGER, ALLOCATABLE :: ipiv(:)      ! vector of pivot indices for GESV
  INTEGER :: gesv_info                 ! control flag for GESV
  INTEGER, SAVE :: mythread, nthreads  ! Thread variables for OpenMP
  INTEGER :: screen_dummy              ! Dummy variable to avoid compiler warning

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

  ! Initialize variable to prevent compiler warning
  screen_dummy = screen

  ! Control screen output
  IF (lastdomain<domain_p .AND. lastdomain>-1) THEN
     screenout = .false.
  ELSE
     screenout = .true.

     ! In case of OpenMP, let only thread 0 write output to the screen
     IF (mythread>0) screenout = .false.

     ! Output, only in case of OpenMP parallelization
#if defined (_OPENMP)
     IF (mype == 0 .AND. screen > 0 .AND. screenout) THEN
        WRITE (*,'(a, 5x,a,i5,a)') &
             'PDAF', '--- Use OpenMP parallelization with ', nthreads, ' threads'
     END IF
#endif
  END IF

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lseik_analysis -- START'

  CALL PDAF_timeit(51, 'old')


! **************************
! *** Compute innovation ***
! ***   d = y - H x      ***
! **************************

  CALL PDAF_timeit(16, 'new')
  CALL PDAF_timeit(20, 'new')
  CALL PDAF_timeit(51, 'new')

  ALLOCATE(innov_l(dim_obs_l))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l)

  innov_l = obs_l - HXbar_l

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lseik_analysis:', debug, '  innovation d_l', innov_l

  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(20, 'old')


! *************************************************
! ***   Compute analyzed matrix Uinv            ***
! ***                                           ***
! ***  -1              T        T  -1           ***
! *** U  = forget*fac T T + (HL)  R  (HL)       ***
! ***                                           ***
! *** Here FAC is a scaling factor according    ***
! *** to the definition of the ensemble         ***
! *** covariance scaled by N^-1 (original SEIK) ***
! *** (N-1)^-1 (SEIK as ensemble KF)            ***
! *************************************************

  CALL PDAF_timeit(21, 'new')
  CALL PDAF_timeit(51, 'new')

  ! HL = [Hx_1 ... Hx_(r+1)] T
  CALL PDAF_seik_matrixT(dim_obs_l, dim_ens, HL_l)

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lseik_analysis:', debug, '  HXT_l', HL_l(:, 1:dim_ens-1)

  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(21, 'old')


  ! ***                RiHL = Rinv HL                 ***
  ! *** this is implemented as a subroutine thus that ***
  ! *** Rinv does not need to be allocated explicitly ***
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lseik_analysis -- call prodRinvA_l'

  ALLOCATE(RiHL_l(dim_obs_l, rank))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l * rank)

  CALL PDAF_timeit(48, 'new')
  CALL U_prodRinvA_l(domain_p, step, dim_obs_l, rank, obs_l, HL_l, RiHL_l)
  CALL PDAF_timeit(48, 'old')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lseik_analysis:', debug, '  R^-1(HXT_l)', RiHL_l
 
  CALL PDAF_timeit(21, 'new')
  CALL PDAF_timeit(51, 'new')

  ! *** Initialize Uinv = (r+1) T^T T ***
  CALL PDAF_seik_Uinv(rank, Uinv_l)

  ! *** Finish computation of Uinv  ***
  ! ***   -1          -1    T       ***
  ! ***  U  = forget U  + HL RiHL   ***

  ALLOCATE(Uinv_inc(rank, rank))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank**2)
  CALL gemmTYPE('t', 'n', rank, rank, dim_obs_l, &
       1.0, HL_l, dim_obs_l, RiHL_l, dim_obs_l, &
       0.0, Uinv_inc, rank)

  Uinv_l = forget * Uinv_l + Uinv_inc

  DEALLOCATE(Uinv_inc)

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lseik_analysis:', debug, '  U^-1_l', Uinv_l

  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(21, 'old')
  CALL PDAF_timeit(16, 'old')


! ************************************
! ***      update model state      ***
! ***                              ***
! ***  a   f          T         f  ***
! *** x = x + L U RiHV  (y - H x ) ***
! ***                              ***
! ************************************

  CALL PDAF_timeit(22, 'new')
  CALL PDAF_timeit(51, 'new')

  ! ************************
  ! *** RiHLd = RiHV^T d ***
  ! ************************
  ALLOCATE(RiHLd_l(rank))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * rank)

  CALL gemvTYPE('t', dim_obs_l, rank, 1.0, RiHL_l, &
       dim_obs_l, innov_l, 1, 0.0, RiHLd_l, 1)

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lseik_analysis:', debug, '  (HXT_l R^-1)^T d_l', RiHLd_l

  DEALLOCATE(RiHL_l, innov_l)


  ! ****************************************
  ! *** Compute  w = U RiHLd  by solving ***
  ! ***           -1                     ***
  ! ***          U  w = RiHLd            ***
  ! *** for w. We use the LAPACK         ***
  ! *** routine GESV.                    ***
  ! ****************************************

  ALLOCATE(Uinv_l_tmp(rank, rank))
  ALLOCATE(ipiv(rank))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank**2)
  IF (allocflag == 0) CALL PDAF_memcount(3, 'i', rank)

  ! save matrix Uinv
  Uinv_l_tmp = Uinv_l

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lseik_analysis:', debug, &
       '  Invert U^-1_l using solver GESV'

  ! call solver (GESV - LU solver)
  CALL gesvTYPE(rank, 1, Uinv_l_tmp, rank, ipiv, &
       RiHLd_l, rank, gesv_info)
  DEALLOCATE(Uinv_l_tmp, ipiv)

  ! *** check if solve was successful
  update: IF (gesv_info /= 0) THEN
     WRITE (*, '(/5x, a, i10, 39a/)') &
          'PDAF-ERROR(1): Domain ', domain_p, 'Problem in solve for state analysis !!!'
     flag = 1

     CALL PDAF_timeit(22, 'old')
  ELSE

     ! **************************
     ! *** Compute vector T w ***
     ! **************************

     ALLOCATE(TRiHLd_l(dim_ens, 1))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)

     CALL PDAF_seik_TtimesA(rank, 1, RiHLd_l, TRiHLd_l)
     DEALLOCATE(RiHLd_l)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lseik_analysis:', debug, '  wbar_l', TRiHLd_l

     CALL PDAF_timeit(22, 'old')


     ! **************************
     ! *** Update model state ***
     ! ***    a   f           ***
     ! ***   x = x + L RiHLd  ***
     ! **************************

     CALL PDAF_timeit(23, 'new')

     CALL gemvTYPE('n', dim_l, dim_ens, 1.0, ens_l, &
          dim_l, TRiHLd_l, 1, 0.0, state_inc_l, 1)
     DEALLOCATE(TRiHLd_l)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lseik_analysis:', debug, '  state_inc_l', state_inc_l

     IF (incremental == 0) THEN
        ! update state only if incremental updating is not used
        state_l = state_l + state_inc_l
     END IF

     CALL PDAF_timeit(23, 'old')
    
  END IF update

  CALL PDAF_timeit(51, 'old')


! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1

  ! Store domain index
  lastdomain = domain_p

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lseik_analysis -- END'

END SUBROUTINE PDAFlseik_analysis

!-------------------------------------------------------------------------------
!> Perform LSEIK ensemble transformation
!!
!! Routine for ensemble transformation in the 
!! LSEIK filter. The routine generates a local
!! ensemble of states that represents the local
!! analysis state und the local analysis 
!! covariance matrix given in factored 
!! form P = L U L$^T$.
!!
!! Variant for domain decomposition. This variant is 
!! also using the more efficient implementation of XT. 
!! Thus ens\_l contains the real state ensemble not
!! the error modes L=XT.
!! In addition this variant uses a block formulation for 
!! the resampling, which reduces the memory requirements.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2005-09 - Lars Nerger - Initial code
!! * Later revisions - see svn log
!!
SUBROUTINE PDAFlseik_resample(domain_p, subtype, dim_l, dim_ens, &
     rank, Uinv_l, state_l, ens_l, OmegaT_in, type_sqrt, screen, flag)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_lseik, &
       ONLY: Nm1vsN, debug
  USE PDAF_mod_filtermpi, &
       ONLY: mype
  USE PDAF_analysis_utils, &
       ONLY: PDAF_seik_TtimesA
#if defined (_OPENMP)
  USE omp_lib, &
       ONLY: omp_get_num_threads, omp_get_thread_num
#endif

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p  !< Current local analysis domain
  INTEGER, INTENT(in) :: subtype   !< Specification of filter subtype
  INTEGER, INTENT(in) :: dim_l     !< State dimension on local analysis domain
  INTEGER, INTENT(in) :: dim_ens   !< Size of ensemble
  INTEGER, INTENT(in) :: rank      !< Rank of initial covariance matrix
  REAL, INTENT(inout) :: Uinv_l(rank, rank)       !< Inverse of matrix U
  REAL, INTENT(inout) :: state_l(dim_l)           !< Local model state
  REAL, INTENT(inout) :: ens_l(dim_l, dim_ens)    !< Local state ensemble
  REAL, INTENT(inout) :: OmegaT_in(rank, dim_ens) !< Matrix Omega
  INTEGER, INTENT(in) :: type_sqrt !< Type of square-root of A
                                   !< (0): symmetric sqrt; (1): Cholesky decomposition
  INTEGER, INTENT(in) :: screen    !< Verbosity flag
  INTEGER, INTENT(inout) :: flag   !< Status flag

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

  CALL PDAF_timeit(24, 'new')
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

        ! *** T A' (A' stored in OmegaT) ***
        ALLOCATE(TA(dim_ens, dim_ens))
        IF (allocflag == 0) CALL PDAF_memcount(4, 'r', dim_ens**2)

        CALL PDAF_seik_TtimesA(rank, dim_ens, OmegaT, TA)

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_lseik_resample:', debug, '  transform', TA

        CALL PDAF_timeit(24, 'old')

        CALL PDAF_timeit(18, 'new')

        ! *** Block formulation for resampling
        maxblksize = 200
        IF (mype == 0 .AND. screen > 0 .AND. screenout) &
             WRITE (*, '(a, 5x, a, i5)') 'PDAF', '--- use blocking with size ', maxblksize

        ALLOCATE(ens_block(maxblksize, dim_ens))
        IF (allocflag == 0) CALL PDAF_memcount(4, 'r', maxblksize * dim_ens)

        blocking: DO blklower = 1, dim_l, maxblksize

           blkupper = MIN(blklower + maxblksize - 1, dim_l)

           ! Store old state ensemble
           DO col = 1, dim_ens
              ens_block(1 : blkupper - blklower + 1, col) &
                   = ens_l(blklower : blkupper, col)
           END DO

           DO col = 1, dim_ens
              ens_l(blklower : blkupper, col) = state_l(blklower : blkupper)
           END DO

           ! *** X = state+ sqrt(FAC) state_ens T A^T (A^T stored in OmegaT) ***
           ! *** Here FAC depends on the use definition of the covariance    ***
           ! *** matrix P using a factor (r+1)^-1 or r^-1.                   ***

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
        END DO blocking
        
        CALL PDAF_timeit(18, 'old')

        DEALLOCATE(ens_block, TA)

     ELSE SolveOK

        ! Solve for A failed
        WRITE (*, '(/5x, a, i10, a/)') &
             'PDAF-ERROR(2): Problem with solve for A in SEIK_RESAMPLE - domain ', &
             domain_p, ' !!!'
        flag = 2

        CALL PDAF_timeit(24, 'old')

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

END SUBROUTINE PDAFlseik_resample

END MODULE PDAF_lseik_analysis
