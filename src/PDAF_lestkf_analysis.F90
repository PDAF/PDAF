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
!> LESTKF analysis/transformation
!!
!! Analysis step of the LESTKF filter with direct
!! transformation of the forecast into the 
!! analysis ensemble. This variant does not
!! compute the analysis state, but only the
!! analysis ensemble, whose mean is the analysis
!! state.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2011-09 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAF_lestkf_analysis

CONTAINS
SUBROUTINE PDAF_lestkf_ana(domain_p, step, dim_l, dim_obs_l, dim_ens, &
     rank, state_l, Ainv_l, ens_l, HL_l, HXbar_l, &
     obs_l, OmegaT_in, forget, U_prodRinvA_l, &
     envar_mode, type_sqrt, TA, screen, debug, flag)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_parallel, &
       ONLY: mype
  USE PDAF_analysis_utils, &
       ONLY: PDAF_estkf_AOmega, PDAF_estkf_OmegaA
#if defined (_OPENMP)
  USE omp_lib, &
       ONLY: omp_get_num_threads, omp_get_thread_num
#endif

  IMPLICIT NONE

! *** Arguments ***
!  Variable naming scheme:
!    suffix _p: Denotes a full variable on the PE-local domain
!    suffix _l: Denotes a local variable on the current analysis domain
  INTEGER, INTENT(in) :: domain_p      !< Current local analysis domain
  INTEGER, INTENT(in) :: step          !< Current time step
  INTEGER, INTENT(in) :: dim_l         !< State dimension on local analysis domain
  INTEGER, INTENT(in) :: dim_obs_l     !< Size of obs. vector on local ana. domain
  INTEGER, INTENT(in) :: dim_ens       !< Size of ensemble 
  INTEGER, INTENT(in) :: rank          !< Rank of initial covariance matrix
  REAL, INTENT(inout) :: state_l(dim_l)           !<state on local analysis domain
  REAL, INTENT(inout) :: Ainv_l(rank, rank)       !< Inverse of matrix U - temporary use only
  REAL, INTENT(inout) :: ens_l(dim_l, dim_ens)    !< Local state ensemble
  REAL, INTENT(inout) :: HL_l(dim_obs_l, dim_ens) !< Local observed state ensemble (perturbation)
  REAL, INTENT(in) :: HXbar_l(dim_obs_l)          !< Local observed ensemble mean
  REAL, INTENT(in) :: obs_l(dim_obs_l)            !< Local observation vector
  REAL, INTENT(in) :: OmegaT_in(rank, dim_ens)    !< Matrix Omega
  REAL, INTENT(inout) :: forget        !< Forgetting factor
  INTEGER, INTENT(in) :: envar_mode    !< Flag whether routine is called from 3DVar for special functionality
  INTEGER, INTENT(in) :: type_sqrt     !< Type of square-root of A
                                       !< (0): symmetric sqrt; (1): Cholesky decomposition
  REAL, INTENT(inout) :: TA(dim_ens, dim_ens)    !< Ensemble transformation matrix
  INTEGER, INTENT(in) :: screen        !< Verbosity flag
  INTEGER, INTENT(in) :: debug         !< Flag for writing debug output
  INTEGER, INTENT(inout) :: flag       !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA_l            !< Provide product R^-1 A for local analysis domain
       
! *** local variables ***
  INTEGER :: i, j, col, row            ! Counters
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  INTEGER :: lib_info                  ! Status flag for LAPACK calls
  INTEGER :: ldwork                    ! Size of work array for SYEVTYPE
  INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
  REAL    :: fac                       ! Temporary variable sqrt(dim_ens) or sqrt(rank)
  INTEGER, SAVE :: lastdomain = -1     ! store domain index
  LOGICAL, SAVE :: screenout = .true.  ! Whether to print information to stdout
  REAL, ALLOCATABLE :: RiHL_l(:,:)     ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: innov_l(:)      ! observation innovation
  REAL, ALLOCATABLE :: RiHLd_l(:)      ! local RiHLd
  REAL, ALLOCATABLE :: VRiHLd_l(:)     ! Temporary vector for analysis
  REAL, ALLOCATABLE :: tmp_Ainv_l(:,:) ! Temporary storage of Ainv
  REAL, ALLOCATABLE :: OmegaT(:,:)     ! Transpose of Omega
  REAL, ALLOCATABLE :: ens_blk(:,:)    ! Temporary blocked state ensemble
  REAL, ALLOCATABLE :: svals(:)        ! Singular values of Ainv
  REAL, ALLOCATABLE :: work(:)         ! Work array for syevTYPE
  INTEGER, ALLOCATABLE :: ipiv(:)      ! vector of pivot indices for GESVTYPE
  INTEGER, SAVE :: mythread, nthreads  ! Thread variables for OpenMP

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

     ! Output, only in case of OpenMP parallelization
#if defined (_OPENMP)
     IF (mype == 0 .AND. screen > 0 .AND. screenout) THEN
        WRITE (*,'(a, 5x, a, i5, a)') &
             'PDAF', '--- Use OpenMP parallelization with ', nthreads, ' threads'
     END IF
#endif
  END IF

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lestkf_analysis -- START'

  CALL PDAF_timeit(51, 'old')


! **************************
! *** Compute innovation ***
! ***   d = y - H x      ***
! **************************

  CALL PDAF_timeit(16, 'new')
  CALL PDAF_timeit(20, 'new')

  haveobsB: IF (dim_obs_l > 0) THEN
     ! *** The residual only exists for domains with observations ***

     CALL PDAF_timeit(51, 'new')

     ALLOCATE(innov_l(dim_obs_l))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l)

     innov_l = obs_l - HXbar_l

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, '  innovation d_l', innov_l

     CALL PDAF_timeit(51, 'old')

  END IF haveobsB
  CALL PDAF_timeit(20, 'old')


! **********************************************
! ***   Compute analyzed matrix Ainv         ***
! ***                                        ***
! ***  -1                       T  -1        ***
! *** A  = forget*(N-1) I + (HL)  R  (HL)    ***
! **********************************************

  CALL PDAF_timeit(21, 'new')

  haveobsA: IF (dim_obs_l > 0) THEN
     ! *** The contribution of observation matrix ist only ***
     ! *** computed for domains with observations          ***

     CALL PDAF_timeit(51, 'new')

     ! Projection HL = [Hx_1 ... Hx_N] Omega
     CALL PDAF_estkf_AOmega(dim_obs_l, dim_ens, HL_l)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, '  HXT_l', HL_l(:, 1:dim_ens-1)

     CALL PDAF_timeit(51, 'old')


     ! ***                RiHL = Rinv HL                 ***
     ! *** this is implemented as a subroutine thus that ***
     ! *** Rinv does not need to be allocated explicitly ***
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lestkf_analysis -- call prodRinvA_l'

     ALLOCATE(RiHL_l(dim_obs_l, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l * rank)

     CALL PDAF_timeit(48, 'new')
     CALL U_prodRinvA_l(domain_p, step, dim_obs_l, rank, obs_l, HL_l, RiHL_l)
     CALL PDAF_timeit(48, 'old')

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, '  R^-1(HXT_l)', RiHL_l
 
     CALL PDAF_timeit(51, 'new')

     ! *** Initialize Ainv = (N-1) I ***
     Ainv_l = 0.0
     DO i = 1, rank
        Ainv_l(i,i) = REAL(dim_ens - 1)
     END DO

     ! ***             T        ***
     ! ***  Compute  HL  RiHL   ***
     ALLOCATE(tmp_Ainv_l(rank, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank**2)

     CALL gemmTYPE('t', 'n', rank, rank, dim_obs_l, &
          1.0, HL_l, dim_obs_l, RiHL_l, dim_obs_l, &
          0.0, tmp_Ainv_l, rank)

     CALL PDAF_timeit(51, 'old')

  ELSE haveobsA
     ! *** For domains with dim_obs_l=0 there is no ***
     ! *** direct observation-contribution to Ainv  ***

     CALL PDAF_timeit(51, 'new')

     ! *** Initialize Ainv = (N-1) I ***
     Ainv_l = 0.0
     DO i = 1, rank
        Ainv_l(i,i) = REAL(dim_ens - 1)
     END DO

     ! No observation-contribution to Ainv from this domain
     ALLOCATE(tmp_Ainv_l(rank, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank**2)

     tmp_Ainv_l = 0.0

     CALL PDAF_timeit(51, 'old')

  END IF haveobsA

  ! *** Complete computation of Ainv  ***
  ! ***   -1                T         ***
  ! ***  A  = forget I  + HL RiHL     ***
  CALL PDAF_timeit(51, 'new')
  Ainv_l = forget * Ainv_l + tmp_Ainv_l
  CALL PDAF_timeit(51, 'old')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, '  A^-1_l', Ainv_l

  CALL PDAF_timeit(21, 'old')
  CALL PDAF_timeit(16, 'old')


! ***********************************************
! *** Compute weight for model state update   ***
! ***                                         ***
! ***              T                    f     ***
! ***    w = A RiHL d  with d = (y - H x )    ***
! ***********************************************

  CALL PDAF_timeit(17, 'new')
  CALL PDAF_timeit(22, 'new')
  CALL PDAF_timeit(51, 'new')

  ! *** Compute RiHLd = RiHL^T d ***
  ALLOCATE(RiHLd_l(rank))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank)

  haveobsC: IF (dim_obs_l > 0) THEN
     ! *** RiHLd_p/=0 only with observations ***

     CALL gemvTYPE('t', dim_obs_l, rank, 1.0, RiHL_l, &
          dim_obs_l, innov_l, 1, 0.0, RiHLd_l, 1)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, '  (HXT_l R^-1)^T d_l', RiHLd_l

     DEALLOCATE(RiHL_l, innov_l)

  ELSE haveobsC

     RiHLd_l = 0.0

  END IF haveobsC


  ! *** Compute weight vector for state analysis:     ***
  ! ***          w = A RiHLd                          ***
  ! *** For this, two variants are implemented:       ***
  ! *** 1. solve for w in:                            ***
  ! ***           -1                                  ***
  ! ***          A  w = RiHLd                         ***
  ! ***   We use the LAPACK routine gesvTYPE          ***
  ! *** 2. Compute singular value decomposition       ***
  ! ***   of Ainv: Ainv = USV^T                       ***
  ! ***   Then: A = U S^(-1) V^T                      ***
  ! ***   This is combined with a symmetric           ***
  ! ***   square-root for the ensemble transformation ***

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lestkf_analysis -- type_sqrt', type_sqrt
  typeainv1: IF (type_sqrt==1) THEN
     ! *** Variant 1: Solve Ainv w= RiHLd for w

     ALLOCATE(ipiv(rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'i', rank)

     ! save matrix Ainv
     tmp_Ainv_l = Ainv_l

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, &
          '  Invert A^-1_l using solver GESV'

     ! call solver (gesvTYPE - LU solver)
     CALL gesvTYPE(rank, 1, tmp_Ainv_l, rank, ipiv, &
          RiHLd_l, rank, lib_info)
     DEALLOCATE(ipiv)

  ELSE typeainv1
     ! *** Variant 2: Invert Ainv using SVD

     ALLOCATE(svals(rank))
     ALLOCATE(work(3 * rank))
     ldwork = 3 * rank
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 4 * rank)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, &
          '  Compute eigenvalue decomposition of A^-1_l'

     ! Compute SVD of Ainv
     CALL syevTYPE('v', 'l', rank, Ainv_l, rank, svals, work, ldwork, lib_info)

     DEALLOCATE(work)

     ! Compute product A RiHLd
     IF (lib_info==0) THEN
        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, '  eigenvalues', svals

        ALLOCATE(VRiHLd_l(rank))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank)

        CALL gemvTYPE('t', rank, rank, 1.0, Ainv_l, &
             rank, RiHLd_l, 1, 0.0, VRiHLd_l, 1)
     
        DO row = 1,rank
           VRiHLd_l(row) = VRiHLd_l(row) / svals(row)
        END DO
  
        CALL gemvTYPE('n', rank, rank, 1.0, Ainv_l, &
             rank, VRiHLd_l, 1, 0.0, RiHLd_l, 1)

        DEALLOCATE(VRiHLd_l)
     END IF
  END IF typeainv1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, '  wbar_l', RiHLd_l

  ! *** check if SVD was successful
  IF (lib_info == 0) THEN
     flag = 0
  ELSE
     WRITE (*, '(/5x, a, i10, a/)') &
          'PDAF-ERROR(1): Domain ', domain_p, ' Problem in computation of analysis weights !!!'
     flag = 1
  END IF

  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(22, 'old')


! ****************************************************************
! ***     Transform state ensemble                             ***
! ***              a   _f   f                                  ***
! ***             X  = X + X  W                                ***
! *** with                                   -T      T         ***
! ***          W = Omega (RiHLd + sqrt(N-1) C   Omega )        ***
! ****************************************************************

! *** Prepare weight matrix for ensemble transformation ***

  CALL PDAF_timeit(51, 'new')
  check1: IF (flag == 0) THEN

     ! allocate fields
     ALLOCATE(OmegaT(rank, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank * dim_ens)

     IF (mype == 0 .AND. screen > 0 .AND. screenout) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', 'Transform state ensemble'
        IF (type_sqrt == 1) THEN
           WRITE (*, '(a, 5x, a)') 'PDAF', '--- use Cholesky square-root of A'
        ELSE
           WRITE (*, '(a, 5x, a)') 'PDAF', '--- use symmetric square-root of A'
        END IF
     END IF

     CALL PDAF_timeit(23, 'new')

     ! Part 1: square-root of A
     typeainv2: IF (type_sqrt == 1) THEN
        ! Variant, if Ainv has been inverted above by solving

        ! Store Ainv for temporary use
        tmp_Ainv_l(:, :) = Ainv_l(:, :)

        ! Cholesky decomposition of tmp_Ainv_l = C C^T
        CALL potrfTYPE('l', rank, tmp_Ainv_l, rank, lib_info)

        ! check if Cholesky decomposition was successful
        CholeskyOK: IF (lib_info == 0) THEN
           ! Decomposition OK, continue
           flag = 0
        ELSE
           ! Decomposition failed
           WRITE (*, '(/5x, a, i8, a/)') &
                '!!! PDAF-ERROR(3): Problem with Cholesky decomposition of Ainv - domain ', &
                domain_p, ' !!!'
           flag = 3
        ENDIF CholeskyOK

     ELSE typeainv2
        ! Variant, if SVD inversion of Ainv has been performed

        ! Use OmegaT as temporary array
        DO col = 1, rank
           DO row = 1, rank
              OmegaT(row, col) = Ainv_l(row, col) / SQRT(svals(col))
           END DO
        END DO

        CALL gemmTYPE('n', 't', rank, rank, rank, &
             1.0, OmegaT, rank, Ainv_l, rank, &
             0.0, tmp_Ainv_l, rank)
        DEALLOCATE(svals)

        ! Set flag
        flag = 0

     END IF typeainv2

  END IF check1

  check2: IF (flag == 0) THEN
    
     ! *** Part 2: Product A^(1/2) Omega ***
    
     IF (type_sqrt == 1) THEN
        ! Initialize the matrix Omega from argument OmegaT_in
        OmegaT = OmegaT_in

       ! A = (Omega C^(-1)) by solving Ct A = OmegaT for A
        CALL trtrsTYPE('L', 'T', 'N', rank, dim_ens, &
             tmp_Ainv_l, rank, OmegaT, rank, lib_info)
     ELSE
        ! TEMP_AINV already contains matrix C (no more inversion)

        CALL gemmTYPE('n', 'n', rank, dim_ens, rank, &
             1.0, tmp_Ainv_l, rank, OmegaT_in, rank, &
             0.0, OmegaT, rank)

        lib_info = 0

     END IF

     ! check if solve was successful
     solveOK: IF (lib_info == 0) THEN
        ! Solve for A OK, continue
        flag = 0
     ELSE
        ! Solve for A failed
        WRITE (*, '(/5x, a, i8, a/)') &
             '!!! PDAF-ERROR(2): Problem in computation of transformation matrix - domain ', &
             domain_p, ' !!!'
        flag = 2
     END IF solveOK

  END IF check2

  check3: IF (flag == 0) THEN

     ! *** Part 4: Add RiHLd and multiply by scaling factor

     ! *** Add RiHLd to A^T stored in OmegaT
     fac = SQRT(REAL(dim_ens - 1))

     IF (envar_mode == 0) THEN
        DO j = 1, dim_ens
           DO i = 1, rank
              OmegaT(i, j) = fac * OmegaT(i, j) + RiHLd_l(i)
           END DO
        END DO
     ELSE
        ! For ensemble 3D-Var update only ensemble perturbations
        DO j = 1, dim_ens
           DO i = 1, rank
              OmegaT(i, j) = fac * OmegaT(i, j)
           END DO
        END DO
     END IF

     DEALLOCATE(RiHLd_l)
      
     ! *** Omega A^T (A^T stored in OmegaT_l) ***
     CALL PDAF_estkf_OmegaA(rank, dim_ens, OmegaT, TA)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, '  transform', TA

     CALL PDAF_timeit(23, 'old')
     CALL PDAF_timeit(17, 'old')


! *** Perform ensemble transformation ***

     CALL PDAF_timeit(18, 'new')

     ! Use block formulation for transformation
     maxblksize = 200
     IF (mype == 0 .AND. screen > 0 .AND. screenout) &
          WRITE (*, '(a, 5x, a, i5)') &
          'PDAF', '--- use blocking with size ', maxblksize

     ALLOCATE(ens_blk(maxblksize, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', maxblksize * dim_ens)

     blocking: DO blklower = 1, dim_l, maxblksize

        blkupper = MIN(blklower + maxblksize - 1, dim_l)

        ! Store old state ensemble
        DO col = 1, dim_ens
           ens_blk(1 : blkupper - blklower + 1, col) &
                = ens_l(blklower : blkupper, col)
        END DO

        DO col = 1, dim_ens
           ens_l(blklower : blkupper, col) = state_l(blklower : blkupper)
        END DO

        !                        a  _f   f    T
        ! Transform ensemble:   X = X + X  T(A )
        CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
             1.0, ens_blk, maxblksize, TA, dim_ens, &
             1.0, ens_l(blklower, 1), dim_l)

     END DO blocking

     CALL PDAF_timeit(18, 'old')
        
     DEALLOCATE(ens_blk)
     DEALLOCATE(OmegaT)

  ELSE
     CALL PDAF_timeit(17, 'old')
  END IF check3

  CALL PDAF_timeit(51, 'old')


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(tmp_Ainv_l)

  IF (allocflag == 0) allocflag = 1

  ! Store domain index
  lastdomain = domain_p

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lestkf_analysis -- END'

END SUBROUTINE PDAF_lestkf_ana

END MODULE PDAF_lestkf_analysis
