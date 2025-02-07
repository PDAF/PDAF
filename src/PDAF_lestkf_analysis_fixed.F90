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
!$Id: PDAF_lestkf_analysis_fixed.F90 1147 2023-03-12 16:14:34Z lnerger $
!BOP
!
! !ROUTINE: PDAF_lestkf_analysis_fixed --- LESTKF analysis without ens transformation
!
! !INTERFACE:
SUBROUTINE PDAF_lestkf_analysis_fixed(domain_p, step, dim_l, dim_obs_l, dim_ens, &
     rank, state_l, Ainv_l, ens_l, HL_l, HXbar_l, &
     obs_l, state_inc_l, forget, &
     U_prodRinvA_l, &
     incremental, type_sqrt, screen, debug, flag)

! !DESCRIPTION:
! Analysis step of the LESTKF filter with direct
! update of the state estimate, but no ensemble 
! transformation. The ensemble is only shifted
! to represent the analysis state. This variant
! is used for the filter variant with a fixed
! covariance matrix.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! __Revision history:__
! 2012-03 - Lars Nerger - Initial code
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
#if defined (_OPENMP)
  USE omp_lib, &
       ONLY: omp_get_num_threads, omp_get_thread_num
#endif

  IMPLICIT NONE

! !ARGUMENTS:
! ! Variable naming scheme:
! !   suffix _p: Denotes a full variable on the PE-local domain
! !   suffix _l: Denotes a local variable on the current analysis domain
  INTEGER, INTENT(in) :: domain_p    ! Current local analysis domain
  INTEGER, INTENT(in) :: step        ! Current time step
  INTEGER, INTENT(in) :: dim_l       ! State dimension on local analysis domain
  INTEGER, INTENT(in) :: dim_obs_l   ! Size of obs. vector on local ana. domain
  INTEGER, INTENT(in) :: dim_ens     ! Size of ensemble 
  INTEGER, INTENT(in) :: rank        ! Rank of initial covariance matrix
  REAL, INTENT(inout) :: state_l(dim_l)        ! on exit: state on local analysis domain
  REAL, INTENT(inout) :: Ainv_l(rank, rank)    ! Inverse of matrix U - temporary use only
  REAL, INTENT(inout) :: ens_l(dim_l, dim_ens) ! Local state ensemble
  REAL, INTENT(in) :: HL_l(dim_obs_l, dim_ens) ! Local observed state ensemble (perturbation)
  REAL, INTENT(in) :: HXbar_l(dim_obs_l)       ! Local observed ensemble mean
  REAL, INTENT(in) :: obs_l(dim_obs_l)         ! Local observation vector
  REAL, INTENT(in) :: state_inc_l(dim_l)       ! Local state increment
  REAL, INTENT(inout) :: forget      ! Forgetting factor
  INTEGER, INTENT(in) :: incremental ! Control incremental updating
  INTEGER, INTENT(in) :: type_sqrt   ! Type of square-root of A
                                     ! (0): symmetric sqrt; (1): Cholesky decomposition
  INTEGER, INTENT(in) :: screen      ! Verbosity flag
  INTEGER, INTENT(in) :: debug       ! Flag for writing debug output
  INTEGER, INTENT(inout) :: flag     ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA_l          ! Provide product R^-1 A for local analysis domain

! !CALLING SEQUENCE:
! Called by: PDAF_lestkf_update
! Calls: U_prodRinvA_l
! Calls: PDAF_estkf_AOmega
! Calls: PDAF_estkf_OmegaA
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: gemvTYPE (BLAS; dgemv or sgemv dependent on precision)
! Calls: gesvTYPE (LAPACK; dgesv or sgesv dependent on precision)
! Calls: syevTYPE (LAPACK; dsyev or ssyev dependent on precision)
!EOP
       
! *** local variables ***
  INTEGER :: i, col, row               ! Counters
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  INTEGER :: lib_info                  ! Status flag for LAPACK calls
  INTEGER :: ldwork                    ! Size of work array for SYEVTYPE
  INTEGER, SAVE :: lastdomain = -1     ! store domain index
  LOGICAL, SAVE :: screenout = .true.  ! Whether to print information to stdout
  REAL, ALLOCATABLE :: RiHL_l(:,:)     ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: innov_l(:)      ! observation innovation
  REAL, ALLOCATABLE :: RiHLd_l(:)      ! local RiHLd
  REAL, ALLOCATABLE :: VRiHLd_l(:)     ! Temporary vector for analysis
  REAL, ALLOCATABLE :: tmp_Ainv_l(:,:) ! Temporary storage of Ainv
  REAL, ALLOCATABLE :: TRiHLd_l(:,:) ! Temporary vector for analysis 
  REAL, ALLOCATABLE :: svals(:)      ! Singular values of Ainv
  REAL, ALLOCATABLE :: work(:)       ! Work array for syevTYPE
  INTEGER, ALLOCATABLE :: ipiv(:)    ! vector of pivot indices for GESVTYPE
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
     IF (screenout .AND. screen > 0) THEN
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

     ! save matrix Ainv
     tmp_Ainv_l = Ainv_l

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, &
          '  Compute eigenvalue decomposition of A^-1_l'

     ! Compute SVD of Ainv
     CALL syevTYPE('v', 'l', rank, tmp_Ainv_l, rank, svals, work, ldwork, lib_info)

     DEALLOCATE(work)

     ! Compute product A RiHLd
     IF (lib_info==0) THEN
        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, '  eigenvalues', svals

        ALLOCATE(VRiHLd_l(rank))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank)

        CALL gemvTYPE('t', rank, rank, 1.0, tmp_Ainv_l, &
             rank, RiHLd_l, 1, 0.0, VRiHLd_l, 1)
     
        DO row = 1,rank
           VRiHLd_l(row) = VRiHLd_l(row) / svals(row)
        END DO
  
        CALL gemvTYPE('n', rank, rank, 1.0, tmp_Ainv_l, &
             rank, VRiHLd_l, 1, 0.0, RiHLd_l, 1)

        DEALLOCATE(svals, VRiHLd_l)
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


! ************************************
! ***      update model state      ***
! ***                              ***
! ***     a   f   f                ***
! ***    x = x + X  Omega RiHLd    ***
! ***                              ***
! ************************************

  check1: IF (flag == 0) THEN

     ! ******************************
     ! *** Compute vector Omega w ***
     ! ******************************

     ALLOCATE(TRiHLd_l(dim_ens, 1))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)

     CALL PDAF_estkf_OmegaA(rank, 1, RiHLd_l, TRiHLd_l)
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, '  transform', TRiHLd_l

     DEALLOCATE(RiHLd_l)

     CALL PDAF_timeit(17, 'old')


     ! *****************************
     ! *** Update state estimate ***
     ! *****************************

     CALL PDAF_timeit(18, 'new')

     CALL gemvTYPE('n', dim_l, dim_ens, 1.0, ens_l, &
          dim_l, TRiHLd_l, 1, 0.0, state_inc_l, 1)
     DEALLOCATE(TRiHLd_l)
     
     ! Shift ensemble
     DO col = 1, dim_ens
        DO row = 1, dim_l
           ens_l(row, col) = ens_l(row, col) + state_inc_l(row)
        END DO
     END DO
     
     IF (incremental == 0) THEN
        ! update state here if incremental updating is not used
        state_l = state_l + state_inc_l
     END IF

     CALL PDAF_timeit(18, 'old')

  ELSE check1

     CALL PDAF_timeit(17, 'old')

  END IF check1
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

END SUBROUTINE PDAF_lestkf_analysis_fixed
