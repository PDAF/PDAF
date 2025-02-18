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
!> LKNETF analysis using T-matrix
!!
!! Synchronous (1-step) analysis update of the LKNETF. The analysis
!! forumulation used a matrix T analogous to the LESTKF and LETKF.
!!
!! The implementation also supports an adaptive forgetting factor.
!!
!! Variant for domain decomposed states.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2017-08 - Lars Nerger - Initial code based on LETKF
!! * Other revisions - see repository log
!!
MODULE PDAF_lknetf_analysis_sync

CONTAINS
SUBROUTINE PDAF_lknetf_analysis_T(domain_p, step, dim_l, dim_obs_l, dim_ens, &
     state_l, Ainv_l, ens_l, HX_l, HXbar_l, obs_l, &
     rndmat, forget, U_prodRinvA_l, U_init_obsvar_l, U_likelihood_l, &
     screen, type_forget, eff_dimens, type_hyb, hyb_g, hyb_k, &
     gamma, skew_mabs, kurt_mabs, flag)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: mype
  USE PDAF_lknetf, &
       ONLY: type_trans, debug
  USE PDAF_mod_filter, &
       ONLY: obs_member
  USE PDAF_analysis_utils, &
       ONLY: PDAF_subtract_rowmean, PDAF_subtract_colmean, PDAF_set_forget, &
       PDAF_set_forget_local
  USE PDAF_lknetf, &
       ONLY: PDAF_lknetf_set_gamma
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
  REAL, INTENT(inout) :: state_l(dim_l)           !< local forecast state
  REAL, INTENT(out)   :: Ainv_l(dim_ens, dim_ens) !< on exit: local weight matrix for ensemble transformation
  REAL, INTENT(inout) :: ens_l(dim_l, dim_ens)    !< Local state ensemble
  REAL, INTENT(in) :: HX_l(dim_obs_l, dim_ens)    !< local observed state ens.
  REAL, INTENT(in) :: HXbar_l(dim_obs_l)          !< local observed ens. mean
  REAL, INTENT(inout) :: rndmat(dim_ens, dim_ens) !< Global random rotation matrix
  REAL, INTENT(in) :: obs_l(dim_obs_l) !< Local observation vector
  REAL, INTENT(inout) :: forget        !< Forgetting factor
  REAL, INTENT(inout) :: eff_dimens(1) !< Effective ensemble size
  INTEGER, INTENT(in) :: screen        !< Verbosity flag
  INTEGER, INTENT(in) :: type_forget   !< Type of forgetting factor
  INTEGER, INTENT(in) :: type_hyb      !< Type of hybrid weight
  REAL, INTENT(in) :: hyb_g            !< Prescribed hybrid weight for state transformation
  REAL, INTENT(in) :: hyb_k            !< Scale factor kappa (for type_hyb 3 and 4)
  REAL, INTENT(inout) :: gamma(1)      !< Hybrid weight for state transformation
  REAL, INTENT(inout) :: skew_mabs(1)  !< Mean absolute skewness
  REAL, INTENT(inout) :: kurt_mabs(1)  !< Mean absolute kurtosis
  INTEGER, INTENT(inout) :: flag       !< Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: &
       U_init_obsvar_l, &              !< Initialize local mean observation error variance
       U_prodRinvA_l, &                !< Provide product R^-1 A for local analysis domain
       U_likelihood_l                  !< Provide likelihood of an ensemble state
       
! *** local variables ***
  INTEGER :: i, j, member, col, row  ! Counters
  INTEGER, SAVE :: allocflag = 0     ! Flag whether first time allocation is done
  INTEGER :: syev_info               ! Status flag for SYEV
  INTEGER :: ldwork                  ! Size of work array for SYEV
  INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
  REAL    :: sqrtNm1                 ! Temporary variable: sqrt(dim_ens-1)
  REAL    :: aforget                 ! adaptive forgetting factor value
  REAL :: weight                     ! Ensemble weight (likelihood)
  REAL :: fac                        ! Multiplication factor
  INTEGER, SAVE :: lastdomain = -1   ! store domain index
  LOGICAL :: screenout = .TRUE.      ! Whether to print information to stdout
  REAL, ALLOCATABLE :: HZ_l(:,:)     ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: RiHZ_l(:,:)   ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: resid_l(:)    ! local observation residual
  REAL, ALLOCATABLE :: w_etkf_l(:)   ! local RiHZd
  REAL, ALLOCATABLE :: VRiHZd_l(:)   ! Temporary vector for analysis
  REAL, ALLOCATABLE :: tmp_Ainv_l(:,:) ! Temporary storage of Ainv
  REAL, ALLOCATABLE :: Asqrt_l(:,:)  ! Temporary for square-root of U
  REAL, ALLOCATABLE :: ens_blk(:,:)  ! Temporary block of state ensemble
  REAL, ALLOCATABLE :: svals(:)      ! Singular values of Ainv
  REAL, ALLOCATABLE :: work(:)       ! Work array for SYEV
  REAL, ALLOCATABLE :: resid_i(:)    ! Residual
  REAL, ALLOCATABLE :: weights(:)    ! weight vector
  REAL, ALLOCATABLE :: A_tmp(:,:)    ! Temporary matrix for NETF transformation matrix
  REAL, ALLOCATABLE :: Anetf(:,:)    ! Temporary matrix for NETF transformation matrix
  REAL, ALLOCATABLE :: Asqrt(:,:)    ! Temporary matrix for NETF transformation matrix
  REAL :: total_weight               ! Sum of weight

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
     screenout = .FALSE.
  ELSE
     screenout = .TRUE.

     ! In case of OpenMP, let only thread 0 write output to the screen
     IF (mythread>0) screenout = .FALSE.

     ! Output, only in case of OpenMP parallelization
#if defined (_OPENMP)
     IF (screenout) THEN
        WRITE (*,'(a, 5x, a, i5, a)') &
             'PDAF', '--- Use OpenMP parallelization with ', nthreads, ' threads'
     END IF
#endif
  END IF

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_analysis -- START'

  CALL PDAF_timeit(51, 'old')


! ****************************************************
! *** 1. Weight vector for ensemble mean from ETKF ***
! ****************************************************


! ************************
! *** Compute residual ***
! ***   d = y - H x    ***
! ************************

  CALL PDAF_timeit(12, 'new')
  CALL PDAF_timeit(51, 'new')

  IF (mype == 0 .AND. screen > 0 .AND. screenout) THEN
     WRITE (*, '(a, 5x, a)') 'PDAF', 'Compute ETKF transform matrix'
  END IF

  haveobsB: IF (dim_obs_l > 0) THEN
     ! *** The residual only exists for domains with observations ***

     ALLOCATE(resid_l(dim_obs_l))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l)

     ! Get residual as difference of observation and observed state
     resid_l = obs_l - HXbar_l

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lknetf_analysis:', debug, '  innovation d_l', resid_l

  END IF haveobsB

  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(12, 'old')


! **********************************************
! ***   Compute analyzed matrix Ainv         ***
! ***                                        ***
! ***     -1                 T  -1           ***
! ***    A  = forget I + (HZ)  R   HZ        ***
! **********************************************

  CALL PDAF_timeit(10, 'new')

  haveobsA: IF (dim_obs_l > 0) THEN
     ! *** The contribution of observation matrix ist only ***
     ! *** computed for domains with observations          ***

     CALL PDAF_timeit(51, 'new')

     ALLOCATE(HZ_l(dim_obs_l, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l*dim_ens)
     HZ_l = HX_l

     ! *** Set the value of the forgetting factor  ***
     ! *** Inserted here, because HZ_l is required ***
     aforget = forget
     IF (type_forget == 6) THEN
        CALL PDAF_set_forget_local(domain_p, step, dim_obs_l, dim_ens, &
             HZ_l, HXbar_l, obs_l, U_init_obsvar_l, forget, aforget)
     ENDIF

     CALL PDAF_timeit(30, 'new')

     ! Subtract ensemble mean: HZ = [Hx_1 ... Hx_N] T
     CALL PDAF_subtract_rowmean(dim_obs_l, dim_ens, HZ_l)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lknetf_analysis:', debug, '  ETKF HXT_l', HZ_l(:, 1:dim_ens-1)

     CALL PDAF_timeit(30, 'old')
     CALL PDAF_timeit(51, 'old')
     CALL PDAF_timeit(31, 'new')


     ! ***                RiHZ = Rinv HZ                
     ! *** This is implemented as a subroutine thus that
     ! *** Rinv does not need to be allocated explicitly.
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_analysis -- call prodRinvA_l'

     ALLOCATE(RiHZ_l(dim_obs_l, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l * dim_ens)

     CALL PDAF_timeit(48, 'new')
     CALL U_prodRinvA_l(domain_p, step, dim_obs_l, dim_ens, obs_l, HZ_l, RiHZ_l)
     CALL PDAF_timeit(48, 'old')

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lknetf_analysis:', debug, '  ETKF R^-1(HXT_l)', RiHZ_l

     CALL PDAF_timeit(51, 'new')

     ! *** Initialize Ainv = (N-1) I ***
     Ainv_l = 0.0
     DO i = 1, dim_ens
        Ainv_l(i, i) = REAL(dim_ens - 1)
     END DO

     ! ***             T        ***
     ! ***  Compute  HZ  RiHZ   ***

     ALLOCATE(tmp_Ainv_l(dim_ens, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

     CALL gemmTYPE('t', 'n', dim_ens, dim_ens, dim_obs_l, &
          1.0, HZ_l, dim_obs_l, RiHZ_l, dim_obs_l, &
          0.0, tmp_Ainv_l, dim_ens)

     DEALLOCATE(HZ_l)
     CALL PDAF_timeit(51, 'old')

  ELSE haveobsA
     ! *** For domains with dim_obs_l=0 there is no ***
     ! *** direct observation-contribution to Ainv  ***

     CALL PDAF_timeit(31, 'new')
     CALL PDAF_timeit(51, 'new')

     ! *** Initialize Ainv = (N-1) I ***
     Ainv_l = 0.0
     DO i = 1, dim_ens
        Ainv_l(i, i) = REAL(dim_ens - 1)
     END DO

     ! No observation-contribution to Ainv from this domain
     ALLOCATE(tmp_Ainv_l(dim_ens, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

     tmp_Ainv_l = 0.0

     CALL PDAF_timeit(51, 'old')

  END IF haveobsA


  ! *** Complete computation of Ainv  ***
  ! ***   -1          -1    T         ***
  ! ***  A  =        A  + HZ RiHZ     ***

  CALL PDAF_timeit(51, 'new')
  IF (type_forget<5) THEN
     ! Usually the forgetting factor is not applied here
     Ainv_l = Ainv_l + tmp_Ainv_l
  ELSE
     Ainv_l = aforget*Ainv_l + tmp_Ainv_l
  END IF
  CALL PDAF_timeit(51, 'new')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lknetf_analysis:', debug, '  ETKF A^-1_l', Ainv_l

  CALL PDAF_timeit(31, 'old')
  CALL PDAF_timeit(10, 'old')


! ***********************************************
! *** Compute weight for model state update   ***
! ***                                         ***
! ***              T                    f     ***
! ***    w = A RiHZ d  with d = (y - H x )    ***
! ***                                         ***
! ***********************************************

  CALL PDAF_timeit(13, 'new')
  CALL PDAF_timeit(51, 'new')

  ! *** Compute RiHZd = RiHZ^T d ***
  ALLOCATE(w_etkf_l(dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)

  haveobsC: IF (dim_obs_l > 0) THEN
     ! *** RiHLd_p/=0 only with observations ***

     CALL gemvTYPE('t', dim_obs_l, dim_ens, 1.0, RiHZ_l, &
          dim_obs_l, resid_l, 1, 0.0, w_etkf_l, 1)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lknetf_analysis:', debug, '  ETKF (HXT_l R^-1)^T d_l', w_etkf_l

     DEALLOCATE(RiHZ_l, resid_l)

  ELSE haveobsC

     w_etkf_l = 0.0

  END IF haveobsC


  ! *** Compute weight vector for state analysis:        ***
  ! ***          w = A RiHZd                             ***
  ! *** Use singular value decomposition of Ainv         ***
  ! ***        Ainv = ASB^T                              ***
  ! *** Then: A = A S^(-1) B                             ***
  ! *** The decomposition is also used for the symmetric ***
  ! *** square-root for the ensemble transformation.     ***

  ! *** Invert Ainv using SVD
  ALLOCATE(svals(dim_ens))
  ALLOCATE(work(3 * dim_ens))
  ldwork = 3 * dim_ens
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3 * dim_ens)

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lknetf_analysis:', debug, &
       '  Compute eigenvalue decomposition of ETKF A^-1_l'

  ! Compute SVD of Ainv
  CALL syevTYPE('v', 'l', dim_ens, Ainv_l, dim_ens, svals, work, ldwork, syev_info)

  ! *** check if SVD was successful
  IF (syev_info == 0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lknetf_analysis:', debug, '  eigenvalues', svals

     flag = 0
  ELSE
     WRITE (*, '(/5x, a, i10, a/)') &
          'PDAF-ERROR(1): Domain ', domain_p, ' Problem in SVD of inverse of A_etkf !!!'
     flag = 1
  END IF

  ! *** Compute w = U RiHZd stored in RiHZd
  check0: IF (flag == 0) THEN

     ALLOCATE(VRiHZd_l(dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)

     CALL gemvTYPE('t', dim_ens, dim_ens, 1.0, Ainv_l, &
          dim_ens, w_etkf_l, 1, 0.0, VRiHZd_l, 1)
     
     DO row = 1, dim_ens
        VRiHZd_l(row) = VRiHZd_l(row) / svals(row)
     END DO
  
     CALL gemvTYPE('n', dim_ens, dim_ens, 1.0, Ainv_l, &
          dim_ens, VRiHZd_l, 1, 0.0, w_etkf_l, 1)

     DEALLOCATE(VRiHZd_l)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lknetf_analysis:', debug, '  ETKF A(HXT_l R^-1)^T d_l', w_etkf_l

  END IF check0

  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(13, 'old')


! **************************************************************
! *** 2. Weight matrix for ensemble transformation from ETKF ***
! **************************************************************

  CALL PDAF_timeit(51, 'new')

  CALL PDAF_timeit(20, 'new')

  ! Part 1: square-root of U
  DO col = 1, dim_ens
     DO row = 1, dim_ens
        tmp_Ainv_l(row, col) = Ainv_l(row, col) / SQRT(svals(col))
     END DO
  END DO

  ALLOCATE(Asqrt_l(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

  sqrtNm1 = SQRT(REAL(dim_ens-1))
  CALL gemmTYPE('n', 't', dim_ens, dim_ens, dim_ens, &
       sqrtNm1, tmp_Ainv_l, dim_ens, Ainv_l, dim_ens, &
       0.0, Asqrt_l, dim_ens)

  ! Optional 
  ! Multiply by orthogonal random matrix with eigenvector (1,...,1)^T
  multrnd: IF (type_trans == 0) THEN
     CALL gemmTYPE('n', 'n', dim_ens, dim_ens, dim_ens, &
          1.0, Asqrt_l, dim_ens, rndmat, dim_ens, &
          0.0, tmp_Ainv_l, dim_ens)
  ELSE
     ! Non-random case
     tmp_Ainv_l = Asqrt_l
  END IF multrnd

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lknetf_analysis:', debug, '  ETKF A^1/2', tmp_Ainv_l

  CALL PDAF_timeit(20, 'old')


! ****************************************************
! *** 3. Weight vector for ensemble mean from NETF ***
! ****************************************************

  ! **********************************************
  ! *** Compute particle weights as likelihood ***
  ! **********************************************

  IF (mype == 0 .AND. screen > 0 .AND. screenout) THEN
     WRITE (*, '(a, 5x, a)') 'PDAF', 'Compute NETF transform matrix'
  END IF

  ! Allocate weight vector
  ALLOCATE(weights(dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)

  ! Allocate temporal array for obs-ens_i
  ALLOCATE(resid_i(dim_obs_l))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l)
     
  CALL PDAF_timeit(22, 'new')
  ! Get residual as difference of observation and observed state for 
  ! each ensemble member only on domains where observations are availible

  CALC_w: DO member = 1,dim_ens

     ! Store current member index
     obs_member = member

     ! Calculate local residual  
     resid_i = obs_l - HX_l(:,member)

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug: ', debug, &
             'PDAF_lknetf_analysis -- member', member
        WRITE (*,*) '++ PDAF-debug PDAF_lknetf_analysis:', debug, '  innovation d_l', resid_i
        WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_analysis -- call likelihood_l'
     end IF

     ! Compute likelihood
     CALL PDAF_timeit(49, 'new')
     CALL U_likelihood_l(domain_p, step, dim_obs_l, obs_l, resid_i, weight)
     CALL PDAF_timeit(49, 'old')
     weights(member) = weight

  END DO CALC_w

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lknetf_analysis:', debug, '  raw weights', weights

  CALL PDAF_timeit(51, 'new')

  ! Normalize weights
  total_weight = 0.0
  DO i = 1, dim_ens
     total_weight = total_weight + weights(i)
  END DO
  IF (total_weight /= 0.0) THEN
     weights = weights / total_weight

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lknetf_analysis:', debug, '  normalized weights', weights
  ELSE
     ! ERROR: weights are zero
     WRITE(*,'(/5x,a/)') 'WARNING: Zero weights in LNETF analysis step - reset to 1/dim_ens'
     weights = 1.0 / REAL(dim_ens)
  END IF

  DEALLOCATE(resid_i)

  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(22, 'old')


! **************************************************************
! *** 4. Weight matrix for ensemble transformation from NETF ***
! **************************************************************

  check1: IF (flag == 0) THEN

     ! ****************************************
     ! *** Calculate the transform matrix   ***
     ! ***      A= (diag(w)-w*w^t)          ***
     ! *** with the weights w               ***
     ! ****************************************

     CALL PDAF_timeit(23, 'new')

     ALLOCATE(Anetf(dim_ens,dim_ens)) 
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens*dim_ens)

     DO j = 1, dim_ens
        DO i = 1, dim_ens
           Anetf(i,j) = -weights(i) * weights(j)
        ENDDO
     ENDDO
     DO i = 1, dim_ens
        Anetf(i,i) = Anetf(i,i) + weights(i)
     END DO

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lknetf_analysis:', debug, '  NETF A', Anetf

     CALL PDAF_timeit(23, 'old')


     ! ***************************************
     ! *** Calculate square root of matrix ***
     ! ***************************************

     CALL PDAF_timeit(24, 'new')

     ! Compute symmetric square-root by SVD
     ALLOCATE(A_tmp(dim_ens,dim_ens))
     ALLOCATE(Asqrt(dim_ens,dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2*dim_ens*dim_ens)
     ldwork = 3*dim_ens
     flag = 0

     CALL PDAF_timeit(31, 'new')

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lknetf_analysis:', debug, &
          '  Compute eigenvalue decomposition of A'

     ! EVD
     CALL syevTYPE('v', 'l', dim_ens, Anetf, dim_ens, svals, work, ldwork, syev_info)

     CALL PDAF_timeit(31, 'old')

     IF (syev_info == 0) THEN
        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_lknetf_analysis:', debug, '  eigenvalues', svals
     ELSE
        WRITE(*,'(/5x,a,i7/)') 'Problem computing svd of W-ww^T in domain', domain_p
        flag = 1
     ENDIF
  
     CALL PDAF_timeit(32,'new')  

     ! Ensure to only use positive singular values - negative ones are numerical error
     DO i = 1, dim_ens
        IF (svals(i)>0.0) THEN
           svals(i) = SQRT(svals(i))
        ELSE
           svals(i) = 0.0
        END IF
     END DO

     DO i = 1, dim_ens
        DO j = 1, dim_ens
           Asqrt(j,i) = Anetf(j,i) * svals(i)
        END DO
     END DO

     DEALLOCATE(svals, work)

     CALL PDAF_timeit(32,'old')

     CALL PDAF_timeit(33,'new')

     ! Multiply with singular vectors
     CALL gemmTYPE('n', 't', dim_ens, dim_ens, dim_ens, 1.0, &
          Asqrt, dim_ens, Anetf, dim_ens, 0.0, A_tmp, dim_ens)

     CALL PDAF_timeit(33,'old')

     ! Multiply Asqrt by sqrt(m) to get unbiased ensemble 
     fac = SQRT(REAL(dim_ens))

     CALL PDAF_timeit(35,'new')

     ! Multiply Asqrt with random matrix and the factor 
     CALL gemmTYPE('n', 'n', dim_ens, dim_ens, dim_ens, &
          fac, A_tmp, dim_ens, rndmat, dim_ens, &
          0.0, Asqrt, dim_ens)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lknetf_analysis:', debug, '  NETF sqrt(N) A^1/2', Asqrt

     CALL PDAF_timeit(35,'old')     

  END IF check1

  
! **********************************
! *** 5. Determine hybrid weight ***
! **********************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_analysis -- Compute hybrid gamma'

  CALL PDAF_timeit(53,'new')
  CALL PDAF_lknetf_set_gamma(domain_p, dim_obs_l, dim_ens, &
       HX_l, HXbar_l, weights, type_hyb, hyb_g, hyb_k, &
       gamma, eff_dimens, skew_mabs, kurt_mabs, &
       screen, flag)
  CALL PDAF_timeit(53,'old')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lknetf_analysis:', debug, '  gamma', gamma


! ************************************************
! ***     6. Transform state ensemble          ***
! ***              a   _f   f                  ***
! ***             X  = X + X  W                ***
! *** The weight matrix W is stored in Ainv_l. ***
! ************************************************


  check2: IF (flag == 0) THEN

     ! **********************************************************
     ! *** Compute transform matrix W including hybridization ***
     ! **********************************************************

     CALL PDAF_timeit(20, 'new')

     IF (mype == 0 .AND. screen > 0 .AND. screenout) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', 'Perform ensemble transformation'
     END IF

     ! Compute hybrid transformation matrix for ensemble perturbations
     tmp_Ainv_l = gamma(1)*tmp_Ainv_l + (1.0-gamma(1))*Asqrt

     ! Total transformation matrix W = sqrt(U) + w (with hybridization)
     DO col = 1, dim_ens
        DO row = 1, dim_ens
           Ainv_l(row, col) = tmp_Ainv_l(row, col) &
                + gamma(1)*w_etkf_l(row) + (1.0-gamma(1))*weights(row)
        END DO
     END DO

     DEALLOCATE(tmp_Ainv_l, Asqrt, A_tmp)
     DEALLOCATE(w_etkf_l, Asqrt_l, weights)
      
     ! Part 4: T W
     CALL PDAF_subtract_colmean(dim_ens, dim_ens, Ainv_l)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lknetf_analysis:', debug, '  Hybrid transform', Ainv_l

     CALL PDAF_timeit(20, 'old')


     ! ***************************************
     ! *** Perform ensemble transformation ***
     ! ***************************************

     ! Use block formulation for transformation
     maxblksize = 200
     IF (mype == 0 .AND. screen > 0 .AND. screenout) &
          WRITE (*, '(a, 5x, a, i5)') &
          'PDAF', '--- use blocking with size ', maxblksize

     ALLOCATE(ens_blk(maxblksize, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(4, 'r', maxblksize * dim_ens)

     blocking: DO blklower = 1, dim_l, maxblksize

        blkupper = MIN(blklower + maxblksize - 1, dim_l)

        ! Store forecast ensemble
        CALL PDAF_timeit(21, 'new')
        DO col = 1, dim_ens
           ens_blk(1 : blkupper - blklower + 1, col) &
                = ens_l(blklower : blkupper, col)
        END DO

        ! Store mean forecast in ensemble matrix
        DO col = 1, dim_ens
           ens_l(blklower : blkupper, col) = state_l(blklower : blkupper)
        END DO
        CALL PDAF_timeit(21, 'old')

        !                        a  _f   f
        ! Transform ensemble:   X = X + X  TW
        CALL PDAF_timeit(22, 'new')

        CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
             1.0, ens_blk, maxblksize, Ainv_l, dim_ens, &
             1.0, ens_l(blklower, 1), dim_l)

        CALL PDAF_timeit(22, 'old')

     END DO blocking
        
     DEALLOCATE(ens_blk)

  END IF check2

  CALL PDAF_timeit(51, 'old')


! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1

  ! Store domain index
  lastdomain = domain_p

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_analysis -- END'

END SUBROUTINE PDAF_lknetf_analysis_T

END MODULE PDAF_lknetf_analysis_sync
