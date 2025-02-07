! Copyright (c) 2014-2025 Paul Kirchgessner
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
!BOP
!
! !ROUTINE: PDAF_lnetf_analysis --- local analysis step of LNETF
!
! !INTERFACE:
SUBROUTINE PDAF_lnetf_analysis(domain_p, step, dim_l, dim_obs_l, &
     dim_ens, ens_l, HX_l, obs_l, rndmat, &
     U_likelihood_l, type_forget, forget, &
     type_winf, limit_winf, cnt_small_svals, eff_dimens, T, &
     screen, debug, flag)

! !DESCRIPTION:
! Analysis step of the LNETF.
!
! Inflation has to be done BEFORE calling this routine !!!
!
! Variant for domain decomposed states.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! __Revision history:__
! 2014-05 - Paul Kirchgessner - Initial code based on LETKF
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
       ONLY: mype
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
  REAL, INTENT(inout) :: ens_l(dim_l, dim_ens)  ! Local state ensemble
  REAL, INTENT(in) :: HX_l(dim_obs_l, dim_ens)    ! Local observed state ensemble (perturbation)
  REAL, INTENT(in) :: obs_l(dim_obs_l)            ! Local observation vector
  REAL, INTENT(in) :: rndmat(dim_ens, dim_ens)  ! Global random rotation matrix
  INTEGER, INTENT(in) :: screen      ! Verbosity flag
  INTEGER, INTENT(in) :: type_forget ! Typ eof forgetting factor
  REAL, INTENT(in) :: forget         ! Forgetting factor
  INTEGER, INTENT(in) :: type_winf   ! Type of weights inflation
  REAL, INTENT(in) :: limit_winf     ! Limit for weights inflation
  INTEGER, INTENT(inout) :: cnt_small_svals   ! Number of small eigen values
  REAL, INTENT(inout) :: eff_dimens(1)        ! Effective ensemble size
  REAL, INTENT(inout) :: T(dim_ens, dim_ens)  ! local ensemble transformation matrix
  INTEGER, INTENT(in) :: debug       ! Flag for writing debug output
  INTEGER, INTENT(inout) :: flag     ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_likelihood_l         ! Compute observation likelihood for an ensemble member

! !CALLING SEQUENCE:
! Called by: PDAF_lnetf_update
! Calls: U_g2l_obs
! Calls: U_init_obs_l
! Calls: PDAF_timeit
! Calls: PDAF_memcount
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: syevTYPE (LAPACK; dsyev or ssyev dependent on precision)
!EOP
       
! *** local variables ***
  INTEGER :: i, j, member, col, row    ! Counters
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  INTEGER :: syev_info                 ! Status flag for SYEV
  INTEGER :: ldwork                    ! Size of work array for SYEV
  INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
  LOGICAL, SAVE :: screenout = .true.  ! Whether to print information to stdout
  REAL :: n_eff                        ! Effective sample size
  REAL :: fac                          ! Multiplication factor
  REAL :: weight                       ! Ensemble weight (likelihood)
  REAL, ALLOCATABLE :: ens_blk(:,:)    ! Temporary block of state ensemble
  REAL, ALLOCATABLE :: svals(:)        ! Singular values of Ainv
  REAL, ALLOCATABLE :: work(:)         ! Work array for SYEV
  REAL, ALLOCATABLE :: A(:,:)          ! Weight transform matrix
  REAL, ALLOCATABLE :: innov_i(:)      ! Innovation
  REAL, ALLOCATABLE :: T_tmp(:,:)      ! Temporary matrix
  REAL, ALLOCATABLE :: weights(:)      ! weight vector
  INTEGER , SAVE :: lastdomain = -1    ! save domain index
  REAL :: total_weight                 ! Sum of weights
  INTEGER :: screen2 = 0               ! Screen flag accounting for screenout
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
     IF (screenout) THEN
        WRITE (*,'(a, 5x, a, i5, a)') &
             'PDAF','--- Use OpenMP parallelization with ', nthreads, ' threads'
     END IF
#endif
  END IF

  IF (screen>0 .AND. screenout) THEN
     screen2 = screen
  ELSE
     screen2 = 0
  END IF

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lnetf_analysis -- START'

  CALL PDAF_timeit(51, 'old')


  ! **********************************************
  ! *** Compute particle weights as likelihood ***
  ! **********************************************

  ! Allocate weight vector
  ALLOCATE(weights(dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)

  ! Allocate temporal array for obs-ens_i
  ALLOCATE(innov_i(dim_obs_l))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l)

  CALL PDAF_timeit(16, 'new')
  CALL PDAF_timeit(20, 'new')

  ! Get residual as difference of observation and observed state for 
  ! each ensemble member only on domains where observations are availible

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, &
       'PDAF_lnetf_analysis -- call g2l_obs and likelihood_l', dim_ens, 'times'

  CALC_w: DO member = 1, dim_ens

     ! Calculate local residual  
     innov_i = obs_l - HX_l(:,member)

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug: ', debug, &
             'PDAF_lnetf_analysis -- member', member
        WRITE (*,*) '++ PDAF-debug PDAF_lnetf_analysis:', debug, '  innovation d_l', innov_i
        WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lnetf_analysis -- call likelihood_l'
     end IF

     ! Compute likelihood
     CALL PDAF_timeit(48, 'new')
     CALL U_likelihood_l(domain_p, step, dim_obs_l, obs_l, innov_i, weight)
     CALL PDAF_timeit(48, 'old')
     weights(member) = weight

  END DO CALC_w

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lnetf_analysis:', debug, '  raw weights', weights

  ! Compute inflation of weights according to N_eff/N>limit_winf
  IF (type_winf == 1) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, &
             'PDAF_lnetf_analysis -- inflate weights '
     CALL PDAF_inflate_weights(screen2, dim_ens, limit_winf, weights)
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lnetf_analysis:', debug, '  inflated weights', weights
  END IF

  CALL PDAF_timeit(51, 'new')

  ! Normalize weights
  total_weight = 0.0
  DO i = 1, dim_ens
     total_weight = total_weight + weights(i)
  END DO

  IF (total_weight /= 0.0) THEN
     weights = weights / total_weight

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lnetf_analysis:', debug, '  normalized weights', weights
  ELSE
     ! ERROR: weights are zero
     WRITE(*,'(/5x,a/)') 'WARNING: Zero weights - reset to 1/dim_ens'
     weights = 1.0 / REAL(dim_ens)
  END IF

  DEALLOCATE(innov_i)

  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(20, 'old')


  ! ****************************************
  ! *** Calculate the transform matrix   ***
  ! ***      A= (diag(w)-w*w^t)          ***
  ! *** with the weights w               ***
  ! ****************************************

  CALL PDAF_timeit(51, 'new')
  CALL PDAF_timeit(21, 'new')

  ALLOCATE(A(dim_ens,dim_ens)) 
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens*dim_ens)

  DO j = 1, dim_ens
     DO i = 1, dim_ens
        A(i,j) = -weights(i) * weights(j)
     ENDDO
  ENDDO
  DO i = 1, dim_ens
     A(i,i) = A(i,i) + weights(i)
  END DO

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lnetf_analysis:', debug, '  A_l', A

  CALL PDAF_timeit(21, 'old')

  ! Compute effective ensemble size
  CALL PDAF_diag_effsample(dim_ens, weights, n_eff)
  eff_dimens(1) = n_eff
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lnetf_analysis:', debug, '  effective sample size', n_eff

  CALL PDAF_timeit(16, 'old')


! ***************************************
! *** Calculate square root of matrix ***
! ***************************************

  CALL PDAF_timeit(17, 'new')

  ! Compute symmetric square-root by SVD
  ALLOCATE(T_tmp(dim_ens,dim_ens))
  ALLOCATE(svals(dim_ens))
  ALLOCATE(work(3*dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3*dim_ens * dim_ens*dim_ens)
  ldwork = 3*dim_ens
  flag = 0

  CALL PDAF_timeit(22, 'new')

  !EVD
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lnetf_analysis:', debug, &
       '  Compute eigenvalue decomposition of A_l'

  CALL syevTYPE('v', 'l', dim_ens, A, dim_ens, svals, work, ldwork, syev_info)

  CALL PDAF_timeit(22, 'old')

  IF (syev_info == 0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lnetf_analysis:', debug, '  eigenvalues', svals
  ELSE
     WRITE(*,'(/5x,a,i7/)') 'Problem computing svd of W-ww^T in domain', domain_p
  ENDIF
 
  ! Check for small eigenvalues
!$OMP CRITICAL
  DO i = 1,dim_ens
     IF (svals(i) < 1.0E-15) THEN
        svals(i) = 0.0
        cnt_small_svals = cnt_small_svals + 1
     ENDIF
  ENDDO
  ! subtract one, because A is rank dim_ens-1
  cnt_small_svals = cnt_small_svals - 1
!$OMP END CRITICAL

  CALL PDAF_timeit(23,'new')  

  ! Ensure to only use positive singular values - negative ones are numerical error
  DO i = 1, dim_ens
     IF (svals(i)>0.0) THEN
        svals(i) = SQRT(svals(i))
     ELSE
        svals(i) = 0.0
     END IF
  END DO

  DO j = 1,dim_ens
     DO i = 1, dim_ens
        T(j,i) = A(j,i) * svals(i)
     END DO
  END DO

  DEALLOCATE(svals, work)

  CALL PDAF_timeit(32,'old')

  CALL PDAF_timeit(33,'new')
  ! Multiply with singular vectors
  CALL gemmTYPE('n', 't', dim_ens, dim_ens, dim_ens, 1.0, &
       T, dim_ens, A, dim_ens, 0.0, T_tmp, dim_ens)
  CALL PDAF_timeit(33,'old')

  ! Multiply T by sqrt(m) to get unbiased ensemble 
  fac = SQRT(REAL(dim_ens))

  CALL PDAF_timeit(34, 'new') 
  IF (type_forget==2 .OR. type_forget==3) fac = fac / SQRT(forget) !analysis inflation
  CALL PDAF_timeit(34, 'old') 

  CALL PDAF_timeit(35,'new')

  ! Multiply T with random matrix and the factor 
  CALL gemmTYPE('n', 'n', dim_ens, dim_ens, dim_ens, &
       fac, T_tmp, dim_ens, rndmat, dim_ens, &
       0.0, T, dim_ens)
  CALL PDAF_timeit(35,'old')

  ! Part 3: W = sqrt(T) + w for efficient ensemble update
  DO col = 1, dim_ens
     DO row = 1, dim_ens
        T(row, col) = T(row, col) + weights(row)
     END DO
  END DO

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lnetf_analysis:', debug, '  transform', T

  DEALLOCATE(weights, A, T_tmp)

  CALL PDAF_timeit(23, 'old')
  CALL PDAF_timeit(17, 'old')


  ! ************************************************
  ! ***     Transform state ensemble             ***
  ! ***              a    f                      ***
  ! ***             X  = X  W                    ***
  ! *** The weight matrix W is stored in T       ***
  ! ************************************************

  CALL PDAF_timeit(18, 'new')

  ! *** Perform ensemble transformation ***
  ! Use block formulation for transformation
  maxblksize = 200
  IF (mype == 0 .AND. screen2 > 0) &
       WRITE (*, '(a, 5x, a, i5)') &
       'PDAF', '--- use blocking with size ', maxblksize

  ALLOCATE(ens_blk(maxblksize, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', maxblksize * dim_ens)

  blocking: DO blklower = 1, dim_l, maxblksize
           
     blkupper = MIN(blklower + maxblksize - 1, dim_l)

     ! Store forecast ensemble
     DO col = 1, dim_ens
        ens_blk(1 : blkupper - blklower + 1, col) &
             = ens_l(blklower : blkupper, col)
     END DO

     !                        a    f
     ! Transform ensemble:   X =  X  W
     CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
          1.0, ens_blk, maxblksize, T, dim_ens, &
          0.0, ens_l(blklower:blkupper, 1), dim_l)

  END DO blocking

  CALL PDAF_timeit(18, 'old')

  DEALLOCATE(ens_blk)

  CALL PDAF_timeit(51, 'old')


! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1

  lastdomain = domain_p

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lnetf_analysis -- END'

END SUBROUTINE PDAF_lnetf_analysis 
