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
!
!> local analysis step of LNETF
!!
!! Analysis step of the LNETF.
!!
!! Inflation has to be done BEFORE calling this routine
!! except for inflation of the posterior ensemble!
!!
!! Variant for domain decomposed states.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2014-05 - Paul Kirchgessner - Initial code based on LETKF
!! * Later revisions - see svn log
!!
MODULE PDAF_lnetf_analysis

CONTAINS
SUBROUTINE PDAF_lnetf_ana(domain_p, step, dim_l, dim_obs_l, &
     dim_ens, ens_l, HX_l, obs_l, rndmat, &
     U_likelihood_l, type_forget, forget, &
     type_winf, limit_winf, cnt_small_svals, eff_dimens, T, &
     screen, debug, flag)

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

! *** Argumeents ***
!  Variable naming scheme:
!    suffix _p: Denotes a full variable on the PE-local domain
!    suffix _l: Denotes a local variable on the current analysis domain
  INTEGER, INTENT(in) :: domain_p               !< Current local analysis domain
  INTEGER, INTENT(in) :: step                   !< Current time step
  INTEGER, INTENT(in) :: dim_l                  !< State dimension on local analysis domain
  INTEGER, INTENT(in) :: dim_obs_l              !< Size of obs. vector on local ana. domain
  INTEGER, INTENT(in) :: dim_ens                !< Size of ensemble 
  REAL, INTENT(inout) :: ens_l(dim_l, dim_ens)  !< Local state ensemble
  REAL, INTENT(in) :: HX_l(dim_obs_l, dim_ens)  !< Local observed state ensemble (perturbation)
  REAL, INTENT(in) :: obs_l(dim_obs_l)          !< Local observation vector
  REAL, INTENT(in) :: rndmat(dim_ens, dim_ens)  !< Global random rotation matrix
  INTEGER, INTENT(in) :: screen                 !< Verbosity flag
  INTEGER, INTENT(in) :: type_forget            !< Typ eof forgetting factor
  REAL, INTENT(in) :: forget                    !< Forgetting factor
  INTEGER, INTENT(in) :: type_winf              !< Type of weights inflation
  REAL, INTENT(in) :: limit_winf                !< Limit for weights inflation
  INTEGER, INTENT(inout) :: cnt_small_svals     !< Number of small eigen values
  REAL, INTENT(inout) :: eff_dimens(1)          !< Effective ensemble size
  REAL, INTENT(inout) :: T(dim_ens, dim_ens)    !< local ensemble transformation matrix
  INTEGER, INTENT(in) :: debug                  !< Flag for writing debug output
  INTEGER, INTENT(inout) :: flag                !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_likelihood_l           !< Compute observation likelihood for an ensemble member
       
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

END SUBROUTINE PDAF_lnetf_ana

!-------------------------------------------------------------------------------
!> Compute ensemble transform of LNETS for smoothing
!!
!! Computation of transform matrix for smoother extension of
!! NETF, which is computed without inflation.
!!
!! Variant for domain decomposed states.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! 2016-11 - Lars Nerger - Initial code based on LNETF_analysis
!! Later revisions - see svn log
!!
SUBROUTINE PDAF_lnetf_smootherT(domain_p, step, dim_obs_f, dim_obs_l, &
     dim_ens, HX_f, rndmat, U_g2l_obs, U_init_obs_l, U_likelihood_l, &
     screen, T, flag)

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

! *** Arguments ***
!  Variable naming scheme:
!    suffix _p: Denotes a full variable on the PE-local domain
!    suffix _l: Denotes a local variable on the current analysis domain
  INTEGER, INTENT(in) :: domain_p    !< Current local analysis domain
  INTEGER, INTENT(in) :: step        !< Current time step
  INTEGER, INTENT(in) :: dim_obs_f   !< PE-local dimension of full observation vector
  INTEGER, INTENT(in) :: dim_obs_l   !< Size of obs. vector on local ana. domain
  INTEGER, INTENT(in) :: dim_ens     !< Size of ensemble 
  REAL, INTENT(in) :: HX_f(dim_obs_f, dim_ens) !< PE-local full observed state ens.
    !< Are equal for subtype==0, may differ for subtype==1. Are stored as -log(wi)
  REAL, INTENT(in) :: rndmat(dim_ens, dim_ens) !< Global random rotation matrix
  INTEGER, INTENT(in) :: screen      !< Verbosity flag
  REAL, INTENT(inout) :: T(dim_ens, dim_ens)  !< local ensemble transformation matrix
  INTEGER, INTENT(inout) :: flag     !< Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_g2l_obs, &   !< Restrict full obs. vector to local analysis domain
       U_init_obs_l, &       !< Init. observation vector on local analysis domain
       U_likelihood_l        !< Compute observation likelihood for an ensemble member

       
! *** local variables ***
  INTEGER :: i, j, member, col, row  ! Counters
  INTEGER, SAVE :: allocflag = 0     ! Flag whether first time allocation is done
  INTEGER, SAVE :: first = 1         ! Flag for very first call to routine
  INTEGER, SAVE :: domain_save = 1   ! Index of domain from last call to routine
  INTEGER :: syev_info               ! Status flag for SYEV
  INTEGER :: ldwork                  ! Size of work array for SYEV
  REAL :: fac                        ! Multiplication factor
  REAL :: weight                     ! Ensemble weight (likelihood)
  REAL, ALLOCATABLE :: obs_l(:)      ! local observation vector
  REAL, ALLOCATABLE :: svals(:)      ! Singular values of Ainv
  REAL, ALLOCATABLE :: work(:)       ! Work array for SYEV
  REAL, ALLOCATABLE :: A(:,:)        ! Weight transform matrix
  REAL, ALLOCATABLE :: resid_i(:)    ! Residual
  REAL, ALLOCATABLE :: Rinvresid(:)  ! R^-1 times residual
  REAL, ALLOCATABLE :: T_tmp(:,:)    ! Temporary matrix
  REAL, ALLOCATABLE :: weights(:)    ! weight vector
  REAL :: total_weight               ! Sum of weights
  INTEGER, SAVE :: mythread, nthreads  ! Thread variables for OpenMP

!$OMP THREADPRIVATE(mythread, nthreads, allocflag, first, domain_save)


! **********************
! *** INITIALIZATION ***
! **********************

  CALL PDAF_timeit(51, 'new')

#if defined (_OPENMP)
  nthreads = omp_get_num_threads()
  mythread = omp_get_thread_num()
#else
  nthreads = 1
  mythread = 0
#endif

  IF ((domain_p <= domain_save) .OR. (first == 1)) THEN
     IF (mype == 0 .AND. screen > 0.AND. mythread==0) THEN
        WRITE (*, '(a, 5x, a)') &
             'PDAF', 'Compute transform matrix for smoother without inflation'
     END IF
  END IF


  ! **********************************************
  ! *** Compute particle weights as likelihood ***
  ! **********************************************

  ! Allocate weight vector
  ALLOCATE(weights(dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)

  ! Allocate temporal array for obs-ens_i
  ALLOCATE(resid_i(dim_obs_l))
  ALLOCATE(Rinvresid(dim_obs_l))
  ALLOCATE(obs_l(dim_obs_l))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3*dim_obs_l)

  CALL PDAF_timeit(51, 'old')
     
  !get local observation vector
  CALL PDAF_timeit(47, 'new')
  CALL U_init_obs_l(domain_p, step, dim_obs_l, obs_l)
  CALL PDAF_timeit(47, 'old')

  ! Get residual as difference of observation and observed state for 
  ! each ensemble member only on domains where observations are availible

  CALC_w: DO member = 1,dim_ens

     ! Restrict global state to local state
     CALL PDAF_timeit(46, 'new')
     CALL U_g2l_obs(domain_p, step, dim_obs_f, dim_obs_l, HX_f(:,member), resid_i)
     CALL PDAF_timeit(46, 'old')

     ! Calculate local residual  
     CALL PDAF_timeit(51, 'new')
     resid_i = obs_l - resid_i
     CALL PDAF_timeit(51, 'old')

     ! Compute likelihood
     CALL PDAF_timeit(47, 'new')
     CALL U_likelihood_l(domain_p, step, dim_obs_l, obs_l, resid_i, weight)
     CALL PDAF_timeit(47, 'old')
     weights(member) = weight

  END DO CALC_w

  CALL PDAF_timeit(51, 'new')

  ! Normalize weights
  total_weight = 0.0
  DO i = 1, dim_ens
     total_weight = total_weight + weights(i)
  END DO
  IF (total_weight /= 0.0) THEN
     weights = weights / total_weight
  ELSE
     ! ERROR: weights are zero
     flag = 3
     WRITE(*,'(/5x,a/)') 'PDAF-ERROR (3): Zero weights in LNETF smoother'
  END IF

  DEALLOCATE(obs_l, resid_i, Rinvresid)


  ! ****************************************
  ! *** Calculate the transform matrix   ***
  ! ***      A= (diag(w)-w*w^t)          ***
  ! *** with the weights w               ***
  ! ****************************************

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


! ***************************************
! *** Calculate square root of Matrix ***
! ***************************************

  ! Compute symmetric square-root by SVD
  ALLOCATE(T_tmp(dim_ens,dim_ens))
  ALLOCATE(svals(dim_ens))
  ALLOCATE(work(3*dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3*dim_ens * dim_ens*dim_ens)
  ldwork = 3*dim_ens
  flag = 0

  !EVD
  CALL syevTYPE('v', 'l', dim_ens, A, dim_ens, svals, work, ldwork, syev_info)

  IF (syev_info /= 0 ) THEN
     WRITE(*,'(/5x,a,i7/)') 'Problem computing svd of W-ww^T in domain', domain_p
  ENDIF
   
  ! Check for small eigenvalues
  DO i = 1,dim_ens
     IF (svals(i) < 1.0E-15) THEN
        svals(i) = 0.0
     ENDIF
  ENDDO

  DO j = 1,dim_ens
     DO i = 1, dim_ens
        T(j,i) = A(j,i) * SQRT(svals(i))
     END DO
  END DO

  DEALLOCATE(svals, work)

  ! Multiply with singular vectors
  CALL gemmTYPE('n', 't', dim_ens, dim_ens, dim_ens, 1.0, &
       T, dim_ens, A, dim_ens, 0.0, T_tmp, dim_ens)

  ! Multiply T by sqrt(m) to get unbiased ensemble 
  fac = SQRT(REAL(dim_ens))

  ! Multiply T with random matrix and the factor 
  CALL gemmTYPE('n', 'n', dim_ens, dim_ens, dim_ens, &
       fac, T_tmp, dim_ens, rndmat, dim_ens, &
       0.0, T, dim_ens)

  ! Part 3: W = sqrt(T) + w for efficient ensemble update
  DO col = 1, dim_ens
     DO row = 1, dim_ens
        T(row, col) = T(row, col) + weights(row)
     END DO
  END DO

  DEALLOCATE(weights, A, T_tmp)

  CALL PDAF_timeit(51, 'old')


! ********************
! *** Finishing up ***
! ********************
  
  ! Increment maxlag counter
  IF ((domain_p <= domain_save) .OR. (first == 1)) THEN
     ! Set flag
     first = 0
  END IF
  domain_save = domain_p

  ! Set flag for memory counting
  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_lnetf_smootherT


!-------------------------------------------------------------------------------
!> Smoother extension for local NETF
!!
!! Smoother extension for the ensemble square-root filters (ETKF, ESTKF). 
!! The routine uses the matrix Ainv computed by the filter analysis
!! to perform the smoothing on past ensembles.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2012-05 - Lars Nerger - Initial code
!! * Later revisions - see svn log
!!
SUBROUTINE PDAF_smoother_lnetf(domain_p, step, dim_p, dim_l, dim_ens, &
     dim_lag, Ainv, ens_l, sens_p, cnt_maxlag, &
     U_g2l_state, U_l2g_state, screen)

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

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p      !< Current local analysis domain
  INTEGER, INTENT(in) :: step          !< Current time step
  INTEGER, INTENT(in) :: dim_p         !< PE-local dimension of model state
  INTEGER, INTENT(in) :: dim_l         !< State dimension on local analysis domain
  INTEGER, INTENT(in) :: dim_ens       !< Size of ensemble
  INTEGER, INTENT(in) :: dim_lag       !< Number of past time instances for smoother
  REAL, INTENT(in)   :: Ainv(dim_ens, dim_ens)  !< Weight matrix for ensemble transformation
  REAL, INTENT(inout) :: ens_l(dim_l, dim_ens)  !< local past ensemble (temporary)
  REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag)   !< PE-local smoother ensemble
  INTEGER, INTENT(inout) :: cnt_maxlag !< Count available number of time steps for smoothing
  INTEGER, INTENT(in) :: screen        !< Verbosity flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_g2l_state, & !< Get state on local ana. domain from global state
       U_l2g_state           !< Init full state from state on local analysis domain

! *** local variables ***
  INTEGER :: member, col, row, lagcol  ! Counters
  INTEGER :: n_lags                    ! Available number of time instances for smoothing
  INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  INTEGER, SAVE :: first = 1           ! Flag for very first call to routine
  INTEGER, SAVE :: domain_save = 1     ! Index of domain from last call to routine
  REAL :: invdimens                    ! Inverse of global ensemble size
  REAL, ALLOCATABLE :: ens_blk(:,:)    ! Temporary block of state ensemble
  REAL, ALLOCATABLE :: W_smooth(:,:)   ! Weight matrix for smoothing
  INTEGER, SAVE :: mythread, nthreads  ! Thread variables for OpenMP

!$OMP THREADPRIVATE(mythread, nthreads, allocflag, first, domain_save)

  
! **********************
! *** INITIALIZATION ***
! **********************

#if defined (_OPENMP)
  nthreads = omp_get_num_threads()
  mythread = omp_get_thread_num()
#else
  nthreads = 1
  mythread = 0
#endif

  ! Determine number of time instances for smoothing
  IF (cnt_maxlag >= dim_lag) THEN
     ! Already performed enough analysis to smooth over full lag
     n_lags = dim_lag
  ELSE
     ! Not yet enough analysis steps to smoother over full lag
     n_lags = cnt_maxlag
  END IF

  IF ((domain_p <= domain_save) .OR. (first == 1)) THEN
     IF (mype == 0 .AND. screen > 0 .AND. n_lags > 0 .AND. mythread==0) THEN
        WRITE (*, '(a, 5x, a, i8)') 'PDAF', 'Perform smoothing up to lag ', n_lags
     END IF
  END IF


! **********************************************
! *** Compute transform matrix for smoothing ***
! ***                                        ***
! *** W_smooth = A_filter + 1_N              ***
! **********************************************

  havelag: IF (n_lags > 0) THEN

     invdimens = 1.0 / REAL(dim_ens)

     ALLOCATE(W_smooth(dim_ens, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

     DO col = 1, dim_ens
        DO row = 1, dim_ens
           W_smooth(row, col) = Ainv(row, col)
        END DO
     END DO


! **********************************************
! *** Perform smoothing                      ***
! *** Transform state ensemble at past times ***
! ***              a    f                    ***
! ***             X  = X  W_smooth           ***
! **********************************************

     ! Use block formulation for transformation
     maxblksize = 200
     ALLOCATE(ens_blk(maxblksize, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', maxblksize * dim_ens)
     lagcol=1

     ! *** Smooth for all available lags ***
     smoothing: DO lagcol = 1, n_lags

        ! *** Get local ensemble ***
        DO member = 1, dim_ens
           CALL U_g2l_state(step, domain_p, dim_p, sens_p(:, member, lagcol), dim_l, &
                ens_l(:, member))
        END DO

        ! Use block formulation for transformation
        blocking: DO blklower = 1, dim_l, maxblksize
           
           blkupper = MIN(blklower + maxblksize - 1, dim_l)

           ! Store former analysis ensemble
           DO member = 1, dim_ens
              ens_blk(1 : blkupper-blklower+1, member) &
                   = ens_l(blklower : blkupper, member)
           END DO

           !                        a   f
           ! Transform ensemble:   X = X  W_smooth
           CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
                1.0, ens_blk(1, 1), maxblksize, W_smooth, dim_ens, &
                0.0, ens_l(blklower, 1), dim_l)

        END DO blocking

        ! *** Initialize global ensemble ***
        DO member = 1, dim_ens
           CALL U_l2g_state(step, domain_p, dim_l, ens_l(:, member), dim_p, &
                sens_p(:, member, lagcol))
        END DO

     END DO smoothing
     
     DEALLOCATE(ens_blk, W_smooth)

  END IF havelag


! ********************
! *** Finishing up ***
! ********************
  
  ! Increment maxlag counter
  IF ((domain_p <= domain_save) .OR. (first == 1)) THEN

     cnt_maxlag = cnt_maxlag + 1
          
     ! Set flag
     first = 0
  END IF
  domain_save = domain_p

  ! Set flag for memory counting
  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_smoother_lnetf

END MODULE PDAF_lnetf_analysis
