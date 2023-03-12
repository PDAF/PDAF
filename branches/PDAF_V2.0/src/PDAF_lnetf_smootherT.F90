! Copyright (c) 2014-2021 Paul Kirchgessner
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
! !ROUTINE: PDAF_lnetf_smootherT --- compute ensemble transform of LNETS for smoothing
!
! !INTERFACE:
SUBROUTINE PDAF_lnetf_smootherT(domain_p, step, dim_obs_f, dim_obs_l, &
     dim_ens, HX_f, rndmat, U_g2l_obs, U_init_obs_l, U_likelihood_l, &
     screen, T, flag)

! !DESCRIPTION:
! Computation of transform matrix for smoother extension of
! NETF, which is computed without inflation.
!
! Variant for domain decomposed states.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2016-11 - Lars Nerger - Initial code based on LNETF_analysis
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
  INTEGER, INTENT(in) :: dim_obs_f   ! PE-local dimension of full observation vector
  INTEGER, INTENT(in) :: dim_obs_l   ! Size of obs. vector on local ana. domain
  INTEGER, INTENT(in) :: dim_ens     ! Size of ensemble 
  REAL, INTENT(in) :: HX_f(dim_obs_f, dim_ens) ! PE-local full observed state ens.
    ! Are equal for subtype==0, may differ for subtype==1. Are stored as -log(wi)!
  REAL, INTENT(in) :: rndmat(dim_ens, dim_ens) ! Global random rotation matrix
  INTEGER, INTENT(in) :: screen      ! Verbosity flag
  REAL, INTENT(inout) :: T(dim_ens, dim_ens)  ! local ensemble transformation matrix
  INTEGER, INTENT(inout) :: flag     ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_g2l_obs, &   ! Restrict full obs. vector to local analysis domain
       U_init_obs_l, &       ! Init. observation vector on local analysis domain
       U_likelihood_l        ! Compute observation likelihood for an ensemble member

! !CALLING SEQUENCE:
! Called by: PDAF_letkf_update
! Calls: U_g2l_obs
! Calls: U_init_obs_l
! Calls: PDAF_memcount
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: syevTYPE (LAPACK; dsyev or ssyev dependent on precision)
!EOP
       
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
  REAL, ALLOCATABLE :: svals(:)      ! Singular values of Uinv
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
