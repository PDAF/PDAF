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
!$Id$
!BOP
!
! !ROUTINE: PDAF_netf_smmotherT --- compute ensemble transform of NETS for smoothing
!
! !INTERFACE:
SUBROUTINE PDAF_netf_smootherT(step, dim_p, dim_obs_p, dim_ens, &
     ens_p, rndmat, T,  &
     U_init_dim_obs, U_obs_op, U_init_obs, U_likelihood, &
     screen, flag)

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
! 2016-11 - Lars Nerger - Initial code based on NETF_analysis
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_mod_filtermpi, &
       ONLY: mype

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step         ! Current time step
  INTEGER, INTENT(in) :: dim_p        ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p   ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens      ! Size of ensemble
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)    ! PE-local state ensemble
  REAL, INTENT(in)    :: rndmat(dim_ens, dim_ens) ! Orthogonal random matrix
  REAL, INTENT(inout) :: T(dim_ens, dim_ens)      ! Ensemble transform matrix
  INTEGER, INTENT(in) :: screen       ! Verbosity flag
  INTEGER, INTENT(inout) :: flag      ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, & ! Initialize dimension of observation vector
       U_obs_op, &              ! Observation operator
       U_init_obs, &            ! Initialize observation vector
       U_likelihood             ! Compute observation likelihood for an ensemble member

! !CALLING SEQUENCE:
! Called by: PDAF_netf_update
! Calls: U_init_dim_obs
! Calls: U_obs_op
! Calls: U_init_obs
! Calls: U_likelihood
! Calls: PDAF_memcount
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: syevTYPE (LAPACK; dsyev or ssyev dependent on precision)
!EOP
       
! *** local variables ***
  INTEGER :: i, j, member, col, row   ! counters
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  INTEGER :: syev_info                ! Status flag for SYEV
  INTEGER :: ldwork                   ! Size of work array for SYEV
  REAL :: fac                         ! Multiplication factor
  REAL :: effN                        ! Effective sample size
  REAL :: weight                      ! Ensemble weight (likelihood)
  INTEGER :: n_small_svals            ! Number of small eigenvalues
  REAL, ALLOCATABLE :: resid_i(:)     ! PE-local observation residual
  REAL, ALLOCATABLE :: obs_p(:)       ! PE-local observation vector
  REAL, ALLOCATABLE :: svals(:)       ! Singular values of Uinv
  REAL, ALLOCATABLE :: work(:)        ! Work array for SYEV
  REAL, ALLOCATABLE :: T_tmp(:,:)     ! Square root of transform matrix
  REAL, ALLOCATABLE :: A(:,:)         ! Full transform matrix
  REAL, ALLOCATABLE :: Rinvresid(:)   ! R^-1 times residual 
  REAL, ALLOCATABLE :: weights(:)     ! Weight vector
  REAL :: total_weight                ! Sum of weights


! **********************
! *** INITIALIZATION ***
! **********************

  CALL PDAF_timeit(16, 'new')

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a)') &
          'PDAF', 'NETF smoother: Compute transform matrix without inflation'
  END IF


! *********************************
! *** Get observation dimension ***
! *********************************

  CALL PDAF_timeit(15, 'new')
  CALL U_init_dim_obs(step, dim_obs_p)
  CALL PDAF_timeit(15, 'old')

  IF (screen > 0) THEN
     WRITE (*, '(a, 5x, a13, 1x, i3, 1x, a, i8)') &
          'PDAF', '--- PE-domain', mype, 'dimension of observation vector', dim_obs_p
  END IF


  ! ***********************************************
  ! *** Compute particle weights                ***
  ! ***   w_i = exp(-0.5*(y-Hx_i)^TR-1(y-Hx_i)) ***
  ! ***********************************************

  ! Allocate weights
  ALLOCATE(weights(dim_ens))   

  haveobs: IF (dim_obs_p > 0) THEN
     ! *** The weights only exist for domains with observations ***

     ALLOCATE(obs_p(dim_obs_p))

     ! get observation vector
     CALL PDAF_timeit(50, 'new')
     CALL U_init_obs(step, dim_obs_p, obs_p)
     CALL PDAF_timeit(50, 'old')

     ! Allocate tempory arrays for obs-ens_i
     ALLOCATE(resid_i(dim_obs_p))
     ALLOCATE(Rinvresid(dim_obs_p))

     ! Get residual as difference of observation and observed state for each ensemble member
     CALC_w: DO member = 1, dim_ens

        CALL PDAF_timeit(44, 'new')
        CALL U_obs_op(step, dim_p, dim_obs_p, ens_p(:, member), resid_i)
        CALL PDAF_timeit(44, 'old')

        CALL PDAF_timeit(51, 'new')
        resid_i = obs_p - resid_i 
        CALL PDAF_timeit(51, 'old')

        ! Compute likelihood
        CALL PDAF_timeit(47, 'new')
        CALL U_likelihood(step, dim_obs_p, obs_p, resid_i, weight)
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
        ! Normalize weights
        weights = weights / total_weight
     ELSE
        ! ERROR: weights are zero
        flag = 3
        WRITE(*,'(/5x,a/)') 'PDAF-ERROR (3): Zero weights in NETF smoother'
     END IF

     DEALLOCATE(obs_p, resid_i, Rinvresid)

     ! Diagnostic: Compute effective sample size
     CALL PDAF_diag_effsample(dim_ens, weights, effN)
     IF (mype == 0 .AND. screen > 0) &
          WRITE (*, '(a, 5x, a, f10.2)') &
          'PDAF', '--- Effective sample size ', effN

     CALL PDAF_timeit(51, 'old')

  ELSE
     ! Without observations, al ensemble member have the same weight

     weights = 1/dim_ens
     
  END IF haveobs


  CALL PDAF_timeit(51, 'new')

  ! ****************************************
  ! *** Calculate the transform matrix   ***
  ! ***      A= (diag(w)-w*w^t)          ***
  ! *** with the weights w               ***
  ! ****************************************

  ALLOCATE(A(dim_ens,dim_ens))

  DO i = 1, dim_ens
     DO j = 1, dim_ens
        A(i,j) = -weights(i) * weights(j)
     ENDDO
  ENDDO
  DO i = 1, dim_ens
     A(i,i) = A(i,i) + weights(i)
  END DO


  ! ********************************************************************
  ! *** Compute ensemble transformation matrix W as square-root of A ***
  ! ********************************************************************

  ! Compute symmetric square-root of A by EVD
  ALLOCATE(svals(dim_ens))
  ALLOCATE(work(3*dim_ens))

  ldwork = 3*dim_ens

  CALL syevTYPE('v', 'l', dim_ens, A, dim_ens, svals, work, ldwork, syev_info)

  IF (syev_info /= 0 ) THEN
      WRITE(*,'(/5x,a/)') 'PDAF-ERROR(4): Problem in computing the SVD of W-ww^T'
      flag = 4
  END IF
  
  ! Check for too small eigenvalues
  n_small_svals = 0
  DO i = 1, dim_ens
     IF (svals(i)<1.0E-15) THEN
        svals(i) = 0.0
        n_small_svals = n_small_svals + 1
     END IF
  END DO
  ! subtract one, because A is rank dim_ens-1
  n_small_svals = n_small_svals - 1
  IF (mype == 0 .AND. screen > 0) &
       WRITE (*, '(a, 5x, a, i5)') &
       'PDAF', '--- number of small singular values ', n_small_svals

  DO j = 1,dim_ens
     DO i = 1, dim_ens
        T(j,i) = A(j,i) * SQRT(svals(i))
     END DO
  END DO

  DEALLOCATE(svals, work)

  ALLOCATE(T_tmp(dim_ens,dim_ens))

  ! Calculate transform matrix T
  CALL gemmTYPE('n', 't', dim_ens, dim_ens, dim_ens, 1.0, &
       T, dim_ens, A, dim_ens, 0.0, T_tmp, dim_ens)

  ! Multiply T by m/(m-1) to get unbiased ensemble
  fac = SQRT(REAL(dim_ens))

  ! Multiply random matrix with quare root of A (T)
  CALL gemmTYPE('n', 'n', dim_ens, dim_ens, dim_ens, &
         fac, T_tmp, dim_ens, rndmat, dim_ens, &
         0.0, T, dim_ens)

  ! Compute W = sqrt(U) + w for efficient ensemble update
  DO col = 1, dim_ens
     DO row = 1, dim_ens
        T(row, col) = T(row, col) + weights(row)
     END DO
  END DO

  DEALLOCATE(weights, A, T_tmp)


! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1
  
  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(16, 'old')

END SUBROUTINE PDAF_netf_smootherT
