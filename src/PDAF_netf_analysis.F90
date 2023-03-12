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
! !ROUTINE: PDAF_netf_analysis --- NETF analysis cf. Toedter & Ahrens (2015)
!
! !INTERFACE:
SUBROUTINE PDAF_netf_analysis(step, dim_p, dim_obs_p, dim_ens, &
     state_p, ens_p, rndmat, T, type_forget, forget, &
     type_winf, limit_winf, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_likelihood, &
     screen, flag)

! !DESCRIPTION:
! Analysis step of the NETF following Toedter and Ahrens (2015)
! A Second-order Exact Ensemble Square Root Filter for Nonlinear 
! Data Assimlation. 
! The update is slightly modified to avoid computing the forecast
! mean state.
!
! Variant for domain decomposed states. 
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2014-03 - Paul Kirchgessner Changed original ETKF code to NETF
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
  USE PDAF_mod_filter, &
       ONLY: obs_member

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step         ! Current time step
  INTEGER, INTENT(in) :: dim_p        ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p   ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens      ! Size of ensemble
  REAL, INTENT(out)   :: state_p(dim_p)          ! on exit: PE-local forecast state
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  REAL, INTENT(in)    :: rndmat(dim_ens, dim_ens) ! Orthogonal random matrix
  REAL, INTENT(inout) :: T(dim_ens, dim_ens)     ! Ensemble transform matrix
  INTEGER, INTENT(in) :: type_forget  ! Type of forgetting factor
  REAL, INTENT(in)    :: forget       ! Forgetting factor
  INTEGER, INTENT(in) :: type_winf    ! Type of weights inflation
  REAL, INTENT(in) :: limit_winf      ! Limit for weights inflation
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
! Calls: PDAF_timeit
! Calls: PDAF_memcount
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: syevTYPE (LAPACK; dsyev or ssyev dependent on precision)
!EOP
       
! *** local variables ***
  INTEGER :: i, j, member, col, row   ! counters
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  INTEGER :: ldwork                   ! Size of work array for SYEV
  INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
  INTEGER :: syev_info                ! Status flag for Lapack
  REAL :: fac                         ! Multiplication factor
  REAL :: effN                        ! Efective sample size
  REAL :: weight                      ! Ensemble weight (likelihood)
  INTEGER :: n_small_svals            ! Number of small eigenvalues
  REAL, ALLOCATABLE :: resid_i(:)     ! PE-local observation residual
  REAL, ALLOCATABLE :: obs_p(:)       ! PE-local observation vector
  REAL, ALLOCATABLE :: ens_blk(:,:)   ! Temporary block of state ensemble
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

  CALL PDAF_timeit(51, 'new')

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a)') &
          'PDAF', 'Compute NETF filter update'
  END IF


! ************************
! *** Inflate ensemble ***
! ************************

  IF (type_forget==0 ) THEN
     CALL PDAF_timeit(34, 'new') ! Apply forgetting factor
     CALL PDAF_inflate_ens(dim_p, dim_ens, state_p, ens_p, forget)
     CALL PDAF_timeit(34, 'old')
  ENDIF
  CALL PDAF_timeit(51, 'old')


! *********************************
! *** Get observation dimension ***
! *********************************

  CALL PDAF_timeit(15, 'new')
  CALL U_init_dim_obs(step, dim_obs_p)
  CALL PDAF_timeit(15, 'old')

  IF (screen > 2) THEN
     WRITE (*, '(a, 5x, a13, 1x, i3, 1x, a, i8)') &
          'PDAF', '--- PE-domain', mype, 'dimension of observation vector', dim_obs_p
  END IF


  ! ***********************************************
  ! *** Compute particle weights (=likelihood)  ***
  ! ***********************************************

  CALL PDAF_timeit(12, 'new')

  ! Allocate weights
  ALLOCATE(weights(dim_ens))   
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)

  haveobs: IF (dim_obs_p > 0) THEN
     ! *** The weights only exist for domains with observations ***

     ALLOCATE(obs_p(dim_obs_p))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p)

     ! Allocate tempory arrays for obs-ens_i
     ALLOCATE(resid_i(dim_obs_p))
     ALLOCATE(Rinvresid(dim_obs_p))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2*dim_obs_p)

     ! Get residual as difference of observation and observed state for each ensemble member
     CALC_w: DO member = 1, dim_ens

        ! Store member index
        obs_member = member

        CALL PDAF_timeit(44, 'new')
        CALL U_obs_op(step, dim_p, dim_obs_p, ens_p(:, member), resid_i)
        CALL PDAF_timeit(44, 'old')

        IF (member==1) THEN
           ! get observation vector (has to be after U_obs_op for OMI)
           CALL PDAF_timeit(50, 'new')
           CALL U_init_obs(step, dim_obs_p, obs_p)
           CALL PDAF_timeit(50, 'old')
        END IF

        CALL PDAF_timeit(51, 'new')
        resid_i = obs_p - resid_i 
        CALL PDAF_timeit(51, 'old')

        ! Compute likelihood
        CALL PDAF_timeit(47, 'new')
        CALL U_likelihood(step, dim_obs_p, obs_p, resid_i, weight)
        CALL PDAF_timeit(47, 'old')
        weights(member) = weight

     END DO CALC_w

     ! Compute inflation of weights according to N_eff
     IF (type_winf == 1) THEN
        CALL PDAF_inflate_weights(screen, dim_ens, limit_winf, weights)
     END IF

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
        flag = 1
        WRITE(*,'(/5x,a/)') 'PDAF-ERROR (1): Zero weights in NETF analysis step'
     END IF

     DEALLOCATE(obs_p, resid_i, Rinvresid)

     ! Diagnostic: Compute effective sample size
     CALL PDAF_diag_effsample(dim_ens, weights, effN)
     IF (mype == 0 .AND. screen > 0) &
          WRITE (*, '(a, 5x, a, f10.2)') 'PDAF', '--- Effective sample size ', effN

     CALL PDAF_timeit(51, 'old')

  ELSE
     ! Without observations, all ensemble member have the same weight

     CALL PDAF_timeit(51, 'new')
     weights = 1/dim_ens
     CALL PDAF_timeit(51, 'old')

  END IF haveobs

  CALL PDAF_timeit(12, 'old')
  CALL PDAF_timeit(51, 'new')


  ! ****************************************
  ! *** Calculate the transform matrix   ***
  ! ***      A= (diag(w)-w*w^t)          ***
  ! *** with the weights w               ***
  ! ****************************************

  CALL PDAF_timeit(10, 'new')

  ALLOCATE(A(dim_ens,dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens*dim_ens)

  DO i = 1, dim_ens
     DO j = 1, dim_ens
        A(i,j) = -weights(i) * weights(j)
     ENDDO
  ENDDO
  DO i = 1, dim_ens
     A(i,i) = A(i,i) + weights(i)
  END DO

  CALL PDAF_timeit(10, 'old')


  ! ********************************************************************
  ! *** Compute ensemble transformation matrix W as square-root of A ***
  ! ********************************************************************

  CALL PDAF_timeit(13, 'new')

  ! Compute symmetric square-root of A by EVD
  ALLOCATE(svals(dim_ens))
  ALLOCATE(work(3*dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 4*dim_ens + dim_ens*dim_ens)
  ldwork = 3*dim_ens

  CALL syevTYPE('v', 'l', dim_ens, A, dim_ens, svals, work, ldwork, syev_info)

  IF (syev_info /= 0 ) THEN
      WRITE(*,'(/5x,a/)') 'PDAF-ERROR(2): Problem in computing the SVD of W-ww^T'
      flag = 2
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
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens*dim_ens)

  ! Calculate transform matrix T
  CALL gemmTYPE('n', 't', dim_ens, dim_ens, dim_ens, 1.0, &
       T, dim_ens, A, dim_ens, 0.0, T_tmp, dim_ens)

  ! Multiply T by m/(m-1) to get unbiased ensemble
  fac = SQRT(REAL(dim_ens))

  IF (type_forget==2) fac = fac / SQRT(forget) !analysis inflation

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

  CALL PDAF_timeit(13, 'old')


! *******************************************
! ***     Transform state ensemble        ***
! ***              a    f                 ***
! ***             X  = X  W               ***
! *** The weight matrix W is stored in T. ***
! *******************************************

  ! Use block formulation for transformation
  maxblksize = 200
  IF (mype == 0 .AND. screen > 0) &
       WRITE (*, '(a, 5x, a, i5)') &
       'PDAF', '--- use blocking with size ', maxblksize
        
  ALLOCATE(ens_blk(maxblksize, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', maxblksize * dim_ens)

  blocking: DO blklower = 1, dim_p, maxblksize
           
     blkupper = MIN(blklower + maxblksize - 1, dim_p)

     ! Store forecast ensemble
     CALL PDAF_timeit(21, 'new')
     DO col = 1, dim_ens
        ens_blk(1 : blkupper - blklower + 1, col) &
             = ens_p(blklower : blkupper, col)
     END DO

     CALL PDAF_timeit(21, 'old')

     !                        a  _f
     ! Transform ensemble:   X = X  W
     CALL PDAF_timeit(22, 'new')

     CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
          1.0, ens_blk, maxblksize, T, dim_ens, &
          0.0, ens_p(blklower:blkupper, 1), dim_p)

     CALL PDAF_timeit(22, 'old')

  END DO blocking

  DEALLOCATE(ens_blk)

  CALL PDAF_timeit(51, 'old')


! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_netf_analysis
