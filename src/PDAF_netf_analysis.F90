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
!> NETF analysis cf. Toedter & Ahrens (2015)
!!
!! Analysis step of the NETF following Toedter and Ahrens (2015)
!! A Second-order Exact Ensemble Square Root Filter for Nonlinear 
!! Data Assimlation. 
!! The update is slightly modified to avoid computing the forecast
!! mean state.
!!
!! Variant for domain decomposed states. 
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2014-03 - Paul Kirchgessner Changed original ETKF code to NETF
!! * Other revisions - see repository log
!!
MODULE PDAF_netf_analysis

CONTAINS
SUBROUTINE PDAF_netf_ana(step, dim_p, dim_obs_p, dim_ens, &
     state_p, ens_p, rndmat, T, type_forget, forget, &
     type_winf, limit_winf, type_noise, noise_amp, &
     HZ_p, obs_p, U_likelihood, screen, debug, flag)

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
       ONLY: PDAF_add_particle_noise, PDAF_inflate_weights
  USE PDAF_diag, &
       ONLY: PDAF_diag_effsample

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step         !< Current time step
  INTEGER, INTENT(in) :: dim_p        !< PE-local dimension of model state
  INTEGER, INTENT(in) :: dim_obs_p    !< PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens      !< Size of ensemble
  REAL, INTENT(out)   :: state_p(dim_p)           !< PE-local forecast state
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)    !< PE-local state ensemble
  REAL, INTENT(in)    :: rndmat(dim_ens, dim_ens) !< Orthogonal random matrix
  REAL, INTENT(inout) :: T(dim_ens, dim_ens)      !< Ensemble transform matrix
  INTEGER, INTENT(in) :: type_forget  !< Type of forgetting factor
  REAL, INTENT(in)    :: forget       !< Forgetting factor
  INTEGER, INTENT(in) :: type_winf    !< Type of weights inflation
  REAL, INTENT(in) :: limit_winf      !< Limit for weights inflation
  INTEGER, INTENT(in) :: type_noise   !< Type of pertubing noise
  REAL, INTENT(in) :: noise_amp       !< Amplitude of noise
  REAL, INTENT(in) :: HZ_p(dim_obs_p, dim_ens)    !< Temporary matrices for analysis
  REAL, INTENT(in) :: obs_p(dim_obs_p)            !< PE-local observation vector
  INTEGER, INTENT(in) :: screen       !< Verbosity flag
  INTEGER, INTENT(in) :: debug        !< Flag for writing debug output
  INTEGER, INTENT(inout) :: flag      !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_likelihood            !< Compute observation likelihood for an ensemble member
       
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
  REAL, ALLOCATABLE :: innov_i(:)     ! PE-local observation innovation
  REAL, ALLOCATABLE :: ens_blk(:,:)   ! Temporary block of state ensemble
  REAL, ALLOCATABLE :: svals(:)       ! Singular values of Ainv
  REAL, ALLOCATABLE :: work(:)        ! Work array for SYEV
  REAL, ALLOCATABLE :: T_tmp(:,:)     ! Square root of transform matrix
  REAL, ALLOCATABLE :: A(:,:)         ! Full transform matrix
  REAL, ALLOCATABLE :: Rinvinnov(:)   ! R^-1 times innovation 
  REAL, ALLOCATABLE :: weights(:)     ! Weight vector
  REAL :: total_weight                ! Sum of weights


! **********************
! *** INITIALIZATION ***
! **********************

  CALL PDAF_timeit(51, 'new')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_netf_analysis -- START'

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a)') &
          'PDAF', 'Compute NETF filter update'
  END IF

  CALL PDAF_timeit(51, 'old')


  ! ***********************************************
  ! *** Compute particle weights (=likelihood)  ***
  ! ***********************************************

  CALL PDAF_timeit(10, 'new')

  ! Allocate weights
  ALLOCATE(weights(dim_ens))   
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)

  haveobs: IF (dim_obs_p > 0) THEN
     ! *** The weights only exist for domains with observations ***

     ! Allocate tempory arrays for obs-ens_i
     ALLOCATE(innov_i(dim_obs_p))
     ALLOCATE(Rinvinnov(dim_obs_p))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2*dim_obs_p)

     ! Get innovation as observation minus observed state for each ensemble member
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, &
          'PDAF_netf_analysis -- call likelihood', dim_ens, 'times'

     CALC_w: DO member = 1, dim_ens

        CALL PDAF_timeit(51, 'new')
        innov_i = obs_p - HZ_p(:, member) 
        CALL PDAF_timeit(51, 'old')

        IF (debug>0) THEN
           WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_netf_analysis -- member', member
           WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_netf_analysis -- innovation d', innov_i
           WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_netf_analysis -- call likelihood'
        end IF

        ! Compute likelihood
        CALL PDAF_timeit(48, 'new')
        CALL U_likelihood(step, dim_obs_p, obs_p, innov_i, weight)
        CALL PDAF_timeit(48, 'old')
        weights(member) = weight

     END DO CALC_w

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_netf_analysis:', debug, '  raw weights', weights

     ! Compute inflation of weights according to N_eff
     IF (type_winf == 1) THEN
        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, &
             'PDAF_netf_analysis -- inflate weights '
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

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_netf_analysis:', debug, '  normalized weights', weights
     ELSE
        ! ERROR: weights are zero
        WRITE(*,'(/5x,a/)') 'WARNING: Zero weights - reset to 1/dim_ens'
        weights = 1.0 / REAL(dim_ens)
     END IF

     DEALLOCATE(innov_i, Rinvinnov)

     ! Diagnostic: Compute effective sample size
     CALL PDAF_diag_effsample(dim_ens, weights, effN)
     IF (mype == 0 .AND. screen > 0) &
          WRITE (*, '(a, 5x, a, f10.2)') 'PDAF', '--- Effective sample size ', effN

     CALL PDAF_timeit(51, 'old')

  ELSE
     ! Without observations, all ensemble members have the same weight

     CALL PDAF_timeit(51, 'new')
     weights = 1/dim_ens
     CALL PDAF_timeit(51, 'old')

  END IF haveobs

  CALL PDAF_timeit(10, 'old')
  CALL PDAF_timeit(51, 'new')


  ! ****************************************
  ! *** Calculate the transform matrix   ***
  ! ***      A= (diag(w)-w*w^t)          ***
  ! *** with the weights w               ***
  ! ****************************************

  CALL PDAF_timeit(11, 'new')

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

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_netf_analysis:', debug, '  A', A

  CALL PDAF_timeit(11, 'old')


  ! ********************************************************************
  ! *** Compute ensemble transformation matrix W as square-root of A ***
  ! ********************************************************************

  CALL PDAF_timeit(20, 'new')

  ! Compute symmetric square-root of A by EVD
  ALLOCATE(svals(dim_ens))
  ALLOCATE(work(3*dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 4*dim_ens + dim_ens*dim_ens)
  ldwork = 3*dim_ens

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_netf_analysis:', debug, &
       '  Compute eigenvalue decomposition of A'

  CALL syevTYPE('v', 'l', dim_ens, A, dim_ens, svals, work, ldwork, syev_info)

  IF (syev_info == 0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_netf_analysis:', debug, '  eigenvalues', svals
  ELSE
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

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_netf_analysis:', debug, '  transform', T

  DEALLOCATE(weights, A, T_tmp)

  CALL PDAF_timeit(20, 'old')


! *******************************************
! ***     Transform state ensemble        ***
! ***              a    f                 ***
! ***             X  = X  W               ***
! *** The weight matrix W is stored in T. ***
! *******************************************

  CALL PDAF_timeit(21, 'new')

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
     DO col = 1, dim_ens
        ens_blk(1 : blkupper - blklower + 1, col) &
             = ens_p(blklower : blkupper, col)
     END DO

     !                        a  _f
     ! Transform ensemble:   X = X  W
     CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
          1.0, ens_blk, maxblksize, T, dim_ens, &
          0.0, ens_p(blklower:blkupper, 1), dim_p)

  END DO blocking

  DEALLOCATE(ens_blk)

  CALL PDAF_timeit(21, 'old')


  ! *****************************************
  ! *** Perturb particles by adding noise ***
  ! *****************************************

  CALL PDAF_timeit(23, 'new')

  IF (type_noise>0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_netf_analysis:', debug, '  add noise to particles'

     CALL PDAF_add_particle_noise(dim_p, dim_ens, state_p, ens_p, type_noise, noise_amp, screen)
  END IF

  CALL PDAF_timeit(23, 'old')

  CALL PDAF_timeit(51, 'old')


! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_netf_analysis -- END'

END SUBROUTINE PDAF_netf_ana

!-------------------------------------------------------------------------------
!> Compute ensemble transform of NETS for smoothing
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
!! * 2016-11 - Lars Nerger - Initial code based on NETF_analysis
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_netf_smootherT(step, dim_p, dim_obs_p, dim_ens, &
     ens_p, rndmat, T,  &
     U_init_dim_obs, U_obs_op, U_init_obs, U_likelihood, &
     screen, flag)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_mod_filtermpi, &
       ONLY: mype
  USE PDAF_diag, &
       ONLY: PDAF_diag_effsample

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step         !< Current time step
  INTEGER, INTENT(in) :: dim_p        !< PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p   !< PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens      !< Size of ensemble
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)    !< PE-local state ensemble
  REAL, INTENT(in)    :: rndmat(dim_ens, dim_ens) !< Orthogonal random matrix
  REAL, INTENT(inout) :: T(dim_ens, dim_ens)      !< Ensemble transform matrix
  INTEGER, INTENT(in) :: screen       !< Verbosity flag
  INTEGER, INTENT(inout) :: flag      !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, &       !< Initialize dimension of observation vector
       U_obs_op, &                    !< Observation operator
       U_init_obs, &                  !< Initialize observation vector
       U_likelihood                   !< Compute observation likelihood for an ensemble member
       
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
  REAL, ALLOCATABLE :: svals(:)       ! Singular values of Ainv
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


!-------------------------------------------------------------------------------
!> Smoother extension for NETF
!!
!! Smoother extension for the ensemble square-root filters (ETKF, ESTKF). 
!! The routine uses the matrix Ainv computed by the filter analysis
!! to perform the smoothing on past ensembles.
!!
!! Variant for domain decomposed states.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2016-08 - Lars Nerger - Initial code adating PDAF_smoother
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_smoother_netf(dim_p, dim_ens, dim_lag, Ainv, sens_p, &
     cnt_maxlag, screen)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: mype

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p         !< PE-local dimension of model state
  INTEGER, INTENT(in) :: dim_ens       !< Size of ensemble
  INTEGER, INTENT(in) :: dim_lag       !< Number of past time instances for smoother
  REAL, INTENT(in)   :: Ainv(dim_ens, dim_ens)             !< Weight matrix for ensemble transformation
  REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag)   !< PE-local smoother ensemble
  INTEGER, INTENT(inout) :: cnt_maxlag !< Count available number of time steps for smoothing
  INTEGER, INTENT(in) :: screen        !< Verbosity flag


! *** local variables ***
  INTEGER :: member, col, row, lagcol  ! Counters
  INTEGER :: n_lags                    ! Available number of time instances for smoothing
  INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  REAL :: invdimens                    ! Inverse of global ensemble size
  REAL, ALLOCATABLE :: ens_blk(:,:)    ! Temporary block of state ensemble
  REAL, ALLOCATABLE :: W_smooth(:,:)   ! Weight matrix for smoothing

  
! **********************
! *** INITIALIZATION ***
! **********************

  ! Determine number of time instances for smoothing
  IF (cnt_maxlag >= dim_lag) THEN
     ! Already performed enough analysis to smooth over full lag
     n_lags = dim_lag
  ELSE
     ! Not yet enough analysis steps to smoother over full lag
     n_lags = cnt_maxlag
  END IF

  IF (mype == 0 .AND. screen > 0 .AND. n_lags > 0) THEN
     WRITE (*, '(a, 5x, a, i8)') 'PDAF', 'Perform smoothing up to lag ', n_lags
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

        ! Use block formulation for transformation
        blocking: DO blklower = 1, dim_p, maxblksize
           
           blkupper = MIN(blklower + maxblksize - 1, dim_p)

           ! Store former analysis ensemble
           DO member = 1, dim_ens
              ens_blk(1 : blkupper-blklower+1, member) &
                   = sens_p(blklower : blkupper, member, lagcol)
           END DO

           !                        a   f
           ! Transform ensemble:   X = X  W_smooth
           CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
                1.0, ens_blk(1, 1), maxblksize, W_smooth, dim_ens, &
                0.0, sens_p(blklower, 1, lagcol), dim_p)

        END DO blocking

     END DO smoothing

     DEALLOCATE(ens_blk, W_smooth)

  END IF havelag


! ********************
! *** Finishing up ***
! ********************
  
  ! Increment maxlag counter
  cnt_maxlag = cnt_maxlag + 1

  ! Set flag for memory counting
  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_smoother_netf

END MODULE PDAF_netf_analysis
