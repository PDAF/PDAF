! Copyright (c) 2014-2025 Lars Nerger
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
!> PF analysis with resampling
!!
!! Analysis step of the PF
!!
!! Variant for domain decomposed states. 
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2019-05 - Lars Nerger
!! * Later revisions - see svn log
!!
MODULE PDAF_pf_analysis

CONTAINS
SUBROUTINE PDAF_pf_ana(step, dim_p, dim_obs_p, dim_ens, &
     state_p, ens_p, type_resample, &
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
       ONLY: PDAF_add_particle_noise

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                  !< Current time step
  INTEGER, INTENT(in) :: dim_p                 !< PE-local dimension of model state
  INTEGER, INTENT(in) :: dim_obs_p             !< PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens               !< Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)        !< PE-local forecast mean state
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) !< PE-local state ensemble
  INTEGER, INTENT(in) :: type_resample         !< Type of resampling scheme
  INTEGER, INTENT(in) :: type_winf             !< Type of weights inflation
  REAL, INTENT(in) :: limit_winf               !< Limit for weights inflation
  INTEGER, INTENT(in) :: type_noise            !< Type of pertubing noise
  REAL, INTENT(in) :: noise_amp                !< Amplitude of noise
  REAL, INTENT(in) :: HZ_p(dim_obs_p, dim_ens) !< Temporary matrices for analysis
  REAL, INTENT(in) :: obs_p(dim_obs_p)         !< PE-local observation vector
  INTEGER, INTENT(in) :: screen                !< Verbosity flag
  INTEGER, INTENT(in) :: debug                 !< Flag for writing debug output
  INTEGER, INTENT(inout) :: flag               !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_likelihood                     !< Compute observation likelihood for an ensemble member

! *** local variables ***
  INTEGER :: i, col, member           ! counters
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  REAL :: effN                        ! Efective sample size
  REAL :: weight                      ! Ensemble weight (likelihood)
  REAL :: total_weight                ! Sum of weights
  INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
  REAL, ALLOCATABLE :: innov_i(:)     ! PE-local observation residual
  REAL, ALLOCATABLE :: Rinvinnov(:)   ! R^-1 times residual 
  REAL, ALLOCATABLE :: weights(:)     ! Weight vector
  REAL, ALLOCATABLE :: ens_blk(:,:)   ! Temporary block of state ensemble
  INTEGER, ALLOCATABLE :: IDs(:)      ! Indices for resampled particles


! **********************
! *** INITIALIZATION ***
! **********************

  CALL PDAF_timeit(51, 'new')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_pf_analysis -- START'

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a)') &
          'PDAF', 'Compute particle filter update'
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

     ! Get residual as difference of observation and observed state for each ensemble member
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, &
          'PDAF_netf_analysis -- call likelihood', dim_ens, 'times'

     CALC_w: DO member = 1, dim_ens

        CALL PDAF_timeit(51, 'new')
        innov_i = obs_p - HZ_p(:, member) 
        CALL PDAF_timeit(51, 'old')

        IF (debug>0) THEN
           WRITE (*,*) '++ PDAF-debug: ', debug, &
                'PDAF_pf_analysis -- member', member
           WRITE (*,*) '++ PDAF-debug PDAF_pf_analysis:', debug, '  innovation d', innov_i
           WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_pf_analysis -- call likelihood'
        end IF

        ! Compute likelihood
        CALL PDAF_timeit(48, 'new')
        CALL U_likelihood(step, dim_obs_p, obs_p, innov_i, weight)
        CALL PDAF_timeit(48, 'old')
        weights(member) = weight

     END DO CALC_w

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_pf_analysis:', debug, '  raw weights', weights

     ! Compute inflation of weights according to N_eff
     IF (type_winf == 1) THEN
        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, &
             'PDAF_pf_analysis -- inflate weights'
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
             WRITE (*,*) '++ PDAF-debug PDAF_pf_analysis:', debug, '  normalized weights', weights
     ELSE
        ! weights are zero - reset to uniform weights
        WRITE(*,'(/5x,a/)') 'WARNING: Zero weights - reset to 1/dim_ens'
        weights = 1.0/REAL(dim_ens)
     END IF

     DEALLOCATE(innov_i, Rinvinnov)

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

  CALL PDAF_timeit(10, 'old')
  CALL PDAF_timeit(51, 'new')


  ! ****************************************
  ! *** Resample particles               ***
  ! ****************************************

  CALL PDAF_timeit(12, 'new')

  ! Determine sample IDs for resampling

  ALLOCATE(IDs(dim_ens))

  CALL PDAF_timeit(21, 'new')

  CALL PDAF_pf_resampling(type_resample, dim_ens, dim_ens, weights, IDs, screen)

  CALL PDAF_timeit(21, 'old')

  ! Perform the resampling (in blocked form to save memory)

  maxblksize = 200
  IF (mype == 0 .AND. screen > 0) &
       WRITE (*, '(a, 5x, a, i5)') &
       'PDAF', '--- use blocking with size ', maxblksize

  ALLOCATE(ens_blk(maxblksize, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', maxblksize * dim_ens)

  blocking: DO blklower = 1, dim_p, maxblksize

     blkupper = MIN(blklower + maxblksize - 1, dim_p)

     ! Store old state ensemble
     CALL PDAF_timeit(22, 'new')
     DO col = 1, dim_ens
        ens_blk(1 : blkupper - blklower + 1, col) &
             = ens_p(blklower : blkupper, col)
     END DO

     ! Perform the resampling 
     DO col = 1, dim_ens
        ens_p(blklower : blkupper, col) = ens_blk(1 : blkupper-blklower+1, IDs(col))
     END DO
     CALL PDAF_timeit(22, 'old')

  END DO blocking


  ! *****************************************
  ! *** Perturb particles by adding noise ***
  ! *****************************************

  CALL PDAF_timeit(23, 'new')

  IF (type_noise>0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_pf_analysis:', debug, '  add noise to particles'
     CALL PDAF_add_particle_noise(dim_p, dim_ens, state_p, ens_p, type_noise, noise_amp, screen)
  END IF

  CALL PDAF_timeit(23, 'old')

  CALL PDAF_timeit(12, 'old')
  CALL PDAF_timeit(51, 'old')


! ********************
! *** Finishing up ***
! ********************
        
  ! Set exit flag
  flag = 0

  DEALLOCATE(ens_blk, IDs)

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_pf_ana

!-------------------------------------------------------------------------------
!> Get particle indices for resampling
!!
!! Determine particle indices for resampling. Implemented
!! are three sampling schemes:
!! (1) probabilistic resampling
!! (2) stochastic unversal resampling
!! (3) residual resampling
!! The schemes follow the algorithms in 
!! Vetra-Carvalho et al., State-of-the-art stochastic data
!! assimilation methods for high-dimensional non-Gaussian problems.
!! Tellus A, 70:1, 1445364, doi:10.1080/16000870.2018.1445364
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2019-05 - Lars Nerger initial code
!! * Later revisions - see svn log
!!
SUBROUTINE PDAF_pf_resampling(method, Nin, Nout, weights, IDs, screen)

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
  INTEGER, INTENT(in) :: method       !< Choose resampling method
                                      !< (1) Probabilistic resampling
                                      !< (2) Stochastic universal resampling
                                      !< (3) Residual resampling
  INTEGER, INTENT(in) :: Nin          !< number of particles
  INTEGER, INTENT(in) :: Nout         !< number of particles to be resampled
  REAL, INTENT(in)    :: weights(Nin) !< Weights
  INTEGER, INTENT(out) :: IDs(Nout)   !< Indices of resampled ensmeble states
  INTEGER, INTENT(in) :: screen       !< Verbosity flag
       
! *** local variables ***
  INTEGER :: i, j                    !< Loop counters
  INTEGER :: c, Nr                   !< Counter for resampling
  INTEGER, SAVE :: first = 1         !< flag for init of random number seed
  INTEGER, SAVE :: iseed(4)          !< seed array for random number routine
  REAL :: rndval                     !< Random value
  REAL, ALLOCATABLE :: w_acc(:)      !< accumulated weights
  INTEGER, ALLOCATABLE :: w_i(:)     !< Integer weights
  REAL, ALLOCATABLE :: w_r(:)        !< residual weights


! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a)') &
          'PDAF', 'Resample particles'
     IF (method == 1) THEN
        WRITE (*, '(a, 5x, a)') &
             'PDAF', '--- Probabilistic resampling'
     ELSE IF (method == 2) THEN
        WRITE (*, '(a, 5x, a)') &
             'PDAF', '--- Stochastic universal resampling'
     ELSE IF (method == 3) THEN
        WRITE (*, '(a, 5x, a)') &
             'PDAF', '--- Residual resampling'
     END IF
  END IF

  IF (first==1) THEN
     ! Set random seed
     iseed(1)=15
     iseed(2)=25
     iseed(3)=47
     iseed(4)=17
     first = 0
  END IF

  ! Initialize accumulated weights
  ALLOCATE(w_acc(Nin))
  w_acc = 0.0

  ! Initialize rasampling IDs
  IDs = 0

   
  Rtype: IF (method == 1) THEN

! ********************************
! *** Probabilistic resampling ***
! ********************************

     ! Get accumulated weights
     w_acc(1) = weights(1)

     DO i = 2, Nin
        W_acc(i) = w_acc(i-1) + weights(i)
     END DO

     c = 1
     DO j = 1, Nout

        ! Init random number
        CALL larnvTYPE(1, iseed, 1, rndval)

        ! Determine index
        checkacc: DO i = 1, Nin
           IF (rndval > w_acc(i)) THEN
              c = c + 1
           ELSE
              EXIT checkacc
           END IF
        END DO checkacc

        IDs(j) = c
        c = 1

     END DO


  ELSE IF (method == 2) THEN

! ***************************************
! *** Stochastic universal resampling ***
! ***************************************

     ! Get accumulated weights
     w_acc(1) = weights(1)

     DO i = 2, Nin
        w_acc(i) = w_acc(i-1) + weights(i)
     END DO

     ! Init random number
     CALL larnvTYPE(1, iseed, 1, rndval)
     rndval = rndval / Nin


     c = 1

     DO j = 1, Nout

        ! Determine index
        checkaccB: DO i = 1, Nin
           IF (rndval > w_acc(i)) THEN
              c = c + 1
           ELSE
              EXIT checkaccB
           END if
        END DO checkaccB

        IDs(j) = c
        
        rndval = rndval + 1.0/REAL(Nin)

        c = 1

     END DO

  ELSE IF (method == 3) THEN

! ***************************
! *** Residual resampling ***
! ***************************

     ALLOCATE(w_i(Nin))
     ALLOCATE(w_r(Nin))

     DO i = 1, Nin
        w_i(i) = FLOOR(weights(i) * Nout)
     END DO

     ! First round setting resamping indices
     Nr = Nout
     c = 1
     IDs = 0
     DO j = 1, Nout
        IF (w_i(j) > 0) THEN
           IDs(c:c+w_i(j)-1) = j
           c = c + w_i(j)
           Nr = Nr - w_i(j)
        END IF
     END DO

     ! Now perform probabilistic resampling for remaining Nr indices

     ! Determine residual weights 

     IF (Nr > 0) THEN
        DO i = 1, Nin
           w_r(i) = (weights(i)*Nout - REAL(w_i(i))) / REAL (Nr)
        END DO
     ELSE
        w_r(:) = 0
     END IF

     ! Get accumulated weights
     w_acc(1) = w_r(1)

     DO i = 2, Nin
        w_acc(i) = w_acc(i-1) + w_r(i)
     END DO

     c = 1
     DO j = Nout-Nr+1, Nout

        ! Init random number
        CALL larnvTYPE(1, iseed, 1, rndval)
        rndval = rndval * REAL(Nout) / REAL(Nin)

        ! Determine index
        checkaccC: DO i = 1, Nin
           IF (rndval > w_acc(i)) THEN
              c = c + 1
           ELSE
              EXIT checkaccC
           END IF
        END DO checkaccC

        IDs(j) = c
        c = 1

     END DO

     DEALLOCATE(w_i, w_r)

    END IF Rtype


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(w_acc)
END SUBROUTINE PDAF_pf_resampling

END MODULE PDAF_pf_analysis
