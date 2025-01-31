! Copyright (c) 2014-2024 Lars Nerger
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
! !ROUTINE: PDAF_pf_analysis --- PF analysis with resampling
!
! !INTERFACE:
SUBROUTINE PDAF_pf_analysis(step, dim_p, dim_obs_p, dim_ens, &
     state_p, ens_p, type_resample, &
     type_winf, limit_winf, noise_type, noise_amp, &
     HZ_p, obs_p, U_likelihood, screen, debug, flag)

! !DESCRIPTION:
! Analysis step of the PF
!
! Variant for domain decomposed states. 
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2019-05 - Lars Nerger
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

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step         ! Current time step
  INTEGER, INTENT(in) :: dim_p        ! PE-local dimension of model state
  INTEGER, INTENT(in) :: dim_obs_p    ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens      ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local forecast mean state
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(in) :: type_resample           ! Type of resampling scheme
  INTEGER, INTENT(in) :: type_winf    ! Type of weights inflation
  REAL, INTENT(in) :: limit_winf      ! Limit for weights inflation
  INTEGER, INTENT(in) :: noise_type   ! Type of pertubing noise
  REAL, INTENT(in) :: noise_amp       ! Amplitude of noise
  REAL, INTENT(in) :: HZ_p(dim_obs_p, dim_ens)    ! Temporary matrices for analysis
  REAL, INTENT(in) :: obs_p(dim_obs_p)            ! PE-local observation vector
  INTEGER, INTENT(in) :: screen       ! Verbosity flag
  INTEGER, INTENT(in) :: debug        ! Flag for writing debug output
  INTEGER, INTENT(inout) :: flag      ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_likelihood            ! Compute observation likelihood for an ensemble member

! !CALLING SEQUENCE:
! Called by: PDAF_pf_update
! Calls: U_likelihood
! Calls: PDAF_timeit
! Calls: PDAF_memcount
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: syevTYPE (LAPACK; dsyev or ssyev dependent on precision)
!EOP

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
        CALL PDAF_timeit(47, 'new')
        CALL U_likelihood(step, dim_obs_p, obs_p, innov_i, weight)
        CALL PDAF_timeit(47, 'old')
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

  IF (noise_type>0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_pf_analysis:', debug, '  add noise to particles'
     CALL PDAF_pf_add_noise(dim_p, dim_ens, state_p, ens_p, noise_type, noise_amp, screen)
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

END SUBROUTINE PDAF_pf_analysis
