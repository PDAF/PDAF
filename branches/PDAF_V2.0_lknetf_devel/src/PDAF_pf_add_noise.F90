! Copyright (c) 2019-2021 Lars Nerger
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
! !ROUTINE: PDAF_pf_add_noise --- Add noise to particles after resampling
!
! !INTERFACE:
SUBROUTINE PDAF_pf_add_noise(dim_p, dim_ens, state_p, ens_p, noise_type, noise_amp, screen)

! !DESCRIPTION:
! Adding noise to particles to avoid identical particles
! after resampling.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2019-07 - Lars Nerger initial code
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
  INTEGER, INTENT(in) :: dim_p          ! State dimension
  INTEGER, INTENT(in) :: dim_ens        ! Number of particles
  REAL, INTENT(inout) :: state_p(dim_p) ! State vector (not filled)
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) ! Ensemble array
  INTEGER, INTENT(in) :: noise_type     ! Type of noise
  REAL, INTENT(in) :: noise_amp         ! Noise amplitude
  INTEGER, INTENT(in) :: screen         ! Verbosity flag

! !CALLING SEQUENCE:
!EOP
       
! *** local variables ***
  INTEGER :: i, member                ! Loop counters
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  INTEGER, SAVE :: first = 1          ! flag for init of random number seed
  INTEGER, SAVE :: iseed(4)           ! seed array for random number routine
  REAL, ALLOCATABLE :: ens_noise(:,:) ! Noise to be added for PF
  REAL :: noisenorm                   ! output argumemt of PDAF_enkf_omega (not used)
  REAL :: invdim_ens                  ! Inverse ensemble size
  REAL :: invdim_ensm1                ! Inverse of ensemble size minus 1
  REAL :: variance                    ! Ensmeble variance

! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a)') &
          'PDAF', 'Perturb particles:'
     IF (noise_type == 1) THEN
        WRITE (*, '(a, 5x, a, f10.4)') &
             'PDAF', '--- Gaussian noise with constant standard deviation', noise_amp
     ELSEIF (noise_type == 2) THEN
        WRITE (*, '(a, 5x, a, es10.3, a)') &
             'PDAF', '--- Gaussian noise with amplitude ', noise_amp,' of ensemble standard deviation'
     END IF
  END IF

  ! Initialized seed for random number routine
  IF (first == 1) THEN
     iseed(1) = 1000
     iseed(2) = 2045+mype
     iseed(3) = 10
     iseed(4) = 3
     first = 2
  END IF

  ! Initialize numbers
  invdim_ens    = 1.0 / REAL(dim_ens)  
  invdim_ensm1  = 1.0 / REAL(dim_ens - 1)


! ******************************
! *** Add noise to particles ***
! ******************************

  ALLOCATE(ens_noise(1, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)


  IF (noise_type == 1) THEN

     ! *** Noise with constant standard deviation ***

     DO i = 1, dim_p
        CALL PDAF_enkf_omega(iseed, dim_ens, 1, ens_noise(1,:), noisenorm, 8, 0)

        DO member = 1, dim_ens
           ens_p(i, member) = ens_p(i, member) + noise_amp * ens_noise(1,member)
        END DO
     END DO

  ELSEIF (noise_type == 2) THEN

     ! *** Noise with fraction of ensemble standard deviation ***

     ! Compute mean state
     state_p = 0.0
     DO member = 1, dim_ens
        DO i = 1, dim_p
           state_p(i) = state_p(i) + ens_p(i, member)
        END DO
     END DO
     state_p(:) = invdim_ens * state_p(:)

     ! Add noise
     DO i = 1, dim_p

        ! Initialize noise vector with zero mean and unit variance
        CALL PDAF_enkf_omega(iseed, dim_ens, 1, ens_noise(1,:), noisenorm, 8, 0)

       ! Compute sampled variance
        variance = 0.0
        DO member = 1, dim_ens
           variance = variance + (ens_p(i, member) - state_p(i))*(ens_p(i, member) - state_p(i))
        END DO
        variance = invdim_ensm1 * variance

        ! Add noise to particles
        DO member = 1, dim_ens
           ens_p(i, member) = ens_p(i, member) + noise_amp * sqrt(variance) * ens_noise(1, member)
        END DO
     END DO

  END IF


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(ens_noise)

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_pf_add_noise
