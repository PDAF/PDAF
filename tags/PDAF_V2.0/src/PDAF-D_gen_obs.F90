! Copyright (c) 2004-2021 Lars Nerger
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
! !ROUTINE: PDAF_gen_obs --- Generate synthetic observations
!
! !INTERFACE:
SUBROUTINE PDAF_gen_obs(step, dim_p, dim_obs_f, dim_ens, &
     state_p, Ainv, ens_p, &
     U_init_dim_obs_f, U_obs_op_f, U_get_obs_f, U_init_obserr_f, &
     U_prepoststep, screen, flag)
  
! !DESCRIPTION:
! This routine generates observations from a model state.
! Used are the general functionality of PDAF provided with
! the call-back routines init_dim_obs, obs_op.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2019-01 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filter, &
       ONLY: obs_member
  USE PDAF_mod_filtermpi, &
       ONLY: mype, dim_ens_l

  IMPLICIT NONE

! !ARGUMENTS:
! ! Variable naming scheme:
! !   suffix _p: Denotes a full variable on the PE-local domain
! !   suffix _l: Denotes a local variable on the current analysis domain
! !   suffix _f: Denotes a full variable of all observations required for the
! !              analysis loop on the PE-local domain
  INTEGER, INTENT(in) :: step        ! Current time step
  INTEGER, INTENT(in) :: dim_p       ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_f  ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens     ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)        ! PE-local model state
  REAL, INTENT(inout) :: Ainv(dim_ens, dim_ens)      ! Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) ! PE-local ensemble matrix
  INTEGER, INTENT(in) :: screen      ! Verbosity flag
  INTEGER, INTENT(inout) :: flag     ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_obs_op_f, &    ! Observation operator
       U_init_dim_obs_f, &     ! Initialize dimension of observation vector
       U_get_obs_f, &          ! Provide observation vector to user
       U_init_obserr_f, &      ! Initialize vector of observation error standard deviations
       U_prepoststep           ! User supplied pre/poststep routine

! !CALLING SEQUENCE:
! Called by:               PDAF_put_state_letkf
! Calls: U_prepoststep
! Calls: U_init_dim_obs
! Calls: U_obs_op
! Calls: U_init_obs
! Calls: PDAF_generate_rndmat
! Calls: PDAF_timeit
! Calls: PDAF_memcount
! Calls: MPI_reduce
!EOP

! *** local variables ***
  INTEGER :: i, member, row        ! Counters
  INTEGER, SAVE :: first = 0       ! Flag whether routine is called first time
  REAL    :: invdimens             ! Inverse global ensemble size
  INTEGER :: minusStep             ! Time step counter
  REAL, ALLOCATABLE :: Hx_f(:)     ! HX for PE-local model state
  REAL, ALLOCATABLE :: rndvec(:)   ! random vector to perturb observations
  REAL, ALLOCATABLE :: rms_obs(:)  ! vector of observation error standard deviations
  INTEGER, SAVE :: iseed(4)        ! seed array for random number routine
  INTEGER :: seedset = 1           ! Choice of seed set for random numbers


! ******************************
! *** Set random number seed ***
! ******************************

  IF (first == 0) THEN
     IF (seedset == 2) THEN
        iseed(1)=1
        iseed(2)=5
        iseed(3)=7
        iseed(4)=9
     ELSE IF (seedset == 3) THEN
        iseed(1)=2
        iseed(2)=5
        iseed(3)=7
        iseed(4)=9
     ELSE IF (seedset == 4) THEN
        iseed(1)=1
        iseed(2)=6
        iseed(3)=7
        iseed(4)=9
     ELSE IF (seedset == 5) THEN
        iseed(1)=1
        iseed(2)=5
        iseed(3)=8
        iseed(4)=9
     ELSE IF (seedset == 6) THEN
        iseed(1)=2
        iseed(2)=5
        iseed(3)=8
        iseed(4)=9
     ELSE IF (seedset == 7) THEN
        iseed(1)=2
        iseed(2)=6
        iseed(3)=8
        iseed(4)=9
     ELSE IF (seedset == 8) THEN
        iseed(1)=2
        iseed(2)=6
        iseed(3)=8
        iseed(4)=11
     ELSE IF (seedset == 9) THEN
        iseed(1)=3
        iseed(2)=6
        iseed(3)=8
        iseed(4)=11
     ELSE IF (seedset == 10) THEN
        iseed(1)=3
        iseed(2)=7
        iseed(3)=8
        iseed(4)=11
     ELSE IF (seedset == 11) THEN
        iseed(1)=13
        iseed(2)=7
        iseed(3)=8
        iseed(4)=11
     ELSE IF (seedset == 12) THEN
        iseed(1)=13
        iseed(2)=11
        iseed(3)=8
        iseed(4)=11
     ELSE IF (seedset == 13) THEN
        iseed(1)=13
        iseed(2)=13
        iseed(3)=8
        iseed(4)=11
     ELSE IF (seedset == 14) THEN
        iseed(1)=13
        iseed(2)=13
        iseed(3)=17
        iseed(4)=11
     ELSE IF (seedset == 15) THEN
        iseed(1)=13
        iseed(2)=13
        iseed(3)=19
        iseed(4)=11
     ELSE IF (seedset == 16) THEN
        iseed(1)=15
        iseed(2)=13
        iseed(3)=19
        iseed(4)=11
     ELSE IF (seedset == 17) THEN
        iseed(1)=15
        iseed(2)=135
        iseed(3)=19
        iseed(4)=11
     ELSE IF (seedset == 18) THEN
        iseed(1)=19
        iseed(2)=135
        iseed(3)=19
        iseed(4)=11
     ELSE IF (seedset == 19) THEN
        iseed(1)=19
        iseed(2)=135
        iseed(3)=19
        iseed(4)=17
     ELSE IF (seedset == 20) THEN
        iseed(1)=15
        iseed(2)=15
        iseed(3)=47
        iseed(4)=17
     ELSE
        ! Standard seed
        iseed(1) = 1000
        iseed(2) = 2034
        iseed(3) = 0
        iseed(4) = 3
     END IF
     first = 2
  END IF

! *************************************
! *** Prestep for forecast ensemble ***
! *************************************


  CALL PDAF_timeit(5, 'new')
  minusStep = - step  ! Indicate forecast by negative time step number
  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a, i7)') 'PDAF', 'Call pre-post routine after forecast; step ', step
  ENDIF
  CALL U_prepoststep(minusStep, dim_p, dim_ens, dim_ens_l, dim_obs_f, &
       state_p, Ainv, ens_p, flag)
  CALL PDAF_timeit(5, 'old')

  IF (mype == 0 .AND. screen > 0) THEN
     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF ', '--- duration of prestep:', PDAF_time_temp(5), 's'
     END IF
     WRITE (*, '(a, 55a)') 'PDAF generate observations ', ('-', i = 1, 42)
  END IF


! ********************************************************************
! *** Generate synthetic observations                              ***
! *** Steps:                                                       ***
! ***   1. Get observed state by applying the observation operator ***
! ***   2. Generate Gaussian random number                         ***
! ***   3. Get vector of observation error standard deviations     ***
! ***   4. Add random noise scaled by obs. error to observed state ***
! ***   5. Provide observation vector to user                      ***
! ********************************************************************

  CALL PDAF_timeit(3, 'new')

! *** 1. Get full observed state vector ***

  ! Get full observation dimension
  CALL U_init_dim_obs_f(step, dim_obs_f)

  IF (screen > 0) THEN
     IF (screen<=2 .AND. mype == 0) THEN
        WRITE (*, '(a, 5x, a, i6, a, i10)') &
             'PDAF', '--- PE-Domain:', mype, &
             ' dimension of PE-local full obs. vector', dim_obs_f
     ELSE IF (screen>2) THEN
        WRITE (*, '(a, 5x, a, i6, a, i10)') &
             'PDAF', '--- PE-Domain:', mype, &
             ' dimension of PE-local full obs. vector', dim_obs_f
     END IF
  END IF

  CALL PDAF_timeit(11, 'new')

  ! Compute mean forecast state
  state_p = 0.0
  invdimens = 1.0 / REAL(dim_ens)
  DO member = 1, dim_ens
     DO row = 1, dim_p
        state_p(row) = state_p(row) + invdimens * ens_p(row, member)
     END DO
  END DO

  CALL PDAF_timeit(11, 'old')

  CALL PDAF_timeit(33, 'new')

  ! Apply observation operator to get Hx for full DIM_OBS_F region on PE-local domain
  ALLOCATE(Hx_f(dim_obs_f))
  IF (first == 0) CALL PDAF_memcount(3, 'r', dim_obs_f)

  obs_member = 1  ! Set this for completeness

  CALL U_obs_op_f(step, dim_p, dim_obs_f, state_p, Hx_f)


! *** 2. Initialize random vector ***

  ALLOCATE(rndvec(dim_obs_f))
  IF (first == 0) CALL PDAF_memcount(3, 'r', dim_ens)

  CALL larnvTYPE(3, iseed, dim_obs_f, rndvec)

! *** 3. Get vector of observation errors ***

  ALLOCATE(rms_obs(dim_obs_f))
  IF (first == 0) CALL PDAF_memcount(3, 'r', dim_ens)

  CALL U_init_obserr_f(step, dim_obs_f, Hx_f, rms_obs)

! *** 4. Generate perturbed observations ***
  DO i = 1, dim_obs_f
     Hx_f(i) = Hx_f(i) + rms_obs(i)*rndvec(i)
  END DO

! *** Provide observation vector to user ***
  CALL U_get_obs_f(step, dim_obs_f, Hx_f)

  CALL PDAF_timeit(33, 'old')

  CALL PDAF_timeit(3, 'old')

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 55a)') 'PDAF Forecast ', ('-', i = 1, 55)
  END IF


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(Hx_f, rndvec, rms_obs)

  IF (first == 0) first = 1

END SUBROUTINE PDAF_gen_obs
