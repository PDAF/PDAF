! Copyright (c) 2004-2024 Lars Nerger
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
! !ROUTINE: PDAF_enkf_update --- Control analysis update of the EnKF
!
! !INTERFACE:
SUBROUTINE  PDAF_enkf_update(step, dim_p, dim_obs_p, dim_ens, state_p, &
     ens_p, rank_ana, U_init_dim_obs, U_obs_op, &
     U_add_obs_err, U_init_obs, U_init_obs_covar, U_prepoststep, screen, &
     subtype, dim_lag, sens_p, cnt_maxlag, flag)

! !DESCRIPTION:
! Routine to control the analysis update of the EnKF.
! 
! The analysis is performed by calling PDAF\_enkf\_analysis.
! In addition, the routine U\_prepoststep is called before
! and after the analysis to allow the user to access the 
! ensemble information.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2003-07 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: mype, dim_ens_l
  USE PDAF_mod_filter, &
       ONLY: forget, debug
  USE PDAFobs, &
       ONLY: PDAFobs_initialize, PDAFobs_dealloc, HX_p, HXbar_p, obs_p

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step       ! Current time step
  INTEGER, INTENT(in) :: dim_p      ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens    ! Size of state ensemble
  INTEGER, INTENT(in) :: rank_ana   ! Rank to be considered for inversion of HPH
  REAL, INTENT(inout) :: state_p(dim_p)        ! PE-local model state
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) ! PE-local state ensemble
  INTEGER, INTENT(in) :: screen     ! Verbosity flag
  INTEGER, INTENT(in) :: subtype    ! Specification of filter subtype
  INTEGER, INTENT(in) :: dim_lag     ! Number of past time instances for smoother
  REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag) ! PE-local smoother ensemble
  INTEGER, INTENT(inout) :: cnt_maxlag ! Count number of past time steps for smoothing
  INTEGER, INTENT(inout) :: flag    ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, & ! Initialize dimension of observation vector
       U_obs_op, &              ! Observation operator
       U_init_obs, &            ! Initialize observation vector
       U_init_obs_covar, &      ! Initialize observation error covariance matrix
       U_prepoststep, &         ! User supplied pre/poststep routine
       U_add_obs_err            ! Add observation error covariance matrix

! !CALLING SEQUENCE:
! Called by: PDAF_put_state_enkf
! Calls: U_prepoststep
! Calls: PDAF_enkf_analysis_rlm
! Calls: PDAF_enkf_analysis_rsm
! Calls: PDAF_timeit
! Calls: PDAF_time_temp
!EOP

! *** local variables ***
  INTEGER :: i, member, row       ! Counters
  INTEGER :: minusStep            ! Time step counter
  INTEGER, SAVE :: allocflag = 0  ! Flag whether first time allocation is done
  REAL :: invdimens               ! Inverse global ensemble size
  REAL :: sqrtinvforget           ! square root of inverse forgetting factor
  REAL :: Uinv(1, 1)              ! Unused array, but required in call to U_prepoststep
  REAL, ALLOCATABLE :: HXB(:,:)   ! Ensemble tranformation matrix


! **********************
! ***  Update phase  ***
! **********************

  IF (debug>0) THEN
  WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_enkf_update -- START'
     DO i = 1, dim_ens
        WRITE (*,*) '++ PDAF-debug PDAF_enkf_update:', debug, 'ensemble member', i, &
             ' forecast values (1:min(dim_p,6)):', ens_p(1:min(dim_p,6),i)
     END DO
  END IF


! ***********************************
! *** Compute mean forecast state ***
! ***********************************

  CALL PDAF_timeit(51, 'new')
  CALL PDAF_timeit(11, 'new')

  state_p = 0.0
  invdimens = 1.0 / REAL(dim_ens)
  DO member = 1, dim_ens
     DO row = 1, dim_p
        state_p(row) = state_p(row) + invdimens * ens_p(row, member)
     END DO
  END DO
  
  CALL PDAF_timeit(11, 'old')
  CALL PDAF_timeit(51, 'old')


! *************************************
! *** Prestep for forecast ensemble ***
! *************************************

  CALL PDAF_timeit(5, 'new')
  minusStep = -step  ! Indicate forecast by negative time step number
  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a, i7)') 'PDAF', 'Call pre-post routine after forecast; step ', step
  ENDIF
  CALL U_prepoststep(minusStep, dim_p, dim_ens, dim_ens_l, dim_obs_p, &
       state_p, Uinv, ens_p, flag)
  CALL PDAF_timeit(5, 'old')

  IF (mype == 0 .AND. screen > 0) THEN
     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- duration of prestep:', PDAF_time_temp(5), 's'
     END IF
  END IF


! ************************
! *** Inflate ensemble ***
! ************************

  CALL PDAF_timeit(51, 'new')
  sqrtinvforget = SQRT(1.0 / forget)
  ENSa: DO member = 1, dim_ens
     ! spread out state ensemble according to forgetting factor
     IF (forget /= 1.0) THEN
        ens_p(:, member) = state_p(:) &
             + (ens_p(:, member) - state_p(:)) * sqrtinvforget
     END IF
  END DO ENSa
  CALL PDAF_timeit(51, 'old')


! *****************************************************
! *** Initialize observations and observed ensemble ***
! *****************************************************

  ! Initialize HX_p, HXbar_p, obs_p in module PDAFomi
  CALL PDAFobs_initialize(step, dim_p, dim_ens, dim_obs_p, &
       state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
       screen, debug, .true., .true., .true., .true.)


! ***********************
! ***  Analysis step  ***
! ***********************

  IF (mype == 0 .AND. screen > 0) &
       WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)

#ifndef PDAF_NO_UPDATE
  ALLOCATE(HXB(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens * dim_ens)

  CALL PDAF_timeit(3, 'new')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_enkf_update', debug, &
          'Configuration: param_int(3) rank_ana_enkf', rank_ana
     WRITE (*,*) '++ PDAF-debug PDAF_enkf_update', debug, &
          'Configuration: param_int(4) -not used-  '
     WRITE (*,*) '++ PDAF-debug PDAF_enkf_update', debug, &
          'Configuration: param_int(5) dim_lag      ', dim_lag

     WRITE (*,*) '++ PDAF-debug PDAF_enkf_update', debug, &
          'Configuration: param_real(1) forget     ', forget
  END IF

  IF (subtype == 0) THEN

     ! *** analysis with representer method - with 2m>n ***
     CALL PDAF_enkf_analysis_rlm(step, dim_p, dim_obs_p, dim_ens, rank_ana, &
          state_p, ens_p, HXB, HX_p, HXbar_p, obs_p, &
          U_add_obs_err, U_init_obs_covar, screen, debug, flag)

     ! *** Perform smoothing of past ensembles ***
     CALL PDAF_timeit(51, 'new')
     CALL PDAF_smoother_enkf(dim_p, dim_ens, dim_lag, HXB, sens_p, &
          cnt_maxlag, forget, screen)
     CALL PDAF_timeit(51, 'old')

  ELSE IF (subtype == 1) THEN

     ! *** analysis with representer method with 2m<n ***
     CALL PDAF_enkf_analysis_rsm(step, dim_p, dim_obs_p, dim_ens, rank_ana, &
          state_p, ens_p, HX_p, HXbar_p, obs_p, &
          U_add_obs_err, U_init_obs_covar, screen, debug, flag)

  END IF

  IF (debug>0) THEN
     DO i = 1, dim_ens
        WRITE (*,*) '++ PDAF-debug PDAF_enkf_update:', debug, 'ensemble member', i, &
             ' analysis values (1:min(dim_p,6)):', ens_p(1:min(dim_p,6),i)
     END DO
  END IF

  CALL PDAF_timeit(3, 'old')

  IF (mype == 0 .AND. screen > 1) THEN
     WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
          'PDAF', '--- update duration:', PDAF_time_temp(3), 's'
  END IF

  DEALLOCATE(HXB)

#else
  WRITE (*,'(/5x,a/)') &
       '!!! PDAF WARNING: ANALYSIS STEP IS DEACTIVATED BY PDAF_NO_UPDATE !!!'
#endif

! *** poststep for analysis ensemble ***
  CALL PDAF_timeit(5, 'new')
  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a)') 'PDAF', 'Call pre-post routine after analysis step'
  ENDIF
  CALL U_prepoststep(step, dim_p, dim_ens, dim_ens_l, dim_obs_p, &
       state_p, Uinv, ens_p, flag)
  CALL PDAF_timeit(5, 'old')
  
  IF (mype == 0 .AND. screen > 0) THEN
     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- duration of poststep:', PDAF_time_temp(5), 's'
     END IF
     WRITE (*, '(a, 55a)') 'PDAF Forecast ', ('-', i = 1, 55)
  END IF


! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1

  ! Deallocate observation arrays
  CALL PDAFobs_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_enkf_update -- END'

END SUBROUTINE PDAF_enkf_update
