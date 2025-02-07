! Copyright (c) 2004-2025 Lars Nerger
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
!> Control analysis update of the SEIK filter
!!
!! Routine to control the analysis update of the SEIK filter.
!! 
!! The analysis is performed by calling PDAF\_seik\_analysis
!! and the resampling is performed in PDAF\_seik\_resample.
!! In addition, the routine U\_prepoststep is called prior
!! to the analysis and after the resampling to allow the user
!! to access the ensemble information.
!!
!! Variant for SEIK with domain decompostion.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2003-07 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
SUBROUTINE  PDAF_seik_update(step, dim_p, dim_obs_p, dim_ens, rank, &
     state_p, Uinv, ens_p, state_inc_p, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, U_init_obsvar, &
     U_prepoststep, screen, subtype, incremental, flag)

  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_mod_filtermpi, &
       ONLY: mype, dim_ens_l
  USE PDAF_seik, &
       ONLY: filterstr, debug, forget, type_forget, &
       type_trans, Nm1vsN, type_sqrt
  USE PDAFobs, &
       ONLY: PDAFobs_init, PDAFobs_dealloc, type_obs_init, &
       HX_p, HXbar_p, obs_p, observe_ens

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step        !< Current time step
  INTEGER, INTENT(in) :: dim_p       !< PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p  !< PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens     !< Size of ensemble
  INTEGER, INTENT(in) :: rank        !< Rank of initial covariance matrix
  REAL, INTENT(inout) :: state_p(dim_p)        !< PE-local model state
  REAL, INTENT(inout) :: Uinv(rank, rank)      !< Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) !< PE-local ensemble matrix
  REAL, INTENT(inout) :: state_inc_p(dim_p)    !< PE-local state analysis increment
  INTEGER, INTENT(in) :: screen      !< Verbosity flag
  INTEGER, INTENT(in) :: subtype     !< Filter subtype
  INTEGER, INTENT(in) :: incremental !< Control incremental updating
  INTEGER, INTENT(inout) :: flag     !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, & !< Initialize dimension of observation vector
       U_obs_op, &              !< Observation operator
       U_init_obs, &            !< Initialize observation vector
       U_init_obsvar, &         !< Initialize mean observation error variance
       U_prepoststep, &         !< User supplied pre/poststep routine
       U_prodRinvA              !< Provide product R^-1 A for SEIK analysis

! *** local variables ***
  INTEGER :: i, j               ! Counters
  INTEGER :: minusStep          ! Time step counter
  REAL :: forget_ana            ! Forgetting factor actually used in analysis
  LOGICAL :: do_init_dim_obs    ! Flag for initializing dim_obs_p in PDAFobs_init


! ***********************************************************
! *** For fixed error space basis compute ensemble states ***
! ***********************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_seik_update -- START'

  CALL PDAF_timeit(3, 'new')
  CALL PDAF_timeit(51, 'new')

  fixed_basis: IF (subtype == 2 .OR. subtype == 3) THEN
     ! *** Add mean/central state to ensemble members ***
     DO j = 1, dim_ens
        DO i = 1, dim_p
           ens_p(i, j) = ens_p(i, j) + state_p(i)
        END DO
     END DO
  END IF fixed_basis

  IF (debug>0) THEN
     DO i = 1, dim_ens
        WRITE (*,*) '++ PDAF-debug PDAF_seik_update:', debug, 'ensemble member', i, &
             ' forecast values (1:min(dim_p,6)):', ens_p(1:min(dim_p,6),i)
     END DO
  END IF
  CALL PDAF_timeit(51, 'old')


! *****************************************************
! *** Initialize observations and observed ensemble ***
! *****************************************************

  IF (type_obs_init==0 .OR. type_obs_init==2) THEN
     ! This call initializes dim_obs_p, HX_p, HXbar_p, obs_p in the module PDAFobs
     ! It also compute the ensemble mean and stores it in state_p
     CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_p, &
          state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
          screen, debug, .true., .true., .true., .true., .true.)
  END IF

  CALL PDAF_timeit(3, 'old')


! *************************************
! *** Prestep for forecast ensemble ***
! *************************************

  CALL PDAF_timeit(5, 'new')
  minusStep = -step  ! Indicate forecast by negative time step number
  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 52a)') 'PDAF Prepoststep ', ('-', i = 1, 52)
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


! *****************************************************
! *** Initialize observations and observed ensemble ***
! *****************************************************

  IF (type_obs_init>0) THEN
     CALL PDAF_timeit(3, 'new')

     IF (type_obs_init==1) THEN
        do_init_dim_obs=.true.
     ELSE
        ! Skip call to U_init_dim_obs when also called before prepoststep
        do_init_dim_obs=.false.   
     END IF

     ! This call initializes dim_obs_p, HX_p, HXbar_p, obs_p in the module PDAFobs
     ! It also compute the ensemble mean and stores it in state_p
     CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_p, &
          state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
          screen, debug, .true., do_init_dim_obs, .true., .true., .true.)

     CALL PDAF_timeit(3, 'old')
  END IF


! ***********************
! ***  Analysis step  ***
! ***********************

  IF (mype == 0 .AND. screen > 0) &
       WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)

#ifndef PDAF_NO_UPDATE
  CALL PDAF_timeit(3, 'new')

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_seik_update', debug, &
          'Configuration: param_int(3) -not used-  '
     WRITE (*,*) '++ PDAF-debug PDAF_seik_update', debug, &
          'Configuration: param_int(4) incremental ', incremental
     WRITE (*,*) '++ PDAF-debug PDAF_seik_update', debug, &
          'Configuration: param_int(5) type_forget ', type_forget
     WRITE (*,*) '++ PDAF-debug PDAF_seik_update', debug, &
          'Configuration: param_int(6) type_trans  ', type_trans
     WRITE (*,*) '++ PDAF-debug PDAF_seik_update', debug, &
          'Configuration: param_int(7) type_sqrt   ', type_sqrt
     WRITE (*,*) '++ PDAF-debug PDAF_seik_update', debug, &
          'Configuration: param_int(8) observe_ens           ', observe_ens

     WRITE (*,*) '++ PDAF-debug PDAF_seik_update', debug, &
          'Configuration: param_real(1) forget     ', forget
  END IF


! *** Compute adaptive forgetting factor ***

  forget_ana = forget
  IF (dim_obs_p > 0 .AND. type_forget == 1) THEN
     CALL PDAF_set_forget(step, filterstr, dim_obs_p, dim_ens, HX_p, &
          HXbar_p, obs_p, U_init_obsvar, forget, forget_ana)
  END IF


! ***  Execute Analysis step  ***

  IF (subtype == 0 .OR. subtype == 2 .OR. subtype == 3) THEN
! *** SEIK analysis with forgetting factor better implementation for T ***
     CALL PDAF_seik_analysis_newT(step, dim_p, dim_obs_p, dim_ens, rank, &
          state_p, Uinv, ens_p, state_inc_p, &
          HX_p, HXbar_p, obs_p, forget_ana, U_prodRinvA, &
          screen, incremental, debug, flag)
  ELSE IF (subtype == 1) THEN
! *** SEIK analysis with forgetting factor ***
     CALL PDAF_seik_analysis(step, dim_p, dim_obs_p, dim_ens, rank, &
          state_p, Uinv, ens_p, state_inc_p, &
          HX_p, HXbar_p, obs_p, forget_ana, U_prodRinvA, &
          screen, incremental, debug, flag)
  ELSE IF (subtype == 4) THEN
! *** SEIK analysis with ensemble transformation ***
     CALL PDAF_seik_analysis_trans(step, dim_p, dim_obs_p, dim_ens, rank, &
          state_p, Uinv, ens_p, state_inc_p, &
          HX_p, HXbar_p, obs_p, forget_ana, U_prodRinvA, &
          screen, incremental, type_sqrt, type_trans, Nm1vsN, debug, flag)
  END IF

  CALL PDAF_timeit(3, 'old')

  IF (mype == 0 .AND. screen > 1) THEN
     WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
          'PDAF', '--- update duration:', PDAF_time_temp(3), 's'
  END IF

! *** Resample the state ensemble
  CALL PDAF_timeit(51, 'new')
  CALL PDAF_timeit(7, 'new')

  IF (subtype == 0 .OR. subtype == 2 .OR. subtype == 3) THEN
     CALL PDAF_seik_resample_newT(subtype, dim_p, dim_ens, rank, &
          Uinv, state_p, ens_p, type_sqrt, type_trans, &
          Nm1vsN, screen, flag)
  ELSE IF (subtype == 1) THEN
     CALL PDAF_seik_resample(subtype, dim_p, dim_ens, rank, &
          Uinv, state_p, ens_p, type_sqrt, type_trans, Nm1vsN, &
          screen, flag)
  END IF

  IF (debug>0) THEN
     DO i = 1, dim_ens
        WRITE (*,*) '++ PDAF-debug PDAF_seik_update:', debug, 'ensemble member', i, &
             ' analysis values (1:min(dim_p,6)):', ens_p(1:min(dim_p,6),i)
     END DO
  END IF

  CALL PDAF_timeit(7, 'old')
  CALL PDAF_timeit(51, 'old')
  IF (mype == 0 .AND. screen > 1) THEN
     WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
          'PDAF', '--- resample duration:', PDAF_time_temp(4), 's'
  END IF
#else
  WRITE (*,'(/5x,a/)') &
       '!!! PDAF WARNING: ANALYSIS STEP IS DEACTIVATED BY PDAF_NO_UPDATE !!!'
#endif

! *** Poststep for analysis ensemble ***
  CALL PDAF_timeit(5, 'new')
  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 52a)') 'PDAF Prepoststep ', ('-', i = 1, 52)
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

  ! Deallocate observation arrays
  CALL PDAFobs_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_seik_update -- END'

END SUBROUTINE PDAF_seik_update
