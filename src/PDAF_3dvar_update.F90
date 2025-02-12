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
!> Control analysis update of 3DVAR
!!
!! The analysis is performend PDAF_3dvar_analysis_cvt.
!! In addition, the routine U\_prepoststep is called prior
!! to the analysis and after the resampling to allow the user
!! to access the ensemble information.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2021-03 - Lars Nerger - Initial code
!! * Later revisions - see svn log
!!
MODULE PDAF_3dvar_update

CONTAINS
SUBROUTINE  PDAF3dvar_update(step, dim_p, dim_obs_p, dim_ens, &
     dim_cvec, state_p, Ainv, ens_p, state_inc_p, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, U_prepoststep, &
     U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
     screen, subtype, incremental, flag)

  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_mod_filtermpi, &
       ONLY: mype, dim_ens_l
  USE PDAF_mod_filter, &
       ONLY: debug
  USE PDAF_3dvar, &
       ONLY: type_opt
  USE PDAFobs, &
       ONLY: PDAFobs_init, PDAFobs_dealloc, type_obs_init, &
       HXbar_p, obs_p
  USE PDAF_3dvar_analysis_cvt, &
       ONLY: PDAF3dvar_analysis_cvt

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step        !< Current time step
  INTEGER, INTENT(in) :: dim_p       !< PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p  !< PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens     !< Size of ensemble
  INTEGER, INTENT(in) :: dim_cvec    !< Size of control vector (parameterized part)
  REAL, INTENT(inout) :: state_p(dim_p)        !< PE-local model state
  REAL, INTENT(inout) :: Ainv(1, 1)            !< Not used in 3D-Var
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) !< PE-local ensemble matrix
  REAL, INTENT(inout) :: state_inc_p(dim_p)    !< PE-local state analysis increment
  INTEGER, INTENT(in) :: screen      !< Verbosity flag
  INTEGER, INTENT(in) :: subtype     !< Filter subtype
  INTEGER, INTENT(in) :: incremental !< Control incremental updating
  INTEGER, INTENT(inout) :: flag     !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, &      !< Initialize dimension of observation vector
       U_obs_op, &                   !< Observation operator
       U_init_obs, &                 !< Initialize observation vector
       U_prepoststep, &              !< User supplied pre/poststep routine
       U_prodRinvA, &                !< Provide product R^-1 A for 3DVAR analysis
       U_cvt, &                      !< Apply control vector transform matrix 
       U_cvt_adj, &                  !< Apply adjoint control vector transform matrix
       U_obs_op_lin, &               !< Linearized observation operator
       U_obs_op_adj                  !< Adjoint observation operator

! *** local variables ***
  INTEGER :: i, j                    ! Counters
  INTEGER :: minusStep               ! Time step counter
  LOGICAL :: do_init_dim_obs         ! Flag for initializing dim_obs_p in PDAFobs_init


! ***********************************************************
! *** For fixed error space basis compute ensemble states ***
! ***********************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_update -- START'

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

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_3dvar_update:', debug, &
       ' forecast state (1:min(dim_p,6)):', ens_p(1:min(dim_p,6),1)


! **********************************************
! *** Initialize state_p from ensemble array ***
! **********************************************

  state_p(:) = ens_p(:, 1)

  CALL PDAF_timeit(51, 'old')


! *****************************************************
! *** Initialize observations and observed ensemble ***
! *****************************************************

  IF (type_obs_init==0 .OR. type_obs_init==2) THEN
     ! This call initializes dim_obs_p, HX_p, HXbar_p, obs_p in the module PDAFobs
     ! It also compute the ensemble mean and stores it in state_p
     CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_p, &
          state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
          screen, debug, .false., .true., .false., .true., .true.)
  END IF
  CALL PDAF_timeit(3, 'new')


! ****************************
! *** Prestep for forecast ***
! ****************************

  CALL PDAF_timeit(5, 'new')
  minusStep = -step  ! Indicate forecast by negative time step number
  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 52a)') 'PDAF Prepoststep ', ('-', i = 1, 52)
     WRITE (*, '(a, 5x, a, i7)') 'PDAF', 'Call pre-post routine after forecast; step ', step
  ENDIF
  CALL U_prepoststep(minusStep, dim_p, dim_ens, dim_ens_l, dim_obs_p, &
       state_p, Ainv, ens_p, flag)
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
          screen, debug, .false., do_init_dim_obs, .false., .true., .true.)

     CALL PDAF_timeit(3, 'new')
  END IF


! ***********************
! ***  Analysis step  ***
! ***********************

  IF (mype == 0 .AND. screen > 0) &
       WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)

#ifndef PDAF_NO_UPDATE
  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_update', debug, &
          'Configuration: param_int(3) solver      ', type_opt
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_update', debug, &
          'Configuration: param_int(4) dim_cvec    ', dim_cvec
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_update', debug, &
          'Configuration: param_int(5) -not used-'
  END IF

  CALL PDAF_timeit(3, 'new')

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 1x, i7, 3x, a)') &
          'PDAF', step, 'Assimilating observations - 3DVAR incremental, transformed'
  END IF

  ! *** 3DVAR analysis ***
  CALL PDAF3dvar_analysis_cvt(step, dim_p, dim_obs_p, dim_cvec, &
       state_p, state_inc_p, HXbar_p, obs_p, U_prodRinvA, &
       U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
       screen, incremental, type_opt, debug, flag)

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_3dvar_update:', debug, &
          ' analysis state (1:min(dim_p,6)):', state_p(1:min(dim_p,6))
  END IF

  ! Initialize ens_p for intregration with PDAF
  ens_p(:, 1) = state_p(:)


  CALL PDAF_timeit(3, 'old')

  IF (mype == 0 .AND. screen > 1) THEN
     WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
          'PDAF', '--- update duration:', PDAF_time_temp(3), 's'
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
       state_p, Ainv, ens_p, flag)
  CALL PDAF_timeit(5, 'old')
  
  IF (mype == 0 .AND. screen > 0) THEN
     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- duration of poststep:', PDAF_time_temp(5), 's'
     END IF
     WRITE (*, '(a, 55a)') 'PDAF Forecast ', ('-', i = 1, 55)
  END IF

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_3dvar_update -- END'


! ********************
! *** Finishing up ***
! ********************

  ! Deallocate observation arrays
  CALL PDAFobs_dealloc()

END SUBROUTINE PDAF3dvar_update

END MODULE PDAF_3dvar_update
