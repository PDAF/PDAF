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
!> Control analysis update of the ESTKF
!!
!! Routine to control the analysis update of the ESTKF.
!! 
!! The analysis and ensemble tranformation are performed by
!! calling PDAF\_estkf\_analysis. In addition, the routine
!! U\_prepoststep is called prior to the analysis and after
!! the resampling to allow the user to access the ensemble
!! information.
!!
!! Variant for ESTKF with domain decompostion.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2011-09 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAF_estkf_update

CONTAINS
SUBROUTINE  PDAFestkf_update(step, dim_p, dim_obs_p, dim_ens, &
     state_p, Ainv, ens_p, U_init_dim_obs, U_obs_op, &
     U_init_obs, U_prodRinvA, U_init_obsvar, U_prepoststep, &
     screen, subtype, envar_mode, &
     dim_lag, sens_p, cnt_maxlag, flag)

  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_parallel, &
       ONLY: mype, dim_ens_l
  USE PDAF_estkf, &
       ONLY: localfilter, debug, forget, type_forget, &
       type_trans, type_sqrt
  USE PDAF_analysis_utils, &
       ONLY: PDAF_set_forget
  USE PDAFobs, &
       ONLY: PDAFobs_init, PDAFobs_dealloc, type_obs_init, &
       observe_ens, HX_p, HXbar_p, obs_p
  USE PDAF_estkf_analysis, &
       ONLY: PDAF_estkf_ana
  USE PDAF_estkf_analysis_fixed, &
       ONLY: PDAF_estkf_ana_fixed
  USE PDAF_smoother, &
       ONLY: PDAF_smoothing

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step        !< Current time step
  INTEGER, INTENT(in) :: dim_p       !< PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p  !< PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens     !< Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)        !< PE-local model state
  REAL, INTENT(inout) :: Ainv(dim_ens-1, dim_ens-1)      !< Inverse of transform matrix A
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) !< PE-local ensemble matrix
  INTEGER, INTENT(in) :: screen      !< Verbosity flag
  INTEGER, INTENT(in) :: subtype     !< Filter subtype
  INTEGER, INTENT(in) :: envar_mode  !< Flag whether routine is called from 3DVar for special functionality
  INTEGER, INTENT(in) :: dim_lag     !< Number of past time instances for smoother
  REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag) !< PE-local smoother ensemble
  INTEGER, INTENT(inout) :: cnt_maxlag !< Count number of past time steps for smoothing
  INTEGER, INTENT(inout) :: flag     !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, & !< Initialize dimension of observation vector
       U_obs_op, &              !< Observation operator
       U_init_obs, &            !< Initialize observation vector
       U_init_obsvar, &         !< Initialize mean observation error variance
       U_prepoststep, &         !< User supplied pre/poststep routine
       U_prodRinvA              !< Provide product R^-1 A for ESTKF analysis

! *** local variables ***
  INTEGER :: i, j                   ! Counters
  INTEGER :: minusStep              ! Time step counter
  REAL :: forget_ana                ! Forgetting factor actually used in analysis
  INTEGER, SAVE :: allocflag = 0    ! Flag whether first time allocation is done
  REAL, ALLOCATABLE :: TA(:,:)      ! Ensemble transform matrix
  LOGICAL :: do_init_dim_obs        ! Flag for initializing dim_obs_p in PDAFobs_init


! ***********************************************************
! *** For fixed error space basis compute ensemble states ***
! ***********************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_estkf_update -- START'

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

  IF (debug>0 .AND. envar_mode<2) THEN
     DO i = 1, dim_ens
        WRITE (*,*) '++ PDAF-debug PDAF_estkf_update:', debug, 'ensemble member', i, &
             ' forecast values (1:min(dim_p,6)):', ens_p(1:min(dim_p,6),i)
     END DO
  END IF
  CALL PDAF_timeit(51, 'old')


! *****************************************************
! *** Initialize observations and observed ensemble ***
! *****************************************************

  IF ((type_obs_init==0 .OR. type_obs_init==2) .AND. envar_mode<2) THEN
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

! *** Prestep for forecast ensemble ***
  IF (envar_mode < 2) THEN
     ! Do prepoststep only if ESTKF is not used in hybrid 3D-Var (envar_mode==2)

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
  END IF


! *****************************************************
! *** Initialize observations and observed ensemble ***
! *****************************************************

  CALL PDAF_timeit(3, 'new')

  IF (envar_mode < 2) THEN
     ! Normal case of direct use of LESTKF
     IF (type_obs_init>0) THEN
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
     END IF
  ELSE
     ! When ESTKF is used in En3DVar or Hyb3DVar
     CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_p, &
          state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
          screen, debug, .true., .false., .true., .true., .false.)
  END IF

  CALL PDAF_timeit(3, 'old')


! ***********************
! ***  Analysis step  ***
! ***********************

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)
     IF (envar_mode<2) THEN
        IF (subtype == 0 .OR. subtype == 2) THEN
           WRITE (*, '(a, i7, 3x, a)') &
                'PDAF ', step, 'Assimilating observations - ESTKF'
        ELSE
           WRITE (*, '(a, i7, 3x, a)') &
                'PDAF ', step, 'Assimilating observations - ESTKF with fixed ensemble'
        END IF
     ELSE
        WRITE (*, '(a, 5x, a)') &
             'PDAF', 'Step 2: Update ensemble perturbations - ESTKF analysis'
     END IF
  END IF

#ifndef PDAF_NO_UPDATE
  CALL PDAF_timeit(3, 'new')

  IF (debug>0) THEN
     IF (envar_mode<2) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_estkf_update', debug, &
             'Configuration: param_int(3) dim_lag     ', dim_lag
        WRITE (*,*) '++ PDAF-debug PDAF_estkf_update', debug, &
             'Configuration: param_int(4) -not used-  '
        WRITE (*,*) '++ PDAF-debug PDAF_estkf_update', debug, &
             'Configuration: param_int(5) type_forget ', type_forget
        WRITE (*,*) '++ PDAF-debug PDAF_estkf_update', debug, &
             'Configuration: param_int(6) type_trans  ', type_trans
        WRITE (*,*) '++ PDAF-debug PDAF_estkf_update', debug, &
             'Configuration: param_int(7) type_sqrt   ', type_sqrt
        WRITE (*,*) '++ PDAF-debug PDAF_estkf_update', debug, &
             'Configuration: param_int(8) observe_ens           ', observe_ens

        WRITE (*,*) '++ PDAF-debug PDAF_estkf_update', debug, &
             'Configuration: param_real(1) forget     ', forget
     ELSE
        WRITE (*,*) '++ PDAF-debug PDAF_estkf_update', debug, &
             'execute LESTKF analysis with default parameters'
     END IF
  END IF


! *** Compute adaptive forgetting factor ***

  forget_ana = forget
  IF (type_forget == 1) THEN
     CALL PDAF_set_forget(step, localfilter, dim_obs_p, dim_ens, HX_p, &
          HXbar_p, obs_p, U_init_obsvar, forget, forget_ana, &
          screen)
  ELSE
     IF (mype == 0 .AND. screen > 0) THEN
        WRITE (*, '(a, 5x, a, F7.2)') &
             'PDAF', '--- apply multiplicative inflation with fixed forget', forget
     END IF
  END IF

! *** Allocate ensemble transform matrix ***
  ALLOCATE(TA(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)
  TA = 0.0

! *** ESTKF analysis ***
  IF (subtype == 0 .OR. subtype == 2) THEN
     ! Analysis with ensemble transformation
     CALL PDAF_estkf_ana(step, dim_p, dim_obs_p, dim_ens, dim_ens-1, &
          state_p, Ainv, ens_p, HX_p, HXbar_p, obs_p, &
          forget_ana, U_prodRinvA, screen, envar_mode, &
          type_sqrt, type_trans, TA, debug, flag)
  ELSE
     ! Analysis with state update but no ensemble transformation
     CALL PDAF_estkf_ana_fixed(step, dim_p, dim_obs_p, dim_ens, dim_ens-1, &
          state_p, Ainv, ens_p, HX_p, HXbar_p, obs_p, &
          forget_ana, U_prodRinvA, &
          screen, type_sqrt, debug, flag)
  END IF

  IF (debug>0) THEN
     IF (envar_mode<2) THEN
        DO i = 1, dim_ens
           WRITE (*,*) '++ PDAF-debug PDAF_estkf_update:', debug, 'ensemble member', i, &
                ' analysis values (1:min(dim_p,6)):', ens_p(1:min(dim_p,6),i)
        END DO
     END IF
  END IF


  ! *** Perform smoothing of past ensembles ***
  IF (dim_lag>0) THEN
     CALL PDAF_timeit(15, 'new')
     CALL PDAF_timeit(51, 'new')
     CALL PDAF_smoothing(dim_p, dim_ens, dim_lag, TA, sens_p, &
          cnt_maxlag, forget_ana, screen)
     CALL PDAF_timeit(51, 'old')
     CALL PDAF_timeit(15, 'old')
  END IF

  CALL PDAF_timeit(3, 'old')

  IF (mype == 0 .AND. screen > 1 .AND. envar_mode < 2) THEN
     WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
          'PDAF', '--- update duration:', PDAF_time_temp(3), 's'
  END IF

  DEALLOCATE(TA)

#else
  WRITE (*,'(/5x,a/)') &
       '!!! PDAF WARNING: ANALYSIS STEP IS DEACTIVATED BY PDAF_NO_UPDATE !!!'
#endif

! *** Poststep for analysis ensemble ***
  IF (envar_mode < 2) THEN
     ! Do prepoststep only if ESTKF is not used in hybrid 3D-Var (envar_mode==2)

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
  END IF


! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1

  ! Deallocate observation arrays
  CALL PDAFobs_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_estkf_update -- END'

END SUBROUTINE PDAFestkf_update

END MODULE PDAF_estkf_update
