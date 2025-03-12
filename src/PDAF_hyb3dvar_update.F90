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
!> Control analysis update of hyb3DVAR/ESTKF and hyb3DVAR/LESTKF
!!
!! Routines to control the analysis update of hybrid 3DVAR
!! using the ESTKF or LESTKF to update the ensemble pertrubations. 
!!
!! __Revision history:__
!! * 2021-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAF_hyb3dvar_update

CONTAINS
!> Control analysis update of Hyb3DVAR/ESTKF
!!
!! Routine to control the analysis update of hybrid 3DVAR
!! using the ESTKF to update the ensemble pertrubations. 
!!
!! The analysis is performend in the analysis routines for
!! hyb3D-Var and ESTKF called by this routine.
!! In addition, the routine U\_prepoststep is called prior
!! to the analysis and after the resampling to allow the user
!! to access the ensemble information.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2021-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAFhyb3dvar_update_estkf(step, dim_p, dim_obs_p, dim_ens, &
     dim_cvec, dim_cvec_ens, state_p, Ainv, ens_p, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, U_prepoststep, &
     U_cvt_ens, U_cvt_adj_ens, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
     U_init_obsvar, screen, subtype, flag)

  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_parallel, &
       ONLY: mype, dim_ens_l
  USE PDAF_mod_core, &
       ONLY: cnt_maxlag, dim_lag, sens
  USE PDAF_3dvar, &
       ONLY: forget, type_opt, beta_3dvar, debug
  USE PDAFobs, &
       ONLY: PDAFobs_init, PDAFobs_dealloc, type_obs_init, &
       HXbar_p, obs_p
  USE PDAFomi_obs_f, &
       ONLY: omi_n_obstypes => n_obstypes, omi_obs_diag => obs_diag
  USE PDAF_hyb3dvar_analysis_cvt, &
       ONLY: PDAFhyb3dvar_analysis_cvt
  USE PDAF_estkf_update, &
       ONLY: PDAFestkf_update

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step         !< Current time step
  INTEGER, INTENT(in) :: dim_p        !< PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p   !< PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens      !< Size of ensemble
  INTEGER, INTENT(in) :: dim_cvec     !< Size of control vector (parameterized part)
  INTEGER, INTENT(in) :: dim_cvec_ens !< Size of control vector (ensemble part)
  REAL, INTENT(inout) :: state_p(dim_p)              !< PE-local model state
  REAL, INTENT(inout) :: Ainv(dim_ens-1, dim_ens-1)  !< Transform matrix
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)       !< PE-local ensemble matrix
  INTEGER, INTENT(in) :: screen       !< Verbosity flag
  INTEGER, INTENT(in) :: subtype      !< Filter subtype
  INTEGER, INTENT(inout) :: flag      !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, & !< Initialize dimension of observation vector
       U_obs_op, &              !< Observation operator
       U_init_obs, &            !< Initialize observation vector
       U_prepoststep, &         !< User supplied pre/poststep routine
       U_prodRinvA, &           !< Provide product R^-1 A for 3DVAR analysis
       U_cvt_ens, &             !< Apply control vector transform matrix (ensemble)
       U_cvt_adj_ens, &         !< Apply adjoint control vector transform matrix (ensemble var)
       U_cvt, &                 !< Apply control vector transform matrix 
       U_cvt_adj, &             !< Apply adjoint control vector transform matrix
       U_obs_op_lin, &          !< Linearized observation operator
       U_obs_op_adj             !< Adjoint observation operator
  EXTERNAL :: U_init_obsvar     !< Initialize mean observation error variance

! *** local variables ***
  INTEGER :: i, j               ! Counters
  INTEGER :: minusStep          ! Time step counter
  INTEGER :: envar_mode         ! Flag to control step executed in analysis routines
  LOGICAL :: do_init_dim_obs    ! Flag for initializing dim_obs_p in PDAFobs_init
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  REAL, ALLOCATABLE :: state_inc_p(:) ! PE-local analysis increment


! **********************
! *** INITIALIZATION ***
! **********************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_en3dvar_update -- START'

  CALL PDAF_timeit(3, 'new')
  CALL PDAF_timeit(51, 'new')

  ALLOCATE(state_inc_p(dim_p))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_p)


! ***********************************************************
! *** For fixed error space basis compute ensemble states ***
! ***********************************************************

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
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_update:', debug, 'ensemble member', i, &
             ' forecast values (1:min(dim_p,6)):', ens_p(1:min(dim_p,6),i)
     END DO
  END IF


! *****************************************************
! *** Initialize observations and observed ensemble ***
! *****************************************************

  IF (type_obs_init==0 .OR. type_obs_init==2) THEN
     ! This call initializes dim_obs_p, HX_p, HXbar_p, obs_p in the module PDAFobs
     ! It also compute the ensemble mean and stores it in state_p
     CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_p, &
          state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
          screen, debug, .true., .true., .false., .true., .true.)
  END IF

  CALL PDAF_timeit(3, 'old')


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
          screen, debug, .true., do_init_dim_obs, .false., .true., .true.)

     CALL PDAF_timeit(3, 'old')
  END IF


! ***********************
! ***  Analysis step  ***
! ***********************

  IF (mype == 0 .AND. screen > 0) &
       WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)

#ifndef PDAF_NO_UPDATE
  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_update', debug, &
          'Configuration: param_int(3) solver      ', type_opt
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_update', debug, &
          'Configuration: param_int(4) dim_cvec    ', dim_cvec
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_update', debug, &
          'Configuration: param_int(5) dim_cvec_ens', dim_cvec_ens

     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_update', debug, &
          'Configuration: param_real(1) forget     ', forget
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_update', debug, &
          'Configuration: param_real(2) beta_3dvar ', beta_3dvar
  END IF

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 1x, i7, 3x, a)') &
          'PDAF', step, 'Assimilating observations - hybrid 3DVAR with ESTKF'
  END IF

  ! *** Step 1: Hybrid 3DVAR analysis - update state estimate ***

  CALL PDAF_timeit(3, 'new')

  envar_mode = 1
  CALL PDAFhyb3dvar_analysis_cvt(step, dim_p, dim_obs_p, dim_ens, &
       dim_cvec, dim_cvec_ens, beta_3dvar, &
       state_p, ens_p, state_inc_p, HXbar_p, obs_p, &
       U_prodRinvA, U_cvt, U_cvt_adj, &
       U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
       screen, type_opt, debug, flag)

  CALL PDAF_timeit(3, 'old')

  ! *** Step 2: ESTKF - update of ensemble perturbations ***

  envar_mode = 2
  CALL PDAFestkf_update(step, dim_p, dim_obs_p, dim_ens, &
       state_p, Ainv, ens_p, U_init_dim_obs, U_obs_op, &
       U_init_obs, U_prodRinvA, U_init_obsvar, U_prepoststep, &
       screen, 0, envar_mode, &
       dim_lag, sens, cnt_maxlag, flag)

  ! *** Step 3: Add state increment from 3D-Var to ensemble *** 

  CALL PDAF_timeit(3, 'new')

  CALL PDAF_timeit(51, 'new')
  DO j = 1, dim_ens
     DO i = 1, dim_p
        ens_p(i, j) = ens_p(i, j) + state_inc_p(i)
     END DO
  END DO
  CALL PDAF_timeit(51, 'old')

  IF (debug>0) THEN
     DO i = 1, dim_ens
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_update:', debug, 'ensemble member', i, &
             ' analysis values (1:min(dim_p,6)):', ens_p(1:min(dim_p,6),i)
     END DO
  END IF

  CALL PDAF_timeit(3, 'old')

  IF (mype == 0 .AND. screen > 1) THEN
     WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
          'PDAF', '--- duration of hyb3D-Var update:', PDAF_time_temp(3), 's'
  END IF


! ******************************************************
! *** Initialize analysis observed ensemble and mean ***
! ******************************************************

  IF (omi_n_obstypes>0 .AND. omi_obs_diag>0) THEN
     ! This call initializes HX_p, HXbar_p in the module PDAFobs
     ! for the analysis ensemble
     CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_p, &
          state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
          screen, debug, .true., .false., .true., .true., .false.)
  END IF

#else
  WRITE (*,'(/5x,a/)') &
       '!!! PDAF WARNING: ANALYSIS STEP IS DEACTIVATED BY PDAF_NO_UPDATE !!!'
#endif
    

! **************************************
! *** Poststep for analysis ensemble ***
! **************************************

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


! ********************
! *** Finishing up ***
! ********************

  ! Deallocate observation arrays
  CALL PDAFobs_dealloc()

  DEALLOCATE(state_inc_p)

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_update -- END'

END SUBROUTINE PDAFhyb3dvar_update_estkf

!-------------------------------------------------------------------------------
!> Control analysis update of Hyb3DVAR/LESTK
!!
!! Routine to control the analysis update of hybrid 3DVAR
!! using the LESTKF to update the ensemble pertrubations.
!! 
!! The analysis is performend in the analysis routines for
!! Hyb3D-Var and LESTKF called by this routines.
!! In addition, the routine U\_prepoststep is called prior
!! to the analysis and after the resampling to allow the user
!! to access the ensemble information.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2021-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE  PDAFhyb3dvar_update_lestkf(step, dim_p, dim_obs_p, dim_ens, &
     dim_cvec, dim_cvec_ens, state_p, Ainv, ens_p, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, U_prepoststep, &
     U_cvt_ens, U_cvt_adj_ens, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
     U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
     U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, &
     U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
     screen, subtype, flag)

  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_parallel, &
       ONLY: mype, dim_ens_l
  USE PDAF_mod_core, &
       ONLY: cnt_maxlag, dim_lag, sens
  USE PDAF_3dvar, &
       ONLY: type_opt, beta_3dvar, localfilter, forget, debug
  USE PDAFomi, &
       ONLY: PDAFomi_dealloc
  USE PDAFobs, &
       ONLY: PDAFobs_init, PDAFobs_dealloc, type_obs_init, &
       HXbar_p, obs_p
  USE PDAFomi_obs_f, &
       ONLY: omi_n_obstypes => n_obstypes, omi_obs_diag => obs_diag
  USE PDAF_hyb3dvar_analysis_cvt, &
       ONLY: PDAFhyb3dvar_analysis_cvt
  USE PDAF_lestkf_update, &
       ONLY: PDAFlestkf_update

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step         !< Current time step
  INTEGER, INTENT(in) :: dim_p        !< PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p   !< PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens      !< Size of ensemble
  INTEGER, INTENT(in) :: dim_cvec     !< Size of control vector (parameterized part)
  INTEGER, INTENT(in) :: dim_cvec_ens !< Size of control vector (ensemble part)
  REAL, INTENT(inout) :: state_p(dim_p)               !< PE-local model state
  REAL, INTENT(inout) :: Ainv(dim_ens-1, dim_ens-1)   !< Transform matrix for LESKTF
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)        !< PE-local ensemble matrix
  INTEGER, INTENT(in) :: screen       !< Verbosity flag
  INTEGER, INTENT(in) :: subtype      !< Filter subtype
  INTEGER, INTENT(inout) :: flag      !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, & !< Initialize dimension of observation vector
       U_obs_op, &              !< Observation operator
       U_init_obs, &            !< Initialize observation vector
       U_prepoststep, &         !< User supplied pre/poststep routine
       U_prodRinvA, &           !< Provide product R^-1 A for 3DVAR analysis
       U_cvt_ens, &             !< Apply control vector transform matrix (ensemble)
       U_cvt_adj_ens, &         !< Apply adjoint control vector transform matrix (ensemble var)
       U_cvt, &                 !< Apply control vector transform matrix 
       U_cvt_adj, &             !< Apply adjoint control vector transform matrix
       U_obs_op_lin, &          !< Linearized observation operator
       U_obs_op_adj             !< Adjoint observation operator
  EXTERNAL :: U_obs_op_f, &     !< Observation operator
       U_init_n_domains_p, &    !< Provide number of local analysis domains
       U_init_dim_l, &          !< Init state dimension for local ana. domain
       U_init_dim_obs_f, &      !< Initialize dimension of observation vector
       U_init_dim_obs_l, &      !< Initialize dim. of obs. vector for local ana. domain
       U_init_obs_f, &          !< Initialize PE-local observation vector
       U_init_obs_l, &          !< Init. observation vector on local analysis domain
       U_init_obsvar, &         !< Initialize mean observation error variance
       U_init_obsvar_l, &       !< Initialize local mean observation error variance
       U_g2l_state, &           !< Get state on local ana. domain from full state
       U_l2g_state, &           !< Init full state from state on local analysis domain
       U_g2l_obs, &             !< Restrict full obs. vector to local analysis domain
       U_prodRinvA_l            !< Provide product R^-1 A on local analysis domain

! *** local variables ***
  INTEGER :: i, j               ! Counters
  INTEGER :: minusStep          ! Time step counter
  INTEGER :: envar_mode         ! Flag to control step executed in analysis routines
  LOGICAL :: do_init_dim_obs    ! Flag for initializing dim_obs_p in PDAFobs_init
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  REAL, ALLOCATABLE :: state_inc_p(:) ! PE-local analysis increment


! **********************
! *** INITIALIZATION ***
! **********************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_en3dvar_update -- START'

  CALL PDAF_timeit(3, 'new')
  CALL PDAF_timeit(51, 'new')

  ALLOCATE(state_inc_p(dim_p))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_p)


! ***********************************************************
! *** For fixed error space basis compute ensemble states ***
! ***********************************************************

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
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_update:', debug, 'ensemble member', i, &
             ' forecast values (1:min(dim_p,6)):', ens_p(1:min(dim_p,6),i)
     END DO
  END IF


! *****************************************************
! *** Initialize observations and observed ensemble ***
! *****************************************************

  IF (type_obs_init==0 .OR. type_obs_init==2) THEN
     ! This call initializes dim_obs_p, HX_p, HXbar_p, obs_p in the module PDAFobs
     ! It also compute the ensemble mean and stores it in state_p
     CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_p, &
          state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
          screen, debug, .true., .true., .false., .true., .true.)
  END IF
  CALL PDAF_timeit(3, 'old')


! ****************************
! *** Prestep for forecast ***
! ****************************

! *** Prestep for forecast ensemble ***
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
          screen, debug, .true., do_init_dim_obs, .false., .true., .true.)

     CALL PDAF_timeit(3, 'old')
  END IF


! ***********************
! ***  Analysis step  ***
! ***********************

  IF (mype == 0 .AND. screen > 0) &
       WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)

#ifndef PDAF_NO_UPDATE
  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_update', debug, &
          'Configuration: param_int(3) solver      ', type_opt
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_update', debug, &
          'Configuration: param_int(4) dim_cvec    ', dim_cvec
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_update', debug, &
          'Configuration: param_int(5) dim_cvec_ens', dim_cvec_ens

     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_update', debug, &
          'Configuration: param_real(1) forget     ', forget
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_update', debug, &
          'Configuration: param_real(2) beta_3dvar ', beta_3dvar
  END IF

  CALL PDAF_timeit(3, 'new')

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 1x, i7, 3x, a)') &
          'PDAF', step, 'Assimilating observations - hybrid 3DVAR with LESTKF'
  END IF

  ! *** Step 1: Hybrid 3DVAR analysis - update state estimate ***

  envar_mode = 1
  localfilter = 0
  CALL PDAFhyb3dvar_analysis_cvt(step, dim_p, dim_obs_p, dim_ens, &
       dim_cvec, dim_cvec_ens, beta_3dvar, &
       state_p, ens_p, state_inc_p, HXbar_p, obs_p, &
       U_prodRinvA, U_cvt, U_cvt_adj, &
       U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
       screen, type_opt, debug, flag)

  CALL PDAF_timeit(3, 'old')

  ! *** Step 2: LESTKF - update of ensemble perturbations ***

  ! Deallocate observations
  CALL PDAFomi_dealloc()
  CALL PDAFobs_dealloc()

  envar_mode = 2
  localfilter = 1
  CALL PDAFlestkf_update(step, dim_p, dim_obs_p, dim_ens, dim_ens-1, &
       state_p, Ainv, ens_p, &
       U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
       U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, &
       U_g2l_obs, U_init_obsvar, U_init_obsvar_l, U_prepoststep, screen, &
       0, envar_mode, dim_lag, sens, cnt_maxlag, flag)
  localfilter = 0

  ! *** Step 3: Add state increment from 3D-Var to ensemble *** 

  CALL PDAF_timeit(3, 'new')

  CALL PDAF_timeit(51, 'new')
  DO j = 1, dim_ens
     DO i = 1, dim_p
        ens_p(i, j) = ens_p(i, j) + state_inc_p(i)
     END DO
  END DO
  CALL PDAF_timeit(51, 'old')

  IF (debug>0) THEN
     DO i = 1, dim_ens
        WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_update:', debug, 'ensemble member', i, &
             ' analysis values (1:min(dim_p,6)):', ens_p(1:min(dim_p,6),i)
     END DO
  END IF

  CALL PDAF_timeit(3, 'old')

  IF (mype == 0 .AND. screen > 1) THEN
     WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
          'PDAF', '--- duration of hyb3D-Var update:', PDAF_time_temp(3), 's'
  END IF


! ******************************************************
! *** Initialize analysis observed ensemble and mean ***
! ******************************************************

  IF (omi_n_obstypes>0 .AND. omi_obs_diag>0) THEN
     ! This call initializes HX_p, HXbar_p in the module PDAFobs
     ! for the analysis ensemble
     CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_p, &
          state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
          screen, debug, .true., .false., .true., .true., .false.)
  END IF

#else
  WRITE (*,'(/5x,a/)') &
       '!!! PDAF WARNING: ANALYSIS STEP IS DEACTIVATED BY PDAF_NO_UPDATE !!!'
#endif
    

! **************************************
! *** Poststep for analysis ensemble ***
! **************************************

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


! ********************
! *** Finishing up ***
! ********************

  ! Deallocate observation arrays
  CALL PDAFobs_dealloc()

  DEALLOCATE(state_inc_p)

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_update -- END'

END SUBROUTINE PDAFhyb3dvar_update_lestkf

END MODULE PDAF_hyb3dvar_update
