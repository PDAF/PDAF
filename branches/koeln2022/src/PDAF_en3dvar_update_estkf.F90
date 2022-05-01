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
! !ROUTINE: PDAF_en3dvar_update_estkf --- Control analysis update of En3DVAR/ESTKF
!
! !INTERFACE:
SUBROUTINE  PDAF_en3dvar_update_estkf(step, dim_p, dim_obs_p, dim_ens, &
     dim_cvec_ens, state_p, Uinv, ens_p, state_inc_p, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, U_prepoststep, &
     U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
     U_init_obsvar, &
     screen, subtype, incremental, type_forget, type_opt, flag)

! !DESCRIPTION:
! Routine to control the analysis update of ensemble 3DVAR
! using the ESTKF to update the ensemble pertrubations.
! 
! The analysis is performend in the analysis routines for
! En3D-Var and ESTKF called by this routines.
! In addition, the routine U\_prepoststep is called prior
! to the analysis and after the resampling to allow the user
! to access the ensemble information.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2021-03 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_mod_filtermpi, &
       ONLY: mype, dim_ens_l
  USE PDAF_mod_filter, &
       ONLY: cnt_maxlag, dim_lag, sens, type_sqrt, forget

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
  INTEGER, INTENT(in) :: dim_p       ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p  ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens     ! Size of ensemble
  INTEGER, INTENT(in) :: dim_cvec_ens ! Size of control vector (ensemble part)
  REAL, INTENT(inout) :: state_p(dim_p)        ! PE-local model state
  REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1)  ! Transform matrix
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) ! PE-local ensemble matrix
  REAL, INTENT(inout) :: state_inc_p(dim_p)    ! PE-local state analysis increment
  INTEGER, INTENT(in) :: screen      ! Verbosity flag
  INTEGER, INTENT(in) :: subtype     ! Filter subtype
  INTEGER, INTENT(in) :: incremental ! Control incremental updating
  INTEGER, INTENT(in) :: type_forget ! Type of forgetting factor
  INTEGER, INTENT(in) :: type_opt    ! Type of minimizer for 3DVar
  INTEGER, INTENT(inout) :: flag     ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, & ! Initialize dimension of observation vector
       U_obs_op, &              ! Observation operator
       U_init_obs, &            ! Initialize observation vector
       U_prepoststep, &         ! User supplied pre/poststep routine
       U_prodRinvA, &           ! Provide product R^-1 A for 3DVAR analysis
       U_cvt_ens, &             ! Apply control vector transform matrix (ensemble)
       U_cvt_adj_ens, &         ! Apply adjoint control vector transform matrix (ensemble var)
       U_obs_op_lin, &          ! Linearized observation operator
       U_obs_op_adj             ! Adjoint observation operator
  EXTERNAL :: U_init_obsvar     ! Initialize mean observation error variance

! !CALLING SEQUENCE:
! Called by: PDAF_put_state_3dvar
! Calls: U_prepoststep
! Calls: PDAF_3dvar_analysis
! Calls: PDAF_timeit
! Calls: PDAF_time_temp
!EOP

! *** local variables ***
  INTEGER :: i, j      ! Counters
  INTEGER :: minusStep ! Time step counter
  INTEGER :: incremental_tmp

! ***********************************************************
! *** For fixed error space basis compute ensemble states ***
! ***********************************************************

  CALL PDAF_timeit(51, 'new')

  fixed_basis: IF (subtype == 2 .OR. subtype == 3) THEN
     ! *** Add mean/central state to ensemble members ***
     DO j = 1, dim_ens
        DO i = 1, dim_p
           ens_p(i, j) = ens_p(i, j) + state_p(i)
        END DO
     END DO
  END IF fixed_basis

  CALL PDAF_timeit(51, 'old')


! **********************
! ***  Update phase  ***
! **********************

! *** Prestep for forecast ensemble ***
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
     WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)
  END IF

#ifndef PDAF_NO_UPDATE
  CALL PDAF_timeit(3, 'new')

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 1x, i7, 3x, a)') &
          'PDAF', step, 'Assimilating observations - ensemble 3DVAR with ESTKF'
  END IF

  ! *** Step 1: Ensemble 3DVAR analysis - update state estimate ***

  incremental_tmp = 1
  CALL PDAF_en3dvar_analysis_cvt(step, dim_p, dim_obs_p, dim_ens, dim_cvec_ens, &
       state_p, ens_p, state_inc_p, &
       U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
       U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
       screen, incremental_tmp, type_opt, flag)

  ! *** Step 2: ESTKF - update of ensemble perturbations ***

  incremental_tmp = 2
  CALL PDAF_estkf_update(step, dim_p, dim_obs_p, dim_ens, dim_ens-1, &
       state_p, Uinv, ens_p, state_inc_p, &
       U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, U_init_obsvar, &
       U_prepoststep, screen, 0, incremental_tmp, type_forget, &
       type_sqrt, dim_lag, sens, cnt_maxlag, flag)

  ! *** Step 3: Add state increment from 3D-Var to ensemble *** 

  IF (incremental==0) THEN
     CALL PDAF_timeit(51, 'new')
     DO j = 1, dim_ens
        DO i = 1, dim_p
           ens_p(i, j) = ens_p(i, j) + state_inc_p(i)
        END DO
     END DO
     CALL PDAF_timeit(51, 'old')
  END IF


  CALL PDAF_timeit(3, 'old')

  IF (mype == 0 .AND. screen > 1) THEN
     WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
          'PDAF', '--- duration of En3D-Var update:', PDAF_time_temp(3), 's'
  END IF

#else
  WRITE (*,'(/5x,a/)') &
       '!!! PDAF WARNING: ANALYSIS STEP IS DEACTIVATED BY PDAF_NO_UPDATE !!!'
#endif
    
! *** Poststep for analysis ensemble ***
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

END SUBROUTINE PDAF_en3dvar_update_estkf
