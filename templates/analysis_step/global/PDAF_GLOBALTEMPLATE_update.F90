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
!
!>Control analysis update of GLOBALTEMPLATE
!!
!! This routine prepares the actual analysis update which
!! is computed the in PDAF_analysis_GLOBALTEMPLATE.
!!
!! ADAPTING THE TEMPLATE
!! The template has the typical structure of an ensemble DA method.
!! The only part that is specific to the GLOBALTEMPLATE DA method
!! is the call to PDAF_GLOBALTEMPLATE_analysis and the particular
!! arguments of this subroutine.
!!
!! __Revision history:__
!! * 2024-12 - Lars Nerger - Initial code for template based on ETKF
!! * Later revisions - see repository log
!!
SUBROUTINE  PDAF_GLOBALTEMPLATE_update(step, dim_p, dim_obs_p, dim_ens, &
     state_p, Ainv, ens_p, U_init_dim_obs, U_obs_op, &
     U_init_obs, U_prodRinvA, U_prepoststep, &
     screen, subtype, incremental, &
     dim_lag, sens_p, cnt_maxlag, flag)

  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_mod_filtermpi, &
       ONLY: mype, dim_ens_l
  USE PDAF_mod_filter, &
       ONLY: forget, type_trans, debug

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step          ! Current time step
  INTEGER, INTENT(in) :: dim_p         ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p    ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens       ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)         ! PE-local model state
  REAL, INTENT(inout) :: Ainv(dim_ens, dim_ens) ! Transform matrix
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)  ! PE-local ensemble matrix
  INTEGER, INTENT(in) :: screen        ! Verbosity flag
  INTEGER, INTENT(in) :: subtype       ! Filter subtype
  INTEGER, INTENT(in) :: incremental   ! Control incremental updating
  INTEGER, INTENT(in) :: dim_lag       ! Number of past time instances for smoother
  REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag) ! PE-local smoother ensemble
  INTEGER, INTENT(inout) :: cnt_maxlag ! Count number of past time steps for smoothing
  INTEGER, INTENT(inout) :: flag       ! Status flag

! ** External subroutines ***
! (PDAF-internal names, real names are defined in the call to PDAF)
  ! Routine for ensemble framework
  EXTERNAL :: U_prepoststep        ! User supplied pre/poststep routine
  ! Observation-related routines for analysis step
  EXTERNAL :: U_init_dim_obs, &    ! Initialize dimension of observation vector
       U_obs_op, &                 ! Observation operator
       U_init_obs, &               ! Initialize observation vector
       U_prodRinvA                 ! Provide product R^-1 A

! *** local variables ***
  INTEGER :: i, j                  ! Counters
  INTEGER :: minusStep             ! Time step counter


! ***********************************************************
! *** For fixed error space basis compute ensemble states ***
! ***********************************************************

! +++ TEMPLATE:
! +++ For fixed-ensemble cases (like Ensemble OI) only the
! +++ central state is integrated by the model. Here, we
! +++ then need to add the ensemble perturbations

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


! *************************************
! *** Prestep for forecast ensemble ***
! *************************************

! +++ TEMPLATE:
! +++ The call to the pre/poststep routine for the forecast
! +++ ensemble is standard and should be kept

  CALL PDAF_timeit(5, 'new')
  minusStep = -step  ! Indicate forecast by negative time step number

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a, i7)') 'PDAF', 'Call pre-post routine after forecast; step ', step
  ENDIF
  CALL U_prepoststep(minusStep, dim_p, dim_ens, dim_ens_l, dim_obs_p, &
       state_p, Ainv, ens_p, flag)
  CALL PDAF_timeit(5, 'old')

  IF (mype == 0 .AND. screen > 0) THEN
     IF (screen > 1) &
          WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
          'PDAF', '--- duration of prestep:', PDAF_time_temp(5), 's'
     WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)
  END IF


! ************************
! *** Perform analysis ***
! ************************

#ifndef PDAF_NO_UPDATE
  CALL PDAF_timeit(3, 'new')

! +++ TEMPLATE:
! +++ The call to PDAF_GLOBALTEMPLATE_analysis should be adapted for
! +++ the call-back routines and variables used by the method

  IF (subtype == 0 .OR. subtype == 2) THEN
     ! *** GLOBALTEMPLATE analysis step
     CALL PDAF_GLOBALTEMPLATE_analysis(step, dim_p, dim_obs_p, dim_ens, &
          state_p, Ainv, ens_p, forget, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
          screen, flag)
  END IF

! +++ TEMPLATE:
! +++ The smoother routine is generic. It can be used
! +++ as long as sens_l(:,lag) * (forget*Ainv + diag(invdimens))
! +++ yields the smoother update for the ensemble at lag 'lag'. 

  ! *** Perform smoothing of past ensembles ***
  CALL PDAF_timeit(51, 'new')
  CALL PDAF_smoother(dim_p, dim_ens, dim_lag, Ainv, sens_p, &
       cnt_maxlag, forget, screen)
  CALL PDAF_timeit(51, 'old')

  CALL PDAF_timeit(3, 'old')

  IF (mype == 0 .AND. screen > 1) THEN
     WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
          'PDAF', '--- update duration:', PDAF_time_temp(3), 's'
  END IF

#else
  WRITE (*,'(/5x,a/)') &
       '!!! PDAF WARNING: ANALYSIS STEP IS DEACTIVATED BY PDAF_NO_UPDATE !!!'
#endif

    
! **************************************
! *** Poststep for analysis ensemble ***
! **************************************

! +++ TEMPLATE:
! +++ The call to the pre/poststep routine for the analysis
! +++ ensemble is standard and should be kept

  CALL PDAF_timeit(5, 'new')
  IF (mype == 0 .AND. screen > 0) THEN
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

END SUBROUTINE PDAF_GLOBALTEMPLATE_update
