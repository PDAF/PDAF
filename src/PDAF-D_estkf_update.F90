! Copyright (c) 2004-2016 Lars Nerger
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
!$Id: PDAF-D_estkf_update.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_estkf_update --- Control analysis update of the ESTKF
!
! !INTERFACE:
SUBROUTINE  PDAF_estkf_update(step, dim_p, dim_obs_p, dim_ens, rank, &
     state_p, Ainv, ens_p, state_inc_p, forget, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, U_init_obsvar, &
     U_prepoststep, screen, subtype, incremental, type_forget, &
     type_sqrt, dim_lag, sens_p, cnt_maxlag, flag)

! !DESCRIPTION:
! Routine to control the analysis update of the ESTKF.
! 
! The analysis and ensemble tranformation are performed by
! calling PDAF\_estkf\_analysis. In addition, the routine
! U\_prepoststep is called prior to the analysis and after
! the resampling to allow the user to access the ensemble
! information.
!
! Variant for ESTKF with domain decompostion.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2011-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: mype, dim_ens_l

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
  INTEGER, INTENT(in) :: dim_p       ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p  ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens     ! Size of ensemble
  INTEGER, INTENT(in) :: rank        ! Rank of initial covariance matrix
  REAL, INTENT(in)    :: forget      ! Forgetting factor
  REAL, INTENT(inout) :: state_p(dim_p)        ! PE-local model state
  REAL, INTENT(inout) :: Ainv(rank, rank)      ! Inverse of transform matrix A
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) ! PE-local ensemble matrix
  REAL, INTENT(inout) :: state_inc_p(dim_p)    ! PE-local state analysis increment
  INTEGER, INTENT(in) :: screen      ! Verbosity flag
  INTEGER, INTENT(in) :: subtype     ! Filter subtype
  INTEGER, INTENT(in) :: incremental ! Control incremental updating
  INTEGER, INTENT(in) :: type_forget ! Type of forgetting factor
  INTEGER, INTENT(in) :: type_sqrt   ! Type of square-root of A
  INTEGER, INTENT(in) :: dim_lag     ! Number of past time instances for smoother
  REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag) ! PE-local smoother ensemble
  INTEGER, INTENT(inout) :: cnt_maxlag ! Count number of past time steps for smoothing
  INTEGER, INTENT(inout) :: flag     ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, & ! Initialize dimension of observation vector
       U_obs_op, &              ! Observation operator
       U_init_obs, &            ! Initialize observation vector
       U_init_obsvar, &         ! Initialize mean observation error variance
       U_prepoststep, &         ! User supplied pre/poststep routine
       U_prodRinvA              ! Provide product R^-1 A for ESTKF analysis

! !CALLING SEQUENCE:
! Called by: PDAF_put_state_estkf
! Calls: U_prepoststep
! Calls: PDAF_estkf_analysis
! Calls: PDAF_timeit
! Calls: PDAF_time_temp
!EOP

! *** local variables ***
  INTEGER :: i, j      ! Counters
  INTEGER :: minusStep ! Time step counter
  INTEGER, SAVE :: allocflag = 0     ! Flag whether first time allocation is done
  REAL :: forget_ana   ! Forgetting factor actually used in analysis
  REAL, ALLOCATABLE :: TA(:,:) ! Ensemble transform matrix


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
       state_p, Ainv, ens_p, flag)
  CALL PDAF_timeit(5, 'old')

  IF (mype == 0 .AND. screen > 0) THEN
     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- duration of prestep:', PDAF_time_temp(5), 's'
     END IF
     WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)
  END IF

#ifndef PDAF_NO_UPDATE
  ! Allocate ensemble transform matrix
  ALLOCATE(TA(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)
  TA = 0.0

  CALL PDAF_timeit(3, 'new')
! *** ESTKF analysis ***
  IF (subtype == 0 .OR. subtype == 2 .OR. subtype == 5) THEN
     ! Analysis with ensemble transformation
     CALL PDAF_estkf_analysis(step, dim_p, dim_obs_p, dim_ens, rank, &
          state_p, Ainv, ens_p, state_inc_p, forget, forget_ana, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_init_obsvar, U_prodRinvA, &
          screen, incremental, type_forget, type_sqrt, TA, flag)
  ELSE
     ! Analysis with state update but no ensemble transformation
     CALL PDAF_estkf_analysis_fixed(step, dim_p, dim_obs_p, dim_ens, rank, &
          state_p, Ainv, ens_p, state_inc_p, forget, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_init_obsvar, U_prodRinvA, &
          screen, incremental, type_forget, type_sqrt, flag)
  END IF

  ! *** Perform smoothing of past ensembles ***
  CALL PDAF_smoother(dim_p, dim_ens, dim_lag, TA, sens_p, &
       cnt_maxlag, forget_ana, screen)

  CALL PDAF_timeit(3, 'old')

  IF (mype == 0 .AND. screen > 1) THEN
     WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
          'PDAF', '--- update duration:', PDAF_time_temp(3), 's'
  END IF

  DEALLOCATE(TA)

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

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_estkf_update
