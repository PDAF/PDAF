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
!$Id: PDAF-D_netf_update.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_netf_update --- Control analysis update of the NETF
!
! !INTERFACE:
SUBROUTINE  PDAF_netf_update(step, dim_p, dim_obs_p, dim_ens, &
     state_p, Uinv, ens_p, type_forget, forget, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_likelihood, &
     U_prepoststep, screen, subtype, &
     dim_lag, sens_p, cnt_maxlag, flag)

! !DESCRIPTION:
! Routine to control the analysis update of the NETF.
! 
! The analysis with ensemble transofrmation is performed by 
! calling PDAF\_netf\_analysis.
! In addition, the routine U\_prepoststep is called prior
! to the analysis and after the resampling to allow the user
! to access the ensemble information.
!
! Variant for NETF with domain decompostion.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2009-07 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filter, &
       ONLY: type_trans
  USE PDAF_mod_filtermpi, &
       ONLY: mype, dim_ens_l

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
  INTEGER, INTENT(in) :: dim_p       ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p  ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens     ! Size of ensemble
  INTEGER, INTENT(in) :: type_forget ! Type of forgetting factor
  REAL, INTENT(in)    :: forget      ! Forgetting factor
  REAL, INTENT(inout) :: state_p(dim_p)        ! PE-local model state
  REAL, INTENT(inout) :: Uinv(dim_ens, dim_ens)! Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) ! PE-local ensemble matrix
  INTEGER, INTENT(in) :: screen      ! Verbosity flag
  INTEGER, INTENT(in) :: subtype     ! Filter subtype
  INTEGER, INTENT(in) :: dim_lag     ! Number of past time instances for smoother
  REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag) ! PE-local smoother ensemble
  INTEGER, INTENT(inout) :: cnt_maxlag ! Count number of past time steps for smoothing
  INTEGER, INTENT(inout) :: flag     ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, & ! Initialize dimension of observation vector
       U_obs_op, &              ! Observation operator
       U_init_obs, &            ! Initialize observation vector
       U_prepoststep, &         ! User supplied pre/poststep routine
       U_likelihood             ! Compute observation likelihood for an ensemble member

! !CALLING SEQUENCE:
! Called by: PDAF_put_state_netf
! Calls: U_prepoststep
! Calls: PDAF_generate_rndmat
! Calls: PDAF_netf_smootherT
! Calls: PDAF_netf_analysis
! Calls: PDAF_smoother_netf
! Calls: PDAF_timeit
! Calls: PDAF_time_temp
!EOP

! *** local variables ***
  INTEGER :: i, j      ! Counters
  INTEGER :: minusStep ! Time step counter
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  REAL, ALLOCATABLE :: TA_noinfl(:,:) ! TA for smoother (without inflation)
  REAL, ALLOCATABLE :: rndmat(:,:)    ! Orthogonal random matrix


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

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 1x, i7, 3x, a)') &
          'PDAF', step, 'Assimilating observations - NETF'
  END IF

  CALL PDAF_timeit(3, 'new')

  ! Generate orthogonal random matrix with eigenvector (1,...,1)^T
  ALLOCATE(rndmat(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens*dim_ens)

  IF (type_trans==0) THEN
     IF (screen > 0 .AND. mype == 0) &
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- Initialize random transformation'
     CALL PDAF_generate_rndmat(dim_ens, rndmat, 2)
  ELSE
     IF (screen > 0 .AND. mype == 0) &
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- Initialize deterministic transformation'
     rndmat = 0.0
     DO i = 1, dim_ens
        rndmat(i,i) = 1.0
     END DO
  END IF

  IF (dim_lag > 0) THEN
     ! Compute transform matrix for smoother

     ALLOCATE(TA_noinfl(dim_ens, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens*dim_ens)

     CALL PDAF_netf_smootherT(step, dim_p, dim_obs_p, dim_ens, &
          ens_p, rndmat, TA_noinfl,  &
          U_init_dim_obs, U_obs_op, U_init_obs, U_likelihood, &
          screen, flag)

  END IF

  ! *** NETF analysis ***
  CALL PDAF_netf_analysis(step, dim_p, dim_obs_p, dim_ens, &
       state_p, ens_p, rndmat, Uinv, forget, &
       U_init_dim_obs, U_obs_op, U_init_obs, U_likelihood, &
       screen, type_forget, flag)

  ! *** Perform smoothing of past ensembles ***
  IF (dim_lag > 0) THEN
     CALL PDAF_smoother_netf(dim_p, dim_ens, dim_lag, TA_noinfl, sens_p, &
          cnt_maxlag, forget, screen)

     DEALLOCATE(TA_noinfl)
  END IF

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

END SUBROUTINE PDAF_netf_update
