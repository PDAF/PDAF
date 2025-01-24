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
! !ROUTINE: PDAF_netf_update --- Control analysis update of the NETF
!
! !INTERFACE:
SUBROUTINE  PDAF_netf_update(step, dim_p, dim_obs_p, dim_ens, &
     state_p, Uinv, ens_p, type_forget, noise_type, noise_amp, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_likelihood, &
     U_prepoststep, screen, subtype, &
     dim_lag, sens_p, cnt_maxlag, flag)

! !DESCRIPTION:
! Routine to control the analysis update of the NETF.
! 
! The analysis with ensemble transformation is performed by 
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
       ONLY: type_trans, type_winf, limit_winf, forget, debug
  USE PDAF_mod_filtermpi, &
       ONLY: mype, dim_ens_l
  USE PDAFobs, &
       ONLY: PDAFobs_initialize, PDAFobs_dealloc, HX_p, obs_p

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
  INTEGER, INTENT(in) :: dim_p       ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p  ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens     ! Size of ensemble
  INTEGER, INTENT(in) :: type_forget ! Type of forgetting factor
  INTEGER, INTENT(in) :: noise_type  ! Type of pertubing noise
  REAL, INTENT(in) :: noise_amp      ! Amplitude of noise
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
  INTEGER :: i, j, member, row        ! Counters
  INTEGER :: minusStep                ! Time step counter
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  REAL :: invdimens                   ! Inverse global ensemble size
  REAL, ALLOCATABLE :: TA_noinfl(:,:) ! TA for smoother (without inflation)
  REAL, ALLOCATABLE :: rndmat(:,:)    ! Orthogonal random matrix


! ***********************************************************
! *** For fixed error space basis compute ensemble states ***
! ***********************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_netf_update -- START'

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
        WRITE (*,*) '++ PDAF-debug PDAF_netf_update:', debug, 'ensemble member', i, &
             ' forecast values (1:min(dim_p,6)):', ens_p(1:min(dim_p,6),i)
     END DO
  END IF


! ***********************************
! *** Compute mean forecast state ***
! ***********************************

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
  END IF


! ************************
! *** Inflate ensemble ***
! ************************

  IF (type_forget==0 ) THEN
     CALL PDAF_timeit(51, 'new')
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_netf_update', debug, &
          'Inflate forecast ensemble'

     CALL PDAF_timeit(34, 'new') ! Apply forgetting factor
     CALL PDAF_inflate_ens(dim_p, dim_ens, state_p, ens_p, forget)
     CALL PDAF_timeit(34, 'old')
     CALL PDAF_timeit(51, 'old')
  ENDIF
  

! *****************************************************
! *** Initialize observations and observed ensemble ***
! *****************************************************

  CALL PDAFobs_initialize(step, dim_p, dim_ens, dim_obs_p, &
       state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
       screen, debug, .true., .true., .false., .true.)


! ***********************
! ***  Analysis step  ***
! ***********************

  IF (mype == 0 .AND. screen > 0) &
       WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)

#ifndef PDAF_NO_UPDATE

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 1x, i7, 3x, a)') &
          'PDAF', step, 'Assimilating observations - NETF'
  END IF
  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_netf_update', debug, &
          'Configuration: param_int(3) dim_lag     ', dim_lag
     WRITE (*,*) '++ PDAF-debug PDAF_netf_update', debug, &
          'Configuration: param_int(4) noise_type  ', noise_type
     WRITE (*,*) '++ PDAF-debug PDAF_netf_update', debug, &
          'Configuration: param_int(5) type_forget ', type_forget
     WRITE (*,*) '++ PDAF-debug PDAF_netf_update', debug, &
          'Configuration: param_int(6) type_trans  ', type_trans
     WRITE (*,*) '++ PDAF-debug PDAF_netf_update', debug, &
          'Configuration: param_int(7) type_winf   ', type_winf

     WRITE (*,*) '++ PDAF-debug PDAF_netf_update', debug, &
          'Configuration: param_real(1) forget     ', forget
     WRITE (*,*) '++ PDAF-debug PDAF_netf_update', debug, &
          'Configuration: param_real(2) limit_winf ', limit_winf
     WRITE (*,*) '++ PDAF-debug PDAF_netf_update', debug, &
          'Configuration: param_real(3) noise amp. ', noise_amp
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
       state_p, ens_p, rndmat, Uinv, type_forget, forget, &
       type_winf, limit_winf, noise_type, noise_amp, &
       HX_p, obs_p, U_likelihood, screen, debug, flag)

  IF (debug>0) THEN
     DO i = 1, dim_ens
        WRITE (*,*) '++ PDAF-debug PDAF_netf_update:', debug, 'ensemble member', i, &
             ' analysis values (1:min(dim_p,6)):', ens_p(1:min(dim_p,6),i)
     END DO
  END IF

  ! *** Perform smoothing of past ensembles ***
  IF (dim_lag > 0) THEN
     CALL PDAF_timeit(51, 'new')

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_netf_update -- call smoother function'
     CALL PDAF_smoother_netf(dim_p, dim_ens, dim_lag, TA_noinfl, sens_p, &
          cnt_maxlag, screen)
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_netf_update -- exit smoother function'

     CALL PDAF_timeit(51, 'old')

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

  ! Deallocate observation arrays
  CALL PDAFobs_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_netf_update -- END'

END SUBROUTINE PDAF_netf_update
