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
!> Control analysis update of the NETF
!!
!! Routine to control the analysis update of the NETF.
!! 
!! The analysis and ensemble tranformation are performed by
!! calling PDAF\_netf\_analysis. In addition, the routine
!! U\_prepoststep is called prior to the analysis and after
!! the resampling to allow the user to access the ensemble
!! information.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2014-05 - Paul Kirchgessner - Initial code based on LETKF
!! * Other revisions - see repository log
!!
MODULE PDAF_netf_update

CONTAINS
SUBROUTINE  PDAFnetf_update(step, dim_p, dim_obs_p, dim_ens, &
     state_p, Ainv, ens_p, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_likelihood, U_prepoststep, &
     screen, subtype, dim_lag, sens_p, cnt_maxlag, flag)

  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_netf, &
       ONLY: type_trans, type_winf, limit_winf, forget, debug
  USE PDAF_netf, &
       ONLY: debug, forget, type_forget, type_trans, type_winf, &
       limit_winf, type_noise, noise_amp
  USE PDAF_mod_parallel, &
       ONLY: mype, dim_ens_l
  USE PDAFobs, &
       ONLY: PDAFobs_init, PDAFobs_dealloc, type_obs_init, &
       HX_p, obs_p
  USE PDAF_analysis_utils, &
       ONLY: PDAF_inflate_ens, PDAF_generate_rndmat
  USE PDAF_netf_analysis, &
       ONLY: PDAF_netf_ana, PDAF_netf_smootherT, PDAF_smoother_netf

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step         !< Current time step
  INTEGER, INTENT(in) :: dim_p        !< PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p   !< PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens      !< Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)        !< PE-local model state
  REAL, INTENT(inout) :: Ainv(dim_ens, dim_ens)!< Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) !< PE-local ensemble matrix
  INTEGER, INTENT(in) :: screen       !< Verbosity flag
  INTEGER, INTENT(in) :: subtype      !< Filter subtype
  INTEGER, INTENT(in) :: dim_lag      !< Number of past time instances for smoother
  REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag) !< PE-local smoother ensemble
  INTEGER, INTENT(inout) :: cnt_maxlag !< Count number of past time steps for smoothing
  INTEGER, INTENT(inout) :: flag      !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, &       !< Initialize dimension of observation vector
       U_obs_op, &                    !< Observation operator
       U_init_obs, &                  !< Initialize observation vector
       U_prepoststep, &               !< User supplied pre/poststep routine
       U_likelihood                   !< Compute observation likelihood for an ensemble member

! *** local variables ***
  INTEGER :: i, j                     ! Counters
  INTEGER :: minusStep                ! Time step counter
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  LOGICAL :: do_init_dim_obs          ! Flag for initializing dim_obs_p in PDAFobs_init
  LOGICAL :: do_ensmean               ! Flag for computing ensemble mean state
  REAL, ALLOCATABLE :: TA_noinfl(:,:) ! TA for smoother (without inflation)
  REAL, ALLOCATABLE :: rndmat(:,:)    ! Orthogonal random matrix


! ***********************************************************
! *** For fixed error space basis compute ensemble states ***
! ***********************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_netf_update -- START'

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
        WRITE (*,*) '++ PDAF-debug PDAF_netf_update:', debug, 'ensemble member', i, &
             ' forecast values (1:min(dim_p,6)):', ens_p(1:min(dim_p,6),i)
     END DO
  END IF
  CALL PDAF_timeit(51, 'old')


! ************************
! *** Inflate ensemble ***
! ************************

  do_ensmean = .true.
  IF (type_obs_init==0 .OR. type_obs_init==2) THEN
     ! We need to call the inflation of the forecast ensemble before
     ! the observed ensemble is initialized in PDAFobs_init

     IF ((type_forget==0 .OR. type_forget==1) .AND. (forget /= 1.0)) THEN
        CALL PDAF_timeit(51, 'new')

        IF (mype == 0 .AND. screen > 0) WRITE (*, '(a, 5x, a, i2, a, f10.3)') &
             'PDAF', 'Inflate forecast ensemble, type_forget=',type_forget,', forget=', forget

        ! Apply forgetting factor
        CALL PDAF_inflate_ens(dim_p, dim_ens, state_p, ens_p, forget, do_ensmean)

        ! PDAF_inflate_ens compute the ensmeble mean; thus don't do this in PDAFobs_init
        do_ensmean = .false.

        CALL PDAF_timeit(51, 'old')
     ENDIF


! *****************************************************
! *** Initialize observations and observed ensemble ***
! ***    optionally before call to U_prepoststep    ***
! *****************************************************

     ! This call initializes dim_obs_p, HX_p, HXbar_p, obs_p in the module PDAFobs
     ! It can also compute the ensemble mean and store it in state_p
     CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_p, &
          state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
          screen, debug, do_ensmean, .true., .true., .true., .true.)
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
       state_p, Ainv, ens_p, flag)
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

  CALL PDAF_timeit(3, 'new')

  do_ensmean = .true.
  IF (type_obs_init==1) THEN
     ! We need to call the inflation of the forecast ensemble before
     ! the observed ensemble is initialized in PDAFobs_init

     IF ((type_forget==0 .OR. type_forget==1) .AND. (forget /= 1.0)) THEN
        CALL PDAF_timeit(51, 'new')

        IF (mype == 0 .AND. screen > 0) WRITE (*, '(a, 5x, a, i2, a, f10.3)') &
             'PDAF', 'Inflate forecast ensemble, type_forget=',type_forget,', forget=', forget

        ! Apply forgetting factor
        CALL PDAF_inflate_ens(dim_p, dim_ens, state_p, ens_p, forget, do_ensmean)

        ! PDAF_inflate_ens compute the ensmeble mean; thus don't do this in PDAFobs_init
        do_ensmean = .false.

        CALL PDAF_timeit(51, 'old')
     END IF
  END IF


! *****************************************************
! *** Initialize observations and observed ensemble ***
! *****************************************************

  IF (type_obs_init>0) THEN

     IF (type_obs_init==1) THEN
        do_init_dim_obs=.true.
     ELSE
        ! Skip call to U_init_dim_obs when also called before prepoststep
        do_init_dim_obs=.false.   
     END IF

     ! This call initializes dim_obs_p, HX_p, HXbar_p, obs_p in the module PDAFobs
     ! It can also compute the ensemble mean and store it in state_p
     CALL PDAFobs_init(step, dim_p, dim_ens, dim_obs_p, &
          state_p, ens_p, U_init_dim_obs, U_obs_op, U_init_obs, &
          screen, debug, do_ensmean, do_init_dim_obs, .true., .true., .true.)

  END IF

  CALL PDAF_timeit(3, 'old')


! ***********************
! ***  Analysis step  ***
! ***********************

  IF (mype == 0 .AND. screen > 0) &
       WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)

#ifndef PDAF_NO_UPDATE
  CALL PDAF_timeit(3, 'new')

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 1x, i7, 3x, a)') &
          'PDAF', step, 'Assimilating observations - NETF'
  END IF
  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_netf_update', debug, &
          'Configuration: param_int(3) dim_lag     ', dim_lag
     WRITE (*,*) '++ PDAF-debug PDAF_netf_update', debug, &
          'Configuration: param_int(4) type_noise  ', type_noise
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

     CALL PDAF_timeit(15, 'new')

     ALLOCATE(TA_noinfl(dim_ens, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens*dim_ens)

     CALL PDAF_netf_smootherT(step, dim_p, dim_obs_p, dim_ens, &
          ens_p, rndmat, TA_noinfl,  &
          U_init_dim_obs, U_obs_op, U_init_obs, U_likelihood, &
          screen, flag)

     CALL PDAF_timeit(15, 'old')
  END IF

  ! *** NETF analysis ***

  CALL PDAF_netf_ana(step, dim_p, dim_obs_p, dim_ens, &
       state_p, ens_p, rndmat, Ainv, type_forget, forget, &
       type_winf, limit_winf, type_noise, noise_amp, &
       HX_p, obs_p, U_likelihood, screen, debug, flag)

  IF (debug>0) THEN
     DO i = 1, dim_ens
        WRITE (*,*) '++ PDAF-debug PDAF_netf_update:', debug, 'ensemble member', i, &
             ' analysis values (1:min(dim_p,6)):', ens_p(1:min(dim_p,6),i)
     END DO
  END IF

  ! *** Perform smoothing of past ensembles ***
  IF (dim_lag > 0) THEN
     CALL PDAF_timeit(15, 'new')
     CALL PDAF_timeit(51, 'new')

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_netf_update -- call smoother function'
     CALL PDAF_smoother_netf(dim_p, dim_ens, dim_lag, TA_noinfl, sens_p, &
          cnt_maxlag, screen)
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_netf_update -- exit smoother function'

     CALL PDAF_timeit(51, 'old')
     CALL PDAF_timeit(15, 'old')

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

  IF (allocflag == 0) allocflag = 1

  ! Deallocate observation arrays
  CALL PDAFobs_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_netf_update -- END'

END SUBROUTINE PDAFnetf_update

END MODULE PDAF_netf_update
