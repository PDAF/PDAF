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
! !ROUTINE: PDAF_get_state --- Interface to control ensemble integration
!
! !INTERFACE:
SUBROUTINE PDAF_get_state(steps, time, doexit, U_next_observation, U_distribute_state, &
     U_prepoststep, outflag)

! !DESCRIPTION:
! Interface routine called from the model before the 
! forecast of each ensemble state to transfer data
! from PDAF to the model.  For the parallelization 
! this involves transfer from filter PEs to model PEs.\\
! At the beginning of a forecast phase sub-ensembles
! are distributed to the model tasks. During the 
! forecast phase each state vector of a sub-ensemble
! is transferred to the model fields by U\_dist\_state.
!
! This version is for PDAF with domain-decomposition.
!
! !  This is a core routine of PDAF and
! should not be changed by the user   !
!
! !REVISION HISTORY:
! 2003-07 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE mpi
  USE PDAF_communicate_ens, &
       ONLY: PDAF_scatter_ens
  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_mod_filter, &
       ONLY: dim_p, dim_eof, dim_ens, local_dim_ens, dim_obs, nsteps, &
       step_obs, step, member_get, member_put=>member, member_save, subtype_filter, &
       ensemblefilter, initevol, epsilon, state, eofV, eofU, &
       firsttime, end_forecast, screen, flag, dim_lag, sens, &
       cnt_maxlag, cnt_steps
  USE PDAF_mod_filtermpi, &
       ONLY: mype_world, mype_filter, mype_couple, npes_couple, task_id, &
       statetask, filterpe, dim_eof_l, dim_ens_l, all_dis_ens_l, &
       all_dim_ens_l, MPIerr, MPIstatus, COMM_couple, filter_no_model, &
       modelpe, mype_model

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(inout) :: steps   ! Flag and number of time steps
  REAL, INTENT(out)      :: time    ! current model time
  INTEGER, INTENT(inout) :: doexit  ! Whether to exit from forecasts
  INTEGER, INTENT(inout) :: outflag ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_next_observation, &  ! Routine to provide time step, time and dimension
                             !   of next observation
       U_distribute_state, &       ! Routine to distribute a state vector
       U_prepoststep         ! User supplied pre/poststep routine

! !CALLING SEQUENCE:
! Called by: model code
! Calls: PDAF_timeit
! Calls: PDAF_time_temp
! Calls: U_prepoststep
! Calls: U_next_observation
! Calls: U_distribute_state
! Calls: MPI_bcast (MPI)
!EOP

! local variables
  INTEGER :: i, j, rank, col_frst, col_last  ! Counters
  LOGICAL :: central_state    ! Perform evolution of only central state in SEEK
  INTEGER, ALLOCATABLE :: MPIreqs(:)      ! Array of MPI requests
  INTEGER, ALLOCATABLE :: MPIstats(:,:)   ! Array of MPI statuses


! **********************
! *** INITIALIZATION ***
! **********************

  firstt: IF (firsttime == 1 .AND. subtype_filter /= 5) THEN

     ! Screen output
     IF (mype_world == 0 .AND. screen > 0) THEN
        WRITE (*, '(//a5, 64a)') 'PDAF ',('-', i = 1, 64)
        WRITE (*, '(a, 20x, a)') 'PDAF', '+++++ ASSIMILATION +++++'
        WRITE (*, '(a5, 64a)') 'PDAF ', ('-', i = 1, 64)
     END IF

     ! *** Initial prestep
     IF (filterpe) THEN
        CALL PDAF_timeit(5, 'new')

        IF (mype_world == 0 .AND. screen > 0) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF','Call pre-post routine at initial time'
        ENDIF

        typef: IF (ensemblefilter) THEN
           ! *** Ensemble-based filter      ***
           CALL U_prepoststep(step_obs, dim_p, dim_ens, dim_ens_l, dim_obs, &
                state, eofU, eofV, flag)
        ELSE
           ! *** Mode-based filter (SEEK)   ***
           CALL U_prepoststep(step_obs, dim_p, dim_eof, dim_eof_l, dim_obs, &
                state, eofU, eofV, flag)
        END IF typef

        CALL PDAF_timeit(5, 'old')
     END IF
    
     IF (mype_world == 0 .AND. screen > 0) THEN
        IF (screen >= 2) THEN
           WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
                'PDAF','--- duration of prestep:', PDAF_time_temp(5), 's'
        END IF
        WRITE (*, '(a, 55a)') 'PDAF Forecast ',('-', i = 1, 55)
     END IF

     ! *** Unset first time flag
     firsttime = 0

  END IF firstt


  ! *** Initialize local ensemble sizes ***
  IF (ensemblefilter) THEN
     local_dim_ens = dim_ens_l
  ELSE
     ! For SEEK take care of central state
     IF (task_id /= statetask) THEN
        local_dim_ens = dim_eof_l
     ELSE
        local_dim_ens = dim_eof_l + 1
     END IF
  END IF

  IF (subtype_filter == 5) THEN
     ! For offline mode of PDAF skip initialization
     ! of forecast and loop through ensemble states
     initevol = 0
     nsteps   = 0
     step_obs = 0
     step     = 0
  END IF

! *********************************
! *** Initialize forecast phase ***
! *********************************
  evolinit: IF (initevol == 1) THEN

     ! *** call timer
     CALL PDAF_timeit(2, 'new')

     ! ********************************************
     ! *** Get information on next observations ***
     ! ********************************************

     CALL U_next_observation(step_obs, nsteps, end_forecast, time)
     ! *** Set time step of next observation ***
     step_obs = step + nsteps - 1
     ! *** Initialize ensemble of error states for SEEK ***
     SEEK1: IF ((.NOT.ensemblefilter) .AND. filterpe &
         .AND. (subtype_filter == 0 .OR. subtype_filter ==1)) THEN
        DO  j = 1, dim_eof
           DO i = 1, dim_p
              eofV(i, j) = state(i) + epsilon * eofV(i, j)
           END DO
        END DO
     END IF SEEK1


     ! **********************************************************
     ! *** Shift states in ensemble array in case of smoother ***
     ! **********************************************************

     IF (nsteps > 0 .AND. dim_lag > 0) THEN
        CALL PDAF_smoother_shift(dim_p, dim_ens, dim_lag, eofV, sens, &
             cnt_maxlag, screen)
     END IF


     ! *****************************************************
     ! *** Distribute ensembles and forecast information ***
     ! *****************************************************

     doevol: IF (nsteps > 0 .AND. end_forecast /= 1) THEN
        ! *** More forecasts to be done             ***
        ! *** Send sub-ensembles to each task which ***
        ! *** not includes the filter PEs           ***

        IF ((mype_world == 0) .AND. (screen > 0)) THEN
           IF (subtype_filter == 2 .OR. subtype_filter == 3) THEN
              ! Output for fixed-basis (only SEEK/SEIK/LSEIK/ESTKF/LESTKF)
              WRITE (*, '(a, 5x, a)') 'PDAF', 'Fixed basis - evolve only state estimate'
           ELSE
              IF (ensemblefilter) THEN
                 ! Ensemble-based filters
                 WRITE (*, '(a, 5x, a)') 'PDAF', 'Evolve state ensemble'
              ELSE
                 ! Mode-based: SEEK with evolution of basis
                 WRITE (*, '(a, 5x, a)') 'PDAF', 'Evolve state and error space basis'
              END IF
           END IF
        END IF


        ! *** Distribute sub-ensembles to all model PEs 0    ***
        ! ***                                                ***
        ! *** This is used if multiple model tasks are used. ***

        ! call timer
        CALL PDAF_timeit(19, 'new')
        CALL PDAF_timeit(49, 'new')

        CALL PDAF_scatter_ens(dim_p, dim_ens_l, eofV, state, screen)
     
        ! call timer
        CALL PDAF_timeit(19, 'old')
        CALL PDAF_timeit(49, 'old')

     END IF doevol

     ENSF1: IF (ensemblefilter .AND. filterpe &
          .AND. (subtype_filter == 2 .OR. subtype_filter == 3)) THEN

        ! *** Compute ensemble mean state ***
        state = 0.0
        DO j = 1, dim_ens
           DO i = 1, dim_p
              state(i) = state(i) + eofV(i, j)
           END DO
        END DO
        state = state / REAL(dim_ens)

        ! *** Remove mean from ensemble members ***
        DO j = 1, dim_ens
           DO i = 1, dim_p
              eofV(i, j) = eofV(i, j) - state(i)
           END DO
        END DO

     END IF ENSF1

     ! *** Set INIT flag ***
     initevol = 2

  ELSE IF (initevol == 2) THEN
     ! Routine is called just after the first ensemble member is evolved
     initevol=0
     
  END IF evolinit


! ********************************************************
! *** Initialize state variables for ensemble forecast ***
! ********************************************************
  doevol1: IF (nsteps > 0) THEN
     IF (ensemblefilter) THEN
        IF (subtype_filter/=2 .AND. subtype_filter/=3) THEN
           IF ((screen > 2) .AND. modelpe .AND. mype_model==0) &
                WRITE (*,*) 'PDAF: get_state - Evolve member ', member_get, &
                'in task ', task_id
        ELSE
           IF ((screen > 2) .AND. modelpe .AND. mype_model==0) &
                WRITE (*,*) 'PDAF: get_state - Evolve ensemble mean state ', &
                'in task ', task_id
        END IF
     ELSE
        IF ((task_id == statetask) .AND. (member_get == dim_eof_l + 1)) THEN
           IF ((screen > 2) .AND. modelpe .AND. mype_model==0) &
                WRITE (*,*) 'PDAF: get_state - Evolve central state ', &
                'in task ', task_id
        ELSE
           IF ((screen > 2) .AND. modelpe .AND. mype_model==0) &
                WRITE (*,*) 'PDAF: get_state - Evolve member ', member_get, &
                'in task ', task_id
        END IF
     END IF

     ! *** Distribute state fields within model communicators ***
     IF (ensemblefilter) THEN
        IF (subtype_filter/=2 .AND. subtype_filter/=3) THEN
           IF ((screen > 2) .AND. modelpe .AND. mype_model==0) &
                WRITE (*,*) 'PDAF: get_state - Distribute state fields ', &
                ' in ', task_id, ', member ', member_get
        ELSE
           IF ((screen > 2) .AND. modelpe .AND. mype_model==0) &
                WRITE (*,*) 'PDAF: get_state - Distribute state fields ', &
                ' in ', task_id, ', ensemble mean state '
        END IF
     ELSE
        IF ((task_id == statetask) .AND. (member_get == dim_eof_l + 1)) THEN
           IF ((screen > 2) .AND. modelpe .AND. mype_model==0) &
                WRITE (*,*) 'PDAF: get_state - Distribute state fields ', &
                ' in ', task_id, ', central state '
        ELSE
           IF ((screen > 2) .AND. modelpe .AND. mype_model==0) &
                WRITE (*,*) 'PDAF: get_state - Distribute state fields ', &
                ' in ', task_id, ', member ', member_get
        END IF
     END IF

     ! *** call timer
     CALL PDAF_timeit(40, 'new')

     filtertype: IF (ensemblefilter) THEN

        modelpes: IF (modelpe) THEN

           ! Store member index for PDAF_get_memberid
           member_save = member_get

           IF (subtype_filter/=2 .AND. subtype_filter/=3) THEN
              ! Dynamic ensemble filter with ensemble forecast

              ! distribute ensemble state
              CALL U_distribute_state(dim_p, eofV(1:dim_p, member_get))
              IF ((screen > 2) .AND. modelpe .AND. mype_model==0) &
                   WRITE (*,*) 'PDAF: get_state - task: ', task_id, &
                   ' evolve sub-member ', member_get

              ! Increment member counter
              member_get = member_get + 1
           ELSE
              ! Ensemble filter with fixed error-space basis 
              ! (Option ONLY for SEIK/LSEIK)

              ! set member to maximum
              member_get=dim_ens_l
              member_put=dim_ens_l

              ! distribute and evolve ensemble mean state
              CALL U_distribute_state(dim_p, state)
              IF ((screen > 2) .AND. modelpe .AND. mype_model==0) &
                   WRITE (*,*) 'PDAF: get_state - task: ', task_id, &
                   ' evolve ensemble mean state '
           END IF

        END IF modelpes

        ! Reset member counting
        IF (member_get == local_dim_ens+1) member_get = 1

     ELSE
        ! Mode-based filter (SEEK)
        !!!! Set CENTRAL_STATE = .TRUE. for evolution of only the central
        !!!! state for computation of free evolution.
        !!!! In this case also set FREE_EVOLUTION=.TRUE. in PDAF-D_SEEK_UPDATE
        central_state = .FALSE.
        IF (central_state) THEN
           WRITE (*,*) 'PDAF-NOTICE: EVOLVE ONLY CENTRAL STATE FOR FREE EVOLUTION !!!'
           member_get = dim_eof_l + 1
           member_put = dim_eof_l + 1
        END IF
        ! For fixed basis SFEK set member to maximum
        IF (subtype_filter == 2 .OR. subtype_filter == 3) THEN
           member_get = dim_eof_l + 1
           member_put = dim_eof_l + 1
        END IF

        modelpesB: IF (modelpe) THEN

           ! Store member index for PDAF_get_memberid
           member_save = member_get

           IF ((task_id == statetask) .AND. (member_get == dim_eof_l + 1)) THEN
              ! distribute central state
              CALL U_distribute_state(dim_p, state)
              IF ((screen > 2) .AND. filterpe) &
                   WRITE (*,*) 'PDAF: get_state - task: ',task_id, &
                   ' evolve central state'

              ! Reset member counting
              IF (member_get == dim_eof_l+1) member_get = 1
           ELSE
              ! distribute ensemble state
              CALL U_distribute_state(dim_p, eofV(1:dim_p, member_get))
              IF ((screen > 2) .AND. filterpe) &
                   WRITE (*,*) 'PDAF: get_state - task: ',task_id, &
                   ' evolve sub-member ',member_get

              ! Increment member counter
              member_get = member_get + 1
           END IF
        END IF modelpesB

        ! Reset member counting
        IF (task_id /= statetask .AND. member_get == dim_eof_l+1) member_get = 1

     END IF filtertype

     ! *** call timer
     CALL PDAF_timeit(40, 'old')

  END IF doevol1


! ********************
! *** finishing up ***
! ********************

  ! Set time step counter for next forecast phase
  ! (only used with PDAF_assimilate_X in fully parallel variant)
  cnt_steps = 0

  ! Set variables for exiting the routine
  steps = nsteps
  doexit = end_forecast
  outflag = flag

END SUBROUTINE PDAF_get_state
