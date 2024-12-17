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
!> Interface to PDAF for GLOBALTEMPLATE DA method
!!
!! Interface routine called from the model or PDAF_assimilate after
!! the forecast of each ensemble state. The routine manages
!! the transfer of a model state from the model to PDAF.
!! For the parallelization this involves transfer from model
!! processes (processing elements, PEs) to filter processes.
!!
!! For the flexible parallelization variance the routine manages
!! the ensemble. Thus, during the forecast phase state vectors are 
!! re-initialized from the forecast model fields by U\_collect\_state. 
!! At the end of a forecast phase when all ensemble members have
!! been integrated by the model sub-ensembles are gathered from
!! the model tasks by the routine PDAF\_gather\_ens). Subsequently
!! the filter update is performed by PDAF\_GLOBALTEMPLATE\_update.
!!
!! The code is very generic. The filter-specific part is the call
!! to the update-routine PDAF\_GLOBALTEMPLATE\_update where the analysis
!! is computed. The filter-specific subroutines that are specified
!! as arguments in the call to PDAF\_put\_state\_GLOBALTEMPLATE are
!! passed though to the update routine PDAF\_GLOBALTEMPLATE\_update.
!!
!! ADAPTING THE TEMPLATE:
!! When implementing a filter, the only required changes to this routine
!! should be
!! - replace 'GLOBALTEMPLATE' by the name of the new method
!! - potentially adapt the argument lists in PDAF\_put\_state\_GLOBALTEMPLATE
!!   and PDAF\_GLOBALTEMPLATE\_update
!!
!! __Revision history:__
!! * 2024-12 - Lars Nerger - Initial code for template based on ETKF
!! * Later revisions - see repository log
!!
SUBROUTINE PDAF_put_state_GLOBALTEMPLATE(U_collect_state, U_init_dim_obs, U_obs_op, &
     U_init_obs, U_prepoststep, U_prodRinvA, outflag)

  USE PDAF_communicate_ens, &
       ONLY: PDAF_gather_ens
  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_mod_filter, &
       ONLY: dim_p, dim_obs, dim_ens, local_dim_ens, &
       nsteps, step_obs, step, member, member_save, &
       subtype_filter, incremental, initevol, state, eofV, &
       eofU, sens, dim_lag, cnt_maxlag, offline_mode, &
       screen, flag
  USE PDAF_mod_filtermpi, &
       ONLY: mype_world, filterpe, &
       dim_ens_l, modelpe, filter_no_model

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag  ! Status flag

! *** External subroutines ***
! (PDAF-internal names, real names are defined in the call to PDAF)
  ! Routines for ensemble framework
  EXTERNAL :: U_collect_state, & ! Write model fields into state vector
       U_prepoststep             ! User supplied pre/poststep routine
  ! Observation-related routines for analysis step
  EXTERNAL :: U_init_dim_obs, &  ! Initialize dimension of observation vector
       U_obs_op, &               ! Observation operator
       U_init_obs, &             ! Initialize observation vector
       U_prodRinvA               ! Provide product R^-1 A

! *** Local variables ***
  INTEGER :: i                   ! Counter


! ***************************************************************
! *** Store forecasted state back to the ensemble array eofV  ***
! *** and increment counter `member` for ensemble state index ***
! *** Only done on the filter processes                       ***
! ***************************************************************

  doevol: IF (nsteps > 0 .OR. .NOT.offline_mode) THEN

     CALL PDAF_timeit(41, 'new')

     modelpes: IF (modelpe) THEN

        ! Store member index for PDAF_get_memberid
        member_save = member

        IF (subtype_filter /= 2 .AND. subtype_filter /= 3) THEN
           ! Save evolved state in ensemble matrix
           CALL U_collect_state(dim_p, eofV(1:dim_p, member))
        ELSE
           ! Save evolved ensemble mean state
           CALL U_collect_state(dim_p, state(1:dim_p))
        END IF
     END IF modelpes

     CALL PDAF_timeit(41, 'old')

     member = member + 1
  ELSE
     member = local_dim_ens + 1
  END IF doevol

  IF (filter_no_model .AND. filterpe) THEN
     member = local_dim_ens + 1
  END IF


! ******************************************************
! *** When forecast phase is completed               ***
! *** - gather forecast sub_ensembles on filter PEs  ***
! *** - perform analysis step                        ***
! *** - re-initialize forecast counters/flags        ***
! ******************************************************

  completeforecast: IF (member == local_dim_ens + 1 &
       .OR. offline_mode) THEN

     ! ****************************************************
     ! *** Gather forecast ensemble on filter processes ***
     ! ****************************************************

     doevolB: IF (nsteps > 0) THEN

        IF (.not.filterpe) THEN
           ! Non filter PEs only store a sub-ensemble
           CALL PDAF_gather_ens(dim_p, dim_ens_l, eofV, screen)
        ELSE
           ! On filter PEs, the ensemble array has full size
           CALL PDAF_gather_ens(dim_p, dim_ens, eofV, screen)
        END IF

     END IF doevolB

     ! *** call timer
     CALL PDAF_timeit(2, 'old')

     ! Screen output
     IF (.NOT.offline_mode .AND. mype_world == 0 .AND. screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- duration of forecast phase:', PDAF_time_temp(2), 's'
     END IF


     ! **************************************
     ! *** Perform analysis on filter PEs ***
     ! **************************************

     ! Screen output
     IF (offline_mode .AND. mype_world == 0 .AND. screen > 0) THEN
        WRITE (*, '(//a5, 64a)') 'PDAF ',('-', i = 1, 64)
        WRITE (*, '(a, 20x, a)') 'PDAF', '+++++ ASSIMILATION +++++'
        WRITE (*, '(a5, 64a)') 'PDAF ', ('-', i = 1, 64)
     ENDIF

     OnFilterPE: IF (filterpe) THEN

        CALL PDAF_GLOBALTEMPLATE_update(step_obs, dim_p, dim_obs, dim_ens, &
             state, eofU, eofV, U_init_dim_obs, U_obs_op, &
             U_init_obs, U_prodRinvA, U_prepoststep, &
             screen, subtype_filter, incremental, &
             dim_lag, sens, cnt_maxlag, flag)

     END IF OnFilterPE


     ! ***********************************
     ! *** Set forecast counters/flags ***
     ! ***********************************

     initevol = 1
     member   = 1
     step     = step_obs + 1

  END IF completeforecast


! ********************
! *** finishing up ***
! ********************

  outflag = flag

END SUBROUTINE PDAF_put_state_GLOBALTEMPLATE
