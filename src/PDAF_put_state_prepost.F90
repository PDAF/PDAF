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
!
!> Interface to transfer state to PDAF
!!
!! Interface routine called from the model during or 
!! after the forecast of each ensemble state to transfer
!! data from the model to PDAF. For the parallelization 
!! this involves transfer from model PEs to filter 
!! PEs.\\
!! This routine can be called at any time during an 
!! ensemble forecast. The routine gathers the state 
!! information from the sub-ensembles.
!! Subsequently the pre-poststep routine U_prepoststep
!! is called to allow to anlize the ensemble.
!!
!! The routine should be called as an alternative to
!! PDAF_put_state_X (with X the name of a filter method)
!! when one wants to analyze the ensemble without
!! performing an analysis step. Note, that the routine
!! does not reset the ensemble, i.e. you should not call
!! PDAF_get_state after this routine.
!!
!! Variant for domain decomposition.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2014-11 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAFput_state_prepost

CONTAINS

SUBROUTINE PDAF_put_state_prepost(U_collect_state, U_prepoststep, outflag)

  USE PDAF_communicate_ens, &
       ONLY: PDAF_gather_ens
  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_mod_filter, &
       ONLY: dim_p, dim_ens, local_dim_ens, &
       nsteps, step_obs, step, member, member_save, subtype_filter, &
       state, ens, Ainv, screen, flag, initevol, offline_mode
  USE PDAF_mod_filtermpi, &
       ONLY: mype_world, filterpe, dim_ens_l, filter_no_model
  USE PDAFobs, &
       ONLY: dim_obs
  USE PDAF_iau, &
       ONLY: step_cnt_iau

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag  !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_collect_state, &   !< Routine to collect a state vector
       U_prepoststep               !< User supplied pre/poststep routine

! *** local variables ***
  INTEGER :: i   ! Counter


! **************************************************
! *** Save forecast state back to the ensemble   ***
! *** Only done on the filter Pes                ***
! **************************************************

  doevol: IF (nsteps > 0) THEN

     ! Store member index for PDAF_get_memberid
     member_save = member

     IF (subtype_filter /= 2 .AND. subtype_filter /= 3) THEN
        ! Save evolved state in ensemble matrix
        CALL U_collect_state(dim_p, ens(1 : dim_p, member))
     ELSE
        ! Save evolved ensemble mean state
        CALL U_collect_state(dim_p, state(1 : dim_p))
     END IF

     member = member + 1

     ! Reset step counter for IAU
     step_cnt_iau = 0
  ELSE
     member = local_dim_ens + 1
  END IF doevol

  IF (filter_no_model .AND. filterpe) THEN
     member = local_dim_ens + 1
  END IF


! ********************************************************
! *** Now (at any time during an ensemble forecast)    ***
! ***   - collect forecast sub_ensembles on filter PEs ***
! ***   - execute pre/poststep routine                 ***
! ********************************************************

  ! ***********************************************
  ! *** Collect forecast ensemble on filter PEs ***
  ! ***********************************************

  doevolB: IF (nsteps > 0) THEN

     IF (.not.filterpe) THEN
        ! Non filter PEs only store a sub-ensemble
        CALL PDAF_gather_ens(dim_p, dim_ens_l, ens, screen)
     ELSE
        ! On filter PEs, the ensemble array has full size
        CALL PDAF_gather_ens(dim_p, dim_ens, ens, screen)
     END IF

  end IF doevolB

  ! *** call timer
  CALL PDAF_timeit(2, 'old')

  IF (.NOT.offline_mode .AND. mype_world == 0 .AND. screen > 1) THEN
     WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
          'PDAF', '--- duration of forecast phase:', PDAF_time_temp(2), 's'
  END IF


  ! ************************
  ! *** Analyze ensemble ***
  ! ************************

  ! Screen output
  IF (mype_world == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 52a)') 'PDAF prepoststep ', ('-', i = 1, 52)
  ENDIF
     
  OnFilterPE: IF (filterpe) THEN
     CALL U_prepoststep(-step_obs, dim_p, dim_ens, dim_ens_l, dim_obs, &
          state, Ainv, ens, flag)
  END IF OnFilterPE


  ! ***********************************
  ! *** Set forecast counters/flags ***
  ! ***********************************

  initevol = 1
  member   = 1
  step     = step_obs + 1


! ********************
! *** finishing up ***
! ********************

  outflag = flag

END SUBROUTINE PDAF_put_state_prepost

END MODULE PDAFput_state_prepost
