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
!> Interface to transfer state to PDAF
!!
!! Interface routine called from the model after the 
!! forecast of each ensemble state to transfer data
!! from the model to PDAF.  For the parallelization 
!! this involves transfer from model PEs to filter 
!! PEs.\\
!! During the forecast phase state vectors are 
!! re-initialized from the forecast model fields
!! by U\_collect\_state. 
!! At the end of a forecast phase (i.e. when all 
!! ensemble members have been integrated by the model)
!! sub-ensembles are gathered from the model tasks.
!! Subsequently the data assimilation update is performed.
!!
!! The code is very generic. Basically the only
!! filter-specific part if the call to the
!! update-routine PDAF\_X\_update where the analysis
!! is computed.  The filter-specific subroutines that
!! are specified in the call to PDAF\_put\_state\_X
!! are passed through to the update routine
!!
!! Variant for hybrid 3DVAR using ESTKF to
!! update the ensemble perturbations.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2021-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAFput_state_hyb3dvar_estkf

CONTAINS

SUBROUTINE PDAF_put_state_hyb3dvar_estkf(U_collect_state, U_init_dim_obs, U_obs_op, &
     U_init_obs, U_prodRinvA, &
     U_cvt, U_cvt_adj, U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
     U_init_obsvar, U_prepoststep, outflag)

  USE PDAF_communicate_ens, &
       ONLY: PDAF_gather_ens
  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_core, &
       ONLY: dim_p, dim_ens, local_dim_ens, assim_flag, &
       nsteps, step_obs, step, member, member_save, subtype_filter, &
       initevol, state, ens, &
       Ainv, screen, flag, offline_mode
  USE PDAF_mod_parallel, &
       ONLY: mype_world, filterpe, &
       dim_ens_l, modelpe, filter_no_model
  USE PDAFobs, &
       ONLY: dim_obs
  USE PDAF_iau, &
       ONLY: step_cnt_iau
  USE PDAF_3dvar, &
       ONLY: dim_cvec, dim_cvec_ens
  USE PDAF_hyb3dvar_update, &
       ONLY: PDAFhyb3dvar_update_estkf

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag  !< Status flag
  
! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_collect_state, &   !< Routine to collect a state vector
       U_init_dim_obs, &           !< Initialize dimension of observation vector
       U_obs_op, &                 !< Observation operator
       U_init_obs, &               !< Initialize observation vector
       U_prepoststep, &            !< User supplied pre/poststep routine
       U_prodRinvA, &              !< Provide product R^-1 A
       U_cvt_ens, &                !< Apply control vector transform matrix (ensemble)
       U_cvt_adj_ens, &            !< Apply adjoint control vector transform matrix (ensemble var)
       U_cvt, &                    !< Apply control vector transform matrix to control vector
       U_cvt_adj, &                !< Apply adjoint control vector transform matrix
       U_obs_op_lin, &             !< Linearized observation operator
       U_obs_op_adj                !< Adjoint observation operator
  EXTERNAL :: U_init_obsvar        !< Initialize mean observation error variance

! *** local variables ***
  INTEGER :: i                     ! Counter


! **************************************************
! *** Save forecasted state back to the ensemble ***
! *** Only done on the filter Pes                ***
! **************************************************

  doevol: IF (nsteps > 0 .OR. .NOT.offline_mode) THEN

     CALL PDAF_timeit(41, 'new')

     modelpes: IF (modelpe) THEN

        ! Store member index for PDAF_get_memberid
        member_save = member

        IF (subtype_filter /= 10 .AND. subtype_filter /= 11) THEN
           ! Save evolved state in ensemble matrix
           CALL U_collect_state(dim_p, ens(1 : dim_p, member))
        ELSE
           ! Save evolved ensemble mean state
           CALL U_collect_state(dim_p, state(1:dim_p))
        END IF
     END IF modelpes

     CALL PDAF_timeit(41, 'old')

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
! *** When forecast phase is completed                 ***
! ***   - collect forecast sub_ensembles on filter PEs ***
! ***   - perform analysis step                        ***
! ***   - re-initialize forecast counters/flags        ***
! ********************************************************
  completeforecast: IF (member == local_dim_ens + 1 &
       .OR. offline_mode) THEN

     ! Set flag for assimilation
     assim_flag = 1


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

     END IF doevolB

     ! *** call timer
     CALL PDAF_timeit(2, 'old')

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
     ELSE IF (.NOT.offline_mode .AND. mype_world==0 .AND. screen > 0) THEN
        WRITE(*,'(a, 5x, a)') 'PDAF', 'Perform assimilation with PDAF'
     ENDIF

     OnFilterPE: IF (filterpe) THEN

        CALL PDAFhyb3dvar_update_estkf(step_obs, dim_p, dim_obs, dim_ens, &
             dim_cvec, dim_cvec_ens, state, Ainv, ens, &
             U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, U_prepoststep, &
             U_cvt, U_cvt_adj, U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
             U_init_obsvar, screen, subtype_filter, flag)

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

END SUBROUTINE PDAF_put_state_hyb3dvar_estkf

END MODULE PDAFput_state_hyb3dvar_estkf
