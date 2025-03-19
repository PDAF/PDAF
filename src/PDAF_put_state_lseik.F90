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
! The code is very generic. Basically the only
! filter-specific part if the call to the
! update-routine PDAF\_X\_update where the analysis
! is computed.  The filter-specific subroutines that
! are specified in the call to PDAF\_put\_state\_X
! are passed through to the update routine
!
! You should have received a copy of the GNU Lesser General Public
! License along with PDAF.  If not, see <http://www.gnu.org/licenses/>.
!
!> Interface to transfer state to PDAF
!!
!! Interface routine called from the model after the 
!! forecast of each ensemble state to transfer data
!! from the model to PDAF. For the parallelization 
!! this involves transfer from model PEs to filter 
!! PEs.\\
!! During the forecast phase state vectors are 
!! re-initialized from the forecast model fields
!! by U\_collect\_state. 
!! At the end of a forecast phase (i.e. when all 
!! ensemble members have been integrated by the model)
!! sub-ensembles are gathered from the model tasks.
!! Subsequently the filter update is performed.
!!
!! Variant for LSEIK with domain decomposition.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2003-09 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAFput_state_lseik

CONTAINS

SUBROUTINE PDAF_put_state_lseik(U_collect_state, U_init_dim_obs, U_obs_op, &
     U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
     U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
     U_init_obsvar, U_init_obsvar_l, outflag)

  USE PDAF_communicate_ens, &
       ONLY: PDAF_gather_ens
  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_mod_core, &
       ONLY: dim_p, dim_ens, local_dim_ens, assim_flag, &
       nsteps, step_obs, step, member, member_save, subtype_filter, &
       initevol, state, ens, Ainv, screen, flag, offline_mode
  USE PDAF_mod_parallel, &
       ONLY: mype_world, filterpe, dim_ens_l, modelpe, filter_no_model
  USE PDAFobs, &
       ONLY: dim_obs
  USE PDAF_iau, &
       ONLY: step_cnt_iau
  USE PDAF_lseik_update, &
       ONLY: PDAFlseik_update

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag  !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_collect_state, &  !< Routine to collect a state vector
       U_obs_op, &             !< Observation operator
       U_init_n_domains_p, &   !< Provide number of local analysis domains
       U_init_dim_l, &         !< Init state dimension for local ana. domain
       U_init_dim_obs, &       !< Initialize dimension of observation vector
       U_init_dim_obs_l, &     !< Initialize dim. of obs. vector for local ana. domain
       U_init_obs, &           !< Initialize PE-local observation vector
       U_init_obs_l, &         !< Init. observation vector on local analysis domain
       U_init_obsvar, &        !< Initialize mean observation error variance
       U_init_obsvar_l, &      !< Initialize local mean observation error variance
       U_g2l_state, &          !< Get state on local ana. domain from full state
       U_l2g_state, &          !< Init full state from state on local analysis domain
       U_g2l_obs, &            !< Restrict full obs. vector to local analysis domain
       U_prodRinvA_l, &        !< Provide product R^-1 A on local analysis domain
       U_prepoststep           !< User supplied pre/poststep routine

! *** Local variables ***
  INTEGER :: i                 ! Counter


! **************************************************
! *** Save forecast state back to the ensemble   ***
! *** Only done on the filter Pes                ***
! **************************************************

  doevol: IF (nsteps > 0) THEN

     CALL PDAF_timeit(41, 'new')

     modelpes: IF (modelpe) THEN

        ! Store member index for PDAF_get_memberid
        member_save = member

        IF (subtype_filter /= 10 .AND. subtype_filter /= 11) THEN
           ! Save evolved state in ensemble matrix
           CALL U_collect_state(dim_p, ens(1 : dim_p, member))
        ELSE
           ! Save evolved ensemble mean state
           CALL U_collect_state(dim_p, state(1 : dim_p))
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

        IF (.NOT.filterpe) THEN
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
        CALL PDAFlseik_update(step_obs, dim_p, dim_obs, dim_ens, dim_ens-1, &
             state, Ainv, ens, U_init_dim_obs, U_obs_op, &
             U_init_obs, U_init_obs_l, U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, &
             U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, U_init_obsvar, &
             U_init_obsvar_l, U_prepoststep, screen, subtype_filter, flag)
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

END SUBROUTINE PDAF_put_state_lseik

END MODULE PDAFput_state_lseik
