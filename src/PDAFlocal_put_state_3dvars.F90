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
!> Interfaces to PDAF for flexible parallelization mode for 3D-Vars using PDAFlocal
!!
!! The interface routines provide the advanced compact
!! interfaces for using PDAF-Local. The routines
!! just call of one the PDAF_put_state interface routines
!! with the full interface. In the call the specific PDAF
!! internal subroutines for PDAF-Local are specified.
!!
!! The interface routines provided here are the older PDAF-2
!! routines using the naming PDAFlocal.
!!
!! !  This is a core file of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code by collecting files into a module
!! * Other revisions - see repository log
!!
MODULE PDAFlocal_put_state_3dvars

CONTAINS

!-------------------------------------------------------------------------------
!> Interface to transfer state to PDAF for En3DVAR/LESTKF
!!
!! __Revision history:__
!! * 2021-03 - Lars Nerger - Initial code
!! * 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
!! * Other revisions - see repository log
!!
SUBROUTINE PDAFlocal_put_state_en3dvar_lestkf(U_collect_state, U_init_dim_obs, U_obs_op, &
     U_init_obs, U_prodRinvA, &
     U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
     U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
     U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
     U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
     U_prepoststep, outflag)

  USE PDAF_communicate_ens, &
       ONLY: PDAF_gather_ens
  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_core, &
       ONLY: dim_p, dim_ens, local_dim_ens, nsteps, step_obs, &
       step, member, member_save, state, ens, Ainv, &
       initevol, subtype_filter, screen, flag, offline_mode
  USE PDAF_mod_parallel, &
       ONLY: mype_world, filterpe, &
       dim_ens_l, modelpe, filter_no_model
  USE PDAFlocal, &
       ONLY: PDAFlocal_g2l_cb, &  ! Project global to local state vector
       PDAFlocal_l2g_cb           ! Project local to global state vecto
  USE PDAF_3dvar, &
       ONLY: dim_cvec_ens
  USE PDAFobs, &
       ONLY: dim_obs
  USE PDAF_en3dvar_update, &
       ONLY: PDAFen3dvar_update_lestkf

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
       U_obs_op_lin, &             !< Linearized observation operator
       U_obs_op_adj                !< Adjoint observation operator
  EXTERNAL :: U_obs_op_f, &        !< Observation operator
       U_init_n_domains_p, &       !< Provide number of local analysis domains
       U_init_dim_l, &             !< Init state dimension for local ana. domain
       U_init_dim_obs_f, &         !< Initialize dimension of observation vector
       U_init_dim_obs_l, &         !< Initialize dim. of obs. vector for local ana. domain
       U_init_obs_f, &             !< Initialize PE-local observation vector
       U_init_obs_l, &             !< Init. observation vector on local analysis domain
       U_init_obsvar, &            !< Initialize mean observation error variance
       U_init_obsvar_l, &          !< Initialize local mean observation error variance
       U_g2l_obs, &                !< Restrict full obs. vector to local analysis domain
       U_prodRinvA_l               !< Provide product R^-1 A on local analysis domain

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
     ENDIF

     OnFilterPE: IF (filterpe) THEN

        CALL PDAFen3dvar_update_lestkf(step_obs, dim_p, dim_obs, dim_ens, &
             dim_cvec_ens, state, Ainv, ens, &
             U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, U_prepoststep, &
             U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
             U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
             U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, PDAFlocal_g2l_cb, &
             PDAFlocal_l2g_cb, U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
             screen, subtype_filter, flag)

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

END SUBROUTINE PDAFlocal_put_state_en3dvar_lestkf


!-------------------------------------------------------------------------------
!> Interface to transfer state to PDAF for Hyb3DVAR/LESTKF
!!
!! __Revision history:__
!! * 2021-03 - Lars Nerger - Initial code
!! * 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
!! * Other revisions - see repository log
!!
SUBROUTINE PDAFlocal_put_state_hyb3dvar_lestkf(U_collect_state, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
     U_cvt_ens, U_cvt_adj_ens, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
     U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
     U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
     U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
     U_prepoststep, outflag)

  USE PDAF_communicate_ens, &
       ONLY: PDAF_gather_ens
  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_core, &
       ONLY: dim_p, dim_ens, local_dim_ens, &
       nsteps, step_obs, step, member, member_save, subtype_filter, &
       initevol, state, ens, Ainv, &
       screen, flag, offline_mode
  USE PDAF_mod_parallel, &
       ONLY: mype_world, filterpe, &
       dim_ens_l, modelpe, filter_no_model
  USE PDAFlocal, &
       ONLY: PDAFlocal_g2l_cb, &  ! Project global to local state vector
       PDAFlocal_l2g_cb ! Project local to global state vector
  USE PDAF_3dvar, &
       ONLY: dim_cvec, dim_cvec_ens
  USE PDAFobs, &
       ONLY: dim_obs
  USE PDAF_hyb3dvar_update, &
       ONLY: PDAFhyb3dvar_update_lestkf

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
  EXTERNAL :: U_obs_op_f, &        !< Observation operator
       U_init_n_domains_p, &       !< Provide number of local analysis domains
       U_init_dim_l, &             !< Init state dimension for local ana. domain
       U_init_dim_obs_f, &         !< Initialize dimension of observation vector
       U_init_dim_obs_l, &         !< Initialize dim. of obs. vector for local ana. domain
       U_init_obs_f, &             !< Initialize PE-local observation vector
       U_init_obs_l, &             !< Init. observation vector on local analysis domain
       U_init_obsvar, &            !< Initialize mean observation error variance
       U_init_obsvar_l, &          !< Initialize local mean observation error variance
       U_g2l_obs, &                !< Restrict full obs. vector to local analysis domain
       U_prodRinvA_l               !< Provide product R^-1 A on local analysis domain

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
     ENDIF

     OnFilterPE: IF (filterpe) THEN

        CALL PDAFhyb3dvar_update_lestkf(step_obs, dim_p, dim_obs, dim_ens, &
             dim_cvec, dim_cvec_ens, state, Ainv, ens, &
             U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, U_prepoststep, &
             U_cvt_ens, U_cvt_adj_ens, U_cvt, U_cvt_adj, U_obs_op_lin, U_obs_op_adj, &
             U_init_dim_obs_f, U_obs_op_f, U_init_obs_f, U_init_obs_l, U_prodRinvA_l, &
             U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, PDAFlocal_g2l_cb, &
             PDAFlocal_l2g_cb, U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
             screen, subtype_filter, flag)

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

END SUBROUTINE PDAFlocal_put_state_hyb3dvar_lestkf

END MODULE PDAFlocal_put_state_3dvars
