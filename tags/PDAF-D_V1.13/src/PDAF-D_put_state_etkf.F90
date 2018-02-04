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
!$Id: PDAF-D_put_state_etkf.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_put_state_etkf --- Interface to transfer state to PDAF
!
! !INTERFACE:
SUBROUTINE PDAF_put_state_etkf(U_collect_state, U_init_dim_obs, U_obs_op, &
     U_init_obs, U_prepoststep, U_prodRinvA, U_init_obsvar, outflag)

! !DESCRIPTION:
! Interface routine called from the model after the 
! forecast of each ensemble state to transfer data
! from the model to PDAF.  For the parallelization 
! this involves transfer from model PEs to filter 
! PEs.\\
! During the forecast phase state vectors are 
! re-initialized from the forecast model fields
! by U\_collect\_state. 
! At the end of a forecast phase (i.e. when all 
! ensemble members have been integrated by the model)
! sub-ensembles are gathered from the model tasks.
! Subsequently the filter update is performed.
!
! The code is very generic. Basically the only
! filter-specific part if the call to the
! update-routine PDAF\_X\_update where the analysis
! is computed.  The filter-specific subroutines that
! are specified in the call to PDAF\_put\_state\_X
! are passed through to the update routine
!
! Variant for ETKF with domain decomposition.
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
       ONLY: dim_p, dim_obs, dim_ens, local_dim_ens, &
       nsteps, step_obs, step, member, subtype_filter, &
       type_forget, incremental, initevol, state, eofV, &
       eofU, state_inc, forget, screen, flag, &
       sens, dim_lag, cnt_maxlag
  USE PDAF_mod_filtermpi, &
       ONLY: mype_world, filterpe, &
       dim_ens_l, modelpe, filter_no_model

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(out) :: outflag  ! Status flag
  
! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
       U_init_dim_obs, &      ! Initialize dimension of observation vector
       U_obs_op, &            ! Observation operator
       U_init_obsvar, &       ! Initialize mean observation error variance
       U_init_obs, &          ! Initialize observation vector
       U_prepoststep, &       ! User supplied pre/poststep routine
       U_prodRinvA            ! Provide product R^-1 A

! !CALLING SEQUENCE:
! Called by: model code  
! Calls: U_collect_state
! Calls: PDAF_gather_ens
! Calls: PDAF_etkf_update
! Calls: PDAF_timeit
!EOP

! local variables
  INTEGER :: i                     ! Counter
  INTEGER, SAVE :: allocflag = 0   ! Flag whether first time allocation is done


! **************************************************
! *** Save forecasted state back to the ensemble ***
! *** Only done on the filter Pes                ***
! **************************************************

  doevol: IF (nsteps > 0 .OR. subtype_filter /= 5) THEN
     modelpes: IF (modelpe) THEN
        IF (subtype_filter /= 2 .AND. subtype_filter /= 3) THEN
           ! Save evolved state in ensemble matrix
           CALL U_collect_state(dim_p, eofV(1 : dim_p, member))
        ELSE
           ! Save evolved ensemble mean state
           CALL U_collect_state(dim_p, state(1:dim_p))
        END IF
     END IF modelpes

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
       .OR. subtype_filter == 5) THEN

     ! ***********************************************
     ! *** Collect forecast ensemble on filter PEs ***
     ! ***********************************************

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

     IF (subtype_filter /= 5 .AND. mype_world == 0 .AND. screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- duration of forecast phase:', PDAF_time_temp(2), 's'
     END IF


     ! **************************************
     ! *** Perform analysis on filter PEs ***
     ! **************************************

     ! Screen output
     IF (subtype_filter == 5 .AND. mype_world == 0 .AND. screen > 0) THEN
        WRITE (*, '(//a5, 64a)') 'PDAF ',('-', i = 1, 64)
        WRITE (*, '(a, 20x, a)') 'PDAF', '+++++ ASSIMILATION +++++'
        WRITE (*, '(a5, 64a)') 'PDAF ', ('-', i = 1, 64)
     ENDIF

     OnFilterPE: IF (filterpe) THEN

        IF (incremental == 0) THEN
           ! Allocate only if no incremental updating is used. 
           ! With incremental STATE_INC is allocated in PDAF_filter_init.
           ALLOCATE(state_inc(dim_p))
           IF (allocflag == 0) THEN
              CALL PDAF_memcount(3, 'r', dim_p)
              allocflag = 1
           END IF
        END IF

        CALL PDAF_etkf_update(step_obs, dim_p, dim_obs, dim_ens, &
             state, eofU, eofV, state_inc, forget, &
             U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, U_init_obsvar, &
             U_prepoststep, screen, subtype_filter, incremental, type_forget, &
             dim_lag, sens, cnt_maxlag, flag)

        IF (incremental == 0) DEALLOCATE(state_inc)

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

END SUBROUTINE PDAF_put_state_etkf
