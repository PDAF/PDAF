! Copyright (c) 2014-2018 Paul Kirchgessner
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
! !ROUTINE: PDAF_put_state_ewpf --- Interface to transfer state to PDAF
!
! !INTERFACE:
SUBROUTINE PDAF_put_state_ewpf(U_collect_state, U_obs_op, U_init_dim_obs, &
     U_init_obs, U_prodRinvA, U_solve_invHQHTpR, U_adjoint_obs_op, &
     U_prodQA, U_prepoststep, U_randVec, outflag)

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
! 2014-05 - Paul Kirchgessner - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_ewpf, &
       ONLY: dim_p, dim_obs, dim_ens, local_dim_ens, &
       nsteps, step_obs, step, member, subtype_filter, &
       type_forget, initevol, state,eofU, eofV, &
       forget, screen, flag, observation, weights, weight
  USE PDAF_mod_filtermpi, &
       ONLY: mype_world, mype_filter, mype_couple, filterpe, dim_ens_l

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(out) :: outflag  ! Status flag
  
! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
       U_obs_op, &            ! Observation operator
       U_prodRinvA, &         ! Provide product R^-1 A
       U_init_dim_obs, &      ! Get observation dimension
       U_init_obs, &          ! Init observations
       U_solve_invHQHTpR,&    ! Calculate (HQH^T+R)^-1
       U_adjoint_obs_op , &   ! Provide adjoint of observation operator
       U_prodQA         ,&    ! Provide product Q*A where A is some vector or matrix
       U_randVec, &           ! Generate model error 
       U_prepoststep

! !CALLING SEQUENCE:
! Called by: model code  
! Calls: U_collect_state
! Calls: PDAF_gather_ens
! Calls: PDAF_etkf_update
! Calls: PDAF_timeit
!EOP

! local variables
  INTEGER :: i                     ! Counter
  INTEGER :: minusStep
  INTEGER, SAVE :: allocflag = 0   ! Flag whether first time allocation is done
  REAL, ALLOCATABLE :: Ainv(:,:)

! **************************************************
! *** Save forecasted state back to the ensemble ***
! *** Only done on the filter Pes                ***
! **************************************************

  doevol: IF (nsteps >= 0) THEN
     ! Save evolved state in ensemble matrix
     CALL U_collect_state(dim_p, eofV(1 : dim_p, member))

     member = member + 1
  ELSE
     member = local_dim_ens + 1
  END IF doevol


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

     doevolB: IF (nsteps >= 0) THEN

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
        WRITE (*, '(8x, a, F10.3, 1x, a)') &
             '--- duration of forecast phase:', PDAF_time_temp(2), 's'
     END IF


     ! **************************************
     ! *** Perform analysis on filter PEs ***
     ! **************************************

     ! Screen output
     IF (subtype_filter == 5 .AND. mype_world == 0 .AND. screen > 0) THEN
        WRITE (*, '(//1x, 64a)') ('-', i = 1, 64)
        WRITE (*, '(20x, a)') '+++++ ASSIMILATION +++++'
        WRITE (*, '(1x, 64a)') ('-', i = 1, 64)
     ENDIF


     OnFilterPEA: IF (filterpe) THEN

        CALL PDAF_timeit(5, 'new')
        minusStep = - step_obs  ! Indicate forecast by negative time step number
        CALL U_prepoststep(minusStep, dim_p, dim_ens, dim_ens, dim_obs, &
             state, Ainv, eofV, flag)
        CALL PDAF_timeit(5, 'old')

        CALL PDAF_ewpf_ew_step(step, dim_p, dim_obs, dim_ens, &
             eofV, weights, U_obs_op, U_init_dim_obs, &
             U_init_obs, U_solve_invHQHTpR, &
             U_adjoint_obs_op, U_prodQA, U_randVec, U_prodRinvA)

        ! *** Poststep for analysis ensemble ***
        CALL PDAF_timeit(5, 'new')
        CALL U_prepoststep(step_obs, dim_p, dim_ens, dim_ens, dim_obs, &
             state, Ainv, eofV, flag)
        CALL PDAF_timeit(5, 'old')

     END IF OnFilterPEA

 
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

  outflag = 0

END SUBROUTINE PDAF_put_state_ewpf
