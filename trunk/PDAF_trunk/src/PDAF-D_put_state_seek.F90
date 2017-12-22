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
!$Id: PDAF-D_put_state_seek.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_put_state_seek --- Interface to transfer state to PDAF
!
! !INTERFACE:
SUBROUTINE PDAF_put_state_seek(U_collect_state, U_init_dim_obs, U_obs_op, &
     U_init_obs, U_prepoststep, U_prodRinvA, outflag)

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
! The code is rather generic. However, in contrast
! to the put_state routines for ensemble-based filters,
! this variant has also to handle the central state
! the is separately forecasted.
!
! Variant for SEEK with domain decomposition.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2003-07 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_mod_filter, &
       ONLY: dim_p, dim_obs, dim_eof, local_dim_ens, nsteps, &
       step_obs, step, member, subtype_filter, &
       type_forget, int_rediag, incremental, initevol, epsilon, &
       state, eofV, eofU, forget, screen, flag
  USE PDAF_mod_filtermpi, &
       ONLY: mype_world, mype_filter, mype_couple, npes_couple, task_id, &
       statetask, filterpe, COMM_couple, MPI_REALTYPE, MPIerr, &
        MPIstatus, dim_eof_l, modelpe, filter_no_model

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(out) :: outflag  ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
       U_init_dim_obs, &       ! Initialize dimension of observation vector
       U_obs_op, &             ! Observation operator
       U_init_obs, &           ! Initialize observation vector
       U_prepoststep, &        ! User supplied pre/poststep routine
       U_prodRinvA             ! Provide product R^-1 HV
  
! !CALLING SEQUENCE:
! Called by: model code  
! Calls: U_collect_state
! Calls: PDAF_seek_update
! Calls: PDAF_timeit
!EOP

! local variables
  INTEGER :: i  ! Counter


! **************************************************
! *** Save forecasted state back to the ensemble ***
! *** Only done on the filter Pes                ***
! **************************************************

  doevol1: IF (nsteps > 0) THEN
     modelpes: IF (modelpe) THEN
        IF ((task_id == statetask) .AND. (member == local_dim_ens)) THEN
           ! save evolved state fields in state vector
           CALL U_collect_state(dim_p, state)
        ELSE
           ! save evolved state in ensemble matrix
           CALL U_collect_state(dim_p, eofV(1:dim_p, member))
        END IF
     END IF modelpes

     member = member + 1
  ELSE
     member = local_dim_ens + 1
  END IF doevol1

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
           CALL PDAF_gather_ens(dim_p, dim_eof_l, eofV, screen)
        ELSE
           ! On filter PEs, the ensemble array has full size
           CALL PDAF_gather_ens(dim_p, dim_eof, eofV, screen)
        END IF


        ! *** Send from model PEs that are not filter PEs ***
        subensS: if (.not.filterpe .and. npes_couple > 1) then

           ! Send central state to filter PEs
           IF (mype_couple+1 == statetask .AND. statetask /= 1) THEN
              CALL MPI_SEND(state, dim_p, MPI_REALTYPE, &
                   0, mype_couple, COMM_couple, MPIerr)

              IF (screen > 2) WRITE (*,*) &
                   'PDAF: put_state - send state from statetask ',statetask, &
                   ' in couple task ', mype_filter + 1
           END IF

        end if subensS

        ! *** Receive on filter PEs ***
        subensR: if (filterpe .and. npes_couple > 1) then

           ! Receive central state on filter PEs
           IF (statetask /= 1) THEN
              CALL MPI_RECV(state, dim_p, MPI_REALTYPE, &
                   statetask-1, statetask-1, COMM_couple, MPIstatus, MPIerr)
              IF (screen > 2) WRITE (*,*) &
                   'PDAF: put_state - recv state from statetask ', &
                   statetask, ' in couple task ', mype_filter + 1
           END IF

        end if subensR

     end IF doevolB


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
        CALL PDAF_seek_update(step_obs, dim_p, dim_obs, dim_eof, state, &
             eofU, eofV, epsilon, forget, int_rediag, &
             U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, U_prepoststep, &
             screen, subtype_filter, incremental, type_forget, flag)
     END IF OnFilterPE


     ! ***********************************
     ! *** Set forecast counters/flags ***
     ! ***********************************
     initevol = 1
     member = 1
     step = step_obs + 1

  END IF completeforecast
    

! ********************
! *** finishing up ***
! ********************

  outflag = flag

END SUBROUTINE PDAF_put_state_seek
