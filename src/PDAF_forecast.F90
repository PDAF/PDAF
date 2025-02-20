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
!> Module for operations during forecast phase
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAF_forecast

CONTAINS
!> Generic routine for operations during the forecast phase
!!
!! This routine collects operations that are performed on the
!! model state during the forecast phase. Possible operations are
!! e.g. the application of incremental analysis updating (IAU)
!! or the collection of the observed ensemble for asynchrous DA.
!!
  SUBROUTINE PDAF_fcst_operations(step, U_collect_state, U_distribute_state, &
       U_init_dim_obs, U_obs_op, U_init_obs, U_init_obsvar, outflag)

    USE PDAF_mod_filter, &
         ONLY: dim_p, ens, step_obs
    USE PDAF_mod_filtermpi, &
         ONLY: mype_world, dim_ens_task
    USE PDAF_iau, &
         ONLY: PDAF_iau_apply_inc

    IMPLICIT NONE
  
! *** Arguments ***
    INTEGER, INTENT(in) :: step       !< Time step in current forecast phase
    INTEGER, INTENT(inout) :: outflag !< Status flag
  
! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
    EXTERNAL :: U_collect_state, &    !< Routine to collect a state vector
         U_obs_op, &                  !< Observation operator
         U_init_dim_obs, &            !< Initialize dimension of observation vector
         U_init_obs, &                !< Initialize PE-local observation vector
         U_init_obsvar, &             !< Initialize mean observation error variance
         U_distribute_state           !< Routine to distribute a state vector


! ********************
! ***  Initialize  ***
! ********************

    outflag = 0


! *********************
! ***   Apply IAU   ***
! *********************

    CALL PDAF_iau_apply_inc(step, dim_p, dim_ens_task, ens, &
         U_collect_state, U_distribute_state)

  END SUBROUTINE PDAF_fcst_operations

END MODULE PDAF_forecast
