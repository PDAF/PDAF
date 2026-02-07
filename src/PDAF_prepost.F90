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
!> Interface to PDAF for calling prepoststep
!!
!! Interface routine called from the model at each time
!! step during the forecast of each ensemble state. If
!! the time of the next analysis step is reached the
!! forecast state is transferred to PDAF and the 
!! user-supplied routine U_prepoststep is executed
!!
!! Variant for domain decomposition.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2014-12 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAFprepost

CONTAINS

  SUBROUTINE PDAF_prepost(U_collect_state, U_distribute_state, &
       U_prepoststep, U_next_observation, outflag)

    USE PDAF_mod_core, &
         ONLY: cnt_steps, nsteps, assim_flag, reset_fcst_flag
    USE PDAF_mod_parallel, &
         ONLY: mype_world
    USE PDAFget_state, &
         ONLY: PDAF_get_state
    USE PDAFput_state_prepost, &
         ONLY: PDAF_put_state_prepost

    IMPLICIT NONE
  
! *** Arguments ***
    INTEGER, INTENT(out) :: outflag  !< Status flag
  
! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
    EXTERNAL :: U_collect_state, &  !< Routine to collect a state vector
         U_prepoststep, &           !< User supplied pre/poststep routine
         U_next_observation, &      !< Routine to provide time step, time and dimension
                                    !<   of next observation
         U_distribute_state         !< Routine to distribute a state vector

! *** Local variables ***
    INTEGER :: steps     ! Number of time steps in next forecast phase
    INTEGER :: doexit    ! Exit flag; not used in this variant
    REAL :: time         ! Current model time; not used in this variant


! *****************************
! ***   At each time step   ***
! *****************************

    ! Increment time step counter
    cnt_steps = cnt_steps + 1


! ********************************
! *** At end of forecast phase ***
! ********************************

    IF (cnt_steps == nsteps) THEN

       ! Set flags for assimilation and forecast
       assim_flag = 0
       reset_fcst_flag = 1

       ! *** Call analysis step ***

       CALL PDAF_put_state_prepost(U_collect_state, U_prepoststep, outflag)

       ! *** Prepare start of next ensemble forecast ***

       IF (outflag==0) THEN
          CALL PDAF_get_state(steps, time, doexit, U_next_observation, U_distribute_state, &
               U_prepoststep, outflag)
       END IF

       nsteps = steps

    ELSE
       assim_flag = 0
       reset_fcst_flag = 0
       outflag = 0
    END IF

  END SUBROUTINE PDAF_prepost


!-------------------------------------------------------------------------------
!> Interface to PDAF call to prepost in offline coupling
!!
!! Interface routine called from the main program
!! for the PDAF offline mode.
!!
!! This case is an alternative to any PDAF_assim_offline
!! code as it's only calling U_prepoststep.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code based on put_state routine
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_prepost_offline(U_prepoststep, outflag)

    USE PDAF_mod_core, &
         ONLY: dim_p, dim_ens, assim_flag, step_obs, &
         subtype_filter, screen, flag, offline_mode, &
         state, ens, Ainv
    USE PDAF_mod_parallel, &
         ONLY: mype_world, filterpe, dim_ens_l
    USE PDAF_utils_filters, &
         ONLY: PDAF_configinfo_filters
    USE PDAFobs, &
         ONLY: dim_obs

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(out) :: outflag  !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
    EXTERNAL :: U_prepoststep        !< User supplied pre/poststep routine

! *** local variables ***
    INTEGER :: i                     ! Counter


! *********************************************
! *** Perform analysis step in offline mode ***
! *********************************************

    ! Set flag for assimilation
    assim_flag = 1

    ! Screen output
    IF (mype_world == 0 .AND. screen > 0) THEN
       ! Print configuration info (if not done before in PDAF_set_offline_mode)
       IF (.NOT.offline_mode) CALL PDAF_configinfo_filters(subtype_filter, 1)

       WRITE (*, '(//a5, 64a)') 'PDAF ',('-', i = 1, 64)
       WRITE (*, '(a, 20x, a)') 'PDAF', '+++++ ASSIMILATION +++++'
       WRITE (*, '(a5, 64a)') 'PDAF ', ('-', i = 1, 64)
    ENDIF

    ! Set flag for offline mode
    offline_mode = .true.

    ! Set number of observations to zero
    dim_obs = 0

    OnFilterPE: IF (filterpe) THEN
       CALL U_prepoststep(-step_obs, dim_p, dim_ens, dim_ens_l, dim_obs, &
            state, Ainv, ens, flag)
    END IF OnFilterPE


! ********************
! *** finishing up ***
! ********************

    outflag = flag

  END SUBROUTINE PDAF_prepost_offline

END MODULE PDAFprepost
