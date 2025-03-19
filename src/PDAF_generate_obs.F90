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
!> Interface to PDAF for observation generation
!!
!! Interface routine called from the model at each time
!! step during the forecast of each ensemble state. If
!! the time of the next analysis step is reached the
!! forecast state is transferred to PDAF and observations
!! are generated according to the observation operator
!! and by adding Gaussian random noise. Subsequently, 
!! PDAF_get_state is called to initialize the next forecast
!! phase. 
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2019-01 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAFgenerate_obs

CONTAINS

SUBROUTINE PDAF_generate_obs(U_collect_state, U_distribute_state, &
     U_init_dim_obs_f, U_obs_op_f, U_init_obserr_f, U_get_obs_f, &
     U_prepoststep, U_next_observation, outflag)

  USE PDAF_mod_core, &
       ONLY: cnt_steps, nsteps, assim_flag, reset_fcst_flag
  USE PDAF_mod_parallel, &
       ONLY: mype_world
  USE PDAFput_state_generate_obs, &
       ONLY: PDAF_put_state_generate_obs
  USE PDAFget_state, &
       ONLY: PDAF_get_state

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag  !< Status flag
  
! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_collect_state, &   !< Routine to collect a state vector
       U_obs_op_f, &               !< Observation operator
       U_init_dim_obs_f, &         !< Initialize dimension of observation vector
       U_init_obserr_f, &          !< Initialize vector of observation error standard deviations
       U_get_obs_f, &              !< Provide observation vector to user 
       U_prepoststep, &            !< User supplied pre/poststep routine
       U_next_observation, &       !< Routine to provide time step, time and dimension of next observation
       U_distribute_state          !< Routine to distribute a state vector

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

     CALL PDAF_put_state_generate_obs(U_collect_state, U_init_dim_obs_f, &
          U_obs_op_f, U_init_obserr_f, U_get_obs_f, U_prepoststep, outflag)

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

END SUBROUTINE PDAF_generate_obs

END MODULE PDAFgenerate_obs
