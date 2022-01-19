! Copyright (c) 2004-2021 Lars Nerger
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
! !ROUTINE: PDAFomi_generate_obs --- Interface to transfer state to PDAF
!
! !INTERFACE:
SUBROUTINE PDAFomi_generate_obs(collect_state_pdaf, distribute_state_pdaf, &
       init_dim_obs_f_pdaf, obs_op_f_pdaf, get_obs_f_pdaf, prepoststep_pdaf, &
       next_observation_pdaf, outflag)

! !DESCRIPTION:
! Interface routine called from the model during the 
! forecast of each ensemble state to transfer data
! from the model to PDAF and to perform the analysis
! step.
!
! This routine provides the simplified interface
! where names of user-provided subroutines are
! fixed. It simply calls the routine with the
! full interface using pre-defined routine names.
!
! Variant for gnerating observations with domain
! decomposition.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2019-01 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAFomi, ONLY: PDAFomi_dealloc

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(inout) :: outflag ! Status flag
  
! ! Names of external subroutines 
  EXTERNAL :: collect_state_pdaf, &  ! Routine to collect a state vector
       distribute_state_pdaf, &      ! Routine to distribute a state vector
       prepoststep_pdaf, &           ! User supplied pre/poststep routine
       next_observation_pdaf         ! Provide time step, time and dimension of next observation
  EXTERNAL :: init_dim_obs_f_pdaf, & ! Initialize dimension of observation vector
       obs_op_f_pdaf, &              ! Observation operator
       get_obs_f_pdaf                ! Initialize observation vector
  EXTERNAL :: PDAFomi_init_obserr_f_cb ! Initialize mean observation error variance

! !CALLING SEQUENCE:
! Called by: model code  
! Calls: PDAF_generate_obs
!EOP


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  CALL PDAF_generate_obs(collect_state_pdaf, distribute_state_pdaf, &
       init_dim_obs_f_pdaf, obs_op_f_pdaf, PDAFomi_init_obserr_f_cb, get_obs_f_pdaf, &
       prepoststep_pdaf, next_observation_pdaf, outflag)


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

END SUBROUTINE PDAFomi_generate_obs
