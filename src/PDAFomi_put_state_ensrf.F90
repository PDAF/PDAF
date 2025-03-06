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
!> Interface to transfer state to PDAF
!!
!! Interface routine called from the model during the 
!! forecast of each ensemble state to transfer data
!! from the model to PDAF and to perform the analysis
!! step.
!!
!! This routine provides the simplified interface
!! where names of user-provided subroutines are
!! fixed. It simply calls the routine with the
!! full interface using pre-defined routine names.
!!
!! Variant for ENSRF with domain decomposition.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2020-06 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAFomi_put_state_ensrf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
     localize_covar_serial_pdaf, prepoststep_pdaf, outflag)

  USE PDAFomi, ONLY: PDAFomi_dealloc

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag  !< Status flag
  
! *** External subroutines ***
  EXTERNAL :: collect_state_pdaf, &    !< Routine to collect a state vector
       prepoststep_pdaf                !< User supplied pre/poststep routine
  EXTERNAL :: init_dim_obs_pdaf, &     !< Initialize dimension of observation vector
       obs_op_pdaf, &                  !< Observation operator
       localize_covar_serial_pdaf      !< Apply localization to HP and HXY
  EXTERNAL :: PDAFomi_init_obs_f_cb, & !< Initialize observation vector
       PDAFomi_init_obsvars_f_cb       !< Initialize vector of observation error variances


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  CALL PDAF_put_state_ensrf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
       PDAFomi_init_obs_f_cb, PDAFomi_init_obsvars_f_cb, localize_covar_serial_pdaf, &
       prepoststep_pdaf, outflag)

! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

END SUBROUTINE PDAFomi_put_state_ensrf
