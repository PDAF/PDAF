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
!$Id: PDAF-D_assimilate_lenkf_si.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_assimilate_lenkf_si --- Interface to transfer state to PDAF
!
! !INTERFACE:
SUBROUTINE PDAF_assimilate_lenkf_si(outflag)

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
! Variant for local ENKF with domain decomposition.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2015-12 - Lars Nerger - Initial code by copying PDAF_assimilate_enkf_si
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(inout) :: outflag ! Status flag
  
! ! Names of external subroutines 
  EXTERNAL :: collect_state_pdaf, & ! Routine to collect a state vector
       distribute_state_pdaf, &     ! Routine to distribute a state vector
       init_dim_obs_pdaf, &         ! Initialize dimension of observation vector
       obs_op_pdaf, &               ! Observation operator
       init_obscovar_pdaf, &        ! Initialize obs. error cov. matrix R in EnKF
       init_obs_pdaf, &             ! Initialize observation vector
       prepoststep_pdaf, &          ! User supplied pre/poststep routine
       add_obs_error_pdaf, &        ! Add obs error covariance R to HPH in EnKF
       localize_covar_pdaf, &       ! Apply localization to HP and HPH^T
       next_observation_pdaf        ! Routine to provide time step, time and dimension
                                    !   of next observation

! !CALLING SEQUENCE:
! Called by: model code  
! Calls: PDAF_assimilate_enkf
!EOP


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  CALL PDAF_assimilate_lenkf(collect_state_pdaf, distribute_state_pdaf, &
       init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf, prepoststep_pdaf, &
       localize_covar_pdaf, add_obs_error_pdaf, init_obscovar_pdaf, &
       next_observation_pdaf, outflag)

END SUBROUTINE PDAF_assimilate_lenkf_si
