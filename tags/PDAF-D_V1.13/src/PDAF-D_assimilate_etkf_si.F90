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
!$Id: PDAF-D_assimilate_etkf_si.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_assimilate_etkf_si --- Interface to transfer state to PDAF
!
! !INTERFACE:
SUBROUTINE PDAF_assimilate_etkf_si(outflag)

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
! Variant for ETKF with domain decomposition.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2013-08 - Lars Nerger - Initial code
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
       init_obsvar_pdaf, &          ! Initialize mean observation error variance
       init_obs_pdaf, &             ! Initialize observation vector
       prepoststep_pdaf, &          ! User supplied pre/poststep routine
       prodRinvA_pdaf, &            ! Provide product R^-1 A
       next_observation_pdaf        ! Routine to provide time step, time and dimension
                                    !   of next observation

! !CALLING SEQUENCE:
! Called by: model code  
! Calls: PDAF_assimilate_etkf
!EOP


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  CALL PDAF_assimilate_etkf(collect_state_pdaf, distribute_state_pdaf, &
       init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf, prepoststep_pdaf, &
       prodRinvA_pdaf, init_obsvar_pdaf, next_observation_pdaf, outflag)

END SUBROUTINE PDAF_assimilate_etkf_si
