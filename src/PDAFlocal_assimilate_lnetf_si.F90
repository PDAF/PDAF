! Copyright (c) 2004-2024 Lars Nerger
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
! !ROUTINE: PDAFlocal_assimilate_lnetf_si --- Interface to transfer state to PDAF
!
! !INTERFACE:
SUBROUTINE PDAFlocal_assimilate_lnetf_si(outflag)

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
! Variant for LNETF with domain decomposition.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2016-11 - Lars Nerger - Initial code
! 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(inout) :: outflag ! Status flag
  
! ! Names of external subroutines 
  EXTERNAL :: collect_state_pdaf, & ! Routine to collect a state vector
       distribute_state_pdaf, &     ! Routine to distribute a state vector
       obs_op_f_pdaf, &             ! Full observation operator
       init_n_domains_pdaf, &       ! Provide number of local analysis domains
       init_dim_l_pdaf, &           ! Init state dimension for local ana. domain
       init_dim_obs_f_pdaf, &       ! Initialize dimension of full observation vector
       init_dim_obs_l_pdaf, &       ! Initialize local dimimension of obs. vector
       init_obs_l_pdaf, &           ! Initialize local observation vector
       g2l_obs_pdaf, &              ! Restrict full obs. vector to local analysis domain
       likelihood_l_pdaf, &         ! Compute observation likelihood for an ensemble member
       prepoststep_pdaf, &          ! User supplied pre/poststep routine
       next_observation_pdaf        ! Routine to provide time step, time and dimension
                                    !   of next observation

! !CALLING SEQUENCE:
! Called by: model code  
! Calls: PDAFlocal_assimilate_lestkf
!EOP


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  CALL PDAFlocal_assimilate_lnetf(collect_state_pdaf, distribute_state_pdaf, &
       init_dim_obs_f_pdaf, obs_op_f_pdaf, init_obs_l_pdaf, prepoststep_pdaf, &
       likelihood_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
        g2l_obs_pdaf, next_observation_pdaf, outflag)

END SUBROUTINE PDAFlocal_assimilate_lnetf_si
