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
! !ROUTINE: PDAFlocalomi_put_state_si --- Interface to transfer state to PDAF
!
! !INTERFACE:
SUBROUTINE PDAFlocalomi_put_state_si(outflag)

! !DESCRIPTION:
! Interface routine called from the model after the 
! forecast of each ensemble state to transfer data
! from the model to PDAF. 
!
! This routine provides the simplified interface
! where names of user-provided subroutines are
! fixed. It simply calls the routine with the
! full interface using pre-defined routine names.
!
! The routine supports all domain-localized filters.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2021-10 - Lars Nerger - Initial code
! 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(inout) :: outflag  ! Status flag

! ! Names of external subroutines 
  EXTERNAL :: collect_state_pdaf, &  ! Routine to collect a state vector
       prepoststep_pdaf              ! User supplied pre/poststep routine
  ! Localization of state vector
  EXTERNAL :: init_n_domains_pdaf, & ! Provide number of local analysis domains
       init_dim_l_pdaf               ! Initialize state dimension for local analysis domain
  ! Interface to PDAF-OMI for local and global filters
  EXTERNAL :: &
       init_dim_obs_pdafomi, &       ! Get dimension of full obs. vector for PE-local domain
       obs_op_pdafomi, &             ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_l_pdafomi, &     ! Get dimension of obs. vector for local analysis domain
       localize_covar_pdafomi        ! Apply localization to covariance matrix in LEnKF

! !CALLING SEQUENCE:
! Called by: model code  
! Calls: PDAFlocalomi_put_state
!EOP


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  CALL PDAFlocalomi_put_state(collect_state_pdaf, init_dim_obs_pdafomi, &
       obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
       init_dim_obs_l_pdafomi,  outflag)

END SUBROUTINE PDAFlocalomi_put_state_si
