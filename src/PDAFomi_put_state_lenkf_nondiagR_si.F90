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
! !ROUTINE: PDAFomi_put_state_lenkf_nondiagR_si --- Interface to transfer state to PDAF
!
! !INTERFACE:
SUBROUTINE PDAFomi_put_state_lenkf_nondiagR_si(outflag)

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
! The routine supports all global filters.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2024-08 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(inout) :: outflag ! Status flag
  
! ! Names of external subroutines 
  EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
       prepoststep_pdaf                ! User supplied pre/poststep routine
  ! Interface to PDAF-OMI for local and global filters
  EXTERNAL :: init_dim_obs_pdafomi, &  ! Get dimension of full obs. vector for PE-local domain
       obs_op_pdafomi, &               ! Obs. operator for full obs. vector for PE-local domain
       localize_covar_pdafomi, &       ! Apply localization to HP and HPH^T
       init_obscovar_pdafomi, &        ! Initialize mean observation error variance
       add_obs_error_pdafomi           ! Provide product R^-1 A

! !CALLING SEQUENCE:
! Called by: model code  
! Calls: PDAFomi_put_state_lenkf_nondiagR
!EOP


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  CALL PDAFomi_put_state_lenkf_nondiagR(collect_state_pdaf, &
       init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, localize_covar_pdafomi, &
       add_obs_error_pdafomi, init_obscovar_pdafomi, outflag)

END SUBROUTINE PDAFomi_put_state_lenkf_nondiagR_si
