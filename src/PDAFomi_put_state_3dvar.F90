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
! !ROUTINE: PDAFomi_put_state_3dvar --- Interface to PDAF for 3D-Var
!
! !INTERFACE:
SUBROUTINE PDAFomi_put_state_3dvar(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
     cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf, outflag)
  
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
! 2021-04 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filter, ONLY: filterstr
  USE PDAFomi, ONLY: PDAFomi_dealloc

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(inout) :: outflag ! Status flag
  
! ! Names of external subroutines 
  EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
       prepoststep_pdaf                ! User supplied pre/poststep routine
  EXTERNAL :: init_dim_obs_pdaf, &     ! Initialize dimension of observation vector
       obs_op_pdaf, &                  ! Observation operator
       cvt_pdaf, &                     ! Apply control vector transform matrix to control vector
       cvt_adj_pdaf, &                 ! Apply adjoint control vector transform matrix
       obs_op_lin_pdaf, &              ! Linearized observation operator
       obs_op_adj_pdaf                 ! Adjoint observation operator
  EXTERNAL :: PDAFomi_init_obs_f_cb, & ! Initialize observation vector
       PDAFomi_prodRinvA_cb            ! Provide product R^-1 A

! !CALLING SEQUENCE:
! Called by: model code  
!EOP


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (TRIM(filterstr) == '3DVAR') THEN
     CALL PDAF_put_state_3dvar(collect_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_prodRinvA_cb, &
          cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          prepoststep_pdaf, outflag)
  ELSE
     WRITE (*,*) 'PDAF-ERROR: No valid filter type for PDAFomi_put_state_3dvar'
     outflag = 200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

END SUBROUTINE PDAFomi_put_state_3dvar
