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
! !ROUTINE: PDAFlocalomi_put_state_en3dvar_lestkf --- Interface to PDAF for En3D-Var/LESTKF
!
! !INTERFACE:
SUBROUTINE PDAFlocalomi_put_state_en3dvar_lestkf(collect_state_pdaf, &
     init_dim_obs_f_pdaf, obs_op_f_pdaf, &
     cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
     init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
      prepoststep_pdaf, outflag)

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
! 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
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
  EXTERNAL :: cvt_ens_pdaf, &          ! Apply control vector transform matrix to control vector
       cvt_adj_ens_pdaf, &             ! Apply adjoint control vector transform matrix
       obs_op_lin_pdaf, &              ! Linearized observation operator
       obs_op_adj_pdaf                 ! Adjoint observation operator
  EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
       init_dim_l_pdaf, &              ! Init state dimension for local ana. domain
       init_dim_obs_f_pdaf, &          ! Initialize dimension of full observation vector
       obs_op_f_pdaf, &                ! Full observation operator
       init_dim_obs_l_pdaf             ! Initialize local dimimension of obs. vector
  EXTERNAL :: PDAFomi_init_obs_f_cb, & ! Initialize observation vector
       PDAFomi_init_obs_l_cb, &        ! Initialize local observation vector
       PDAFomi_init_obsvar_cb, &       ! Initialize mean observation error variance
       PDAFomi_init_obsvar_l_cb, &     ! Initialize local mean observation error variance
       PDAFomi_prodRinvA_cb, &         ! Provide product R^-1 A
       PDAFomi_g2l_obs_cb, &           ! Restrict full obs. vector to local analysis domain
       PDAFomi_prodRinvA_l_cb, &       ! Provide product R^-1 A on local analysis domain
       PDAFomi_likelihood_l_cb         ! Compute likelihood and apply localization

! !CALLING SEQUENCE:
! Called by: model code  
!EOP


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (TRIM(filterstr) == '3DVAR') THEN
     CALL PDAFlocal_put_state_en3dvar_lestkf(collect_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_prodRinvA_cb, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
          PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, init_dim_l_pdaf, &
          init_dim_obs_l_pdaf,  PDAFomi_g2l_obs_cb, &
          PDAFomi_init_obsvar_cb, PDAFomi_init_obsvar_l_cb, prepoststep_pdaf, outflag)
  ELSE
     WRITE (*,*) 'PDAF-ERROR: No valid filter type for PDAFlocalomi_put_state_3dvar'
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

END SUBROUTINE PDAFlocalomi_put_state_en3dvar_lestkf
