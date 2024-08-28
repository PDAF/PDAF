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
! !ROUTINE: PDAFlocalomi_put_state_nondiagR --- Interface to transfer state to PDAF
!
! !INTERFACE:
SUBROUTINE PDAFlocalomi_put_state_nondiagR(collect_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdafomi, prodRinvA_l_pdafomi, &
           outflag)

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
! The routine supports all domain-localized filters.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2024-07 - Lars Nerger - Initial code
! 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filter, ONLY: filterstr, debug
  USE PDAFomi, ONLY: PDAFomi_dealloc

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(inout) :: outflag ! Status flag
  
! ! Names of external subroutines 
  EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector
       prepoststep_pdaf                ! User supplied pre/poststep routine
  EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
       init_dim_l_pdaf                 ! Init state dimension for local ana. domain
  EXTERNAL :: init_dim_obs_pdafomi, &  ! Initialize dimension of full observation vector
       obs_op_pdafomi, &               ! Full observation operator
       init_dim_obs_l_pdafomi, &       ! Initialize local dimimension of obs. vector
       prodRinvA_l_pdafomi             ! Provide product of inverse of R with matrix A
  EXTERNAL :: PDAFomi_init_obs_f_cb, & ! Initialize full observation vector
       PDAFomi_init_obs_l_cb, &        ! Initialize local observation vector
       PDAFomi_init_obsvar_cb, &       ! Initialize mean observation error variance
       PDAFomi_init_obsvar_l_cb, &     ! Initialize local mean observation error variance
       PDAFomi_g2l_obs_cb              ! Restrict full obs. vector to local analysis domain

! !CALLING SEQUENCE:
! Called by: model code  
!EOP


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFlocalomi_put_state_nondiagR -- START'

  IF (TRIM(filterstr) == 'LSEIK') THEN
     CALL PDAFlocal_put_state_lseik(collect_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
          prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
           PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
          PDAFomi_init_obsvar_l_cb, outflag)
  ELSE IF (TRIM(filterstr) == 'LETKF') THEN
     CALL PDAFlocal_put_state_letkf(collect_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
          prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
           PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
          PDAFomi_init_obsvar_l_cb, outflag)
  ELSE IF (TRIM(filterstr) == 'LESTKF') THEN
     CALL PDAFlocal_put_state_lestkf(collect_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
          prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
           PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
          PDAFomi_init_obsvar_l_cb, outflag)
  ELSE IF (TRIM(filterstr) == 'LNETF') THEN
     WRITE (*,*) 'PDAF-ERROR: Use PDAFlocalomi_put_state_lnetf_nondiagR for LNETF'
     outflag=200
  ELSE IF (TRIM(filterstr) == 'LKNETF') THEN
     WRITE (*,*) 'PDAF-ERROR: Use PDAFlocalomi_put_state_lknetf_nondiagR for LKNETF'
     outflag=200
  ELSE
     WRITE (*,*) 'PDAF-ERROR: Invalid filter choice for PDAFlocalomi_put_state_nondiagR'
     outflag=200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFlocalomi_put_state_nondiagR -- END'

END SUBROUTINE PDAFlocalomi_put_state_nondiagR
