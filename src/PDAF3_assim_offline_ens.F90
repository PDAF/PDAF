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
!> Interfaces to PDAF for offline mode for ensemble filters
!!
!! The interface routines provide the advanced compact
!! interfaces for using PDAF-OMI and PDAF-Local. The routines
!! just call of one the PDAF_assim_offline interface routines
!! with the full interface. In the call the specific PDAF
!! internal subroutines for PDAF-OMI and PDAF-Local are 
!! specified.
!!
!! !  This is a core file of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-04 - Lars Nerger - Initial code by adapting PDAF3_put_state_ens
!! * Other revisions - see repository log
!!
MODULE PDAF3_assim_offline_ens

CONTAINS

!-------------------------------------------------------------------------------
!!> Universal interface routine to PDAF for all filters
!!
!! This variant of the universal routines uses PDAF-local
!!
!! __Revision history:__
!! * 2025-04 - Lars Nerger - Initial code based on PDAF3_put_state_ens
!! Other revisions - see repository log
!!
  SUBROUTINE PDAF3_assim_offline(init_dim_obs_pdaf, obs_op_pdaf, &
       init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
       prepoststep_pdaf, outflag)

    USE PDAF_mod_core, ONLY: filterstr
    USE PDAF_cb_procedures
    USE PDAFomi, ONLY: PDAFomi_dealloc
    USE PDAFlocal, &
         ONLY: PDAFlocal_g2l_cb, &       !< Project global to local state vector
         PDAFlocal_l2g_cb                !< Project local to global state vecto
    USE PDAFassimilate_lseik, ONLY: PDAF_assim_offline_lseik
    USE PDAFassimilate_letkf, ONLY: PDAF_assim_offline_letkf
    USE PDAFassimilate_lestkf, ONLY: PDAF_assim_offline_lestkf
    USE PDAFassimilate_lnetf, ONLY: PDAF_assim_offline_lnetf
    USE PDAFassimilate_lknetf, ONLY: PDAF_assim_offline_lknetf
    USE PDAFassimilate_ensrf, ONLY: PDAF_assim_offline_ensrf
    USE PDAFassimilate_seik, ONLY: PDAF_assim_offline_seik
    USE PDAFassimilate_enkf, ONLY: PDAF_assim_offline_enkf
    USE PDAFassimilate_lenkf, ONLY: PDAF_assim_offline_lenkf
    USE PDAFassimilate_etkf, ONLY: PDAF_assim_offline_etkf
    USE PDAFassimilate_estkf, ONLY: PDAF_assim_offline_estkf
    USE PDAFassimilate_netf, ONLY: PDAF_assim_offline_netf
    USE PDAFassimilate_pf, ONLY: PDAF_assim_offline_pf

    IMPLICIT NONE
  
! *** Arguments ***
    INTEGER, INTENT(inout) :: outflag    !< Status flag

! *** Argument procedures ***
    PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf     !< Initialize dimension of full observation vector
    PROCEDURE(obs_op_cb) :: obs_op_pdaf                 !< Full observation operator
    PROCEDURE(init_n_domains_cb) :: init_n_domains_pdaf !< Provide number of local analysis domains
    PROCEDURE(init_dim_l_cb) :: init_dim_l_pdaf         !< Init state dimension for local ana. domain
    PROCEDURE(init_dim_obs_l_cb) :: init_dim_obs_l_pdaf !< Initialize local dimimension of obs. vector
    PROCEDURE(prepost_cb) :: prepoststep_pdaf           !< User supplied pre/poststep routine

! *** OMI-provided procedures ***
    PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb         !< Initialize full observation vector
    PROCEDURE(init_obs_l_cb) :: PDAFomi_init_obs_l_cb      ! < Initialize local observation vector
    PROCEDURE(init_obsvar_cb) :: PDAFomi_init_obsvar_cb     !< Initialize mean observation error variance
    PROCEDURE(init_obsvar_l_cb) :: PDAFomi_init_obsvar_l_cb !< Initialize local mean observation error variance
    PROCEDURE(init_obsvars_cb) :: PDAFomi_init_obsvars_f_cb !< Initialize vector of observation error variances
    PROCEDURE(localize_serial_cb) :: PDAFomi_localize_covar_serial_cb !< Apply localization to HP and BXY
    PROCEDURE(g2l_obs_cb) :: PDAFomi_g2l_obs_cb             !< Restrict full obs. vector to local analysis domain
    PROCEDURE(prodRinvA_l_cb) :: PDAFomi_prodRinvA_l_cb     !< Provide product R^-1 A on local analysis domain
    PROCEDURE(likelihood_l_cb) :: PDAFomi_likelihood_l_cb   !< Compute likelihood and apply localization
    PROCEDURE(init_obscovar_cb) :: PDAFomi_init_obscovar_cb !< Initialize mean observation error variance
    PROCEDURE(localize_cb) :: PDAFomi_localize_covar_cb     !< Apply localization to HP and HPH^T
    PROCEDURE(add_obs_error_cb) :: PDAFomi_add_obs_error_cb !< Add observation error covariance matrix
    PROCEDURE(prodRinvA_cb) :: PDAFomi_prodRinvA_cb         !< Provide product R^-1 A
    PROCEDURE(likelihood_cb) :: PDAFomi_likelihood_cb       !< Compute likelihood
    PROCEDURE(prodRinvA_hyb_l_cb) :: PDAFomi_prodRinvA_hyb_l_cb     !< Product R^-1 A on local analysis domain with hybrid weight
    PROCEDURE(likelihood_hyb_l_cb) :: PDAFomi_likelihood_hyb_l_cb   !< Compute likelihood and apply localization with tempering



! ******************************************************
! *** Call the full assim_offline interface routine  ***
! ******************************************************

    IF (TRIM(filterstr) == 'LSEIK') THEN
       CALL PDAF_assim_offline_lseik(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
            PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
            PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
            PDAFomi_init_obsvar_l_cb, outflag)
    ELSE IF (TRIM(filterstr) == 'LETKF') THEN
       CALL PDAF_assim_offline_letkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
            PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
            PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
            PDAFomi_init_obsvar_l_cb, outflag)
    ELSE IF (TRIM(filterstr) == 'LESTKF') THEN
       CALL PDAF_assim_offline_lestkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
            PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
            PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
            PDAFomi_init_obsvar_l_cb, outflag)
    ELSE IF (TRIM(filterstr) == 'LNETF') THEN
       CALL PDAF_assim_offline_lnetf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, PDAFomi_likelihood_l_cb, &
            init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf,  &
            PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, PDAFomi_g2l_obs_cb, outflag)
    ELSE IF (TRIM(filterstr) == 'LKNETF') THEN
       CALL PDAF_assim_offline_lknetf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
            PDAFomi_prodRinvA_l_cb, PDAFomi_prodRinvA_hyb_l_cb, &
            init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
            PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
            PDAFomi_init_obsvar_l_cb, PDAFomi_likelihood_l_cb, PDAFomi_likelihood_hyb_l_cb, outflag)
    ELSE IF (TRIM(filterstr) == 'ENSRF') THEN
       CALL PDAF_assim_offline_ensrf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_init_obsvars_f_cb, PDAFomi_localize_covar_serial_cb, &
            prepoststep_pdaf, outflag)
    ELSE IF (TRIM(filterstr) == 'SEIK') THEN
       CALL PDAF_assim_offline_seik(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
            PDAFomi_prodRinvA_cb, PDAFomi_init_obsvar_cb, outflag)
    ELSEIF (TRIM(filterstr) == 'ENKF') THEN
       CALL PDAF_assim_offline_enkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, PDAFomi_add_obs_error_cb, &
            PDAFomi_init_obscovar_cb, outflag)
    ELSEIF (TRIM(filterstr) == 'LENKF') THEN
       CALL PDAF_assim_offline_lenkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, PDAFomi_localize_covar_cb, &
            PDAFomi_add_obs_error_cb, PDAFomi_init_obscovar_cb, outflag)
    ELSEIF (TRIM(filterstr) == 'ETKF') THEN
       CALL PDAF_assim_offline_etkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
            PDAFomi_prodRinvA_cb, PDAFomi_init_obsvar_cb, outflag)
    ELSEIF (TRIM(filterstr) == 'ESTKF') THEN
       CALL PDAF_assim_offline_estkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
            PDAFomi_prodRinvA_cb, PDAFomi_init_obsvar_cb, outflag)
    ELSEIF (TRIM(filterstr) == 'NETF') THEN
       CALL PDAF_assim_offline_netf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, PDAFomi_likelihood_cb, outflag)
    ELSEIF (TRIM(filterstr) == 'PF') THEN
       CALL PDAF_assim_offline_pf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, PDAFomi_likelihood_cb, outflag)
    END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

    CALL PDAFomi_dealloc()

  END SUBROUTINE PDAF3_assim_offline


!-------------------------------------------------------------------------------
!!> Universal interface routine to PDAF for all filters
!!
!! This variant of the universal routines does not use PDAF-local.
!! Compared to PDAF3_assim_offline, there are two additional arguments
!! g2l_state_pdaf and l2g_state_pdaf.
!!
!! __Revision history:__
!! * 2025-04 - Lars Nerger - Initial code based on PDAF3_put_state_ens
!! Other revisions - see repository log
!!
  SUBROUTINE PDAF3_assim_offline_local(init_dim_obs_pdaf, obs_op_pdaf, &
       init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
       g2l_state_pdaf, l2g_state_pdaf, prepoststep_pdaf, outflag)

    USE PDAF_mod_core, ONLY: filterstr
    USE PDAF_cb_procedures
    USE PDAFomi, ONLY: PDAFomi_dealloc
    USE PDAFassimilate_lseik, ONLY: PDAF_assim_offline_lseik
    USE PDAFassimilate_letkf, ONLY: PDAF_assim_offline_letkf
    USE PDAFassimilate_lestkf, ONLY: PDAF_assim_offline_lestkf
    USE PDAFassimilate_lnetf, ONLY: PDAF_assim_offline_lnetf
    USE PDAFassimilate_lknetf, ONLY: PDAF_assim_offline_lknetf
    USE PDAFassimilate_ensrf, ONLY: PDAF_assim_offline_ensrf
    USE PDAFassimilate_seik, ONLY: PDAF_assim_offline_seik
    USE PDAFassimilate_enkf, ONLY: PDAF_assim_offline_enkf
    USE PDAFassimilate_lenkf, ONLY: PDAF_assim_offline_lenkf
    USE PDAFassimilate_etkf, ONLY: PDAF_assim_offline_etkf
    USE PDAFassimilate_estkf, ONLY: PDAF_assim_offline_estkf
    USE PDAFassimilate_netf, ONLY: PDAF_assim_offline_netf
    USE PDAFassimilate_pf, ONLY: PDAF_assim_offline_pf

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: outflag    !< Status flag

! *** Argument procedures ***
    PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf     !< Initialize dimension of full observation vector
    PROCEDURE(obs_op_cb) :: obs_op_pdaf                 !< Full observation operator
    PROCEDURE(init_n_domains_cb) :: init_n_domains_pdaf !< Provide number of local analysis domains
    PROCEDURE(init_dim_l_cb) :: init_dim_l_pdaf         !< Init state dimension for local ana. domain
    PROCEDURE(init_dim_obs_l_cb) :: init_dim_obs_l_pdaf !< Initialize local dimimension of obs. vector
    PROCEDURE(g2l_state_cb) :: g2l_state_pdaf           !< Get local state from full state
    PROCEDURE(l2g_state_cb) :: l2g_state_pdaf           !< Init full state from local state
    PROCEDURE(prepost_cb) :: prepoststep_pdaf           !< User supplied pre/poststep routine

! *** OMI-provided procedures ***
    PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb         !< Initialize full observation vector
    PROCEDURE(init_obs_l_cb) :: PDAFomi_init_obs_l_cb      ! < Initialize local observation vector
    PROCEDURE(init_obsvar_cb) :: PDAFomi_init_obsvar_cb     !< Initialize mean observation error variance
    PROCEDURE(init_obsvar_l_cb) :: PDAFomi_init_obsvar_l_cb !< Initialize local mean observation error variance
    PROCEDURE(init_obsvars_cb) :: PDAFomi_init_obsvars_f_cb !< Initialize vector of observation error variances
    PROCEDURE(localize_serial_cb) :: PDAFomi_localize_covar_serial_cb !< Apply localization to HP and BXY
    PROCEDURE(g2l_obs_cb) :: PDAFomi_g2l_obs_cb             !< Restrict full obs. vector to local analysis domain
    PROCEDURE(prodRinvA_l_cb) :: PDAFomi_prodRinvA_l_cb     !< Provide product R^-1 A on local analysis domain
    PROCEDURE(likelihood_l_cb) :: PDAFomi_likelihood_l_cb   !< Compute likelihood and apply localization
    PROCEDURE(init_obscovar_cb) :: PDAFomi_init_obscovar_cb !< Initialize mean observation error variance
    PROCEDURE(localize_cb) :: PDAFomi_localize_covar_cb     !< Apply localization to HP and HPH^T
    PROCEDURE(add_obs_error_cb) :: PDAFomi_add_obs_error_cb !< Add observation error covariance matrix
    PROCEDURE(prodRinvA_cb) :: PDAFomi_prodRinvA_cb         !< Provide product R^-1 A
    PROCEDURE(likelihood_cb) :: PDAFomi_likelihood_cb       !< Compute likelihood
    PROCEDURE(prodRinvA_hyb_l_cb) :: PDAFomi_prodRinvA_hyb_l_cb     !< Product R^-1 A on local analysis domain with hybrid weight
    PROCEDURE(likelihood_hyb_l_cb) :: PDAFomi_likelihood_hyb_l_cb   !< Compute likelihood and apply localization with tempering


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

    IF (TRIM(filterstr) == 'LSEIK') THEN
       CALL PDAF_assim_offline_lseik(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
            PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
            g2l_state_pdaf, l2g_state_pdaf, PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
            PDAFomi_init_obsvar_l_cb, outflag)
    ELSE IF (TRIM(filterstr) == 'LETKF') THEN
       CALL PDAF_assim_offline_letkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
            PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
            g2l_state_pdaf, l2g_state_pdaf, PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
            PDAFomi_init_obsvar_l_cb, outflag)
    ELSE IF (TRIM(filterstr) == 'LESTKF') THEN
       CALL PDAF_assim_offline_lestkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
            PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
            g2l_state_pdaf, l2g_state_pdaf, PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
            PDAFomi_init_obsvar_l_cb, outflag)
    ELSE IF (TRIM(filterstr) == 'LNETF') THEN
       CALL PDAF_assim_offline_lnetf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, PDAFomi_likelihood_l_cb, &
            init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf,  &
            g2l_state_pdaf, l2g_state_pdaf, PDAFomi_g2l_obs_cb, outflag)
    ELSE IF (TRIM(filterstr) == 'LKNETF') THEN
       CALL PDAF_assim_offline_lknetf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
            PDAFomi_prodRinvA_l_cb, PDAFomi_prodRinvA_hyb_l_cb, &
            init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
            g2l_state_pdaf, l2g_state_pdaf, PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
            PDAFomi_init_obsvar_l_cb, PDAFomi_likelihood_l_cb, PDAFomi_likelihood_hyb_l_cb, outflag)
    ELSE IF (TRIM(filterstr) == 'ENSRF') THEN
       CALL PDAF_assim_offline_ensrf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_init_obsvars_f_cb, PDAFomi_localize_covar_serial_cb, &
            prepoststep_pdaf, outflag)
    ELSE IF (TRIM(filterstr) == 'SEIK') THEN
       CALL PDAF_assim_offline_seik(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
            PDAFomi_prodRinvA_cb, PDAFomi_init_obsvar_cb, outflag)
    ELSEIF (TRIM(filterstr) == 'ENKF') THEN
       CALL PDAF_assim_offline_enkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, PDAFomi_add_obs_error_cb, &
            PDAFomi_init_obscovar_cb, outflag)
    ELSEIF (TRIM(filterstr) == 'LENKF') THEN
       CALL PDAF_assim_offline_lenkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, PDAFomi_localize_covar_cb, &
            PDAFomi_add_obs_error_cb, PDAFomi_init_obscovar_cb, outflag)
    ELSEIF (TRIM(filterstr) == 'ETKF') THEN
       CALL PDAF_assim_offline_etkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
            PDAFomi_prodRinvA_cb, PDAFomi_init_obsvar_cb, outflag)
    ELSEIF (TRIM(filterstr) == 'ESTKF') THEN
       CALL PDAF_assim_offline_estkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
            PDAFomi_prodRinvA_cb, PDAFomi_init_obsvar_cb, outflag)
    ELSEIF (TRIM(filterstr) == 'NETF') THEN
       CALL PDAF_assim_offline_netf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, PDAFomi_likelihood_cb, outflag)
    ELSEIF (TRIM(filterstr) == 'PF') THEN
       CALL PDAF_assim_offline_pf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, PDAFomi_likelihood_cb, outflag)
    END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

    CALL PDAFomi_dealloc()

  END SUBROUTINE PDAF3_assim_offline_local

!-------------------------------------------------------------------------------
!> Interface to PDAF for global filters
!!
!! __Revision history:__
!! * 2025-04 - Lars Nerger - Initial code based on PDAF3_put_state_ens
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF3_assim_offline_global(init_dim_obs_pdaf, obs_op_pdaf, &
       prepoststep_pdaf, outflag)

    USE PDAF_mod_core, ONLY: filterstr
    USE PDAF_cb_procedures
    USE PDAFomi, ONLY: PDAFomi_dealloc
    USE PDAFassimilate_seik, ONLY: PDAF_assim_offline_seik
    USE PDAFassimilate_enkf, ONLY: PDAF_assim_offline_enkf
    USE PDAFassimilate_lenkf, ONLY: PDAF_assim_offline_lenkf
    USE PDAFassimilate_etkf, ONLY: PDAF_assim_offline_etkf
    USE PDAFassimilate_estkf, ONLY: PDAF_assim_offline_estkf
    USE PDAFassimilate_netf, ONLY: PDAF_assim_offline_netf
    USE PDAFassimilate_pf, ONLY: PDAF_assim_offline_pf

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: outflag    !< Status flag

! *** Argument procedures ***
    PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf     !< Initialize dimension of full observation vector
    PROCEDURE(obs_op_cb) :: obs_op_pdaf                 !< Full observation operator
    PROCEDURE(prepost_cb) :: prepoststep_pdaf           !< User supplied pre/poststep routine

! *** OMI-provided procedures ***
    PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb         !< Initialize full observation vector
    PROCEDURE(init_obsvar_cb) :: PDAFomi_init_obsvar_cb     !< Initialize mean observation error variance
    PROCEDURE(init_obscovar_cb) :: PDAFomi_init_obscovar_cb !< Initialize mean observation error variance
    PROCEDURE(localize_cb) :: PDAFomi_localize_covar_cb     !< Apply localization to HP and HPH^T
    PROCEDURE(add_obs_error_cb) :: PDAFomi_add_obs_error_cb !< Add observation error covariance matrix
    PROCEDURE(prodRinvA_cb) :: PDAFomi_prodRinvA_cb         !< Provide product R^-1 A
    PROCEDURE(likelihood_cb) :: PDAFomi_likelihood_cb       !< Compute likelihood


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

    IF (TRIM(filterstr) == 'SEIK') THEN
       CALL PDAF_assim_offline_seik(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
            PDAFomi_prodRinvA_cb, PDAFomi_init_obsvar_cb, outflag)
    ELSEIF (TRIM(filterstr) == 'ENKF') THEN
       CALL PDAF_assim_offline_enkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, PDAFomi_add_obs_error_cb, &
            PDAFomi_init_obscovar_cb, outflag)
    ELSEIF (TRIM(filterstr) == 'LENKF') THEN
       CALL PDAF_assim_offline_lenkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, PDAFomi_localize_covar_cb, &
            PDAFomi_add_obs_error_cb, PDAFomi_init_obscovar_cb, outflag)
    ELSEIF (TRIM(filterstr) == 'ETKF') THEN
       CALL PDAF_assim_offline_etkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
            PDAFomi_prodRinvA_cb, PDAFomi_init_obsvar_cb, outflag)
    ELSEIF (TRIM(filterstr) == 'ESTKF') THEN
       CALL PDAF_assim_offline_estkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
            PDAFomi_prodRinvA_cb, PDAFomi_init_obsvar_cb, outflag)
    ELSEIF (TRIM(filterstr) == 'NETF') THEN
       CALL PDAF_assim_offline_netf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, PDAFomi_likelihood_cb, outflag)
    ELSEIF (TRIM(filterstr) == 'PF') THEN
       CALL PDAF_assim_offline_pf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, prepoststep_pdaf, PDAFomi_likelihood_cb, outflag)
    END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

    CALL PDAFomi_dealloc()

  END SUBROUTINE PDAF3_assim_offline_global

!-------------------------------------------------------------------------------
!> Interface to transfer state to PDAF for LEnKF
!!
!! __Revision history:__
!! * 2025-04 - Lars Nerger - Initial code based on PDAF3_put_state_ens
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF3_assim_offline_lenkf(init_dim_obs_pdaf, obs_op_pdaf, &
       localize_pdaf, prepoststep_pdaf, outflag)

    USE PDAFomi, ONLY: PDAFomi_dealloc
    USE PDAF_cb_procedures
    USE PDAFassimilate_lenkf, ONLY: PDAF_assim_offline_lenkf

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: outflag    !< Status flag

! *** Argument procedures ***
    PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf     !< Initialize dimension of full observation vector
    PROCEDURE(obs_op_cb) :: obs_op_pdaf                 !< Full observation operator
    PROCEDURE(localize_cb) :: localize_pdaf !< Apply localization to HP and HPH^T
    PROCEDURE(prepost_cb) :: prepoststep_pdaf           !< User supplied pre/poststep routine

! *** OMI-provided procedures ***
    PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb         !< Initialize full observation vector
    PROCEDURE(init_obscovar_cb) :: PDAFomi_init_obscovar_cb !< Initialize mean observation error variance
    PROCEDURE(add_obs_error_cb) :: PDAFomi_add_obs_error_cb !< Add observation error covariance matrix


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

    CALL PDAF_assim_offline_lenkf(init_dim_obs_pdaf, obs_op_pdaf, &
         PDAFomi_init_obs_f_cb, prepoststep_pdaf, localize_pdaf, &
         PDAFomi_add_obs_error_cb, PDAFomi_init_obscovar_cb, outflag)


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

    CALL PDAFomi_dealloc()

  END SUBROUTINE PDAF3_assim_offline_lenkf

!-------------------------------------------------------------------------------
!> Interface to transfer state to PDAF for ENSRF
!!
!! __Revision history:__
!! * 2025-04 - Lars Nerger - Initial code based on PDAF3_put_state_ens
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF3_assim_offline_ensrf(init_dim_obs_pdaf, obs_op_pdaf, &
       localize_serial_pdaf, prepoststep_pdaf, outflag)

    USE PDAF_cb_procedures
    USE PDAFomi, ONLY: PDAFomi_dealloc
    USE PDAFassimilate_ensrf, ONLY: PDAF_assim_offline_ensrf

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: outflag    !< Status flag

! *** Argument procedures ***
    PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf     !< Initialize dimension of full observation vector
    PROCEDURE(obs_op_cb) :: obs_op_pdaf                 !< Full observation operator
    PROCEDURE(localize_serial_cb) :: localize_serial_pdaf  !< Apply localization to HP and BXY for single observation
    PROCEDURE(prepost_cb) :: prepoststep_pdaf           !< User supplied pre/poststep routine

! *** OMI-provided procedures ***
    PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb         !< Initialize full observation vector
    PROCEDURE(init_obsvars_cb) :: PDAFomi_init_obsvars_f_cb !< Initialize vector of observation error variances
    PROCEDURE(add_obs_error_cb) :: PDAFomi_add_obs_error_cb !< Add observation error covariance matrix


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

    CALL PDAF_assim_offline_ensrf(init_dim_obs_pdaf, obs_op_pdaf, &
         PDAFomi_init_obs_f_cb, PDAFomi_init_obsvars_f_cb, localize_serial_pdaf, &
         prepoststep_pdaf, outflag)


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

    CALL PDAFomi_dealloc()

  END SUBROUTINE PDAF3_assim_offline_ensrf


!-------------------------------------------------------------------------------
!> Interface to transfer state to PDAF
!!
!! __Revision history:__
!! * 2025-04 - Lars Nerger - Initial code based on PDAF3_put_state_ens
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF3_generate_obs_offline(init_dim_obs_pdaf, &
       obs_op_pdaf, get_obs_pdaf, prepoststep_pdaf, outflag)

    USE PDAF_cb_procedures
    USE PDAFomi, ONLY: PDAFomi_dealloc
    USE PDAFgenerate_obs, ONLY: PDAF_generate_obs_offline

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(out) :: outflag  !< Status flag

! *** Argument procedures ***
    PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf     !< Initialize dimension of full observation vector
    PROCEDURE(obs_op_cb) :: obs_op_pdaf                 !< Full observation operator
    PROCEDURE(get_obs_cb) :: get_obs_pdaf               !< Initialize observation vector
    PROCEDURE(prepost_cb) :: prepoststep_pdaf           !< User supplied pre/poststep routine

! *** OMI-provided procedures ***
    PROCEDURE(init_obserr_cb) :: PDAFomi_init_obserr_f_cb !< Initialize mean observation error variance


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

    CALL PDAF_generate_obs_offline(init_dim_obs_pdaf, obs_op_pdaf, &
         PDAFomi_init_obserr_f_cb, get_obs_pdaf, prepoststep_pdaf, outflag)


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

    CALL PDAFomi_dealloc()

  END SUBROUTINE PDAF3_generate_obs_offline

END MODULE PDAF3_assim_offline_ens
