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
!> Interfaces to PDAF for fully-parallel mode
!!
!! The interface routines provide the advanced compact
!! interfaces for using PDAF-OMI and PDAF-Local. The routines
!! just call of one the PDAF_assimilate interface routines
!! with the full interface. In the call the specific PDAF
!! internal subroutines for PDAF-OMI and PDAF-Local are 
!! specified.
!!
!! !  This is a core file of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code by collecting files into a module
!! * Other revisions - see repository log
!!
MODULE PDAF_assimilate

CONTAINS
!!> Interface routine to PDAF for local filters
!!
!! __Revision history:__
!! * 2020-11 - Lars Nerger - Initial code
!! * 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
!! Other revisions - see repository log
!!
SUBROUTINE PDAF3_assimilate_local(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, prepoststep_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdaf,  &
          next_observation_pdaf, outflag)

  USE PDAF_mod_filter, ONLY: filterstr, debug
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFlocal, &
       ONLY: PDAFlocal_g2l_cb, &       !< Project global to local state vector
       PDAFlocal_l2g_cb                !< Project local to global state vecto

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(inout) :: outflag    !< Status flag

! *** Names of external subroutines ***
  EXTERNAL :: collect_state_pdaf, &    !< Routine to collect a state vector
       distribute_state_pdaf, &        !< Routine to distribute a state vector
       next_observation_pdaf, &        !< Provide time step, time and dimension of next observation
       prepoststep_pdaf                !< User supplied pre/poststep routine
  EXTERNAL :: init_n_domains_pdaf, &   !< Provide number of local analysis domains
       init_dim_l_pdaf, &              !< Init state dimension for local ana. domain
       init_dim_obs_f_pdaf, &          !< Initialize dimension of full observation vector
       obs_op_f_pdaf, &                !< Full observation operator
       init_dim_obs_l_pdaf             !< Initialize local dimimension of obs. vector
  EXTERNAL :: PDAFomi_init_obs_f_cb, & !< Initialize full observation vector
       PDAFomi_init_obs_l_cb, &        !< Initialize local observation vector
       PDAFomi_init_obsvar_cb, &       !< Initialize mean observation error variance
       PDAFomi_init_obsvar_l_cb, &     !< Initialize local mean observation error variance
       PDAFomi_init_obsvars_f_cb, &    !< Initialize vector of observation error variances
       PDAFomi_g2l_obs_cb, &           !< Restrict full obs. vector to local analysis domain
       PDAFomi_prodRinvA_l_cb, &       !< Provide product R^-1 A on local analysis domain
       PDAFomi_likelihood_l_cb         !< Compute likelihood and apply localization
  EXTERNAL :: PDAFomi_prodRinvA_hyb_l_cb, &  !< Product R^-1 A on local analysis domain with hybrid weight
       PDAFomi_likelihood_hyb_l_cb     !< Compute likelihood and apply localization with tempering


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFlocalomi_assimilate -- START'

  IF (TRIM(filterstr) == 'LSEIK') THEN
     CALL PDAF_assimilate_lseik(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
          prepoststep_pdaf, PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdaf, PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, &
          PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, PDAFomi_init_obsvar_l_cb, &
          next_observation_pdaf, outflag)
  ELSE IF (TRIM(filterstr) == 'LETKF') THEN
     CALL PDAF_assimilate_letkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
          prepoststep_pdaf, PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdaf, PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, &
          PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, PDAFomi_init_obsvar_l_cb, &
          next_observation_pdaf, outflag)
  ELSE IF (TRIM(filterstr) == 'LESTKF') THEN
     CALL PDAF_assimilate_lestkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
          prepoststep_pdaf, PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdaf, PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, &
          PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, PDAFomi_init_obsvar_l_cb, &
          next_observation_pdaf, outflag)
  ELSE IF (TRIM(filterstr) == 'LNETF') THEN
     CALL PDAF_assimilate_lnetf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
          prepoststep_pdaf, PDAFomi_likelihood_l_cb, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdaf, PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, &
          PDAFomi_g2l_obs_cb, next_observation_pdaf, outflag)
  ELSE IF (TRIM(filterstr) == 'LKNETF') THEN
     CALL PDAF_assimilate_lknetf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
          PDAFomi_prodRinvA_l_cb, PDAFomi_prodRinvA_hyb_l_cb, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
          PDAFomi_init_obsvar_l_cb, PDAFomi_likelihood_l_cb, PDAFomi_likelihood_hyb_l_cb, &
          next_observation_pdaf, outflag)
  ELSEIF (TRIM(filterstr) == 'ENSRF') THEN
     CALL PDAF_assimilate_ensrf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obsvars_f_cb, &
          prepoststep_pdaf, next_observation_pdaf, outflag)
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFlocalomi_assimilate -- END'

END SUBROUTINE PDAF3_assimilate_local

!> Interface to PDAF for global filters
!!
!! __Revision history:__
!! * 2020-11 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF3_assimilate_global(collect_state_pdaf, distribute_state_pdaf, &
     init_dim_obs_pdaf, obs_op_pdaf, prepoststep_pdaf, next_observation_pdaf, outflag)

  USE PDAF_mod_filter, ONLY: filterstr, debug
  USE PDAFomi, ONLY: PDAFomi_dealloc

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag      !< Status flag
  
! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: collect_state_pdaf, &    !< Routine to collect a state vector
       distribute_state_pdaf, &        !< Routine to distribute a state vector
       next_observation_pdaf, &        !< Provide time step, time and dimension of next observation
       prepoststep_pdaf                !< User supplied pre/poststep routine
  EXTERNAL :: init_dim_obs_pdaf, &     !< Initialize dimension of observation vector
       obs_op_pdaf                     !< Observation operator
  EXTERNAL :: PDAFomi_init_obs_f_cb, & !< Initialize observation vector
       PDAFomi_init_obsvar_cb, &       !< Initialize mean observation error variance
       PDAFomi_init_obsvars_f_cb, &    !< Initialize vector of observation error variances
       PDAFomi_init_obscovar_cb, &     !< Initialize mean observation error variance
       PDAFomi_add_obs_error_cb, &     !< Add observation error covariance matrix
       PDAFomi_prodRinvA_cb, &         !< Provide product R^-1 A
       PDAFomi_likelihood_cb           !< Compute likelihood


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_assimilate_global -- START'

  IF (TRIM(filterstr) == 'SEIK') THEN
     CALL PDAF_assimilate_seik(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
          PDAFomi_prodRinvA_cb, PDAFomi_init_obsvar_cb, next_observation_pdaf, outflag)
  ELSEIF (TRIM(filterstr) == 'ENKF') THEN
     CALL PDAF_assimilate_enkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
          PDAFomi_add_obs_error_cb, PDAFomi_init_obscovar_cb, next_observation_pdaf, outflag)
  ELSEIF (TRIM(filterstr) == 'ETKF') THEN
     CALL PDAF_assimilate_etkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
          PDAFomi_prodRinvA_cb, PDAFomi_init_obsvar_cb, next_observation_pdaf, outflag)
  ELSEIF (TRIM(filterstr) == 'ESTKF') THEN
     CALL PDAF_assimilate_estkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
          PDAFomi_prodRinvA_cb, PDAFomi_init_obsvar_cb, next_observation_pdaf, outflag)
  ELSEIF (TRIM(filterstr) == 'NETF') THEN
     CALL PDAF_assimilate_netf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
          PDAFomi_likelihood_cb, next_observation_pdaf, outflag)
  ELSEIF (TRIM(filterstr) == 'PF') THEN
     CALL PDAF_assimilate_pf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
          PDAFomi_likelihood_cb, next_observation_pdaf, outflag)
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_assimilate_global -- END'

END SUBROUTINE PDAF3_assimilate_global

!> Interface to transfer state to PDAF
!!
!! __Revision history:__
!! * 2020-06 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF3_assimilate_lenkf(collect_state_pdaf, distribute_state_pdaf, &
     init_dim_obs_pdaf, obs_op_pdaf, prepoststep_pdaf, localize_covar_pdaf, &
     next_observation_pdaf, outflag)

  USE PDAF_mod_filter, ONLY: debug
  USE PDAFomi, ONLY: PDAFomi_dealloc

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag      !< Status flag
  
! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: collect_state_pdaf, &    !< Routine to collect a state vector
       distribute_state_pdaf, &        !< Routine to distribute a state vector
       next_observation_pdaf, &        !< Provide time step, time and dimension of next observation
       prepoststep_pdaf                !< User supplied pre/poststep routine
  EXTERNAL :: init_dim_obs_pdaf, &     !< Initialize dimension of observation vector
       obs_op_pdaf, &                  !< Observation operator
       localize_covar_pdaf             !< Apply localization to HP and HPH^T
  EXTERNAL :: PDAFomi_init_obs_f_cb, & !< Initialize observation vector
       PDAFomi_init_obscovar_cb, &     !< Initialize mean observation error variance
       PDAFomi_add_obs_error_cb        !< Provide product R^-1 A


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_assimilate_lenkf -- START'

  CALL PDAF_assimilate_lenkf(collect_state_pdaf, distribute_state_pdaf, &
       init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
       localize_covar_pdaf, PDAFomi_add_obs_error_cb, PDAFomi_init_obscovar_cb, &
       next_observation_pdaf, outflag)


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_assimilate_lenkf -- END'

END SUBROUTINE PDAF3_assimilate_lenkf

END MODULE PDAF_assimilate
