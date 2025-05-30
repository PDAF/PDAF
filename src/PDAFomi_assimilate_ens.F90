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
!> Interfaces to PDAF for fully-parallel mode using PDAF-OMI
!!
!! The interface routines provide the advanced compact
!! interfaces for using PDAF-OMI. The routines
!! just call of one the PDAF_assimilate interface routines
!! with the full interface. In the call the specific PDAF
!! internal subroutines for PDAF-OMI are specified.
!!
!! !  This is a core file of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code by collecting files into a module
!! * Other revisions - see repository log
!!
MODULE PDAFomi_assimilate_ens

CONTAINS

!-------------------------------------------------------------------------------
!> Interface to transfer state to PDAF
!!
!! __Revision history:__
!! 2020-11 - Lars Nerger - Initial code
!! Other revisions - see repository log
!!
SUBROUTINE PDAFomi_assimilate_local(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, prepoststep_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
          next_observation_pdaf, outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAFomi_obs_l, ONLY: PDAFomi_dealloc
  USE PDAFassimilate_lseik, ONLY: PDAF_assimilate_lseik
  USE PDAFassimilate_letkf, ONLY: PDAF_assimilate_letkf
  USE PDAFassimilate_lestkf, ONLY: PDAF_assimilate_lestkf
  USE PDAFassimilate_lnetf, ONLY: PDAF_assimilate_lnetf
  USE PDAFassimilate_lknetf, ONLY: PDAF_assimilate_lknetf
  USE PDAFassimilate_ensrf, ONLY: PDAF_assimilate_ensrf

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag      !< Status flag
  
! *** External subroutines ***
  EXTERNAL :: collect_state_pdaf, &          !< Routine to collect a state vector
       distribute_state_pdaf, &              !< Routine to distribute a state vector
       next_observation_pdaf, &              !< Provide time step, time and dimension of next observation
       prepoststep_pdaf                      !< User supplied pre/poststep routine
  EXTERNAL :: init_n_domains_pdaf, &         !< Provide number of local analysis domains
       init_dim_l_pdaf, &                    !< Init state dimension for local ana. domain
       g2l_state_pdaf, &                     !< Get state on local ana. domain from full state
       l2g_state_pdaf, &                     !< Init full state from local state
       init_dim_obs_f_pdaf, &                !< Initialize dimension of full observation vector
       obs_op_f_pdaf, &                      !< Full observation operator
       init_dim_obs_l_pdaf                   !< Initialize local dimimension of obs. vector
  EXTERNAL :: PDAFomi_init_obs_f_cb, &       !< Initialize full observation vector
       PDAFomi_init_obs_l_cb, &              !< Initialize local observation vector
       PDAFomi_init_obsvar_cb, &             !< Initialize mean observation error variance
       PDAFomi_init_obsvar_l_cb, &           !< Initialize local mean observation error variance
       PDAFomi_init_obsvars_f_cb, &          !< Initialize vector of observation error variances
       PDAFomi_localize_covar_serial_cb, &   !< Apply localization to HP and BXY
       PDAFomi_g2l_obs_cb, &                 !< Restrict full obs. vector to local analysis domain
       PDAFomi_prodRinvA_l_cb, &             !< Provide product R^-1 A on local analysis domain
       PDAFomi_likelihood_l_cb               !< Compute likelihood and apply localization
  EXTERNAL :: PDAFomi_prodRinvA_hyb_l_cb, &  !< Product R^-1 A on local analysis domain with hybrid weight
       PDAFomi_likelihood_hyb_l_cb           !< Compute likelihood and apply localization with tempering


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_assimilate_local -- START'

  IF (TRIM(filterstr) == 'LSEIK') THEN
     CALL PDAF_assimilate_lseik(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
          prepoststep_pdaf, PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
          PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, PDAFomi_init_obsvar_l_cb, &
          next_observation_pdaf, outflag)
  ELSE IF (TRIM(filterstr) == 'LETKF') THEN
     CALL PDAF_assimilate_letkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
          prepoststep_pdaf, PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
          PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, PDAFomi_init_obsvar_l_cb, &
          next_observation_pdaf, outflag)
  ELSE IF (TRIM(filterstr) == 'LESTKF') THEN
     CALL PDAF_assimilate_lestkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
          prepoststep_pdaf, PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
          PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, PDAFomi_init_obsvar_l_cb, &
          next_observation_pdaf, outflag)
  ELSE IF (TRIM(filterstr) == 'LNETF') THEN
     CALL PDAF_assimilate_lnetf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
          prepoststep_pdaf, PDAFomi_likelihood_l_cb, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
          PDAFomi_g2l_obs_cb, next_observation_pdaf, outflag)
  ELSE IF (TRIM(filterstr) == 'LKNETF') THEN
     CALL PDAF_assimilate_lknetf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
          PDAFomi_prodRinvA_l_cb, PDAFomi_prodRinvA_hyb_l_cb, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          g2l_state_pdaf, l2g_state_pdaf, PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
          PDAFomi_init_obsvar_l_cb, PDAFomi_likelihood_l_cb, PDAFomi_likelihood_hyb_l_cb, &
          next_observation_pdaf, outflag)
  ELSE IF (TRIM(filterstr) == 'ENSRF') THEN
     CALL PDAF_assimilate_ensrf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_f_pdaf, obs_op_f_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obsvars_f_cb, &
          PDAFomi_localize_covar_serial_cb, prepoststep_pdaf, next_observation_pdaf, outflag)
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_assimilate_local -- END'

END SUBROUTINE PDAFomi_assimilate_local


!-------------------------------------------------------------------------------
!> Interface to PDAF for global filters
!!
!! __Revision history:__
!! * 2020-11 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAFomi_assimilate_global(collect_state_pdaf, distribute_state_pdaf, &
     init_dim_obs_pdaf, obs_op_pdaf, prepoststep_pdaf, next_observation_pdaf, outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAFomi_obs_l, ONLY: PDAFomi_dealloc
  USE PDAFassimilate_seik, ONLY: PDAF_assimilate_seik
  USE PDAFassimilate_enkf, ONLY: PDAF_assimilate_enkf
  USE PDAFassimilate_lenkf, ONLY: PDAF_assimilate_lenkf
  USE PDAFassimilate_etkf, ONLY: PDAF_assimilate_etkf
  USE PDAFassimilate_estkf, ONLY: PDAF_assimilate_estkf
  USE PDAFassimilate_netf, ONLY: PDAF_assimilate_netf
  USE PDAFassimilate_pf, ONLY: PDAF_assimilate_pf

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
       PDAFomi_init_obscovar_cb, &     !< Initialize mean observation error variance
       PDAFomi_localize_covar_cb, &    !< Apply localization to HP and HPH^T
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
  ELSEIF (TRIM(filterstr) == 'LENKF') THEN
     CALL PDAF_assimilate_lenkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
          PDAFomi_localize_covar_cb, PDAFomi_add_obs_error_cb, PDAFomi_init_obscovar_cb, &
          next_observation_pdaf, outflag)
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

END SUBROUTINE PDAFomi_assimilate_global


!-------------------------------------------------------------------------------
!> Interface to transfer state to PDAF
!!
!! __Revision history:__
!! * 2020-06 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAFomi_assimilate_lenkf(collect_state_pdaf, distribute_state_pdaf, &
     init_dim_obs_pdaf, obs_op_pdaf, prepoststep_pdaf, localize_covar_pdaf, &
     next_observation_pdaf, outflag)

  USE PDAF_mod_core, ONLY: debug
  USE PDAFomi_obs_l, ONLY: PDAFomi_dealloc
  USE PDAFassimilate_lenkf, ONLY: PDAF_assimilate_lenkf

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

END SUBROUTINE PDAFomi_assimilate_lenkf


!-------------------------------------------------------------------------------
!> Interface to transfer state to PDAF
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAFomi_assimilate_ensrf(collect_state_pdaf, distribute_state_pdaf, &
     init_dim_obs_pdaf, obs_op_pdaf, prepoststep_pdaf, localize_covar_serial_pdaf, &
     next_observation_pdaf, outflag)

  USE PDAF_mod_core, ONLY: debug
  USE PDAFomi_obs_l, ONLY: PDAFomi_dealloc
  USE PDAFassimilate_ensrf, ONLY: PDAF_assimilate_ensrf

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
       localize_covar_serial_pdaf      !< Apply localization to HP and BXY
  EXTERNAL :: PDAFomi_init_obs_f_cb, & !< Initialize observation vector
       PDAFomi_init_obsvars_f_cb       !< Initialize vector of observation error variances


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_assimilate_ensrf -- START'

  CALL PDAF_assimilate_ensrf(collect_state_pdaf, distribute_state_pdaf, &
       init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obsvars_f_cb, &
          localize_covar_serial_pdaf, prepoststep_pdaf, next_observation_pdaf, outflag)


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_assimilate_ensrf -- END'

END SUBROUTINE PDAFomi_assimilate_ensrf

!-------------------------------------------------------------------------------
!> Interface to transfer state to PDAF
!!
!! Interface routine called from the model during the 
!! forecast of each ensemble state to transfer data
!! from the model to PDAF and to perform the analysis
!! step.
!!
!! This routine provides the simplified interface
!! where names of user-provided subroutines are
!! fixed. It simply calls the routine with the
!! full interface using pre-defined routine names.
!!
!! Variant for gnerating observations with domain
!! decomposition.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! 2019-01 - Lars Nerger - Initial code
!! Other revisions - see repository log
!!
SUBROUTINE PDAFomi_generate_obs(collect_state_pdaf, distribute_state_pdaf, &
       init_dim_obs_f_pdaf, obs_op_f_pdaf, get_obs_f_pdaf, prepoststep_pdaf, &
       next_observation_pdaf, outflag)

  USE PDAFomi_obs_l, ONLY: PDAFomi_dealloc
  USE PDAFgenerate_obs, ONLY: PDAF_generate_obs

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(inout) :: outflag !< Status flag
  
! *** Names of external subroutines ***
  EXTERNAL :: collect_state_pdaf, &    !< Routine to collect a state vector
       distribute_state_pdaf, &        !< Routine to distribute a state vector
       prepoststep_pdaf, &             !< User supplied pre/poststep routine
       next_observation_pdaf           !< Provide time step, time and dimension of next observation
  EXTERNAL :: init_dim_obs_f_pdaf, &   !< Initialize dimension of observation vector
       obs_op_f_pdaf, &                !< Observation operator
       get_obs_f_pdaf                  !< Initialize observation vector
  EXTERNAL :: PDAFomi_init_obserr_f_cb !< Initialize mean observation error variance


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  CALL PDAF_generate_obs(collect_state_pdaf, distribute_state_pdaf, &
       init_dim_obs_f_pdaf, obs_op_f_pdaf, PDAFomi_init_obserr_f_cb, get_obs_f_pdaf, &
       prepoststep_pdaf, next_observation_pdaf, outflag)


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

END SUBROUTINE PDAFomi_generate_obs

END MODULE PDAFomi_assimilate_ens
