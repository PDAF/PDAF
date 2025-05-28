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
!! This variant of the interfaces is for non-diagonal R-matrices.
!! To support this non-diagonal matrix the observation-related
!! routine, prodRinvA_pdaf or likelihood_pdaf is includes
!! as an argument.
!!
!! !  This is a core file of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code by collecting files into a module
!! * Other revisions - see repository log
!!
MODULE PDAF3_assimilate_ens_nondiagR

CONTAINS

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
!! The routine supports all domain-localized filters.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2024-07 - Lars Nerger - Initial code
!! * 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF3_assimilate_local_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdaf, prodRinvA_l_pdaf, &
          prepoststep_pdaf, next_observation_pdaf, outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug 
  USE PDAF_cb_procedures
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFlocal, ONLY: PDAFlocal_g2l_cb, PDAFlocal_l2g_cb 
  USE PDAFassimilate_lseik, ONLY: PDAF_assimilate_lseik
  USE PDAFassimilate_letkf, ONLY: PDAF_assimilate_letkf
  USE PDAFassimilate_lestkf, ONLY: PDAF_assimilate_lestkf

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(inout) :: outflag    !< Status flag

! *** Argument procedures ***
  PROCEDURE(collect_cb) :: collect_state_pdaf         !< Routine to collect a state vector
  PROCEDURE(distribute_cb) :: distribute_state_pdaf   !< Routine to distribute a state vector
  PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf     !< Initialize dimension of full observation vector
  PROCEDURE(obs_op_cb) :: obs_op_pdaf                 !< Full observation operator
  PROCEDURE(init_n_domains_cb) :: init_n_domains_pdaf !< Provide number of local analysis domains
  PROCEDURE(init_dim_l_cb) :: init_dim_l_pdaf         !< Init state dimension for local ana. domain
  PROCEDURE(init_dim_obs_l_cb) :: init_dim_obs_l_pdaf !< Initialize local dimimension of obs. vector
  PROCEDURE(prepost_cb) :: prepoststep_pdaf           !< User supplied pre/poststep routine
  PROCEDURE(next_obs_cb) :: next_observation_pdaf     !< Provide information on next forecast
  PROCEDURE(prodRinvA_l_cb) :: prodRinvA_l_pdaf       !< Provide product of inverse of R with matrix A

! *** OMI-provided procedures ***
  PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb         !< Initialize full observation vector
  PROCEDURE(init_obs_l_cb) :: PDAFomi_init_obs_l_cb       !< Initialize local observation vector
  PROCEDURE(init_obsvar_cb) :: PDAFomi_init_obsvar_cb     !< Initialize mean observation error variance
  PROCEDURE(init_obsvar_l_cb) :: PDAFomi_init_obsvar_l_cb !< Initialize local mean observation error variance
  PROCEDURE(g2l_obs_cb) :: PDAFomi_g2l_obs_cb             !< Restrict full obs. vector to local analysis domain


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAF3_assimilate_nondiagR -- START'

  IF (TRIM(filterstr) == 'LSEIK') THEN
     CALL PDAF_assimilate_lseik(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
          prepoststep_pdaf, prodRinvA_l_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdaf, PDAFlocal_g2l_cb, PDAFlocal_l2g_cb,  &
          PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, PDAFomi_init_obsvar_l_cb, &
          next_observation_pdaf, outflag)
  ELSE IF (TRIM(filterstr) == 'LETKF') THEN
     CALL PDAF_assimilate_letkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
          prepoststep_pdaf, prodRinvA_l_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdaf, PDAFlocal_g2l_cb, PDAFlocal_l2g_cb,  &
          PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, PDAFomi_init_obsvar_l_cb, &
          next_observation_pdaf, outflag)
  ELSE IF (TRIM(filterstr) == 'LESTKF') THEN
     CALL PDAF_assimilate_lestkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
          prepoststep_pdaf, prodRinvA_l_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdaf, PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, &
          PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, PDAFomi_init_obsvar_l_cb, &
          next_observation_pdaf, outflag)
  ELSE IF (TRIM(filterstr) == 'LNETF') THEN
     WRITE (*,*) 'PDAF-ERROR: Use PDAF3_assimilate_lnetf_nondiagR for LNETF'
     outflag=200
  ELSE IF (TRIM(filterstr) == 'LKNETF') THEN
     WRITE (*,*) 'PDAF-ERROR: Use PDAF3_assimilate_lknetf_nondiagR for LKNETF'
     outflag=200
  ELSE
     WRITE (*,*) 'PDAF-ERROR: Invalid filter choice for PDAF3_assimilate_nondiagR'
     outflag=200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAF3_assimilate_nondiagR -- END'

END SUBROUTINE PDAF3_assimilate_local_nondiagR


!-------------------------------------------------------------------------------
!> Interface to PDAF for global filters
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF3_assimilate_global_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
     init_dim_obs_pdaf, obs_op_pdaf, prodRinvA_pdaf, &
     prepoststep_pdaf, next_observation_pdaf, outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAF_cb_procedures
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFassimilate_seik, ONLY: PDAF_assimilate_seik
  USE PDAFassimilate_etkf, ONLY: PDAF_assimilate_etkf
  USE PDAFassimilate_estkf, ONLY: PDAF_assimilate_estkf

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag      !< Status flag
  
! *** Argument procedures ***
  PROCEDURE(collect_cb) :: collect_state_pdaf         !< Routine to collect a state vector
  PROCEDURE(distribute_cb) :: distribute_state_pdaf   !< Routine to distribute a state vector
  PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf     !< Initialize dimension of full observation vector
  PROCEDURE(obs_op_cb) :: obs_op_pdaf                 !< Full observation operator
  PROCEDURE(prepost_cb) :: prepoststep_pdaf           !< User supplied pre/poststep routine
  PROCEDURE(next_obs_cb) :: next_observation_pdaf     !< Provide information on next forecast
  PROCEDURE(prodRinvA_cb) :: prodRinvA_pdaf           !< Provide product of inverse of R with matrix A

! *** OMI-provided procedures ***
  PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb         !< Initialize full observation vector
  PROCEDURE(init_obsvar_cb) :: PDAFomi_init_obsvar_cb     !< Initialize mean observation error variance
  PROCEDURE(init_obscovar_cb) :: PDAFomi_init_obscovar_cb !< Initialize mean observation error variance
  PROCEDURE(add_obs_error_cb) :: PDAFomi_add_obs_error_cb !< Add observation error covariance matrix


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAF3_assimilate_global_nondiagR -- START'

  IF (TRIM(filterstr) == 'SEIK') THEN
     CALL PDAF_assimilate_seik(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
          prodRinvA_pdaf, PDAFomi_init_obsvar_cb, next_observation_pdaf, outflag)
  ELSEIF (TRIM(filterstr) == 'ENKF') THEN
     WRITE (*,*) 'PDAF-ERROR: Use PDAF3_assimilate_enkf_nondiagR for EnKF'
     outflag=200
  ELSEIF (TRIM(filterstr) == 'ETKF') THEN
     CALL PDAF_assimilate_etkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
          prodRinvA_pdaf, PDAFomi_init_obsvar_cb, next_observation_pdaf, outflag)
  ELSEIF (TRIM(filterstr) == 'ESTKF') THEN
     CALL PDAF_assimilate_estkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
          prodRinvA_pdaf, PDAFomi_init_obsvar_cb, next_observation_pdaf, outflag)
  ELSEIF (TRIM(filterstr) == 'NETF') THEN
     WRITE (*,*) 'PDAF-ERROR: Use PDAF3_assimilate_nonlin_nondiagR for NETF and PF'
     outflag=200
  ELSEIF (TRIM(filterstr) == 'PF') THEN
     WRITE (*,*) 'PDAF-ERROR: Use PDAF3_assimilate_nonlin_nondiagR for NETF and PF'
     outflag=200
  ELSE
     WRITE (*,*) 'PDAF-ERROR: Invalid filter choice for PDAF3_assimilate_global_nondiagR'
     outflag=200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAF3_assimilate_global_nondiagR -- END'

END SUBROUTINE PDAF3_assimilate_global_nondiagR


!-------------------------------------------------------------------------------
!> Interface to transfer state to PDAF
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF3_assimilate_lnetf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdaf, likelihood_l_pdaf,  &
          prepoststep_pdaf, next_observation_pdaf, outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAF_cb_procedures
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFlocal, &
       ONLY: PDAFlocal_g2l_cb, &       !< Project global to local state vector
       PDAFlocal_l2g_cb                !< Project local to global state vecto
  USE PDAFassimilate_lnetf, ONLY: PDAF_assimilate_lnetf

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(inout) :: outflag    !< Status flag

! *** Argument procedures ***
  PROCEDURE(collect_cb) :: collect_state_pdaf         !< Routine to collect a state vector
  PROCEDURE(distribute_cb) :: distribute_state_pdaf   !< Routine to distribute a state vector
  PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf     !< Initialize dimension of full observation vector
  PROCEDURE(obs_op_cb) :: obs_op_pdaf                 !< Full observation operator
  PROCEDURE(init_n_domains_cb) :: init_n_domains_pdaf !< Provide number of local analysis domains
  PROCEDURE(init_dim_l_cb) :: init_dim_l_pdaf         !< Init state dimension for local ana. domain
  PROCEDURE(init_dim_obs_l_cb) :: init_dim_obs_l_pdaf !< Initialize local dimimension of obs. vector
  PROCEDURE(prepost_cb) :: prepoststep_pdaf           !< User supplied pre/poststep routine
  PROCEDURE(next_obs_cb) :: next_observation_pdaf     !< Provide information on next forecast
  PROCEDURE(likelihood_l_cb) :: likelihood_l_pdaf     !< Compute likelihood and apply localization

! *** OMI-provided procedures ***
  PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb         !< Initialize full observation vector
  PROCEDURE(init_obs_l_cb) :: PDAFomi_init_obs_l_cb       !< Initialize local observation vector
  PROCEDURE(init_obsvar_cb) :: PDAFomi_init_obsvar_cb     !< Initialize mean observation error variance
  PROCEDURE(init_obsvar_l_cb) :: PDAFomi_init_obsvar_l_cb !< Initialize local mean observation error variance
  PROCEDURE(g2l_obs_cb) :: PDAFomi_g2l_obs_cb             !< Restrict full obs. vector to local analysis domain


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAF3_assimilate_lnetf_nondiagR -- START'

  IF (TRIM(filterstr) == 'LNETF') THEN
     CALL PDAF_assimilate_lnetf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
          prepoststep_pdaf, likelihood_l_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdaf, PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, &
          PDAFomi_g2l_obs_cb, next_observation_pdaf, outflag)
  ELSE
     WRITE (*,*) 'PDAF-ERROR: Invalid filter choice for PDAF3_assimilate_lnetf_nondiagR'
     outflag=200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAF3_assimilate_lnetf_nondiagR -- END'

END SUBROUTINE PDAF3_assimilate_lnetf_nondiagR


!-------------------------------------------------------------------------------
!> Interface to transfer state to PDAF
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF3_assimilate_lknetf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdaf, prodRinvA_l_pdaf, prodRinvA_hyb_l_pdaf, &
          likelihood_l_pdaf, likelihood_hyb_l_pdaf,  &
          prepoststep_pdaf, next_observation_pdaf, outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAF_cb_procedures
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFlocal, &
       ONLY: PDAFlocal_g2l_cb, &       !< Project global to local state vector
       PDAFlocal_l2g_cb                !< Project local to global state vecto
  USE PDAFassimilate_lknetf, ONLY: PDAF_assimilate_lknetf

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(inout) :: outflag    !< Status flag

! *** Argument procedures ***
  PROCEDURE(collect_cb) :: collect_state_pdaf         !< Routine to collect a state vector
  PROCEDURE(distribute_cb) :: distribute_state_pdaf   !< Routine to distribute a state vector
  PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf     !< Initialize dimension of full observation vector
  PROCEDURE(obs_op_cb) :: obs_op_pdaf                 !< Full observation operator
  PROCEDURE(init_n_domains_cb) :: init_n_domains_pdaf !< Provide number of local analysis domains
  PROCEDURE(init_dim_l_cb) :: init_dim_l_pdaf         !< Init state dimension for local ana. domain
  PROCEDURE(init_dim_obs_l_cb) :: init_dim_obs_l_pdaf !< Initialize local dimimension of obs. vector
  PROCEDURE(prepost_cb) :: prepoststep_pdaf           !< User supplied pre/poststep routine
  PROCEDURE(next_obs_cb) :: next_observation_pdaf     !< Provide information on next forecast
  PROCEDURE(prodRinvA_l_cb) :: prodRinvA_l_pdaf       !< Provide product of inverse of R with matrix A
  PROCEDURE(likelihood_l_cb) :: likelihood_l_pdaf     !< Compute likelihood and apply localization
  PROCEDURE(prodRinvA_hyb_l_cb) :: prodRinvA_hyb_l_pdaf   !< Product R^-1 A on local analysis domain with hybrid weight
  PROCEDURE(likelihood_hyb_l_cb) :: likelihood_hyb_l_pdaf !< Compute likelihood and apply localization with tempering

! *** OMI-provided procedures ***
  PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb         !< Initialize full observation vector
  PROCEDURE(init_obs_l_cb) :: PDAFomi_init_obs_l_cb       !< Initialize local observation vector
  PROCEDURE(init_obsvar_cb) :: PDAFomi_init_obsvar_cb     !< Initialize mean observation error variance
  PROCEDURE(init_obsvar_l_cb) :: PDAFomi_init_obsvar_l_cb !< Initialize local mean observation error variance
  PROCEDURE(g2l_obs_cb) :: PDAFomi_g2l_obs_cb             !< Restrict full obs. vector to local analysis domain


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAF3_assimilate_lknetf_nondiagR -- START'

  IF (TRIM(filterstr) == 'LKNETF') THEN
     CALL PDAF_assimilate_lknetf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
          prodRinvA_l_pdaf, prodRinvA_hyb_l_pdaf, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          PDAFlocal_g2l_cb, PDAFlocal_l2g_cb,  PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
          PDAFomi_init_obsvar_l_cb, likelihood_l_pdaf, likelihood_hyb_l_pdaf, &
          next_observation_pdaf, outflag)
  ELSE
     WRITE (*,*) 'PDAF-ERROR: Invalid filter choice for PDAF3_assimilate_lknetf_nondiagR'
     outflag=200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAF3_assimilate_lknetf_nondiagR -- END'

END SUBROUTINE PDAF3_assimilate_lknetf_nondiagR


!-------------------------------------------------------------------------------
!> Interface to PDAF for global filters
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF3_assimilate_enkf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
     init_dim_obs_pdaf, obs_op_pdaf, add_obs_error_pdaf, init_obscovar_pdaf, &
     prepoststep_pdaf, next_observation_pdaf, outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAF_cb_procedures
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFassimilate_enkf, ONLY: PDAF_assimilate_enkf
  USE PDAFassimilate_lenkf, ONLY: PDAF_assimilate_lenkf

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag      !< Status flag
  
! *** Argument procedures ***
  PROCEDURE(collect_cb) :: collect_state_pdaf         !< Routine to collect a state vector
  PROCEDURE(distribute_cb) :: distribute_state_pdaf   !< Routine to distribute a state vector
  PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf     !< Initialize dimension of full observation vector
  PROCEDURE(obs_op_cb) :: obs_op_pdaf                 !< Full observation operator
  PROCEDURE(prepost_cb) :: prepoststep_pdaf           !< User supplied pre/poststep routine
  PROCEDURE(next_obs_cb) :: next_observation_pdaf     !< Provide information on next forecast
  PROCEDURE(prodRinvA_l_cb) :: prodRinvA_l_pdaf       !< Provide product of inverse of R with matrix A
  PROCEDURE(init_obscovar_cb) :: init_obscovar_pdaf   !< Initialize mean observation error variance
  PROCEDURE(add_obs_error_cb) :: add_obs_error_pdaf   !< Add observation error covariance matrix

! *** OMI-provided procedures ***
  PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb     !< Initialize full observation vector
  PROCEDURE(localize_cb) :: PDAFomi_localize_covar_cb !< Apply localization to HP and HPH^T


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAF3_assimilate_enkf_nondiagR -- START'

  IF (TRIM(filterstr) == 'ENKF') THEN
     CALL PDAF_assimilate_enkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
          add_obs_error_pdaf, init_obscovar_pdaf, &
          next_observation_pdaf, outflag)
  ELSEIF (TRIM(filterstr) == 'LENKF') THEN
     CALL PDAF_assimilate_lenkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
          PDAFomi_localize_covar_cb, add_obs_error_pdaf, init_obscovar_pdaf, &
          next_observation_pdaf, outflag)
  ELSE
     WRITE (*,*) 'PDAF-ERROR: Invalid filter choice for PDAF3_assimilate_enkf_nondiagR'
     outflag=200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAF3_assimilate_enkf_nondiagR -- END'

END SUBROUTINE PDAF3_assimilate_enkf_nondiagR


!-------------------------------------------------------------------------------
!> Interface to transfer state to PDAF
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF3_assimilate_lenkf_nondiagR(collect_state_pdaf, &
     init_dim_obs_pdaf, obs_op_pdaf, prepoststep_pdaf, localize_pdaf, &
     add_obs_error_pdaf, init_obscovar_pdaf, outflag)

  USE PDAF_mod_core, ONLY: debug
  USE PDAF_cb_procedures
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFassimilate_lenkf, ONLY: PDAF_assimilate_lenkf

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag  !< Status flag
  
! *** Argument procedures ***
  PROCEDURE(collect_cb) :: collect_state_pdaf         !< Routine to collect a state vector
  PROCEDURE(distribute_cb) :: distribute_state_pdaf   !< Routine to distribute a state vector
  PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf     !< Initialize dimension of full observation vector
  PROCEDURE(obs_op_cb) :: obs_op_pdaf                 !< Full observation operator
  PROCEDURE(prepost_cb) :: prepoststep_pdaf           !< User supplied pre/poststep routine
  PROCEDURE(next_obs_cb) :: next_observation_pdaf     !< Provide information on next forecast
  PROCEDURE(prodRinvA_l_cb) :: prodRinvA_l_pdaf       !< Provide product of inverse of R with matrix A
  PROCEDURE(init_obscovar_cb) :: init_obscovar_pdaf   !< Initialize mean observation error variance
  PROCEDURE(add_obs_error_cb) :: add_obs_error_pdaf   !< Add observation error covariance matrix
  PROCEDURE(localize_cb) :: localize_pdaf             !< Apply covariance localization

! *** OMI-provided procedures ***
  PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb     !< Initialize full observation vector


! **************************************************
! *** Call the full assimilate interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAF3_assimilate_lenkf_nondiagR -- START'

  CALL PDAF_assimilate_lenkf(collect_state_pdaf, distribute_state_pdaf, &
       init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
       localize_pdaf, add_obs_error_pdaf, init_obscovar_pdaf, &
       next_observation_pdaf, outflag)


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAF3_assimilate_lenkf_nondiagR -- END'

END SUBROUTINE PDAF3_assimilate_lenkf_nondiagR


!-------------------------------------------------------------------------------
!> Interface to PDAF for global filters
!!
!! __Revision history:__
!! 2024-08 - Lars Nerger - Initial code
!! Other revisions - see repository log
!!
SUBROUTINE PDAF3_assimilate_nonlin_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
     init_dim_obs_pdaf, obs_op_pdaf, likelihood_pdaf, &
     prepoststep_pdaf, next_observation_pdaf, outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAF_cb_procedures
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFassimilate_netf, ONLY: PDAF_assimilate_netf
  USE PDAFassimilate_pf, ONLY: PDAF_assimilate_pf

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag      !< Status flag
  
! *** Argument procedures ***
  PROCEDURE(collect_cb) :: collect_state_pdaf         !< Routine to collect a state vector
  PROCEDURE(distribute_cb) :: distribute_state_pdaf   !< Routine to distribute a state vector
  PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf     !< Initialize dimension of full observation vector
  PROCEDURE(obs_op_cb) :: obs_op_pdaf                 !< Full observation operator
  PROCEDURE(prepost_cb) :: prepoststep_pdaf           !< User supplied pre/poststep routine
  PROCEDURE(next_obs_cb) :: next_observation_pdaf     !< Provide information on next forecast
  PROCEDURE(likelihood_l_cb) :: likelihood_pdaf       !< Compute likelihood

! *** OMI-provided procedures ***
  PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb         !< Initialize full observation vector


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAF3_assimilate_nonlin_nondiagR -- START'

  IF (TRIM(filterstr) == 'NETF') THEN
     CALL PDAF_assimilate_netf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
          likelihood_pdaf, next_observation_pdaf, outflag)
  ELSEIF (TRIM(filterstr) == 'PF') THEN
     CALL PDAF_assimilate_pf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
          likelihood_pdaf, next_observation_pdaf, outflag)
  ELSE
     WRITE (*,*) 'PDAF-ERROR: Invalid filter choice for PDAF3_assimilate_nonlin_nondiagR'
     outflag=200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAF3_assimilate_nonlin_nondiagR -- END'

END SUBROUTINE PDAF3_assimilate_nonlin_nondiagR

END MODULE PDAF3_assimilate_ens_nondiagR
