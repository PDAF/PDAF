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
!> Interfaces to PDAF for flexible parallelization mode for ensemble filters using OMI
!!
!! The interface routines provide the advanced compact
!! interfaces for using PDAF-OMI. The routines
!! just call of one the PDAF_put_state interface routines
!! with the full interface. In the call the specific PDAF
!! internal subroutines for PDAF-OMI are specified.
!!
!! The interfaces in this module are for handling nondiagonal R matrices
!!
!! !  This is a core file of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code by collecting files into a module
!! * Other revisions - see repository log
!!
MODULE PDAFomi_put_state_ens_nondiagR

CONTAINS

!-------------------------------------------------------------------------------
!> Interface to transfer state to PDAF
!!
!! __Revision history:__
!! * 2024-07 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAFomi_put_state_local_nondiagR(collect_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdafomi, prodRinvA_l_pdafomi, &
          g2l_state_pdaf, l2g_state_pdaf, outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAFomi_obs_l, ONLY: PDAFomi_dealloc
  USE PDAFput_state_lseik, ONLY: PDAF_put_state_lseik
  USE PDAFput_state_letkf, ONLY: PDAF_put_state_letkf
  USE PDAFput_state_lestkf, ONLY: PDAF_put_state_lestkf

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag  !< Status flag
  
! *** External subroutines ***
  EXTERNAL :: collect_state_pdaf, &    !< Routine to collect a state vector
       prepoststep_pdaf                !< User supplied pre/poststep routine
  EXTERNAL :: init_n_domains_pdaf, &   !< Provide number of local analysis domains
       init_dim_l_pdaf, &              !< Init state dimension for local ana. domain
       g2l_state_pdaf, &               !< Get state on local ana. domain from full state
       l2g_state_pdaf                  !< Init full state from local state
  EXTERNAL :: init_dim_obs_pdafomi, &  !< Initialize dimension of full observation vector
       obs_op_pdafomi, &               !< Full observation operator
       init_dim_obs_l_pdafomi, &       !< Initialize local dimimension of obs. vector
       prodRinvA_l_pdafomi             !< Provide product of inverse of R with matrix A
  EXTERNAL :: PDAFomi_init_obs_f_cb, & !< Initialize full observation vector
       PDAFomi_init_obs_l_cb, &        !< Initialize local observation vector
       PDAFomi_init_obsvar_cb, &       !< Initialize mean observation error variance
       PDAFomi_init_obsvar_l_cb, &     !< Initialize local mean observation error variance
       PDAFomi_g2l_obs_cb              !< Restrict full obs. vector to local analysis domain


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_put_state_local_nondiagR -- START'

  IF (TRIM(filterstr) == 'LSEIK') THEN
     CALL PDAF_put_state_lseik(collect_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
          prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
          g2l_state_pdaf, l2g_state_pdaf, PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
          PDAFomi_init_obsvar_l_cb, outflag)
  ELSE IF (TRIM(filterstr) == 'LETKF') THEN
     CALL PDAF_put_state_letkf(collect_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
          prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
          g2l_state_pdaf, l2g_state_pdaf, PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
          PDAFomi_init_obsvar_l_cb, outflag)
  ELSE IF (TRIM(filterstr) == 'LESTKF') THEN
     CALL PDAF_put_state_lestkf(collect_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
          prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
          g2l_state_pdaf, l2g_state_pdaf, PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
          PDAFomi_init_obsvar_l_cb, outflag)
  ELSE IF (TRIM(filterstr) == 'LNETF') THEN
     WRITE (*,*) 'PDAF-ERROR: Use PDAFomi_put_state_lnetf_nondiagR for LNETF'
     outflag=200
  ELSE IF (TRIM(filterstr) == 'LKNETF') THEN
     WRITE (*,*) 'PDAF-ERROR: Use PDAFomi_put_state_lknetf_nondiagR for LKNETF'
     outflag=200
  ELSE
     WRITE (*,*) 'PDAF-ERROR: Invalid filter choice for PDAFomi_put_state_local_nondiagR'
     outflag=200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_put_state_local_nondiagR -- END'

END SUBROUTINE PDAFomi_put_state_local_nondiagR


!-------------------------------------------------------------------------------
!> Interface to PDAF for global filters
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAFomi_put_state_global_nondiagR(collect_state_pdaf, &
     init_dim_obs_pdaf, obs_op_pdaf, prodRinvA_pdaf, prepoststep_pdaf, outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAFomi_obs_l, ONLY: PDAFomi_dealloc
  USE PDAFput_state_seik, ONLY: PDAF_put_state_seik
  USE PDAFput_state_enkf, ONLY: PDAF_put_state_enkf
  USE PDAFput_state_estkf, ONLY: PDAF_put_state_estkf

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag  !< Status flag
  
! *** External subroutines ***
  EXTERNAL :: collect_state_pdaf, &    !< Routine to collect a state vector
       prepoststep_pdaf                !< User supplied pre/poststep routine
  EXTERNAL :: init_dim_obs_pdaf, &     !< Initialize dimension of observation vector
       obs_op_pdaf                     !< Observation operator
  EXTERNAL :: prodRinvA_pdaf           !< Provide product R^-1 A
  EXTERNAL :: PDAFomi_init_obs_f_cb, & !< Initialize observation vector
       PDAFomi_init_obsvar_cb, &       !< Initialize mean observation error variance
       PDAFomi_init_obscovar_cb, &     !< Initialize mean observation error variance
       PDAFomi_add_obs_error_cb        !< Add observation error covariance matrix
       

! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_put_state_global_nondiagR -- START'

  IF (TRIM(filterstr) == 'SEIK') THEN
     CALL PDAF_put_state_seik(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
          PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
          prodRinvA_pdaf, PDAFomi_init_obsvar_cb, outflag)
  ELSEIF (TRIM(filterstr) == 'ENKF') THEN
     WRITE (*,*) 'PDAF-ERROR: Use PDAFomi_put_state_enkf_nondiagR for EnKF'
     outflag=200
  ELSEIF (TRIM(filterstr) == 'ETKF') THEN
     CALL PDAF_put_state_etkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
          PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
          prodRinvA_pdaf, PDAFomi_init_obsvar_cb, outflag)
  ELSEIF (TRIM(filterstr) == 'ESTKF') THEN
     CALL PDAF_put_state_estkf(collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
          PDAFomi_init_obs_f_cb, prepoststep_pdaf, &
          prodRinvA_pdaf, PDAFomi_init_obsvar_cb, outflag)
  ELSEIF (TRIM(filterstr) == 'NETF') THEN
     WRITE (*,*) 'PDAF-ERROR: Use PDAFomi_put_state_nonlin_nondiagR for NETF and PF'
     outflag=200
  ELSEIF (TRIM(filterstr) == 'PF') THEN
     WRITE (*,*) 'PDAF-ERROR: Use PDAFomi_put_state_nonlin_nondiagR for NETF and PF'
     outflag=200
  ELSE
     WRITE (*,*) 'PDAF-ERROR: Invalid filter choice for PDAFomi_put_state_global_nondiagR'
     outflag=200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_put_state_global_nondiagR -- END'

END SUBROUTINE PDAFomi_put_state_global_nondiagR


!-------------------------------------------------------------------------------
!> Interface to PDAF for global filters
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAFomi_put_state_enkf_nondiagR(collect_state_pdaf, &
     init_dim_obs_pdafomi, obs_op_pdafomi, add_obs_error_pdafomi, init_obscovar_pdafomi, &
     prepoststep_pdaf, outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAFomi_obs_l, ONLY: PDAFomi_dealloc
  USE PDAFput_state_enkf, ONLY: PDAF_put_state_enkf

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag  !< Status flag
  
! *** External subroutines ***
  EXTERNAL :: collect_state_pdaf, &    !< Routine to collect a state vector
       prepoststep_pdaf                !< User supplied pre/poststep routine
  EXTERNAL :: init_dim_obs_pdafomi, &  !< Initialize dimension of observation vector
       obs_op_pdafomi, &               !< Observation operator
       init_obscovar_pdafomi, &        !< Initialize mean observation error variance
       add_obs_error_pdafomi           !< Add observation error covariance matrix
  EXTERNAL :: PDAFomi_init_obs_f_cb, & !< Initialize observation vector
       PDAFomi_prodRinvA_cb, &         !< Provide product R^-1 A
       PDAFomi_likelihood_cb           !< Compute likelihood


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_put_state_enkf_nondiagR -- START'

  IF (TRIM(filterstr) == 'ENKF') THEN
     CALL PDAF_put_state_enkf(collect_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
          PDAFomi_init_obs_f_cb, prepoststep_pdaf, add_obs_error_pdafomi, &
          init_obscovar_pdafomi, outflag)
  ELSE
     WRITE (*,*) 'PDAF-ERROR: Invalid filter choice for PDAFomi_put_state_enkf_nondiagR'
     outflag=200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_put_state_enkf_nondiagR -- END'

END SUBROUTINE PDAFomi_put_state_enkf_nondiagR


!-------------------------------------------------------------------------------
!> Interface to transfer state to PDAF
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAFomi_put_state_lenkf_nondiagR(collect_state_pdaf, &
     init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, localize_covar_pdafomi, &
     add_obs_error_pdafomi, init_obscovar_pdafomi, outflag)

  USE PDAF_mod_core, ONLY: debug
  USE PDAFomi_obs_l, ONLY: PDAFomi_dealloc
  USE PDAFput_state_lenkf, ONLY: PDAF_put_state_lenkf

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag  !< Status flag
  
! *** External subroutines ***
  EXTERNAL :: collect_state_pdaf, &    !< Routine to collect a state vector
       prepoststep_pdaf                !< User supplied pre/poststep routine
  EXTERNAL :: init_dim_obs_pdafomi, &  !< Initialize dimension of observation vector
       obs_op_pdafomi, &               !< Observation operator
       localize_covar_pdafomi, &       !< Apply localization to HP and HPH^T
       init_obscovar_pdafomi, &        !< Initialize mean observation error variance
       add_obs_error_pdafomi           !< Provide product R^-1 A
  EXTERNAL :: PDAFomi_init_obs_f_cb    !< Initialize observation vector


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_put_state_lenkf_nondiagR -- START'

  CALL PDAF_put_state_lenkf(collect_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
       PDAFomi_init_obs_f_cb, prepoststep_pdaf, localize_covar_pdafomi, &
       add_obs_error_pdafomi, init_obscovar_pdafomi, outflag)


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_put_state_lenkf_nondiagR -- END'

END SUBROUTINE PDAFomi_put_state_lenkf_nondiagR


!-------------------------------------------------------------------------------
!> Interface to PDAF for global filters
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAFomi_put_state_nonlin_nondiagR(collect_state_pdaf, &
     init_dim_obs_pdafomi, obs_op_pdafomi, likelihood_pdafomi, prepoststep_pdaf, &
     outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAFomi_obs_l, ONLY: PDAFomi_dealloc
  USE PDAFput_state_netf, ONLY: PDAF_put_state_netf
  USE PDAFput_state_pf, ONLY: PDAF_put_state_pf

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(out) :: outflag  !< Status flag
  
! *** External subroutines ***
  EXTERNAL :: collect_state_pdaf, &    !< Routine to collect a state vector
       prepoststep_pdaf                !< User supplied pre/poststep routine
  EXTERNAL :: init_dim_obs_pdafomi, &  !< Initialize dimension of observation vector
       obs_op_pdafomi, &               !< Observation operator
       likelihood_pdafomi              !< Compute likelihood
  EXTERNAL :: PDAFomi_init_obs_f_cb, & !< Initialize observation vector
       PDAFomi_init_obsvar_cb, &       !< Initialize mean observation error variance
       PDAFomi_init_obscovar_cb, &     !< Initialize mean observation error variance
       PDAFomi_add_obs_error_cb        !< Add observation error covariance matrix


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_put_state_nonlin_nondiagR -- START'

  IF (TRIM(filterstr) == 'NETF') THEN
     CALL PDAF_put_state_netf(collect_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
          PDAFomi_init_obs_f_cb, prepoststep_pdaf, likelihood_pdafomi, outflag)
  ELSEIF (TRIM(filterstr) == 'PF') THEN
     CALL PDAF_put_state_pf(collect_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
          PDAFomi_init_obs_f_cb, prepoststep_pdaf, likelihood_pdafomi, outflag)
  ELSE
     WRITE (*,*) 'PDAF-ERROR: Invalid filter choice for PDAFomi_put_state_nonlin_nondiagR'
     outflag=200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_put_state_nonlin_nondiagR -- END'

END SUBROUTINE PDAFomi_put_state_nonlin_nondiagR


!-------------------------------------------------------------------------------
!> Interface to transfer state to PDAF
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAFomi_put_state_lnetf_nondiagR(collect_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdafomi, likelihood_l_pdafomi, &
          g2l_state_pdaf, l2g_state_pdaf, outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAFomi_obs_l, ONLY: PDAFomi_dealloc
  USE PDAFput_state_lnetf, ONLY: PDAF_put_state_lnetf

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag  !< Status flag
  
! *** External subroutines ***
  EXTERNAL :: collect_state_pdaf, &    !< Routine to collect a state vector
       prepoststep_pdaf                !< User supplied pre/poststep routine
  EXTERNAL :: init_n_domains_pdaf, &   !< Provide number of local analysis domains
       init_dim_l_pdaf, &              !< Init state dimension for local ana. domain
       g2l_state_pdaf, &               !< Get state on local ana. domain from full state
       l2g_state_pdaf                  !< Init full state from local state
  EXTERNAL :: init_dim_obs_pdafomi, &  !< Initialize dimension of full observation vector
       obs_op_pdafomi, &               !< Full observation operator
       init_dim_obs_l_pdafomi, &       !< Initialize local dimimension of obs. vector
       likelihood_l_pdafomi            !< Compute likelihood and apply localization
  EXTERNAL :: PDAFomi_init_obs_f_cb, & !< Initialize full observation vector
       PDAFomi_init_obs_l_cb, &        !< Initialize local observation vector
       PDAFomi_init_obsvar_cb, &       !< Initialize mean observation error variance
       PDAFomi_init_obsvar_l_cb, &     !< Initialize local mean observation error variance
       PDAFomi_g2l_obs_cb              !< Restrict full obs. vector to local analysis domain


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_put_state_lnetf_nondiagR -- START'

  IF (TRIM(filterstr) == 'LNETF') THEN
     CALL PDAF_put_state_lnetf(collect_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, likelihood_l_pdafomi, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
          g2l_state_pdaf, l2g_state_pdaf, PDAFomi_g2l_obs_cb, outflag)
  ELSE
     WRITE (*,*) 'PDAF-ERROR: Invalid filter choice for PDAFomi_put_state_lnetf_nondiagR'
     outflag=200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_put_state_lnetf_nondiagR -- END'

END SUBROUTINE PDAFomi_put_state_lnetf_nondiagR


!-------------------------------------------------------------------------------
!> Interface to transfer state to PDAF
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAFomi_put_state_lknetf_nondiagR(collect_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdafomi, prodRinvA_l_pdafomi, prodRinvA_hyb_l_pdafomi, &
          likelihood_l_pdafomi, likelihood_hyb_l_pdafomi, g2l_state_pdaf, l2g_state_pdaf, &
          outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAFomi_obs_l, ONLY: PDAFomi_dealloc
  USE PDAFput_state_lknetf, ONLY: PDAF_put_state_lknetf

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag  !< Status flag
  
! *** External subroutines ***
  EXTERNAL :: collect_state_pdaf, &    !< Routine to collect a state vector
       prepoststep_pdaf                !< User supplied pre/poststep routine
  EXTERNAL :: init_n_domains_pdaf, &   !< Provide number of local analysis domains
       init_dim_l_pdaf, &              !< Init state dimension for local ana. domain
       g2l_state_pdaf, &               !< Get state on local ana. domain from full state
       l2g_state_pdaf                  !< Init full state from local state
  EXTERNAL :: init_dim_obs_pdafomi, &  !< Initialize dimension of full observation vector
       obs_op_pdafomi, &               !< Full observation operator
       init_dim_obs_l_pdafomi, &       !< Initialize local dimimension of obs. vector
       prodRinvA_l_pdafomi, &          !< Provide product R^-1 A on local analysis domain
       likelihood_l_pdafomi, &         !< Compute likelihood and apply localization
       prodRinvA_hyb_l_pdafomi, &      !< Product R^-1 A on local analysis domain with hybrid weight
       likelihood_hyb_l_pdafomi        !< Compute likelihood and apply localization with tempering
  EXTERNAL :: PDAFomi_init_obs_f_cb, & !< Initialize full observation vector
       PDAFomi_init_obs_l_cb, &        !< Initialize local observation vector
       PDAFomi_init_obsvar_cb, &       !< Initialize mean observation error variance
       PDAFomi_init_obsvar_l_cb, &     !< Initialize local mean observation error variance
       PDAFomi_g2l_obs_cb              !< Restrict full obs. vector to local analysis domain


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_put_state_lknetf_nondiagR -- START'

  IF (TRIM(filterstr) == 'LKNETF') THEN
     CALL PDAF_put_state_lknetf(collect_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
          prodRinvA_l_pdafomi, prodRinvA_hyb_l_pdafomi, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
          g2l_state_pdaf, l2g_state_pdaf, PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
          PDAFomi_init_obsvar_l_cb, likelihood_l_pdafomi, likelihood_hyb_l_pdafomi, &
          outflag)
  ELSE
     WRITE (*,*) 'PDAF-ERROR: Invalid filter choice for PDAFomi_put_state_lknetf_nondiagR'
     outflag=200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_put_state_lknetf_nondiagR -- END'

END SUBROUTINE PDAFomi_put_state_lknetf_nondiagR

END MODULE PDAFomi_put_state_ens_nondiagR
