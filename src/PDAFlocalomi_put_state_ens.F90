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
!> Interfaces to PDAF for flexible parallelization mode for ensemble filters
!!
!! The interface routines provide the advanced compact
!! interfaces for using PDAF-OMI and PDAF-Local. The routines
!! just call of one the PDAF_put_state interface routines
!! with the full interface. In the call the specific PDAF
!! internal subroutines for PDAF-OMI and PDAF-Local are 
!! specified.
!!
!! The interface routines provided here are the PDAF-2
!! routines using the naming PDAFlocalomi.
!!
!! !  This is a core file of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code by collecting files into a module
!! * Other revisions - see repository log
!!
MODULE PDAFlocalomi_put_state_ens

CONTAINS

!-------------------------------------------------------------------------------
!> Interface to PDAF for domain-local filters with diagonal R
!!
!! __Revision history:__
!! * 2020-11 - Lars Nerger - Initial code
!! * 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
!! * Other revisions - see repository log
!!
SUBROUTINE PDAFlocalomi_put_state(collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
     prepoststep_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
      outflag)

  USE PDAF_mod_filter, ONLY: filterstr
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFlocal, &
       ONLY: PDAFlocal_g2l_cb, &       !< Project global to local state vector
       PDAFlocal_l2g_cb                !< Project local to global state vecto

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(inout) :: outflag ! Status flag
  
! *** Names of external subroutines ***
  EXTERNAL :: collect_state_pdaf, &          !< Routine to collect a state vector
       prepoststep_pdaf                      !< User supplied pre/poststep routine
  EXTERNAL :: init_n_domains_pdaf, &         !< Provide number of local analysis domains
       init_dim_l_pdaf, &                    !< Init state dimension for local ana. domain
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

  IF (TRIM(filterstr) == 'LSEIK') THEN
     CALL PDAF_put_state_lseik(collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
          PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
          PDAFomi_init_obsvar_l_cb, outflag)
  ELSE IF (TRIM(filterstr) == 'LETKF') THEN
     CALL PDAF_put_state_letkf(collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
          PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
          PDAFomi_init_obsvar_l_cb, outflag)
  ELSE IF (TRIM(filterstr) == 'LESTKF') THEN
     CALL PDAF_put_state_lestkf(collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
          PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
          PDAFomi_init_obsvar_l_cb, outflag)
  ELSE IF (TRIM(filterstr) == 'LNETF') THEN
     CALL PDAF_put_state_lnetf(collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, PDAFomi_likelihood_l_cb, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf,  &
          PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, PDAFomi_g2l_obs_cb, outflag)
  ELSE IF (TRIM(filterstr) == 'LKNETF') THEN
     CALL PDAF_put_state_lknetf(collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
          PDAFomi_prodRinvA_l_cb, PDAFomi_prodRinvA_hyb_l_cb, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
          PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
          PDAFomi_init_obsvar_l_cb, PDAFomi_likelihood_l_cb, PDAFomi_likelihood_hyb_l_cb, outflag)
  ELSE IF (TRIM(filterstr) == 'ENSRF') THEN
     CALL PDAF_put_state_ensrf(collect_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obsvars_f_cb, PDAFomi_localize_covar_serial_cb, &
          prepoststep_pdaf, outflag)
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

END SUBROUTINE PDAFlocalomi_put_state


!-------------------------------------------------------------------------------
!> Interface to transfer state to PDAF for local filters with nondiagonal R
!!
!! __Revision history:__
!! * 2024-07 - Lars Nerger - Initial code
!! * 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
!! * Other revisions - see repository log
!!
SUBROUTINE PDAFlocalomi_put_state_nondiagR(collect_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdafomi, prodRinvA_l_pdafomi, &
           outflag)

  USE PDAF_mod_filter, ONLY: filterstr, debug
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFlocal, &
       ONLY: PDAFlocal_g2l_cb, &       !< Project global to local state vector
       PDAFlocal_l2g_cb                !< Project local to global state vecto

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag  !< Status flag
  
! *** External subroutines ***
  EXTERNAL :: collect_state_pdaf, &    !< Routine to collect a state vector
       prepoststep_pdaf                !< User supplied pre/poststep routine
  EXTERNAL :: init_n_domains_pdaf, &   !< Provide number of local analysis domains
       init_dim_l_pdaf                 !< Init state dimension for local ana. domain
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
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFlocalomi_put_state_nondiagR -- START'

  IF (TRIM(filterstr) == 'LSEIK') THEN
     CALL PDAF_put_state_lseik(collect_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
          prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
          PDAFlocal_g2l_cb, PDAFlocal_l2g_cb,  PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
          PDAFomi_init_obsvar_l_cb, outflag)
  ELSE IF (TRIM(filterstr) == 'LETKF') THEN
     CALL PDAF_put_state_letkf(collect_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
          prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
          PDAFlocal_g2l_cb, PDAFlocal_l2g_cb,  PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
          PDAFomi_init_obsvar_l_cb, outflag)
  ELSE IF (TRIM(filterstr) == 'LESTKF') THEN
     CALL PDAF_put_state_lestkf(collect_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
          prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
          PDAFlocal_g2l_cb, PDAFlocal_l2g_cb,  PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
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


!-------------------------------------------------------------------------------
!> Interface to transfer state to PDAF for LNETF with nondiagonal R
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
!! * Other revisions - see repository log
!!
SUBROUTINE PDAFlocalomi_put_state_lnetf_nondiagR(collect_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdafomi, likelihood_l_pdafomi,  &
          outflag)

  USE PDAF_mod_filter, ONLY: filterstr, debug
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFlocal, &
       ONLY: PDAFlocal_g2l_cb, &       !< Project global to local state vector
       PDAFlocal_l2g_cb                !< Project local to global state vecto

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag  !< Status flag
  
! *** External subroutines ***
  EXTERNAL :: collect_state_pdaf, &    !< Routine to collect a state vector
       prepoststep_pdaf                !< User supplied pre/poststep routine
  EXTERNAL :: init_n_domains_pdaf, &   !< Provide number of local analysis domains
       init_dim_l_pdaf                 !< Init state dimension for local ana. domain
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
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFlocalomi_put_state_lnetf_nondiagR -- START'

  IF (TRIM(filterstr) == 'LNETF') THEN
     CALL PDAF_put_state_lnetf(collect_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, likelihood_l_pdafomi, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
          PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, PDAFomi_g2l_obs_cb, outflag)
  ELSE
     WRITE (*,*) 'PDAF-ERROR: Invalid filter choice for PDAFlocalomi_put_state_lnetf_nondiagR'
     outflag=200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFlocalomi_put_state_lnetf_nondiagR -- END'

END SUBROUTINE PDAFlocalomi_put_state_lnetf_nondiagR


!-------------------------------------------------------------------------------
!> Interface to transfer state to PDAF for LKNETF with nondiagonal R
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
!! * Other revisions - see repository log
!!
SUBROUTINE PDAFlocalomi_put_state_lknetf_nondiagR(collect_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdafomi, prodRinvA_l_pdafomi, prodRinvA_hyb_l_pdafomi, &
          likelihood_l_pdafomi, likelihood_hyb_l_pdafomi,  &
          outflag)

  USE PDAF_mod_filter, ONLY: filterstr, debug
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFlocal, &
       ONLY: PDAFlocal_g2l_cb, &       !< Project global to local state vector
       PDAFlocal_l2g_cb                !< Project local to global state vecto

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag  !< Status flag
  
! *** External subroutines ***
  EXTERNAL :: collect_state_pdaf, &    !< Routine to collect a state vector
       prepoststep_pdaf                !< User supplied pre/poststep routine
  EXTERNAL :: init_n_domains_pdaf, &   !< Provide number of local analysis domains
       init_dim_l_pdaf                 !< Init state dimension for local ana. domain
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
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFlocalomi_put_state_lknetf_nondiagR -- START'

  IF (TRIM(filterstr) == 'LKNETF') THEN
     CALL PDAF_put_state_lknetf(collect_state_pdaf, init_dim_obs_pdafomi, obs_op_pdafomi, &
          PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, prepoststep_pdaf, &
          prodRinvA_l_pdafomi, prodRinvA_hyb_l_pdafomi, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
          PDAFlocal_g2l_cb, PDAFlocal_l2g_cb,  PDAFomi_g2l_obs_cb, PDAFomi_init_obsvar_cb, &
          PDAFomi_init_obsvar_l_cb, likelihood_l_pdafomi, likelihood_hyb_l_pdafomi, &
          outflag)
  ELSE
     WRITE (*,*) 'PDAF-ERROR: Invalid filter choice for PDAFlocalomi_put_state_lknetf_nondiagR'
     outflag=200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFlocalomi_put_state_lknetf_nondiagR -- END'

END SUBROUTINE PDAFlocalomi_put_state_lknetf_nondiagR

END MODULE PDAFlocalomi_put_state_ens
