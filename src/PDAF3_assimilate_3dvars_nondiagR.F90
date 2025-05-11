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
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code by collecting files into a module
!! * Other revisions - see repository log
!!
MODULE PDAF3_assimilate_3dvars_nondiagR

CONTAINS

!-------------------------------------------------------------------------------
!> Interface to PDAF for 3D-Var
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF3_assimilate_3dvar_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
     init_dim_obs_pdaf, obs_op_pdaf, prodRinvA_pdaf, &
     cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
     prepoststep_pdaf, next_observation_pdaf, outflag)
  
  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAF_cb_procedures
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFassimilate_3dvar, ONLY: PDAF_assimilate_3dvar

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
  PROCEDURE(cvt_cb) :: cvt_pdaf                       !< Apply control vector transform matrix to control vector
  PROCEDURE(cvt_adj_cb) :: cvt_adj_pdaf               !< Apply adjoint control vector transform matrix
  PROCEDURE(obs_op_lin_cb) :: obs_op_lin_pdaf         !< Linearized observation operator
  PROCEDURE(obs_op_adj_cb) :: obs_op_adj_pdaf         !< Adjoint observation operator
  PROCEDURE(prodRinvA_cb) :: prodRinvA_pdaf           !< Provide product R^-1 A

! *** OMI-provided procedures ***
  PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb     !< Initialize full observation vector


! **************************************************
! *** Call the full assimilate interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_assimilate_3dvar_nondiagR -- START'

  IF (TRIM(filterstr) == '3DVAR') THEN
     CALL PDAF_assimilate_3dvar(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, prodRinvA_pdaf, &
          cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          prepoststep_pdaf, next_observation_pdaf, outflag)
  ELSE
     WRITE (*,*) 'PDAF-ERROR: No valid filter type for PDAFomi_assimilate_3dvar_nondiagR'
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_assimilate_3dvar_nondiagR -- END'

END SUBROUTINE PDAF3_assimilate_3dvar_nondiagR


!-------------------------------------------------------------------------------
!> Interface to PDAF for En3D-Var/ESTKF
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! Other revisions - see repository log
!!
SUBROUTINE PDAF3_assimilate_en3dvar_estkf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
                init_dim_obs_pdaf, obs_op_pdaf, prodRinvA_pdaf, &
                cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
                prepoststep_pdaf, next_observation_pdaf, outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAF_cb_procedures
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFassimilate_en3dvar_estkf, ONLY: PDAF_assimilate_en3dvar_estkf

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
  PROCEDURE(cvt_ens_cb) :: cvt_ens_pdaf               !< Apply control vector transform matrix to control vector
  PROCEDURE(cvt_adj_ens_cb) :: cvt_adj_ens_pdaf       !< Apply adjoint control vector transform matrix
  PROCEDURE(obs_op_lin_cb) :: obs_op_lin_pdaf         !< Linearized observation operator
  PROCEDURE(obs_op_adj_cb) :: obs_op_adj_pdaf         !< Adjoint observation operator
  PROCEDURE(prodRinvA_cb) :: prodRinvA_pdaf           !< Provide product R^-1 A

! *** OMI-provided procedures ***
  PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb     !< Initialize full observation vector
  PROCEDURE(init_obsvar_cb) :: PDAFomi_init_obsvar_cb !< Initialize mean observation error variance


! **************************************************
! *** Call the full assimilate interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_assimilate_en3dvar_estkf_nondiagR -- START'

  IF (TRIM(filterstr) == '3DVAR') THEN
     CALL PDAF_assimilate_en3dvar_estkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, prodRinvA_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          PDAFomi_init_obsvar_cb, prepoststep_pdaf, next_observation_pdaf, outflag)
  ELSE
     WRITE (*,*) 'PDAF-ERROR: No valid filter type for PDAFomi_assimilate_en3dvar_estkf_nondiagR'
     outflag = 200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_assimilate_en3dvar_estkf_nondiagR -- END'

END SUBROUTINE PDAF3_assimilate_en3dvar_estkf_nondiagR

!-------------------------------------------------------------------------------
!> Interface to PDAF for En3D-Var/LESTKF
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
!! Other revisions - see repository log
!!
SUBROUTINE PDAF3_assimilate_en3dvar_lestkf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
     init_dim_obs_pdaf, obs_op_pdaf, prodRinvA_pdaf, &
     cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
     prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
      prepoststep_pdaf, next_observation_pdaf, outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAF_cb_procedures
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFlocal, &
       ONLY: PDAFlocal_g2l_cb, &       !< Project global to local state vector
       PDAFlocal_l2g_cb                !< Project local to global state vecto
  USE PDAFassimilate_en3dvar_lestkf, ONLY: PDAF_assimilate_en3dvar_lestkf

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
  PROCEDURE(cvt_ens_cb) :: cvt_ens_pdaf               !< Apply control vector transform matrix to control vector
  PROCEDURE(cvt_adj_ens_cb) :: cvt_adj_ens_pdaf       !< Apply adjoint control vector transform matrix
  PROCEDURE(cvt_cb) :: cvt_pdaf                       !< Apply control vector transform matrix to control vector
  PROCEDURE(cvt_adj_cb) :: cvt_adj_pdaf               !< Apply adjoint control vector transform matrix
  PROCEDURE(obs_op_lin_cb) :: obs_op_lin_pdaf         !< Linearized observation operator
  PROCEDURE(obs_op_adj_cb) :: obs_op_adj_pdaf         !< Adjoint observation operator
  PROCEDURE(prodRinvA_cb) :: prodRinvA_pdaf           !< Provide product R^-1 A
  PROCEDURE(prodRinvA_l_cb) :: prodRinvA_l_pdaf       !< Provide product R^-1 A and apply localizations

! *** OMI-provided procedures ***
  PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb         !< Initialize full observation vector
  PROCEDURE(init_obs_l_cb) :: PDAFomi_init_obs_l_cb       !< Initialize local observation vector
  PROCEDURE(init_obsvar_cb) :: PDAFomi_init_obsvar_cb     !< Initialize mean observation error variance
  PROCEDURE(init_obsvar_l_cb) :: PDAFomi_init_obsvar_l_cb !< Initialize local mean observation error variance
  PROCEDURE(g2l_obs_cb) :: PDAFomi_g2l_obs_cb             !< Restrict full obs. vector to local analysis domain


! **************************************************
! *** Call the full assimilate interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFlocalomi_assimilate_en3dvar_lestkf_nondiagR -- START'

  IF (TRIM(filterstr) == '3DVAR') THEN
     CALL PDAF_assimilate_en3dvar_lestkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, prodRinvA_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
          prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
          init_dim_obs_l_pdaf, PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, PDAFomi_g2l_obs_cb, &
          PDAFomi_init_obsvar_cb, PDAFomi_init_obsvar_l_cb, &
          prepoststep_pdaf, next_observation_pdaf, outflag)
  ELSE
     WRITE (*,*) 'PDAF-ERROR: No valid filter type for PDAFlocalomi_assimilate_en3dvar_lestkf_nondiagR'
     outflag = 200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFlocalomi_assimilate_en3dvar_lestkf_nondiagR -- END'

END SUBROUTINE PDAF3_assimilate_en3dvar_lestkf_nondiagR


!-------------------------------------------------------------------------------
!> Interface to PDAF for Hyb3D-Var/ESTKF
!!
!! __Revision history:__
!! 2024-08 - Lars Nerger - Initial code
!! Other revisions - see repository log
!!
SUBROUTINE PDAF3_assimilate_hyb3dvar_estkf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
                init_dim_obs_pdaf, obs_op_pdaf, prodRinvA_pdaf, &
                cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
                obs_op_lin_pdaf, obs_op_adj_pdaf, &
                prepoststep_pdaf, next_observation_pdaf, outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAF_cb_procedures
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFassimilate_hyb3dvar_estkf, ONLY: PDAF_assimilate_hyb3dvar_estkf

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
  PROCEDURE(cvt_ens_cb) :: cvt_ens_pdaf               !< Apply control vector transform matrix to control vector
  PROCEDURE(cvt_adj_ens_cb) :: cvt_adj_ens_pdaf       !< Apply adjoint control vector transform matrix
  PROCEDURE(cvt_cb) :: cvt_pdaf                       !< Apply control vector transform matrix to control vector
  PROCEDURE(cvt_adj_cb) :: cvt_adj_pdaf               !< Apply adjoint control vector transform matrix
  PROCEDURE(obs_op_lin_cb) :: obs_op_lin_pdaf         !< Linearized observation operator
  PROCEDURE(obs_op_adj_cb) :: obs_op_adj_pdaf         !< Adjoint observation operator
  PROCEDURE(prodRinvA_cb) :: prodRinvA_pdaf           !< Provide product R^-1 A

! *** OMI-provided procedures ***
  PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb     !< Initialize full observation vector
  PROCEDURE(init_obsvar_cb) :: PDAFomi_init_obsvar_cb !< Initialize mean observation error variance


! **************************************************
! *** Call the full assimilate interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_assimilate_hyb3dvar_estkf_nondiagR -- START'

  IF (TRIM(filterstr) == '3DVAR') THEN
     CALL PDAF_assimilate_hyb3dvar_estkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, prodRinvA_pdaf, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
          obs_op_lin_pdaf, obs_op_adj_pdaf, &
          PDAFomi_init_obsvar_cb, prepoststep_pdaf, next_observation_pdaf, outflag)
  ELSE
     WRITE (*,*) 'PDAF-ERROR: No valid filter type for PDAFomi_assimilate_hyb3dvar_estkf_nondiagR'
     outflag = 200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_assimilate_hyb3dvar_estkf_nondiagR -- END'

END SUBROUTINE PDAF3_assimilate_hyb3dvar_estkf_nondiagR


!-------------------------------------------------------------------------------
!> Interface to PDAF for Hyb3D-Var/LESTKF
!!
!! __Revision history:__
!! 2024-08 - Lars Nerger - Initial code
!! 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
!! Other revisions - see repository log
!!
SUBROUTINE PDAF3_assimilate_hyb3dvar_lestkf_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
     init_dim_obs_pdaf, obs_op_pdaf, prodRinvA_pdaf, &
     cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
     prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
      prepoststep_pdaf, next_observation_pdaf, outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAF_cb_procedures
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFlocal, &
       ONLY: PDAFlocal_g2l_cb, &       !< Project global to local state vector
       PDAFlocal_l2g_cb                !< Project local to global state vecto
  USE PDAFassimilate_hyb3dvar_lestkf, ONLY: PDAF_assimilate_hyb3dvar_lestkf

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
  PROCEDURE(cvt_ens_cb) :: cvt_ens_pdaf               !< Apply control vector transform matrix to control vector
  PROCEDURE(cvt_adj_ens_cb) :: cvt_adj_ens_pdaf       !< Apply adjoint control vector transform matrix
  PROCEDURE(cvt_cb) :: cvt_pdaf                       !< Apply control vector transform matrix to control vector
  PROCEDURE(cvt_adj_cb) :: cvt_adj_pdaf               !< Apply adjoint control vector transform matrix
  PROCEDURE(obs_op_lin_cb) :: obs_op_lin_pdaf         !< Linearized observation operator
  PROCEDURE(obs_op_adj_cb) :: obs_op_adj_pdaf         !< Adjoint observation operator
  PROCEDURE(prodRinvA_cb) :: prodRinvA_pdaf           !< Provide product R^-1 A
  PROCEDURE(prodRinvA_l_cb) :: prodRinvA_l_pdaf       !< Provide product R^-1 A and apply localizations

! *** OMI-provided procedures ***
  PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb         !< Initialize full observation vector
  PROCEDURE(init_obs_l_cb) :: PDAFomi_init_obs_l_cb       !< Initialize local observation vector
  PROCEDURE(init_obsvar_cb) :: PDAFomi_init_obsvar_cb     !< Initialize mean observation error variance
  PROCEDURE(init_obsvar_l_cb) :: PDAFomi_init_obsvar_l_cb !< Initialize local mean observation error variance
  PROCEDURE(g2l_obs_cb) :: PDAFomi_g2l_obs_cb             !< Restrict full obs. vector to local analysis domain


! **************************************************
! *** Call the full assimilate interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFlocalomi_assimilate_hyb3dvar_lestkf_nondiagR -- START'

  IF (TRIM(filterstr) == '3DVAR') THEN
     CALL PDAF_assimilate_hyb3dvar_lestkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, &
          prodRinvA_pdaf, cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
          obs_op_lin_pdaf, obs_op_adj_pdaf, &
          init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
          prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf, &
          init_dim_obs_l_pdaf,  PDAFlocal_g2l_cb, PDAFlocal_g2l_cb, PDAFomi_g2l_obs_cb, &
          PDAFomi_init_obsvar_cb, PDAFomi_init_obsvar_l_cb, &
          prepoststep_pdaf, next_observation_pdaf, outflag)
  ELSE
     WRITE (*,*) 'PDAF-ERROR: No valid filter type for PDAFlocalomi_assimilate_hyb3dvar_lestkf_nondiagR'
     outflag = 200
  END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

  CALL PDAFomi_dealloc()

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFlocalomi_assimilate_hyb3dvar_lestkf_nondiagR -- END'

END SUBROUTINE PDAF3_assimilate_hyb3dvar_lestkf_nondiagR

END MODULE PDAF3_assimilate_3dvars_nondiagR
