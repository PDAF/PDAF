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
!! routine, prodRinvA_pdafomi or likelihood_pdafomi is includes
!! as an argument.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code by collecting files into a module
!! * Other revisions - see repository log
!!
MODULE PDAF_assimilate_3dvars_nondiagR

CONTAINS

!-------------------------------------------------------------------------------
!> Interface to PDAF for 3D-Var
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF3_assimilate_3dvar_nondiagR(collect_state_pdaf, distribute_state_pdaf, &
     init_dim_obs_pdafomi, obs_op_pdafomi, prodRinvA_pdafomi, &
     cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
     prepoststep_pdaf, next_observation_pdaf, outflag)
  
  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFassimilate_3dvar, ONLY: PDAF_assimilate_3dvar

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag      !< Status flag
  
! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: collect_state_pdaf, &    !< Routine to collect a state vector
       distribute_state_pdaf, &        !< Routine to distribute a state vector
       next_observation_pdaf, &        !< Provide time step, time and dimension of next observation
       prepoststep_pdaf                !< User supplied pre/poststep routine
  EXTERNAL :: cvt_pdaf, &              !< Apply control vector transform matrix to control vector
       cvt_adj_pdaf                    !< Apply adjoint control vector transform matrix
  EXTERNAL :: init_dim_obs_pdafomi, &  !< Initialize dimension of observation vector
       obs_op_pdafomi, &               !< Observation operator
       obs_op_lin_pdafomi, &           !< Linearized observation operator
       obs_op_adj_pdafomi, &           !< Adjoint observation operator
       prodRinvA_pdafomi               !< Provide product R^-1 A
  EXTERNAL :: PDAFomi_init_obs_f_cb    !< Initialize observation vector


! **************************************************
! *** Call the full assimilate interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_assimilate_3dvar_nondiagR -- START'

  IF (TRIM(filterstr) == '3DVAR') THEN
     CALL PDAF_assimilate_3dvar(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, PDAFomi_init_obs_f_cb, prodRinvA_pdafomi, &
          cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
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
                init_dim_obs_pdafomi, obs_op_pdafomi, prodRinvA_pdafomi, &
                cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
                prepoststep_pdaf, next_observation_pdaf, outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFassimilate_en3dvar_estkf, ONLY: PDAF_assimilate_en3dvar_estkf

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag      !< Status flag
  
! *** External subroutines ***
  EXTERNAL :: collect_state_pdaf, &    !< Routine to collect a state vector
       distribute_state_pdaf, &        !< Routine to distribute a state vector
       next_observation_pdaf, &        !< Provide time step, time and dimension of next observation
       prepoststep_pdaf                !< User supplied pre/poststep routine
  EXTERNAL :: cvt_ens_pdaf, &          !< Apply control vector transform matrix to control vector
       cvt_adj_ens_pdaf                !< Apply adjoint control vector transform matrix
  EXTERNAL :: init_dim_obs_pdafomi, &  !< Initialize dimension of observation vector
       prodRinvA_pdafomi, &            !< Provide product R^-1 A
       obs_op_pdafomi, &               !< Observation operator
       obs_op_lin_pdafomi, &           !< Linearized observation operator
       obs_op_adj_pdafomi              !< Adjoint observation operator
  EXTERNAL :: PDAFomi_init_obs_f_cb, & !< Initialize observation vector
       PDAFomi_init_obsvar_cb          !< Initialize mean observation error variance


! **************************************************
! *** Call the full assimilate interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_assimilate_en3dvar_estkf_nondiagR -- START'

  IF (TRIM(filterstr) == '3DVAR') THEN
     CALL PDAF_assimilate_en3dvar_estkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, PDAFomi_init_obs_f_cb, prodRinvA_pdafomi, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
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
     init_dim_obs_pdafomi, obs_op_pdafomi, prodRinvA_pdafomi, &
     cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
     prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
      prepoststep_pdaf, next_observation_pdaf, outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFlocal, &
       ONLY: PDAFlocal_g2l_cb, &       !< Project global to local state vector
       PDAFlocal_l2g_cb                !< Project local to global state vecto
  USE PDAFassimilate_en3dvar_lestkf, ONLY: PDAF_assimilate_en3dvar_lestkf

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(inout) :: outflag    !< Status flag

! *** Names of external subroutines ***
  EXTERNAL :: collect_state_pdaf, &    !< Routine to collect a state vector
       distribute_state_pdaf, &        !< Routine to distribute a state vector
       next_observation_pdaf, &        !< Provide time step, time and dimension of next observation
       prepoststep_pdaf                !< User supplied pre/poststep routine
  EXTERNAL :: cvt_ens_pdaf, &          !< Apply control vector transform matrix to control vector
       cvt_adj_ens_pdaf                !< Apply adjoint control vector transform matrix
  EXTERNAL :: init_n_domains_pdaf, &   !< Provide number of local analysis domains
       init_dim_l_pdaf                 !< Init state dimension for local ana. domain
  EXTERNAL :: init_dim_obs_pdafomi, &  !< Initialize dimension of full observation vector
       obs_op_pdafomi, &               !< Full observation operator
       obs_op_lin_pdafomi, &           !< Linearized observation operator
       obs_op_adj_pdafomi, &           !< Adjoint observation operator
       init_dim_obs_l_pdafomi, &       !< Initialize local dimimension of obs. vector
       prodRinvA_pdafomi, &            !< Provide product R^-1 A
       prodRinvA_l_pdafomi             !< Provide product R^-1 A with localization
  EXTERNAL :: PDAFomi_init_obs_f_cb, & !< Initialize observation vector
       PDAFomi_init_obs_l_cb, &        !< Initialize local observation vector
       PDAFomi_init_obsvar_cb, &       !< Initialize mean observation error variance
       PDAFomi_init_obsvar_l_cb, &     !< Initialize local mean observation error variance
       PDAFomi_g2l_obs_cb              !< Restrict full obs. vector to local analysis domain


! **************************************************
! *** Call the full assimilate interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFlocalomi_assimilate_en3dvar_lestkf_nondiagR -- START'

  IF (TRIM(filterstr) == '3DVAR') THEN
     CALL PDAF_assimilate_en3dvar_lestkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, PDAFomi_init_obs_f_cb, prodRinvA_pdafomi, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
          init_dim_obs_pdafomi, obs_op_pdafomi, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
          prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, &
          init_dim_obs_l_pdafomi, PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, PDAFomi_g2l_obs_cb, &
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
                init_dim_obs_pdafomi, obs_op_pdafomi, prodRinvA_pdafomi, &
                cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
                obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
                prepoststep_pdaf, next_observation_pdaf, outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFassimilate_hyb3dvar_estkf, ONLY: PDAF_assimilate_hyb3dvar_estkf

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag      !< Status flag
  
! *** External subroutines ***
  EXTERNAL :: collect_state_pdaf, &    !< Routine to collect a state vector
       distribute_state_pdaf, &        !< Routine to distribute a state vector
       next_observation_pdaf, &        !< Provide time step, time and dimension of next observation
       prepoststep_pdaf                !< User supplied pre/poststep routine
  EXTERNAL :: cvt_pdaf, &              !< Apply control vector transform matrix to control vector
       cvt_adj_pdaf, &                 !< Apply adjoint control vector transform matrix
       cvt_ens_pdaf, &                 !< Apply ensemble control vector transform matrix to control vector
       cvt_adj_ens_pdaf                !< Apply adjoint ensemble control vector transform matrix
  EXTERNAL :: init_dim_obs_pdafomi, &  !< Initialize dimension of observation vector
       obs_op_pdafomi, &               !< Observation operator
       obs_op_lin_pdafomi, &           !< Linearized observation operator
       obs_op_adj_pdafomi, &           !< Adjoint observation operator
       prodRinvA_pdafomi               !< Provide product R^-1 A
  EXTERNAL :: PDAFomi_init_obs_f_cb, & !< Initialize observation vector
       PDAFomi_init_obsvar_cb          !< Initialize mean observation error variance


! **************************************************
! *** Call the full assimilate interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFomi_assimilate_hyb3dvar_estkf_nondiagR -- START'

  IF (TRIM(filterstr) == '3DVAR') THEN
     CALL PDAF_assimilate_hyb3dvar_estkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, PDAFomi_init_obs_f_cb, prodRinvA_pdafomi, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
          obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
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
     init_dim_obs_pdafomi, obs_op_pdafomi, prodRinvA_pdafomi, &
     cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
     prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
      prepoststep_pdaf, next_observation_pdaf, outflag)

  USE PDAF_mod_core, ONLY: filterstr, debug
  USE PDAFomi, ONLY: PDAFomi_dealloc
  USE PDAFlocal, &
       ONLY: PDAFlocal_g2l_cb, &       !< Project global to local state vector
       PDAFlocal_l2g_cb                !< Project local to global state vecto
  USE PDAFassimilate_hyb3dvar_lestkf, ONLY: PDAF_assimilate_hyb3dvar_lestkf

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(inout) :: outflag    !< Status flag

! *** Names of external subroutines ***
 EXTERNAL :: collect_state_pdaf, &     !< Routine to collect a state vector
       distribute_state_pdaf, &        !< Routine to distribute a state vector
       next_observation_pdaf, &        !< Provide time step, time and dimension of next observation
       prepoststep_pdaf                !< User supplied pre/poststep routine
  EXTERNAL :: cvt_ens_pdaf, &          !< Apply control vector transform matrix to control vector
       cvt_adj_ens_pdaf, &             !< Apply adjoint control vector transform matrix
       cvt_pdaf, &                     !< Apply control vector transform matrix to control vector
       cvt_adj_pdaf                    !< Apply adjoint control vector transform matrix
  EXTERNAL :: init_n_domains_pdaf, &   !< Provide number of local analysis domains
       init_dim_l_pdaf                 !< Init state dimension for local ana. domain
  EXTERNAL :: init_dim_obs_pdafomi, &  !< Initialize dimension of full observation vector
       obs_op_pdafomi, &               !< Full observation operator
       obs_op_lin_pdafomi, &           !< Linearized observation operator
       obs_op_adj_pdafomi, &           !< Adjoint observation operator
       init_dim_obs_l_pdafomi, &       !< Initialize local dimimension of obs. vector
       prodRinvA_pdafomi, &            !< Provide product R^-1 A
       prodRinvA_l_pdafomi             !< Provide product R^-1 A
  EXTERNAL :: PDAFomi_init_obs_f_cb, & !< Initialize observation vector
       PDAFomi_init_obs_l_cb, &        !< Initialize local observation vector
       PDAFomi_init_obsvar_cb, &       !< Initialize mean observation error variance
       PDAFomi_init_obsvar_l_cb, &     !< Initialize local mean observation error variance
       PDAFomi_prodRinvA_cb, &         !< Provide product R^-1 A
       PDAFomi_g2l_obs_cb              !< Restrict full obs. vector to local analysis domain


! **************************************************
! *** Call the full assimilate interface routine  ***
! **************************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAFomi-debug: ', debug, 'PDAFlocalomi_assimilate_hyb3dvar_lestkf_nondiagR -- START'

  IF (TRIM(filterstr) == '3DVAR') THEN
     CALL PDAF_assimilate_hyb3dvar_lestkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, PDAFomi_init_obs_f_cb, &
          prodRinvA_pdafomi, cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
          obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
          init_dim_obs_pdafomi, obs_op_pdafomi, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
          prodRinvA_l_pdafomi, init_n_domains_pdaf, init_dim_l_pdaf, &
          init_dim_obs_l_pdafomi,  PDAFlocal_g2l_cb, PDAFlocal_g2l_cb, PDAFomi_g2l_obs_cb, &
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

END MODULE PDAF_assimilate_3dvars_nondiagR
