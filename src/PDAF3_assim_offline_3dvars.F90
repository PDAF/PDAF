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
!> Interfaces to PDAF for offline mode for 3D-Vars
!!
!! The interface routines provide the advanced compact
!! interfaces for using PDAF-OMI and PDAF-Local. The routines
!! just call of one the PDAF_put_state interface routines
!! with the full interface. In the call the specific PDAF
!! internal subroutines for PDAF-OMI and PDAF-Local are 
!! specified.
!!
!! !  This is a core file of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-04 - Lars Nerger - Initial code by adapting PDAF3_put_state_ens_3dvars
!! * Other revisions - see repository log
!!
MODULE PDAF3_assim_offline_3dvars

CONTAINS

!-------------------------------------------------------------------------------
!> Universal interface to PDAF for all 3D-Var methods
!!
!! __Revision history:__
!! * 2025-04 - Lars Nerger - Initial code based on PDAF3_put_state code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF3_assim_offline_3dvar_all(init_dim_obs_pdaf, obs_op_pdaf, &
       cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
       init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
       prepoststep_pdaf, outflag)

    USE PDAF_mod_core, ONLY: filterstr, debug, subtype_filter
    USE PDAF_cb_procedures
    USE PDAFomi, ONLY: PDAFomi_dealloc
    USE PDAFlocal, &
         ONLY: PDAFlocal_g2l_cb, &       !< Project global to local state vector
         PDAFlocal_l2g_cb                !< Project local to global state vecto
    USE PDAFassimilate_hyb3dvar_lestkf, ONLY: PDAF_assim_offline_hyb3dvar_lestkf
    USE PDAFassimilate_hyb3dvar_estkf, ONLY: PDAF_assim_offline_hyb3dvar_estkf
    USE PDAFassimilate_en3dvar_lestkf, ONLY: PDAF_assim_offline_en3dvar_lestkf
    USE PDAFassimilate_en3dvar_estkf, ONLY: PDAF_assim_offline_en3dvar_estkf
    USE PDAFassimilate_3dvar, ONLY: PDAF_assim_offline_3dvar

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(out) :: outflag      !< Status flag

! *** Argument procedures ***
    PROCEDURE(prepost_cb) :: prepoststep_pdaf           !< User supplied pre/poststep routine
    PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf     !< Initialize dimension of full observation vector
    PROCEDURE(obs_op_cb) :: obs_op_pdaf                 !< Full observation operator
    PROCEDURE(init_n_domains_cb) :: init_n_domains_pdaf !< Provide number of local analysis domains
    PROCEDURE(init_dim_l_cb) :: init_dim_l_pdaf         !< Init state dimension for local ana. domain
    PROCEDURE(init_dim_obs_l_cb) :: init_dim_obs_l_pdaf !< Initialize local dimimension of obs. vector
    PROCEDURE(cvt_ens_cb) :: cvt_ens_pdaf               !< Apply control vector transform matrix to control vector
    PROCEDURE(cvt_adj_ens_cb) :: cvt_adj_ens_pdaf       !< Apply adjoint control vector transform matrix
    PROCEDURE(cvt_cb) :: cvt_pdaf                       !< Apply control vector transform matrix to control vector
    PROCEDURE(cvt_adj_cb) :: cvt_adj_pdaf               !< Apply adjoint control vector transform matrix
    PROCEDURE(obs_op_lin_cb) :: obs_op_lin_pdaf         !< Linearized observation operator
    PROCEDURE(obs_op_adj_cb) :: obs_op_adj_pdaf         !< Adjoint observation operator

! *** OMI-provided procedures ***
    PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb         !< Initialize full observation vector
    PROCEDURE(init_obs_l_cb) :: PDAFomi_init_obs_l_cb       !< Initialize local observation vector
    PROCEDURE(init_obsvar_cb) :: PDAFomi_init_obsvar_cb     !< Initialize mean observation error variance
    PROCEDURE(init_obsvar_l_cb) :: PDAFomi_init_obsvar_l_cb !< Initialize local mean observation error variance
    PROCEDURE(prodRinvA_cb) :: PDAFomi_prodRinvA_cb         !< Provide product R^-1 A
    PROCEDURE(g2l_obs_cb) :: PDAFomi_g2l_obs_cb             !< Restrict full obs. vector to local analysis domain
    PROCEDURE(prodRinvA_l_cb) :: PDAFomi_prodRinvA_l_cb     !< Provide product R^-1 A on local analysis domain


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

    IF (debug>0) &
         WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF3_assim_offline_3dvar_all -- START'

    IF (TRIM(filterstr) == '3DVAR' .AND. subtype_filter==0) THEN
       CALL PDAF_assim_offline_3dvar(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_prodRinvA_cb, &
            cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
            prepoststep_pdaf, outflag)
    ELSE IF (TRIM(filterstr) == '3DVAR' .AND. subtype_filter==1) THEN
       CALL PDAF_assim_offline_en3dvar_lestkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_prodRinvA_cb, &
            cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
            init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
            PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, init_dim_l_pdaf, &
            init_dim_obs_l_pdaf,  PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, PDAFomi_g2l_obs_cb, &
            PDAFomi_init_obsvar_cb, PDAFomi_init_obsvar_l_cb, prepoststep_pdaf, outflag)
    ELSEIF (TRIM(filterstr) == '3DVAR' .AND. subtype_filter==2) THEN
       CALL PDAF_assim_offline_en3dvar_estkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_prodRinvA_cb, &
            cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
            PDAFomi_init_obsvar_cb, prepoststep_pdaf, outflag)
    ELSE IF (TRIM(filterstr) == '3DVAR' .AND. subtype_filter==3) THEN
       CALL PDAF_assim_offline_hyb3dvar_lestkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_prodRinvA_cb, &
            cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
            init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
            PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, init_dim_l_pdaf, &
            init_dim_obs_l_pdaf,  PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, PDAFomi_g2l_obs_cb, &
            PDAFomi_init_obsvar_cb, PDAFomi_init_obsvar_l_cb, prepoststep_pdaf, outflag)
    ELSE IF (TRIM(filterstr) == '3DVAR' .AND. subtype_filter==4) THEN
       CALL PDAF_assim_offline_hyb3dvar_estkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_prodRinvA_cb, &
            cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
            obs_op_lin_pdaf, obs_op_adj_pdaf, PDAFomi_init_obsvar_cb, &
            prepoststep_pdaf, outflag)
    ELSE
       WRITE (*,*) 'PDAF-ERROR: No valid filter type for PDAF3_assim_offline_3dvar_all'
    END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

    CALL PDAFomi_dealloc()

    IF (debug>0) &
         WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF3_assimilate_3dvar_all -- END'

  END SUBROUTINE PDAF3_assim_offline_3dvar_all



!-------------------------------------------------------------------------------
!> Interface to PDAF for 3dvar
!!
!! __Revision history:__
!! * 2025-04 - Lars Nerger - Initial code based on PDAF3_put_state code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF3_assim_offline_3dvar(init_dim_obs_pdaf, obs_op_pdaf, &
       cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf, outflag)

    USE PDAF_mod_core, ONLY: filterstr, subtype_filter
    USE PDAF_cb_procedures
    USE PDAFomi, ONLY: PDAFomi_dealloc
    USE PDAFassimilate_3dvar, ONLY: PDAF_assim_offline_3dvar

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(out) :: outflag      !< Status flag

! *** Argument procedures ***
    PROCEDURE(prepost_cb) :: prepoststep_pdaf           !< User supplied pre/poststep routine
    PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf     !< Initialize dimension of full observation vector
    PROCEDURE(obs_op_cb) :: obs_op_pdaf                 !< Full observation operator
    PROCEDURE(cvt_cb) :: cvt_pdaf                       !< Apply control vector transform matrix to control vector
    PROCEDURE(cvt_adj_cb) :: cvt_adj_pdaf               !< Apply adjoint control vector transform matrix
    PROCEDURE(obs_op_lin_cb) :: obs_op_lin_pdaf         !< Linearized observation operator
    PROCEDURE(obs_op_adj_cb) :: obs_op_adj_pdaf         !< Adjoint observation operator

! *** OMI-provided procedures ***
    PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb         !< Initialize full observation vector
    PROCEDURE(prodRinvA_cb) :: PDAFomi_prodRinvA_cb         !< Provide product R^-1 A


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

    IF (TRIM(filterstr) == '3DVAR' .AND. subtype_filter==0) THEN
       CALL PDAF_assim_offline_3dvar(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_prodRinvA_cb, &
            cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
            prepoststep_pdaf, outflag)
    ELSE
       WRITE (*,*) 'PDAF-ERROR: No valid filter type for PDAF3_assim_offline_3dvar'
       outflag = 200
    END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

    CALL PDAFomi_dealloc()

  END SUBROUTINE PDAF3_assim_offline_3dvar

!-------------------------------------------------------------------------------
!> Universal interface to PDAF for En3D-Var
!!
!! __Revision history:__
!! * 2025-04 - Lars Nerger - Initial code based on PDAF3_put_state code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF3_assim_offline_en3dvar(init_dim_obs_pdaf, obs_op_pdaf, &
       cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
       init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
       prepoststep_pdaf, outflag)

    USE PDAF_mod_core, ONLY: filterstr, subtype_filter
    USE PDAF_cb_procedures
    USE PDAFomi, ONLY: PDAFomi_dealloc
    USE PDAFlocal, ONLY: PDAFlocal_g2l_cb, PDAFlocal_l2g_cb
    USE PDAFassimilate_en3dvar_lestkf, ONLY: PDAF_assim_offline_en3dvar_lestkf
    USE PDAFassimilate_en3dvar_estkf, ONLY: PDAF_assim_offline_en3dvar_estkf

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(out) :: outflag      !< Status flag

! *** Argument procedures ***
    PROCEDURE(prepost_cb) :: prepoststep_pdaf           !< User supplied pre/poststep routine
    PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf     !< Initialize dimension of full observation vector
    PROCEDURE(obs_op_cb) :: obs_op_pdaf                 !< Full observation operator
    PROCEDURE(init_n_domains_cb) :: init_n_domains_pdaf !< Provide number of local analysis domains
    PROCEDURE(init_dim_l_cb) :: init_dim_l_pdaf         !< Init state dimension for local ana. domain
    PROCEDURE(init_dim_obs_l_cb) :: init_dim_obs_l_pdaf !< Initialize local dimimension of obs. vector
    PROCEDURE(cvt_ens_cb) :: cvt_ens_pdaf               !< Apply control vector transform matrix to control vector
    PROCEDURE(cvt_adj_ens_cb) :: cvt_adj_ens_pdaf       !< Apply adjoint control vector transform matrix
    PROCEDURE(obs_op_lin_cb) :: obs_op_lin_pdaf         !< Linearized observation operator
    PROCEDURE(obs_op_adj_cb) :: obs_op_adj_pdaf         !< Adjoint observation operator

! *** OMI-provided procedures ***
    PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb         !< Initialize full observation vector
    PROCEDURE(init_obs_l_cb) :: PDAFomi_init_obs_l_cb       !< Initialize local observation vector
    PROCEDURE(init_obsvar_cb) :: PDAFomi_init_obsvar_cb     !< Initialize mean observation error variance
    PROCEDURE(init_obsvar_l_cb) :: PDAFomi_init_obsvar_l_cb !< Initialize local mean observation error variance
    PROCEDURE(prodRinvA_cb) :: PDAFomi_prodRinvA_cb         !< Provide product R^-1 A
    PROCEDURE(g2l_obs_cb) :: PDAFomi_g2l_obs_cb             !< Restrict full obs. vector to local analysis domain
    PROCEDURE(prodRinvA_l_cb) :: PDAFomi_prodRinvA_l_cb     !< Provide product R^-1 A on local analysis domain


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

    IF (TRIM(filterstr) == '3DVAR' .AND. subtype_filter==1) THEN
       CALL PDAF_assim_offline_en3dvar_lestkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_prodRinvA_cb, &
            cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
            init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
            PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, init_dim_l_pdaf, &
            init_dim_obs_l_pdaf,  PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, PDAFomi_g2l_obs_cb, &
            PDAFomi_init_obsvar_cb, PDAFomi_init_obsvar_l_cb, prepoststep_pdaf, outflag)

    ELSEIF (TRIM(filterstr) == '3DVAR' .AND. subtype_filter==2) THEN
       CALL PDAF_assim_offline_en3dvar_estkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_prodRinvA_cb, &
            cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
            PDAFomi_init_obsvar_cb, prepoststep_pdaf, outflag)
    ELSE
       WRITE (*,*) 'PDAF-ERROR: No valid filter type for PDAF3_assim_offline_3dvar'
    END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

    CALL PDAFomi_dealloc()

  END SUBROUTINE PDAF3_assim_offline_en3dvar



!-------------------------------------------------------------------------------
!> Interface to PDAF for en3dvar/ESKTF
!!
!! __Revision history:__
!! * 2025-04 - Lars Nerger - Initial code based on PDAF3_put_state code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF3_assim_offline_en3dvar_estkf(init_dim_obs_pdaf, obs_op_pdaf, &
       cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
       prepoststep_pdaf, outflag)

    USE PDAF_mod_core, ONLY: filterstr, subtype_filter
    USE PDAF_cb_procedures
    USE PDAFomi, ONLY: PDAFomi_dealloc
    USE PDAFassimilate_en3dvar_estkf, ONLY: PDAF_assim_offline_en3dvar_estkf

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(out) :: outflag      !< Status flag

! *** Argument procedures ***
    PROCEDURE(prepost_cb) :: prepoststep_pdaf           !< User supplied pre/poststep routine
    PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf     !< Initialize dimension of full observation vector
    PROCEDURE(obs_op_cb) :: obs_op_pdaf                 !< Full observation operator
    PROCEDURE(cvt_ens_cb) :: cvt_ens_pdaf               !< Apply control vector transform matrix to control vector
    PROCEDURE(cvt_adj_ens_cb) :: cvt_adj_ens_pdaf       !< Apply adjoint control vector transform matrix
    PROCEDURE(obs_op_lin_cb) :: obs_op_lin_pdaf         !< Linearized observation operator
    PROCEDURE(obs_op_adj_cb) :: obs_op_adj_pdaf         !< Adjoint observation operator

! *** OMI-provided procedures ***
    PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb     !< Initialize full observation vector
    PROCEDURE(init_obsvar_cb) :: PDAFomi_init_obsvar_cb !< Initialize mean observation error variance
    PROCEDURE(prodRinvA_cb) :: PDAFomi_prodRinvA_cb     !< Provide product R^-1 A


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

    IF (TRIM(filterstr) == '3DVAR' .AND. subtype_filter==1) THEN
       CALL PDAF_assim_offline_en3dvar_estkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_prodRinvA_cb, &
            cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
            PDAFomi_init_obsvar_cb, prepoststep_pdaf, outflag)
    ELSE
       WRITE (*,*) 'PDAF-ERROR: No valid filter type for PDAF3_assim_offline_3dvar'
    END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

    CALL PDAFomi_dealloc()

  END SUBROUTINE PDAF3_assim_offline_en3dvar_estkf

!-------------------------------------------------------------------------------
!> Interface to PDAF for en3dvar/LESKTF
!!
!! __Revision history:__
!! * 2025-04 - Lars Nerger - Initial code based on PDAF3_put_state code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF3_assim_offline_en3dvar_lestkf(init_dim_obs_pdaf, obs_op_pdaf, &
       cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
       init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
       prepoststep_pdaf, outflag)

    USE PDAF_mod_core, ONLY: filterstr, subtype_filter
    USE PDAF_cb_procedures
    USE PDAFomi, ONLY: PDAFomi_dealloc
    USE PDAFlocal, &
         ONLY: PDAFlocal_g2l_cb, &       !< Project global to local state vector
         PDAFlocal_l2g_cb                !< Project local to global state vecto
    USE PDAFassimilate_en3dvar_lestkf, ONLY: PDAF_assim_offline_en3dvar_lestkf

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(out) :: outflag      !< Status flag

! *** Argument procedures ***
    PROCEDURE(prepost_cb) :: prepoststep_pdaf           !< User supplied pre/poststep routine
    PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf     !< Initialize dimension of full observation vector
    PROCEDURE(obs_op_cb) :: obs_op_pdaf                 !< Full observation operator
    PROCEDURE(init_n_domains_cb) :: init_n_domains_pdaf !< Provide number of local analysis domains
    PROCEDURE(init_dim_l_cb) :: init_dim_l_pdaf         !< Init state dimension for local ana. domain
    PROCEDURE(init_dim_obs_l_cb) :: init_dim_obs_l_pdaf !< Initialize local dimimension of obs. vector
    PROCEDURE(cvt_ens_cb) :: cvt_ens_pdaf               !< Apply control vector transform matrix to control vector
    PROCEDURE(cvt_adj_ens_cb) :: cvt_adj_ens_pdaf       !< Apply adjoint control vector transform matrix
    PROCEDURE(obs_op_lin_cb) :: obs_op_lin_pdaf         !< Linearized observation operator
    PROCEDURE(obs_op_adj_cb) :: obs_op_adj_pdaf         !< Adjoint observation operator

! *** OMI-provided procedures ***
    PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb         !< Initialize full observation vector
    PROCEDURE(init_obs_l_cb) :: PDAFomi_init_obs_l_cb       !< Initialize local observation vector
    PROCEDURE(init_obsvar_cb) :: PDAFomi_init_obsvar_cb     !< Initialize mean observation error variance
    PROCEDURE(init_obsvar_l_cb) :: PDAFomi_init_obsvar_l_cb !< Initialize local mean observation error variance
    PROCEDURE(prodRinvA_cb) :: PDAFomi_prodRinvA_cb         !< Provide product R^-1 A
    PROCEDURE(g2l_obs_cb) :: PDAFomi_g2l_obs_cb             !< Restrict full obs. vector to local analysis domain
    PROCEDURE(prodRinvA_l_cb) :: PDAFomi_prodRinvA_l_cb     !< Provide product R^-1 A on local analysis domain


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

    IF (TRIM(filterstr) == '3DVAR' .AND. subtype_filter==1) THEN
       CALL PDAF_assim_offline_en3dvar_lestkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_prodRinvA_cb, &
            cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
            init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
            PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, init_dim_l_pdaf, &
            init_dim_obs_l_pdaf,  PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, PDAFomi_g2l_obs_cb, &
            PDAFomi_init_obsvar_cb, PDAFomi_init_obsvar_l_cb, prepoststep_pdaf, outflag)
    ELSE
       WRITE (*,*) 'PDAF-ERROR: No valid filter type for PDAF3_assim_offline_3dvar'
    END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

    CALL PDAFomi_dealloc()

  END SUBROUTINE PDAF3_assim_offline_en3dvar_lestkf



!-------------------------------------------------------------------------------
!> Universal interface to PDAF for Hyb3D-Var
!!
!! This routine is just an alias for PDAF3_assim_offline_3dvar_all
!!
!! __Revision history:__
!! * 2025-04 - Lars Nerger - Initial code based on PDAF3_put_state code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF3_assim_offline_hyb3dvar(init_dim_obs_pdaf, obs_op_pdaf, &
       cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
       init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
       prepoststep_pdaf, outflag)

    USE PDAF_mod_core, ONLY: filterstr
    USE PDAF_cb_procedures

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(out) :: outflag      !< Status flag

! *** Argument procedures ***
    PROCEDURE(prepost_cb) :: prepoststep_pdaf           !< User supplied pre/poststep routine
    PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf     !< Initialize dimension of full observation vector
    PROCEDURE(obs_op_cb) :: obs_op_pdaf                 !< Full observation operator
    PROCEDURE(init_n_domains_cb) :: init_n_domains_pdaf !< Provide number of local analysis domains
    PROCEDURE(init_dim_l_cb) :: init_dim_l_pdaf         !< Init state dimension for local ana. domain
    PROCEDURE(init_dim_obs_l_cb) :: init_dim_obs_l_pdaf !< Initialize local dimimension of obs. vector
    PROCEDURE(cvt_ens_cb) :: cvt_ens_pdaf               !< Apply control vector transform matrix to control vector
    PROCEDURE(cvt_adj_ens_cb) :: cvt_adj_ens_pdaf       !< Apply adjoint control vector transform matrix
    PROCEDURE(cvt_cb) :: cvt_pdaf                       !< Apply control vector transform matrix to control vector
    PROCEDURE(cvt_adj_cb) :: cvt_adj_pdaf               !< Apply adjoint control vector transform matrix
    PROCEDURE(obs_op_lin_cb) :: obs_op_lin_pdaf         !< Linearized observation operator
    PROCEDURE(obs_op_adj_cb) :: obs_op_adj_pdaf         !< Adjoint observation operator

! *** OMI-provided procedures ***
    PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb         !< Initialize full observation vector
    PROCEDURE(init_obs_l_cb) :: PDAFomi_init_obs_l_cb       !< Initialize local observation vector
    PROCEDURE(init_obsvar_cb) :: PDAFomi_init_obsvar_cb     !< Initialize mean observation error variance
    PROCEDURE(init_obsvar_l_cb) :: PDAFomi_init_obsvar_l_cb !< Initialize local mean observation error variance
    PROCEDURE(prodRinvA_cb) :: PDAFomi_prodRinvA_cb         !< Provide product R^-1 A
    PROCEDURE(g2l_obs_cb) :: PDAFomi_g2l_obs_cb             !< Restrict full obs. vector to local analysis domain
    PROCEDURE(prodRinvA_l_cb) :: PDAFomi_prodRinvA_l_cb     !< Provide product R^-1 A on local analysis domain


! **************************************************
! *** Call universal put_state interface routine ***
! **************************************************

    IF (TRIM(filterstr) == '3DVAR') THEN
       CALL PDAF3_assim_offline_3dvar_all(init_dim_obs_pdaf, obs_op_pdaf, &
            cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
            init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
            prepoststep_pdaf, outflag)
    ELSE
       WRITE (*,*) 'PDAF-ERROR: No valid filter type for PDAF3_assim_offline_hyb3dvar'
    END IF

  END SUBROUTINE PDAF3_assim_offline_hyb3dvar


!-------------------------------------------------------------------------------
!> Interface to PDAF for hyb3dvar/ESKTF
!!
!! __Revision history:__
!! * 2025-04 - Lars Nerger - Initial code based on PDAF3_put_state code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF3_assim_offline_hyb3dvar_estkf(init_dim_obs_pdaf, obs_op_pdaf, &
       cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
       obs_op_lin_pdaf, obs_op_adj_pdaf, prepoststep_pdaf, outflag)

    USE PDAF_mod_core, ONLY: filterstr, subtype_filter
    USE PDAF_cb_procedures
    USE PDAFomi, ONLY: PDAFomi_dealloc
    USE PDAFassimilate_hyb3dvar_estkf, ONLY: PDAF_assim_offline_hyb3dvar_estkf

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(out) :: outflag      !< Status flag

! *** Argument procedures ***
    PROCEDURE(prepost_cb) :: prepoststep_pdaf       !< User supplied pre/poststep routine
    PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf !< Initialize dimension of full observation vector
    PROCEDURE(obs_op_cb) :: obs_op_pdaf             !< Full observation operator
    PROCEDURE(cvt_ens_cb) :: cvt_ens_pdaf           !< Apply control vector transform matrix to control vector
    PROCEDURE(cvt_adj_ens_cb) :: cvt_adj_ens_pdaf   !< Apply adjoint control vector transform matrix
    PROCEDURE(cvt_cb) :: cvt_pdaf                   !< Apply control vector transform matrix to control vector
    PROCEDURE(cvt_adj_cb) :: cvt_adj_pdaf           !< Apply adjoint control vector transform matrix
    PROCEDURE(obs_op_lin_cb) :: obs_op_lin_pdaf     !< Linearized observation operator
    PROCEDURE(obs_op_adj_cb) :: obs_op_adj_pdaf     !< Adjoint observation operator

! *** OMI-provided procedures ***
    PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb         !< Initialize full observation vector
    PROCEDURE(init_obsvar_cb) :: PDAFomi_init_obsvar_cb     !< Initialize mean observation error variance
    PROCEDURE(prodRinvA_cb) :: PDAFomi_prodRinvA_cb         !< Provide product R^-1 A
    PROCEDURE(prodRinvA_l_cb) :: PDAFomi_prodRinvA_l_cb     !< Provide product R^-1 A on local analysis domain


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

    IF (TRIM(filterstr) == '3DVAR' .AND. subtype_filter==4) THEN
       CALL PDAF_assim_offline_hyb3dvar_estkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_prodRinvA_cb, &
            cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
            obs_op_lin_pdaf, obs_op_adj_pdaf, PDAFomi_init_obsvar_cb, &
            prepoststep_pdaf, outflag)
    ELSE
       WRITE (*,*) 'PDAF-ERROR: No valid filter type for PDAF3_assim_offline_3dvar'
    END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

    CALL PDAFomi_dealloc()

  END SUBROUTINE PDAF3_assim_offline_hyb3dvar_estkf

!-------------------------------------------------------------------------------
!> Interface to PDAF for hyb3dvar/LESKTF
!!
!! __Revision history:__
!! * 2025-04 - Lars Nerger - Initial code based on PDAF3_put_state code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF3_assim_offline_hyb3dvar_lestkf(init_dim_obs_pdaf, obs_op_pdaf, &
       cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
       init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdaf, &
       prepoststep_pdaf, outflag)

    USE PDAF_mod_core, ONLY: filterstr, subtype_filter
    USE PDAF_cb_procedures
    USE PDAFomi, ONLY: PDAFomi_dealloc
    USE PDAFlocal, &
         ONLY: PDAFlocal_g2l_cb, &       !< Project global to local state vector
         PDAFlocal_l2g_cb                !< Project local to global state vecto
    USE PDAFassimilate_hyb3dvar_lestkf, ONLY: PDAF_assim_offline_hyb3dvar_lestkf

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(out) :: outflag      !< Status flag

! *** Argument procedures ***
    PROCEDURE(prepost_cb) :: prepoststep_pdaf           !< User supplied pre/poststep routine
    PROCEDURE(init_dim_obs_cb) :: init_dim_obs_pdaf     !< Initialize dimension of full observation vector
    PROCEDURE(obs_op_cb) :: obs_op_pdaf                 !< Full observation operator
    PROCEDURE(init_n_domains_cb) :: init_n_domains_pdaf !< Provide number of local analysis domains
    PROCEDURE(init_dim_l_cb) :: init_dim_l_pdaf         !< Init state dimension for local ana. domain
    PROCEDURE(init_dim_obs_l_cb) :: init_dim_obs_l_pdaf !< Initialize local dimimension of obs. vector
    PROCEDURE(cvt_ens_cb) :: cvt_ens_pdaf               !< Apply control vector transform matrix to control vector
    PROCEDURE(cvt_adj_ens_cb) :: cvt_adj_ens_pdaf       !< Apply adjoint control vector transform matrix
    PROCEDURE(cvt_cb) :: cvt_pdaf                       !< Apply control vector transform matrix to control vector
    PROCEDURE(cvt_adj_cb) :: cvt_adj_pdaf               !< Apply adjoint control vector transform matrix
    PROCEDURE(obs_op_lin_cb) :: obs_op_lin_pdaf         !< Linearized observation operator
    PROCEDURE(obs_op_adj_cb) :: obs_op_adj_pdaf         !< Adjoint observation operator

! *** OMI-provided procedures ***
    PROCEDURE(init_obs_cb) :: PDAFomi_init_obs_f_cb         !< Initialize full observation vector
    PROCEDURE(init_obs_l_cb) :: PDAFomi_init_obs_l_cb       !< Initialize local observation vector
    PROCEDURE(init_obsvar_cb) :: PDAFomi_init_obsvar_cb     !< Initialize mean observation error variance
    PROCEDURE(init_obsvar_l_cb) :: PDAFomi_init_obsvar_l_cb !< Initialize local mean observation error variance
    PROCEDURE(prodRinvA_cb) :: PDAFomi_prodRinvA_cb         !< Provide product R^-1 A
    PROCEDURE(g2l_obs_cb) :: PDAFomi_g2l_obs_cb             !< Restrict full obs. vector to local analysis domain
    PROCEDURE(prodRinvA_l_cb) :: PDAFomi_prodRinvA_l_cb     !< Provide product R^-1 A on local analysis domain


! **************************************************
! *** Call the full put_state interface routine  ***
! **************************************************

    IF (TRIM(filterstr) == '3DVAR' .AND. subtype_filter==3) THEN
       CALL PDAF_assim_offline_hyb3dvar_lestkf(init_dim_obs_pdaf, obs_op_pdaf, &
            PDAFomi_init_obs_f_cb, PDAFomi_prodRinvA_cb, &
            cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdaf, obs_op_adj_pdaf, &
            init_dim_obs_pdaf, obs_op_pdaf, PDAFomi_init_obs_f_cb, PDAFomi_init_obs_l_cb, &
            PDAFomi_prodRinvA_l_cb, init_n_domains_pdaf, init_dim_l_pdaf, &
            init_dim_obs_l_pdaf,  PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, PDAFomi_g2l_obs_cb, &
            PDAFomi_init_obsvar_cb, PDAFomi_init_obsvar_l_cb, prepoststep_pdaf, outflag)
    ELSE
       WRITE (*,*) 'PDAF-ERROR: No valid filter type for PDAF3_assim_offline_3dvar'
    END IF


! *******************************************
! *** Deallocate and re-init observations ***
! *******************************************

    CALL PDAFomi_dealloc()

  END SUBROUTINE PDAF3_assim_offline_hyb3dvar_lestkf

END MODULE PDAF3_assim_offline_3dvars
