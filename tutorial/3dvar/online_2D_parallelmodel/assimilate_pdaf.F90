!>  Routine to call PDAF for analysis step in fully-parallel mode
!!
!! This routine is called during the model integrations at each time 
!! step. It calls the filter-specific assimilation routine of PDAF 
!! (PDAFomi_assimilate_X), which checks whether the forecast phase
!! is completed. If so, the analysis step is computed inside PDAF.
!!
!! __Revision history:__
!! * 2021-12 - Lars Nerger - Initial code for 3D-Vars
!! * Later revisions - see repository log
!!
SUBROUTINE assimilate_pdaf()

  USE pdaf_interfaces_module, &   ! Interface definitions to PDAF core routines
       ONLY: PDAFomi_assimilate_3dvar, PDAFomi_assimilate_en3dvar_estkf, &
       PDAFomi_assimilate_hyb3dvar_estkf
  USE PDAFlocal, &                ! Interface definitions for PDAFlocal
       ONLY: PDAFlocalomi_assimilate_en3dvar_lestkf, PDAFlocalomi_assimilate_hyb3dvar_lestkf
  USE mod_parallel_model, &       ! Parallelization variables
       ONLY: mype_world, abort_parallel
  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: subtype

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: status_pdaf          ! PDAF status flag


! *** External subroutines ***
! Subroutine names are passed over to PDAF in the calls to 
! PDAF_get_state and PDAF_put_state_X. This allows the user 
! to specify the actual name of a routine.  
! The PDAF-internal name of a subroutine might be different
! from the external name!

  ! Interface between model and PDAF, and prepoststep
  EXTERNAL :: collect_state_pdaf, &   ! Collect a state vector from model fields
       distribute_state_pdaf, &       ! Distribute a state vector to model fields
       next_observation_pdaf, &       ! Provide time step of next observation
       prepoststep_ens_pdaf           ! User supplied pre/poststep routine
  ! Localization of state vector
  EXTERNAL :: init_n_domains_pdaf, &  ! Provide number of local analysis domains
       init_dim_l_pdaf                ! Initialize state dimension for local analysis domain
  ! Interface to PDAF-OMI for local and global filters
  EXTERNAL :: init_dim_obs_pdafomi, & ! Get dimension of full obs. vector for PE-local domain
       obs_op_pdafomi, &              ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_l_pdafomi         ! Get dimension of obs. vector for local analysis domain
  ! Variational methods: 3D-Var/En3D-Var/Hybrid 3D-Var 
  EXTERNAL :: cvt_ens_pdaf, &         ! Transform control vector into state vector (ensemble var)
       cvt_adj_ens_pdaf, &            ! Apply adjoint control vector transform matrix (ensemble var)
       cvt_pdaf, &                    ! Apply control vector transform matrix to control vector
       cvt_adj_pdaf, &                ! Apply adjoint control vector transform matrix
       prepoststep_3dvar_pdaf, &      ! User supplied pre/poststep routine for parameterized 3D-Var
       obs_op_lin_pdafomi, &          ! PDAF-OMI: Linearized observation operator
       obs_op_adj_pdafomi             ! PDAF-OMI: Adjoint observation operator


! *********************************
! *** Call assimilation routine ***
! *********************************

  IF (subtype==0) THEN
     ! parameterized 3D-Var
     CALL PDAFomi_assimilate_3dvar(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, &
          cvt_pdaf, cvt_adj_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
          prepoststep_3dvar_pdaf, next_observation_pdaf, status_pdaf)
  ELSEIF (subtype==1) THEN
     ! Ensemble 3D-Var with local ESTKF update of ensemble perturbations
     CALL PDAFlocalomi_assimilate_en3dvar_lestkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
          prepoststep_ens_pdaf, next_observation_pdaf, status_pdaf)
  ELSEIF (subtype==4) THEN
     ! Ensemble 3D-Var with global ESTKF update of ensemble perturbations
     CALL PDAFomi_assimilate_en3dvar_estkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
          prepoststep_ens_pdaf, next_observation_pdaf, status_pdaf)
  ELSEIF (subtype==6) THEN
     ! Hybrid 3D-Var with local ESTKF update of ensemble perturbations
     CALL PDAFlocalomi_assimilate_hyb3dvar_lestkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
          obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
          prepoststep_ens_pdaf, next_observation_pdaf, status_pdaf)
  ELSEIF (subtype==7) THEN
     ! Hybrid 3D-Var with global ESTKF update of ensemble perturbations
     CALL PDAFomi_assimilate_hyb3dvar_estkf(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, &
          cvt_ens_pdaf, cvt_adj_ens_pdaf, cvt_pdaf, cvt_adj_pdaf, &
          obs_op_lin_pdafomi, obs_op_adj_pdafomi, &
          prepoststep_ens_pdaf, next_observation_pdaf, status_pdaf)
  END IF


! ************************
! *** Check error flag ***
! ************************

  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in PDAFomi_assimilate - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF

END SUBROUTINE assimilate_pdaf
