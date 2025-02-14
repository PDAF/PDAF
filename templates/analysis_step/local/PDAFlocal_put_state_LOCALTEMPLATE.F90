!> Interface to transfer state to PDAF
!!
!! Interface routine called from the model after the 
!! forecast of each ensemble state to transfer data
!! from the model to PDAF. For the parallelization 
!! this involves transfer from model PEs to filter 
!! PEs.
!!
!! This is just an interface for PDAFlocal to call
!! PDAF_put_state_LOCALTEMPLATE.
!!
!! __Revision history:__
!! * 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
!! * 2024-12 - Lars Nerger - Adaption of LETKF code variant for template
!! *  Later revisions - see repository log
!!
SUBROUTINE PDAFlocal_put_state_LOCALTEMPLATE(U_collect_state, U_init_dim_obs, U_obs_op, &
     U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
     U_init_dim_l, U_init_dim_obs_l, U_g2l_obs, &
     U_init_obsvar, U_init_obsvar_l, outflag)

  USE PDAFlocal, &                ! PDAFlocal routines for state localization
       ONLY: PDAFlocal_g2l_cb, &  ! Project global to local state vector
       PDAFlocal_l2g_cb           ! Project local to global state vecto

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(out) :: outflag  !< Status flag
  
! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_collect_state, &    ! Write model fields into state vector
       U_prepoststep                ! User supplied pre/poststep routine
  ! Observation-related routines for analysis step
  EXTERNAL :: U_init_dim_obs, &     !< Initialize dimension of observation vector
       U_obs_op, &                  !< Observation operator
       U_init_obs, &                !< Initialize observation vector
       U_init_obsvar, &             !< Initialize mean observation error variance
       U_init_dim_obs_l, &          !< Initialize dim. of obs. vector for local ana. domain
       U_init_obs_l, &              !< Init. observation vector on local analysis domain
       U_g2l_obs, &                 !< Restrict full obs. vector to local analysis domain
       U_init_obsvar_l, &           !< Initialize local mean observation error variance
       U_prodRinvA_l                !< Provide product R^-1 A on local analysis domain
  ! Routines for state localization
  EXTERNAL :: U_init_n_domains_p, & !< Provide number of local analysis domains
       U_init_dim_l                 !< Init state dimension for local ana. domain


! *******************************************
! *** Call put_state routine of DA method ***
! *******************************************

  CALL PDAF_put_state_LOCALTEMPLATE(U_collect_state, U_init_dim_obs, U_obs_op, &
       U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
       U_init_dim_l, U_init_dim_obs_l, PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, U_g2l_obs, &
       outflag)

END SUBROUTINE PDAFlocal_put_state_LOCALTEMPLATE
