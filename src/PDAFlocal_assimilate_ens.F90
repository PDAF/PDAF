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
!> Interfaces to PDAF for fully-parallel mode for ensemble filters for PDAFlocal
!!
!! The interface routines provide the advanced compact
!! interfaces for using PDAF-Local. The routines
!! just call of one the PDAF_assimilate interface routines
!! with the full interface. In the call the specific PDAF
!! internal subroutines for PDAF-Local are specified.
!!
!! The interface routines provided here are the PDAF-2
!! routines using the naming PDAFlocal_.
!!
!! !  This is a core file of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code by collecting files into a module
!! * Other revisions - see repository log
!!
MODULE PDAFlocal_assimilate_ens

CONTAINS

!-------------------------------------------------------------------------------
!> Interface to PDAF for LSEIK
!!
!! __Revision history:__
!! 2013-08 - Lars Nerger - Initial code
!! 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
!! Other revisions - see repository log
!!
SUBROUTINE PDAFlocal_assimilate_lseik(U_collect_state, U_distribute_state, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
     U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
     U_next_observation, outflag)

  USE PDAF_mod_core, &
       ONLY: cnt_steps, nsteps, assim_flag, use_PDAF_assim
  USE PDAF_mod_parallel, &
       ONLY: mype_world
  USE PDAF_forecast, &
       ONLY: PDAF_fcst_operations
  USE PDAFlocal, &
       ONLY: PDAFlocal_g2l_cb, &       !< Project global to local state vector
       PDAFlocal_l2g_cb                !< Project local to global state vecto
  USE PDAFput_state_lseik, ONLY: PDAF_put_state_lseik

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(out) :: outflag      !< Status flag
  
! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_collect_state, &       !< Routine to collect a state vector
       U_obs_op, &                     !< Observation operator
       U_init_n_domains_p, &           !< Provide number of local analysis domains
       U_init_dim_l, &                 !< Init state dimension for local ana. domain
       U_init_dim_obs, &               !< Initialize dimension of observation vector
       U_init_dim_obs_l, &             !< Initialize dim. of obs. vector for local ana. domain
       U_init_obs, &                   !< Initialize PE-local observation vector
       U_init_obs_l, &                 !< Init. observation vector on local analysis domain
       U_init_obsvar, &                !< Initialize mean observation error variance
       U_init_obsvar_l, &              !< Initialize local mean observation error variance
       U_g2l_obs, &                    !< Restrict full obs. vector to local analysis domain
       U_prodRinvA_l, &                !< Provide product R^-1 A on local analysis domain
       U_prepoststep, &                !< User supplied pre/poststep routine
       U_next_observation, &           !< Routine to provide time step, time and dimension
                                       !<   of next observation
       U_distribute_state              !< Routine to distribute a state vector

! Local variables
  INTEGER :: steps     ! Number of time steps in next forecast phase
  INTEGER :: doexit    ! Exit flag; not used in this variant
  REAL :: time         ! Current model time; not used in this variant


! *****************************
! ***   At each time step   ***
! *****************************

  ! Set flag for using PDAF_assimilate
  use_PDAF_assim = .TRUE.

  ! Increment time step counter
  cnt_steps = cnt_steps + 1

  ! *** Call generic routine for operations during time stepping.          ***
  ! *** Operations are, e.g., IAU or handling of asynchronous observations ***

  CALL PDAF_fcst_operations(cnt_steps, U_collect_state, U_distribute_state, &
     U_init_dim_obs, U_obs_op, U_init_obs, outflag)


! ********************************
! *** At end of forecast phase ***
! ********************************

  IF (cnt_steps == nsteps) THEN

     IF (mype_world==0) WRITE(*,'(a, 5x, a)') 'PDAF', 'Perform assimilation with PDAF'

     ! Set flag for assimilation
     assim_flag = 1

     ! *** Call analysis step ***

     CALL PDAF_put_state_lseik(U_collect_state, U_init_dim_obs, U_obs_op, &
     U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
     U_init_dim_l, U_init_dim_obs_l, PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, U_g2l_obs, &
     U_init_obsvar, U_init_obsvar_l, outflag)

     ! *** Prepare start of next ensemble forecast ***

     IF (outflag==0) THEN
        CALL PDAF_get_state(steps, time, doexit, U_next_observation, U_distribute_state, &
             U_prepoststep, outflag)
     END IF

     nsteps = steps

  ELSE
     assim_flag = 0
     outflag = 0
  END IF

END SUBROUTINE PDAFlocal_assimilate_lseik


!-------------------------------------------------------------------------------
!> Interface to PDAF for LETKF
!!
!! __Revision history:__
!! 2013-08 - Lars Nerger - Initial code
!! 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
!! Other revisions - see repository log
!!
SUBROUTINE PDAFlocal_assimilate_letkf(U_collect_state, U_distribute_state, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
     U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
     U_next_observation, outflag)

  USE PDAF_mod_core, &
       ONLY: cnt_steps, nsteps, assim_flag, use_PDAF_assim
  USE PDAF_mod_parallel, &
       ONLY: mype_world
  USE PDAF_forecast, &
       ONLY: PDAF_fcst_operations
  USE PDAFlocal, &
       ONLY: PDAFlocal_g2l_cb, &       !< Project global to local state vector
       PDAFlocal_l2g_cb                !< Project local to global state vecto
  USE PDAFput_state_letkf, ONLY: PDAF_put_state_letkf

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag      !< Status flag
  
! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_collect_state, &       !< Routine to collect a state vector
       U_obs_op, &                     !< Observation operator
       U_init_n_domains_p, &           !< Provide number of local analysis domains
       U_init_dim_l, &                 !< Init state dimension for local ana. domain
       U_init_dim_obs, &               !< Initialize dimension of observation vector
       U_init_dim_obs_l, &             !< Initialize dim. of obs. vector for local ana. domain
       U_init_obs, &                   !< Initialize PE-local observation vector
       U_init_obs_l, &                 !< Init. observation vector on local analysis domain
       U_init_obsvar, &                !< Initialize mean observation error variance
       U_init_obsvar_l, &              !< Initialize local mean observation error variance
       U_g2l_obs, &                    !< Restrict full obs. vector to local analysis domain
       U_prodRinvA_l, &                !< Provide product R^-1 A on local analysis domain
       U_prepoststep, &                !< User supplied pre/poststep routine
       U_next_observation, &           !< Routine to provide time step, time and dimension
                                       !<   of next observation
       U_distribute_state              !< Routine to distribute a state vector

! Local variables
  INTEGER :: steps     ! Number of time steps in next forecast phase
  INTEGER :: doexit    ! Exit flag; not used in this variant
  REAL :: time         ! Current model time; not used in this variant


! *****************************
! ***   At each time step   ***
! *****************************

  ! Set flag for using PDAF_assimilate
  use_PDAF_assim = .TRUE.

  ! Increment time step counter
  cnt_steps = cnt_steps + 1

  ! *** Call generic routine for operations during time stepping.          ***
  ! *** Operations are, e.g., IAU or handling of asynchronous observations ***

  CALL PDAF_fcst_operations(cnt_steps, U_collect_state, U_distribute_state, &
     U_init_dim_obs, U_obs_op, U_init_obs, outflag)


! ********************************
! *** At end of forecast phase ***
! ********************************

  IF (cnt_steps == nsteps) THEN

     IF (mype_world==0) WRITE(*,'(a, 5x, a)') 'PDAF', 'Perform assimilation with PDAF'

     ! Set flag for assimilation
     assim_flag = 1

     ! *** Call analysis step ***

     CALL PDAF_put_state_letkf(U_collect_state, U_init_dim_obs, U_obs_op, &
     U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
     U_init_dim_l, U_init_dim_obs_l,  PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, U_g2l_obs, &
     U_init_obsvar, U_init_obsvar_l, outflag)

     ! *** Prepare start of next ensemble forecast ***
     IF (outflag==0) THEN
        CALL PDAF_get_state(steps, time, doexit, U_next_observation, U_distribute_state, &
             U_prepoststep, outflag)
     END IF

     nsteps = steps

  ELSE
     assim_flag = 0
     outflag = 0
  END IF

END SUBROUTINE PDAFlocal_assimilate_letkf

!-------------------------------------------------------------------------------
!> Interface to PDAF for LESTKF
!!
!! __Revision history:__
!! 2013-08 - Lars Nerger - Initial code
!! 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
!! Other revisions - see repository log
!!
SUBROUTINE PDAFlocal_assimilate_lestkf(U_collect_state, U_distribute_state, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
     U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
     U_next_observation, outflag)

  USE PDAF_mod_core, &
       ONLY: cnt_steps, nsteps, assim_flag, use_PDAF_assim
  USE PDAF_mod_parallel, &
       ONLY: mype_world
  USE PDAF_forecast, &
       ONLY: PDAF_fcst_operations
  USE PDAFlocal, &
       ONLY: PDAFlocal_g2l_cb, &       !< Project global to local state vector
       PDAFlocal_l2g_cb                !< Project local to global state vecto
  USE PDAFput_state_lestkf, ONLY: PDAF_put_state_lestkf

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag      !< Status flag
  
! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_collect_state, &       !< Routine to collect a state vector
       U_obs_op, &                     !< Observation operator
       U_init_n_domains_p, &           !< Provide number of local analysis domains
       U_init_dim_l, &                 !< Init state dimension for local ana. domain
       U_init_dim_obs, &               !< Initialize dimension of observation vector
       U_init_dim_obs_l, &             !< Initialize dim. of obs. vector for local ana. domain
       U_init_obs, &                   !< Initialize PE-local observation vector
       U_init_obs_l, &                 !< Init. observation vector on local analysis domain
       U_init_obsvar, &                !< Initialize mean observation error variance
       U_init_obsvar_l, &              !< Initialize local mean observation error variance
       U_g2l_obs, &                    !< Restrict full obs. vector to local analysis domain
       U_prodRinvA_l, &                !< Provide product R^-1 A on local analysis domain
       U_prepoststep, &                !< User supplied pre/poststep routine
       U_next_observation, &           !< Routine to provide time step, time and dimension
                                       !<   of next observation
       U_distribute_state              !< Routine to distribute a state vector

! Local variables
  INTEGER :: steps     ! Number of time steps in next forecast phase
  INTEGER :: doexit    ! Exit flag; not used in this variant
  REAL :: time         ! Current model time; not used in this variant


! *****************************
! ***   At each time step   ***
! *****************************

  ! Set flag for using PDAF_assimilate
  use_PDAF_assim = .TRUE.

  ! Increment time step counter
  cnt_steps = cnt_steps + 1

  ! *** Call generic routine for operations during time stepping.          ***
  ! *** Operations are, e.g., IAU or handling of asynchronous observations ***

  CALL PDAF_fcst_operations(cnt_steps, U_collect_state, U_distribute_state, &
     U_init_dim_obs, U_obs_op, U_init_obs, outflag)


! ********************************
! *** At end of forecast phase ***
! ********************************

  IF (cnt_steps == nsteps) THEN

     IF (mype_world==0) WRITE(*,'(a, 5x, a)') 'PDAF', 'Perform assimilation with PDAF'

     ! Set flag for assimilation
     assim_flag = 1

     ! *** Call analysis step ***

     CALL PDAF_put_state_lestkf(U_collect_state, U_init_dim_obs, U_obs_op, &
     U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
     U_init_dim_l, U_init_dim_obs_l,  PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, U_g2l_obs, &
     U_init_obsvar, U_init_obsvar_l, outflag)

     ! *** Prepare start of next ensemble forecast ***

     IF (outflag==0) THEN
        CALL PDAF_get_state(steps, time, doexit, U_next_observation, U_distribute_state, &
             U_prepoststep, outflag)
     END IF

     nsteps = steps

  ELSE
     assim_flag = 0
     outflag = 0
  END IF

END SUBROUTINE PDAFlocal_assimilate_lestkf


!-------------------------------------------------------------------------------
!> Interface to PDAF for LNETF
!!
!! __Revision history:__
!! 2014-05 - Paul Kirchgessner - Initial code based on ETKF
!! 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
!! Other revisions - see repository log
!!
SUBROUTINE PDAFlocal_assimilate_lnetf(U_collect_state, U_distribute_state, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
     U_likelihood_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
      U_g2l_obs, U_next_observation, outflag)

  USE PDAF_mod_core, &
       ONLY: cnt_steps, nsteps, assim_flag, use_PDAF_assim
  USE PDAF_mod_parallel, &
       ONLY: mype_world
  USE PDAF_forecast, &
       ONLY: PDAF_fcst_operations
  USE PDAFlocal, &
       ONLY: PDAFlocal_g2l_cb, &       !< Project global to local state vector
       PDAFlocal_l2g_cb                !< Project local to global state vecto
  USE PDAFput_state_lnetf, ONLY: PDAF_put_state_lnetf

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(out) :: outflag      !< Status flag
  
! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_collect_state, &       !< Routine to collect a state vector
       U_obs_op, &                     !< Observation operator
       U_init_n_domains_p, &           !< Provide number of local analysis domains
       U_init_dim_l, &                 !< Init state dimension for local ana. domain
       U_init_dim_obs, &               !< Initialize dimension of observation vector
       U_init_dim_obs_l, &             !< Initialize dim. of obs. vector for local ana. domain
       U_init_obs, &                   !< Initialize PE-local observation vector
       U_init_obs_l, &                 !< Init. observation vector on local analysis domain
       U_g2l_obs, &                    !< Restrict full obs. vector to local analysis domain
       U_likelihood_l, &               !< Compute observation likelihood for an ensemble member
       U_prepoststep, &                !< User supplied pre/poststep routine
       U_next_observation, &           !< Routine to provide time step, time and dimension
                                       !<   of next observation
       U_distribute_state              !< Routine to distribute a state vector

! Local variables
  INTEGER :: steps     ! Number of time steps in next forecast phase
  INTEGER :: doexit    ! Exit flag; not used in this variant
  REAL :: time         ! Current model time; not used in this variant


! *****************************
! ***   At each time step   ***
! *****************************

  ! Set flag for using PDAF_assimilate
  use_PDAF_assim = .TRUE.

  ! Increment time step counter
  cnt_steps = cnt_steps + 1

  ! *** Call generic routine for operations during time stepping.          ***
  ! *** Operations are, e.g., IAU or handling of asynchronous observations ***

  CALL PDAF_fcst_operations(cnt_steps, U_collect_state, U_distribute_state, &
     U_init_dim_obs, U_obs_op, U_init_obs, outflag)


! **********************************************
! ***   At observation time - analysis step  ***
! **********************************************

  IF (cnt_steps == nsteps) THEN
     IF (mype_world==0) WRITE(*,'(a, 5x, a)') 'PDAF', 'Perform assimilation with PDAF - LNETF'

     ! Set flag for assimilation
     assim_flag = 1

     ! *** Call analysis step ***

     CALL PDAF_put_state_lnetf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_init_obs_l, U_prepoststep, U_likelihood_l, U_init_n_domains_p, &
          U_init_dim_l, U_init_dim_obs_l, PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, U_g2l_obs, &
          outflag)

     ! *** Prepare start of next ensemble forecast ***

     IF (outflag==0) THEN
        CALL PDAF_get_state(steps, time, doexit, U_next_observation, &
             U_distribute_state, U_prepoststep, outflag)
     END IF

     nsteps = steps

  ELSE
     assim_flag = 0
     outflag = 0
  END IF


END SUBROUTINE PDAFlocal_assimilate_lnetf


!-------------------------------------------------------------------------------
!> Interface to PDAF for LKNETF
!!
!! __Revision history:__
!! 2017-08 - Lars Nerger - Initial code based on LETKF
!! 2024-08 - Yumeng Chen - Initial code based on non-PDAFlocal routine
!! Other revisions - see repository log
!!
! !USES:
SUBROUTINE PDAFlocal_assimilate_lknetf(U_collect_state, U_distribute_state, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
     U_prodRinvA_l, U_prodRinvA_hyb_l, U_init_n_domains_p, U_init_dim_l, &
     U_init_dim_obs_l, &
      U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
     U_likelihood_l, U_likelihood_hyb_l, &
     U_next_observation, outflag)

  USE PDAF_mod_core, &
       ONLY: cnt_steps, nsteps, assim_flag, use_PDAF_assim
  USE PDAF_mod_parallel, &
       ONLY: mype_world
  USE PDAF_forecast, &
       ONLY: PDAF_fcst_operations
  USE PDAFlocal, &
       ONLY: PDAFlocal_g2l_cb, &       !< Project global to local state vector
       PDAFlocal_l2g_cb                !< Project local to global state vecto
  USE PDAFput_state_lknetf, ONLY: PDAF_put_state_lknetf

  IMPLICIT NONE
  
  INTEGER, INTENT(out) :: outflag      !< Status flag
  
! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_collect_state, &       !< Routine to collect a state vector
       U_obs_op, &                     !< Observation operator
       U_init_n_domains_p, &           !< Provide number of local analysis domains
       U_init_dim_l, &                 !< Init state dimension for local ana. domain
       U_init_dim_obs, &               !< Initialize dimension of observation vector
       U_init_dim_obs_l, &             !< Initialize dim. of obs. vector for local ana. domain
       U_init_obs, &                   !< Initialize PE-local observation vector
       U_init_obs_l, &                 !< Init. observation vector on local analysis domain
       U_init_obsvar, &                !< Initialize mean observation error variance
       U_init_obsvar_l, &              !< Initialize local mean observation error variance
       U_g2l_obs, &                    !< Restrict full obs. vector to local analysis domain
       U_prodRinvA_l, &                !< Provide product R^-1 A on local analysis domain
       U_prodRinvA_hyb_l, &            !< Provide product R^-1 A on local analysis domain with hybrid weight
       U_likelihood_l, &               !< Compute likelihood
       U_likelihood_hyb_l, &           !< Compute likelihood with hybrid weight
       U_prepoststep, &                !< User supplied pre/poststep routine
       U_next_observation, &           !< Routine to provide time step, time and dimension
                                       !<   of next observation
       U_distribute_state              !< Routine to distribute a state vector

! Local variables
  INTEGER :: steps     ! Number of time steps in next forecast phase
  INTEGER :: doexit    ! Exit flag; not used in this variant
  REAL :: time         ! Current model time; not used in this variant


! *****************************
! ***   At each time step   ***
! *****************************

  ! Set flag for using PDAF_assimilate
  use_PDAF_assim = .TRUE.

  ! Increment time step counter
  cnt_steps = cnt_steps + 1

  ! *** Call generic routine for operations during time stepping.          ***
  ! *** Operations are, e.g., IAU or handling of asynchronous observations ***

  CALL PDAF_fcst_operations(cnt_steps, U_collect_state, U_distribute_state, &
     U_init_dim_obs, U_obs_op, U_init_obs, outflag)


! **********************************************
! ***   At observation time - analysis step  ***
! **********************************************

  IF (cnt_steps == nsteps) THEN

     IF (mype_world==0) WRITE(*,'(a, 5x, a)') 'PDAF', 'Perform assimilation with PDAF - LKNETF'

     ! Set flag for assimilation
     assim_flag = 1

     ! *** Call analysis step ***

     CALL PDAF_put_state_lknetf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_prodRinvA_hyb_l, &
          U_init_n_domains_p, &
          U_init_dim_l, U_init_dim_obs_l,  PDAFlocal_g2l_cb, PDAFlocal_l2g_cb, U_g2l_obs, &
          U_init_obsvar, U_init_obsvar_l, U_likelihood_l, U_likelihood_hyb_l, outflag)

     ! *** Prepare start of next ensemble forecast ***

     IF (outflag==0) THEN
        CALL PDAF_get_state(steps, time, doexit, U_next_observation, &
             U_distribute_state, U_prepoststep, outflag)
     END IF

     nsteps = steps

  ELSE
     assim_flag = 0
     outflag = 0
  END IF

END SUBROUTINE PDAFlocal_assimilate_lknetf

END MODULE PDAFlocal_assimilate_ens
