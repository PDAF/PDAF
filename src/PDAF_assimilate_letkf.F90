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
!> Interface to PDAF for LETKF
!!
!! __Revision history:__
!! * 2013-08 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAFassimilate_letkf

CONTAINS

!> Interface to PDAF for LETKF
!!
!! Interface routine called from the model at each time
!! step during the forecast of each ensemble state. If
!! the time of the next analysis step is reached the
!! forecast state is transferred to PDAF and the analysis
!! is computed by calling PDAF_put_state_letkf. Subsequently, 
!! PDAF_get_state is called to initialize the next forecast
!! phase. 
!!
!! The code is very generic. Basically the only
!! filter-specific part are the calls to the
!! routines PDAF\_put\_state\_X where the analysis
!! is computed and PDAF\_get\_state to initialize the next
!! forecast phase. The filter-specific call-back subroutines 
!! are specified in the calls to the two core routines.
!!
!! Variant for LETKF.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2013-08 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_assimilate_letkf(U_collect_state, U_distribute_state, &
       U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
       U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
       U_g2l_state, U_l2g_state, U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
       U_next_observation, outflag)

    USE PDAF_mod_core, &
         ONLY: cnt_steps, nsteps, assim_flag, reset_fcst_flag, use_PDAF_assim
    USE PDAF_forecast, &
         ONLY: PDAF_fcst_operations
    USE PDAFget_state, &
         ONLY: PDAF_get_state
    USE PDAFput_state_letkf, &
         ONLY: PDAF_put_state_letkf

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(out) :: outflag  !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
    EXTERNAL :: U_collect_state, &   !< Routine to collect a state vector
         U_obs_op, &                 !< Observation operator
         U_init_n_domains_p, &       !< Provide number of local analysis domains
         U_init_dim_l, &             !< Init state dimension for local ana. domain
         U_init_dim_obs, &           !< Initialize dimension of observation vector
         U_init_dim_obs_l, &         !< Initialize dim. of obs. vector for local ana. domain
         U_init_obs, &               !< Initialize PE-local observation vector
         U_init_obs_l, &             !< Init. observation vector on local analysis domain
         U_init_obsvar, &            !< Initialize mean observation error variance
         U_init_obsvar_l, &          !< Initialize local mean observation error variance
         U_g2l_state, &              !< Get state on local ana. domain from full state
         U_l2g_state, &              !< Init full state from state on local analysis domain
         U_g2l_obs, &                !< Restrict full obs. vector to local analysis domain
         U_prodRinvA_l, &            !< Provide product R^-1 A on local analysis domain
         U_prepoststep, &            !< User supplied pre/poststep routine
         U_next_observation, &       !< Routine to provide time step, time and dimension
                                     !<   of next observation
         U_distribute_state          !< Routine to distribute a state vector

! *** Local variables ***
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

       ! Set flags for assimilation and forecast
       assim_flag = 0
       reset_fcst_flag = 1

       ! *** Call analysis step ***

       CALL PDAF_put_state_letkf(U_collect_state, U_init_dim_obs, U_obs_op, &
            U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
            U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
            U_init_obsvar, U_init_obsvar_l, outflag)

       ! *** Prepare start of next ensemble forecast ***

       IF (outflag==0) THEN
          CALL PDAF_get_state(steps, time, doexit, U_next_observation, U_distribute_state, &
               U_prepoststep, outflag)
       END IF

       nsteps = steps

    ELSE
       assim_flag = 0
       reset_fcst_flag = 0
       outflag = 0
    END IF

  END SUBROUTINE PDAF_assimilate_letkf


!-------------------------------------------------------------------------------
!> Interface to PDAF analysis step in offline coupling
!!
!! Interface routine called from the main program
!! for the PDAF offline mode.
!!
!! The code is very generic. Basically the only
!! filter-specific part is the call to the
!! update-routine PDAF\_X\_update where the analysis
!! is computed.  The filter-specific subroutines that
!! are specified in the call to PDAF\_assim\_offline\_X
!! are passed through to the update routine
!!
!! Variant for LETKF.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-04 - Lars Nerger - Initial code based on put_state routine
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_assim_offline_letkf(U_init_dim_obs, U_obs_op, &
       U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
       U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
       U_init_obsvar, U_init_obsvar_l, outflag)

    USE PDAF_mod_core, &
         ONLY: dim_p, dim_ens, assim_flag, step_obs, &
         subtype_filter, screen, flag, offline_mode, &
         state, ens, Ainv, &
         sens, dim_lag, cnt_maxlag
    USE PDAF_mod_parallel, &
         ONLY: mype_world, filterpe
    USE PDAF_utils_filters, &
         ONLY: PDAF_configinfo_filters
    USE PDAFobs, &
         ONLY: dim_obs
    USE PDAF_letkf_update, &
         ONLY: PDAFletkf_update

    IMPLICIT NONE
  
! *** Arguments ***
    INTEGER, INTENT(out) :: outflag  !< Status flag
  
! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
    EXTERNAL :: U_init_dim_obs, &    !< Initialize dimension of observation vector
         U_obs_op, &                 !< Observation operator
         U_init_obs, &               !< Initialize observation vector
         U_init_obs_l, &             !< Init. observation vector on local analysis domain
         U_prodRinvA_l, &            !< Provide product R^-1 A on local analysis domain
         U_init_n_domains_p, &       !< Provide number of local analysis domains
         U_init_dim_l, &             !< Init state dimension for local ana. domain
         U_init_dim_obs_l, &         !< Initialize dim. of obs. vector for local ana. domain
         U_g2l_state, &              !< Get state on local ana. domain from full state
         U_l2g_state, &              !< Init full state from state on local analysis domain
         U_g2l_obs, &                !< Restrict full obs. vector to local analysis domain
         U_init_obsvar, &            !< Initialize mean observation error variance
         U_init_obsvar_l, &          !< Initialize local mean observation error variance
         U_prepoststep               !< User supplied pre/poststep routine

! ** local variables ***
    INTEGER :: i   ! Counter


! *********************************************
! *** Perform analysis step in offline mode ***
! *********************************************

    ! Set flag for assimilation
    assim_flag = 1

    ! Screen output
    IF (mype_world == 0 .AND. screen > 0) THEN
       ! Print configuration info (if not done before in PDAF_set_offline_mode)
       IF (.NOT.offline_mode) CALL PDAF_configinfo_filters(subtype_filter, 1)

       WRITE (*, '(//a5, 64a)') 'PDAF ',('-', i = 1, 64)
       WRITE (*, '(a, 20x, a)') 'PDAF', '+++++ ASSIMILATION +++++'
       WRITE (*, '(a5, 64a)') 'PDAF ', ('-', i = 1, 64)
    ENDIF

    ! Set flag for offline mode
    offline_mode = .true.
 
    OnFilterPE: IF (filterpe) THEN
       CALL PDAFletkf_update(step_obs, dim_p, dim_obs, dim_ens, state, &
            Ainv, ens, U_init_dim_obs, U_obs_op, U_init_obs, &
            U_init_obs_l, U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, &
            U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
            U_init_obsvar, U_init_obsvar_l, U_prepoststep, screen, &
            subtype_filter, dim_lag, sens, cnt_maxlag, flag)
    END IF OnFilterPE


! ********************
! *** finishing up ***
! ********************

    outflag = flag

  END SUBROUTINE PDAF_assim_offline_letkf

END MODULE PDAFassimilate_letkf
