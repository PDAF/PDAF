!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz, Guowei He and Mukund Pondkule (Forschungszentrum Juelich GmbH)
!
!This file is part of TerrSysMP-PDAF
!
!TerrSysMP-PDAF is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TerrSysMP-PDAF is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU LesserGeneral Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public License
!along with TerrSysMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------
!assimilate_pdaf.F90: TerrSysMP-PDAF implementation of routine
!                     'assimilate_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: assimilate_pdaf.F90 1442 2013-10-04 10:35:19Z lnerger $
!BOP
!
! !ROUTINE: assimilate_pdaf - Routine to control perform analysis step
!
! !INTERFACE:
SUBROUTINE assimilate_pdaf()

    ! !DESCRIPTION:
    ! This routine is called during the model integrations at each time
    ! step. It check whether the forecast phase is completed. If so,
    ! PDAF_put_state_X is called to perform the analysis step.
    !
    ! !REVISION HISTORY:
    ! 2013-08 - Lars Nerger - Initial code for NEMO
    ! Later revisions - see svn log
    !
    ! !USES:
    !   USE mod_parallel_model, &     ! Parallelization variables
    !        ONLY: mype_world, abort_parallel
    USE mod_assimilation, &      ! Variables for assimilation
        ONLY: filtertype

    IMPLICIT NONE

    ! !CALLING SEQUENCE:
    ! Called by: step
    ! CAlls: PDAF_assimilate_X
    !EOP

    ! Local variables
    INTEGER :: status_pdaf       ! PDAF status flag


    ! External subroutines
    EXTERNAL :: collect_state_pdaf, & ! Routine to collect a state vector from model fields
        init_dim_obs_pdaf, &         ! Initialize Dimension Of Observation Vector
        obs_op_pdaf, &               ! Implementation of the Observation operator
        init_obs_pdaf, &             ! Routine to provide vector of measurements
        prepoststep_ens_pdaf, &      ! User supplied pre/poststep routine
        prodRinvA_pdaf, &            ! Provide product R^-1 A for some matrix A
        init_obsvar_pdaf, &          ! Initialize mean observation error variance
        next_observation_pdaf, &     ! Provide time step, model time, &
                                     ! and dimension of next observation
        distribute_state_pdaf        ! Routine to distribute a state vector to model fields
    EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
        init_dim_l_pdaf, &              ! Initialize state dimension for local ana. domain
        init_dim_obs_l_pdaf,&           ! Initialize dim. of obs. vector for local ana. domain
        g2l_state_pdaf, &               ! Get state on local ana. domain from global state
        l2g_state_pdaf, &               ! Init global state from state on local analysis domain
        g2l_obs_pdaf, &                 ! Restrict a global obs. vector to local analysis domain
        init_obs_l_pdaf, &              ! Provide vector of measurements for local ana. domain
        prodRinvA_l_pdaf, &             ! Provide product R^-1 A for some local matrix A
        init_obsvar_l_pdaf, &           ! Initialize local mean observation error variance
        init_obs_f_pdaf, &              ! Provide full vector of measurements for PE-local domain
        obs_op_f_pdaf, &                ! Obs. operator for full obs. vector for PE-local domain
        init_dim_obs_f_pdaf, &             ! Get dimension of full obs. vector for PE-local domain
        add_obs_error_pdaf, &
        init_obscovar_pdaf, &
        localize_covar_pdaf

    ! *********************************
    ! *** Call assimilation routine ***
    ! *********************************

    IF (filtertype == 6) THEN
        CALL PDAF_assimilate_estkf(collect_state_pdaf, distribute_state_pdaf, &
            init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf, prepoststep_ens_pdaf, &
            prodRinvA_pdaf, init_obsvar_pdaf, next_observation_pdaf, status_pdaf)
    ELSEIF (filtertype == 7) THEN
        CALL PDAF_assimilate_lestkf(collect_state_pdaf, distribute_state_pdaf, &
            init_dim_obs_f_pdaf, obs_op_f_pdaf, init_obs_f_pdaf, init_obs_l_pdaf, &
            prepoststep_ens_pdaf, prodRinvA_l_pdaf, init_n_domains_pdaf, &
            init_dim_l_pdaf, init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
            g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, next_observation_pdaf, status_pdaf)
    elseif (filtertype == 2) then
        CALL PDAF_assimilate_enkf(collect_state_pdaf, distribute_state_pdaf, &
            init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf, prepoststep_ens_pdaf, &
            add_obs_error_pdaf, init_obscovar_pdaf, next_observation_pdaf, status_pdaf)
    elseif (filtertype == 4) then
        CALL PDAF_assimilate_etkf(collect_state_pdaf, distribute_state_pdaf, &
            init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf, prepoststep_ens_pdaf, &
            prodRinvA_pdaf, init_obsvar_pdaf, next_observation_pdaf, status_pdaf)
    elseif (filtertype == 5) then
        CALL PDAF_assimilate_letkf(collect_state_pdaf, distribute_state_pdaf, &
            init_dim_obs_f_pdaf, obs_op_f_pdaf, init_obs_f_pdaf, init_obs_l_pdaf, &
            prepoststep_ens_pdaf, prodRinvA_l_pdaf, init_n_domains_pdaf, &
            init_dim_l_pdaf, init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf, &
            g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf, next_observation_pdaf, status_pdaf)
    elseif (filtertype == 8) then
        STOP "PDAF_assimilate_lenkf not provided in out version of PDAF1_1"
        !call PDAF_assimilate_lenkf(collect_state_pdaf  , distribute_state_pdaf, &
        !    init_dim_obs_pdaf, obs_op_pdaf, init_obs_pdaf, prepoststep_ens_pdaf, &
        !    localize_covar_pdaf, add_obs_error_pdaf, init_obscovar_pdaf, &
        !    next_observation_pdaf, status_pdaf)
    END IF

  ! Check for errors during execution of PDAF

!   IF (status_pdaf /= 0) THEN
!      WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
!           'ERROR ', status_pdaf, &
!           ' in PDAF_put_state - stopping! (PE ', mype_world,')'
!      CALL  abort_parallel()
!   END IF

END SUBROUTINE assimilate_pdaf
