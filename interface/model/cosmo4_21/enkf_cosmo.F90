subroutine cosmo_init(pdaf_id)  bind(C,name="cosmo_init")

use iso_C_binding
use enkf_cosmo_mod

integer(c_int), intent(in) :: pdaf_id

! Set suffix for INPUT_IO from input
cosmo_input_suffix = pdaf_id

!------------------------------------------------------------------------------
!- End of Header
!------------------------------------------------------------------------------
#ifdef FPEABORT
! Floating point exception trapping
      CALL initialize_fpe_trap(.TRUE., ierrstat)
      IF ( ierrstat /= 0 ) THEN
        STOP 'Error initializing ftp_trap'
      END IF
#endif

  ierrstat = 0_iintegers
  izerror  = 0_iintegers

!------------------------------------------------------------------------------
!- Section 1: Setup of the model and Namelist Input for all components
!------------------------------------------------------------------------------

  ! Section 1 has to be computed by all PEs (Compute and IO-PEs)
  ! After Section 1 the PEs for computing and IO are splitted

  CALL organize_setup

  IF (my_cart_id == 0) THEN
    ! Print the default information to stdout:
    CALL info_define ('lmparbin')          ! Pre-define the program name
    CALL info_readnl ('INPUT_COSMO')       ! Read additional information from namelist file
    CALL info_print ()                     ! Print the information to stdout
  ENDIF

  ! Initialize, whether debug output shall be done
  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF

  ! Input of the dynamics namelist
  CALL organize_dynamics ('input', izerror, yzerrmsg, dt, .FALSE.)
  IF (izerror /= 0_iintegers) THEN
    CALL model_abort (my_world_id, 100+izerror, yzerrmsg,               &
                                   'organize_dynamics: input')
  ENDIF

  ! Input of the physics namelist
  CALL organize_physics ('input', izerror, yzerrmsg)
  IF (izerror /= 0_iintegers) THEN
    CALL model_abort (my_world_id, 100+izerror, yzerrmsg,             &
                                   'organize_physics: input')
  ENDIF

#ifdef COSMOART
  ! Input of the COSMO_ART namelist
  IF (l_cosmo_art) THEN
    CALL organize_cosmo_art ('input', ydate_ini, izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_world_id, 100+izerror, yzerrmsg,             &
                                   'organize_cosmo_art: input')
    ENDIF
  ENDIF
#endif

#ifdef POLLEN
  ! Input of the Pollen namelist
  IF (l_pollen) THEN
    CALL organize_pollen ('input', ydate_ini, izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_world_id, 100+izerror, yzerrmsg,             &
                                   'organize_pollen: input')
    ENDIF
  ENDIF
#endif

  ! Input of the diagnosis namelist
  CALL organize_diagnosis ('input', izerror, yzerrmsg)
  IF (izerror /= 0_iintegers) THEN
    CALL model_abort (my_world_id, 100+izerror, yzerrmsg,             &
                                   'organize_diagnosis: input')
  ENDIF

#ifdef NUDGING
  ! Input of the assimilation namelist
  IF (luseobs) THEN
    CALL organize_assimilation ('input', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_world_id, 100+izerror, yzerrmsg,             &
                                     'organize_assimilation: input')
    ENDIF
  ENDIF
#endif

  ! Input of the EPS namelist
  IF (leps) THEN
    CALL organize_eps ('input', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_world_id, 100+izerror, yzerrmsg,             &
                                      'organize_eps: input')
    ENDIF
  ENDIF

  ! Input of the namelists for the I/O-package
  zgrids_dt(1) = dt
  CALL organize_data ('input-init', 0, 1, 1, zgrids_dt,    &
                       izerror, yzerrmsg)
  IF (izerror /= 0_iintegers) THEN
    CALL model_abort (my_world_id, 100+izerror, yzerrmsg,               &
                                   'organize_data: input-init')
  ENDIF

#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
  IF (luse_rttov) THEN
    ! Input of the namelists for the I/O-package
    CALL organize_satellites ('input', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_world_id, 100+izerror, yzerrmsg,             &
                                     'organize_satellites: input')
    ENDIF
  ENDIF
#endif

!AK (20.03.12)
  ! Input of the namelists for tracer transport
  CALL organize_tracer(izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,              &
                                     'organize_tracer')
    ENDIF
!AK (20.03.12)

  ! Initialize the timings
  CALL get_free_unit (nutiming)
  CALL init_timings (nstart, nstop, dt, itype_timing, ldfi,             &
                     lphys, luseobs, l2tls, lsemi_imp, l_cosmo_art, izerror)
  IF (izerror /= 0) THEN
    ! no system clock present
    ltime = .FALSE.
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Allocation of space and computation of constant fields
!------------------------------------------------------------------------------

! Now comes the part for the compute PEs. The part for the IO-PEs is in
! the ELSE-part at the end of the program.

!comp_pe: IF (lcompute_pe) THEN

  ! allocate space
  IF (izdebug > 0) THEN
    PRINT *,'    ALLOCATE SPACE'
  ENDIF

  ! fields for the meteorological variables
  CALL organize_allocation ('default', izerror)

  IF (izerror /= 0) THEN
    ierrstat = 1004
    yzerrmsg  = ' ERROR    *** Allocation of space for meteofields failed ***'
    CALL model_abort (my_cart_id, ierrstat, yzerrmsg, 'allocation: default')
  ENDIF
  
#ifdef COSMOART
  IF (l_cosmo_art) THEN
    CALL organize_cosmo_art ('allocate', ydate_ini, izerror, yzerrmsg)

    IF (izerror /= 0) THEN
      ierrstat = 1005
      yzerrmsg  = ' ERROR    *** Allocation of space for COSMO_ART failed ***'
      CALL model_abort (my_cart_id, ierrstat, yzerrmsg,             &
                                    'organize_cosmo_art: allocate')
    ENDIF
  ENDIF
#endif

#ifdef POLLEN
  IF (l_pollen) THEN
    CALL organize_pollen ('allocate', ydate_ini, izerror, yzerrmsg)

    IF (izerror /= 0) THEN
      ierrstat = 1005
      yzerrmsg  = ' ERROR    *** Allocation of space for Pollen failed ***'
      CALL model_abort (my_cart_id, ierrstat, yzerrmsg,             &
                                    'organize_pollen: allocate')
    ENDIF
  ENDIF
#endif

  CALL constant_fields

!------------------------------------------------------------------------------
!- Section 3: Input of first data sets
!------------------------------------------------------------------------------

  ! Read or generate initial data and the first boundary data sets
  zgrids_dt(1) = dt
  CALL organize_data ('start', 0, 1, 1, zgrids_dt, izerror, yzerrmsg)
  IF (izerror /= 0_iintegers) THEN
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,                &
                                   'start: input-init')
  ENDIF

!------------------------------------------------------------------------------
!- Section 4: Initializations and allocation of extra space
!------------------------------------------------------------------------------

  IF (izdebug > 0) THEN
    PRINT *, '  INITIALIZATIONS'
  ENDIF

  ! 3.1:  Initialization of different packages
  CALL organize_dynamics ('init', izerror, yzerrmsg, dt, .FALSE.)
  IF (izerror /= 0_iintegers) THEN
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,                &
                                   'organize_dynamics: init')
  ENDIF

!MU (17.10.12)
  ! Initialization of tracer
  CALL organize_tracer_init(izerror,yzerrmsg)
  IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,              &
                                     'organize_tracer_init: init')
  ENDIF
!MU (17.10.12)

  ! Initialization of the physics
  IF (lphys) THEN
    CALL organize_physics ('init', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,              &
                                     'organize_physics: init')
    ENDIF

    ! allocate fields for meteorological variables concerned to the canopy
    CALL organize_allocation ('canopy', izerror)
    IF (izerror /= 0) THEN
      ierrstat = 1004
      yzerrmsg  = ' ERROR    *** Allocation of extra space failed ***'
      CALL model_abort (my_cart_id, ierrstat, yzerrmsg, 'allocation: canopy')
    ENDIF
  ENDIF
  
#ifdef COSMOART
  ! Initialization of COSMO_ART
  IF (l_cosmo_art) THEN
    CALL organize_cosmo_art ('init', ydate_ini, izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,              &
                                    'organize_cosmo_art: init')
    ENDIF

    ! Initial profiles of gas phase species of LM_ART
    IF (lgas) THEN
      CALL organize_cosmo_art ('start_gas', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,              &
                                    'organize_cosmo_art: start_gas')
      ENDIF
    ENDIF
  ENDIF
#endif

#ifdef POLLEN
  ! Initialization of the pollen
  IF (l_pollen) THEN
    CALL organize_pollen ('init', ydate_ini, izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,              &
                                     'organize_pollen: init')
    ENDIF
  ENDIF
#endif

  ! Initialization of the diagnosis
  IF (ldiagnos) THEN
    CALL organize_diagnosis ('init', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,              &
                                     'organize_diagnosis: init')
    ENDIF
  ENDIF

#ifdef NUDGING
  ! Initialization of the assimilation
  IF (luseobs) THEN
    CALL organize_assimilation ('init', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,              &
                                     'organize_assimilation: init')
    ENDIF
  ENDIF
#endif

#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
  IF (luse_rttov) THEN
    ! initialization of variables for the synthetic satellite computations
    CALL organize_satellites ('init', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,              &
                                     'organize_satellites: init')
    ENDIF
  ENDIF
#endif

!------------------------------------------------------------------------------
! Section 5: Model initialization by digital filtering
!------------------------------------------------------------------------------

  IF (ldfi) THEN
    IF (izdebug > 0) THEN
      PRINT *, '    DIGITAL FILTERING'
    ENDIF
    IF(l2tls) THEN
      IF (izdebug > 0) THEN
        PRINT *, ' **** CAUTION **** CAUTION **** CAUTION **** CAUTION ****'
        PRINT *, ' ****   DIGITAL FILTERING not tested for 2TL-scheme  ****'
        PRINT *, ' **** **** ****  Proceed on your own risk  **** **** ****'
      ENDIF
    ENDIF
    CALL dfi_initialization (lana_qi, llb_qi, llb_qr_qs, llb_qg,      &
                             lbd_frame, undef, izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,            &
                                     'dfi_initialization')
    ENDIF
  ENDIF
  
  ! Close file for control output
  IF (my_cart_id == 0) THEN
    CLOSE (nuspecif, STATUS='KEEP')
  ENDIF

  ! halve the time step, if ntstep = 0
  ! (in src_setup: nstep = nstart)
  IF (.NOT. l2tls) THEN
    IF (ntstep == 0) THEN
      zforecasttime = 0.0
      dt = 0.5 * dt
    ELSE
      zforecasttime = 0.5*dt + (ntstep-1)*dt
      nzdiv         = INT (zforecasttime / 3600.0_ireals, iintegers)
      zforecasttime = zforecasttime - nzdiv * 3600_ireals
    ENDIF
  ELSE
    IF (ntstep == 0) THEN
      zforecasttime = 0.0
    ELSE
      zforecasttime = ntstep*dt
      nzdiv         = INT (zforecasttime / 3600.0_ireals, iintegers)
      zforecasttime = zforecasttime - nzdiv * 3600_ireals
    ENDIF
  ENDIF

  IF (lbdclim) THEN
    ynote         = '...... FORECAST TIME IS NOW xxxxxx DAYS ......'
  ELSE
    ynote         = '...... FORECAST TIME IS NOW xxx HOURS   ......'
  ENDIF

  IF (ltime) THEN
    CALL get_timings (i_initializations, ntstep, dt, izerror)
  ENDIF

  !------------------------------------------------------------------------------
  ! 5a. Temperature disturbance(s) (either in the air or at the soil surface
  !    or within the soil) in the initial conditions
  !------------------------------------------------------------------------------

  ! There are different types of possible disturbances, see the documentation
  ! of the corresponding namelist parameters in INPUT_IDEAL, and there is the possibility
  ! to specify more than one disturbance (up to 50 right now).
  
  IF (lartif_data) THEN
    CALL set_tempdist(nnew)
    ! Initial condition on t_so (takes only effect if lsoil=.true.)
    CALL set_tempdist_tso(nnew)
  END IF

#if defined COUP_OAS_COS
  CALL oas_cos_define
!OASIS4 only
!  CALL oas_cos_update_time(0)
#endif

!kuw initialize cos_start
cos_start = nstart
!kuw end

end subroutine cosmo_init


subroutine cosmo_advance(cos_dt)  bind(C,name="cosmo_advance")

use iso_C_binding
use enkf_cosmo_mod
integer(c_int),intent(in) :: cos_dt
!------------------------------------------------------------------------------
!- Section 6: Time stepping
!------------------------------------------------------------------------------

  IF (izdebug > 0) THEN
    PRINT *, '  TIME STEPPING'
  ENDIF

  !timeloop: DO ntstep = nstart , nstop
  write(*,*)'advancing cosmo from ',cos_start,' to ',(cos_start+cos_dt-1)

  timeloop: DO ntstep = cos_start,(cos_start+cos_dt-1)

    IF ( (izdebug > 1) .AND. (.NOT. lbdclim)) THEN
      PRINT *, '    STEP ',ntstep
    ENDIF

    !--------------------------------------------------------------------------
    !- Section 6.1: Initialization of this time step
    !--------------------------------------------------------------------------

!AK (20.03.12)
    tracertens(:,:,:,:) = 0._ireals
    CALL organize_tracer_bound(izerror,yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                     'organize_tracer_bound: compute')
    ENDIF

    CALL organize_tracer_source(izerror,yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_tracer_source: compute')
    ENDIF
!AK (20.03.12)

    CALL initialize_loop (ntstep, nbd1, nbd2, nold, nnow, nnew)
 
    IF (ltime) CALL get_timings (i_add_computations, ntstep, dt, izerror)

    !--------------------------------------------------------------------------
    !- Section 6.2.1: physics 
    !--------------------------------------------------------------------------

#ifdef COSMOART
    ! CK 20101204 unit conversion necessary before physics
    ! more universal approach: get a general injection point for ART
    IF (l_cosmo_art) THEN
      CALL organize_cosmo_art ('prepare_physics', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_cosmo_art: prepare_physics')
      ENDIF
    ENDIF
#endif

#ifdef POLLEN
    ! Preparations for Pollen
    IF (l_pollen) THEN
      CALL organize_pollen ('prepare_physics', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_pollen: prepare_physics')
      ENDIF
    ENDIF
#endif

    IF (lphys) CALL organize_physics ('compute', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                     'organize_physics: compute')
    ENDIF

    IF (ltime) CALL get_timings (i_phy_computations, ntstep, dt, izerror)

    !--------------------------------------------------------------------------
    !- Section 6.2.1a: After the call to the soil model, add extra
    !  artificial heating rate disturbances on the soil temperature if desired
    !  (formulated as Euler forward updates of the temperatures)
    !--------------------------------------------------------------------------

    IF (lartif_data) THEN
      ! Set possible artificial heating rate disturbance(s) in the soil 
      ! (affects t_so or t_s/t_m/t_cl depending on soil model
      ! and takes effect only IF lsoil=.TRUE.). 
      ! Because the soil model has already done the time integration,
      ! the artificial disturbances have to be imposed on
      ! timelevel nnew:
      CALL artif_heatrate_dist_tso(nnew)
    END IF

#ifdef COSMOART
    ! CK 20101204 unit conversion necessary after physics
    ! more universal approach: get a general injection point for ART
    IF (l_cosmo_art) THEN
      CALL organize_cosmo_art ('finalize_physics', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_cosmo_art: finalize_physics')
      ENDIF
    ENDIF
#endif

#ifdef POLLEN
    ! Clean up after Pollen
    IF (l_pollen) THEN
      CALL organize_pollen ('finalize_physics', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_pollen: finalize_physics')
      ENDIF
    ENDIF
#endif

#ifdef COSMOART
    !--------------------------------------------------------------------------
    !- Section 6.2.2: emissions for COSMO_ART
    !--------------------------------------------------------------------------

    IF (l_cosmo_art) THEN
      CALL organize_cosmo_art ('emiss', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_cosmo_art: emiss')
      ENDIF
    ENDIF
#endif

#ifdef POLLEN
    IF (l_pollen) THEN
      CALL organize_pollen ('emiss', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_pollen: emiss')
      ENDIF
    ENDIF
#endif

    !--------------------------------------------------------------------------
    !- Section 6.3.1: dynamics
    !--------------------------------------------------------------------------

#ifdef COSMOART
    ! CK 20101204 unit conversion necessary before dynamics
    ! more universal approach: get a general injection point for ART
    IF (l_cosmo_art) THEN
      CALL organize_cosmo_art ('prepare_dynamics', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_cosmo_art: prepare_dynamics')
      ENDIF
    ENDIF
#endif

#ifdef POLLEN
    ! Preparations for Pollen
    IF (l_pollen) THEN
      CALL organize_pollen ('prepare_dynamics', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_pollen: prepare_dynamics')
      ENDIF
    ENDIF
#endif

    CALL organize_dynamics ('compute', izerror, yzerrmsg, dt, .FALSE.)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                    'organize_dynamics: compute')
    ENDIF

#ifdef COSMOART
    ! CK 20101204 unit conversion necessary after dynamics
    ! more universal approach: get a general injection point for ART
    IF (l_cosmo_art) THEN
      CALL organize_cosmo_art ('finalize_dynamics', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_cosmo_art: finalize_dynamics')
      ENDIF
    ENDIF
#endif

#ifdef POLLEN
    ! Clean up after Pollen and washout
    IF (l_pollen) THEN
      CALL organize_pollen ('finalize_dynamics', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_pollen: finalize_dynamics')
      ENDIF
    ENDIF
#endif

    CALL set_qrqsqg_boundaries

    IF (ltime) CALL get_timings (i_dyn_computations, ntstep, dt, izerror)

#ifdef COSMOART
    !--------------------------------------------------------------------------
    !- Section 6.3a: chemistry and aerosol dynamics, COSMO_ART; Pollen
    !--------------------------------------------------------------------------

    IF (l_cosmo_art) THEN

      IF (laero) THEN
        CALL organize_cosmo_art ('init_aero', ydate_ini,izerror, yzerrmsg)
      ENDIF

      CALL organize_cosmo_art ('compute', ydate_ini, izerror, yzerrmsg)

      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_cosmo_art: chemie')
      ENDIF

    ENDIF
#endif

#ifdef POLLEN
    IF (l_pollen) THEN
      WHERE (cpollen(:,:,:,:,nnew) < 1.0E-15_ireals)
        cpollen(:,:,:,:,nnew) = 1.0E-14_ireals
      ENDWHERE
    ENDIF
#endif

    !--------------------------------------------------------------------------
    !- Section 6.4: nudging
    !--------------------------------------------------------------------------

#ifdef NUDGING
    IF (luseobs) CALL organize_assimilation ('nudge', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                    'organize_assimilation: nudge')
    ENDIF
#endif

    !--------------------------------------------------------------------------
    !- Section 6.4a: latent heat nudging (LHN)
    !--------------------------------------------------------------------------

#ifdef NUDGING
    IF (luseobs) CALL organize_assimilation ('lhn', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                    'organize_assimilation: lhn')
    ENDIF

    IF (ltime) CALL get_timings (i_lhn_computations, ntstep, dt, izerror)
#endif

    !--------------------------------------------------------------------------
    !- Section 6.5: water budget
    !--------------------------------------------------------------------------
 
    IF (ldiagnos .AND. (l2tls .OR. (ntstep > 0))) THEN
      ! for the leapfrog scheme the summations in diagbudget must not be done
      ! in the first intermediate step ntstep==0. These calculations are done
      ! again for ntstep==1.
      CALL organize_diagnosis ('diagbudget', izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,             &
                                     'organize_diagnosis: diagbudget')
      ENDIF
    ENDIF

    IF (ltime) CALL get_timings (i_add_computations, ntstep, dt, izerror)

    !--------------------------------------------------------------------------
    !- Section 6.6: spectral nudging and relaxation
    !--------------------------------------------------------------------------

    IF (lspecnudge) THEN
      CALL organize_dynamics ('specnudge', izerror, yzerrmsg, dt, .FALSE.)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                      'organize_dynamics: specnudge')
      ENDIF
      IF (ltime) CALL get_timings (i_spectr_nudging, ntstep, dt, izerror)
    ENDIF

#ifdef COSMOART
    ! CK 20101204 unit conversion necessary before relaxation
    ! more universal approach: get a general injection point for ART
    IF (l_cosmo_art) THEN
      CALL organize_cosmo_art ('prepare_relaxation', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_cosmo_art: prepare_relaxation')
      ENDIF
    ENDIF
#endif

#ifdef POLLEN
    ! Preparations for Pollen
    IF (l_pollen) THEN
      CALL organize_pollen ('prepare_relaxation', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_pollen: prepare_relaxation')
      ENDIF
    ENDIF
#endif

    CALL organize_dynamics ('relaxation', izerror, yzerrmsg, dt, .FALSE.)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                    'organize_dynamics: relaxation')
    ENDIF
    IF (ltime) CALL get_timings (i_relaxation, ntstep, dt, izerror)

#ifdef COSMOART
    ! CK 20101204 unit conversion necessary after dynamics
    ! more universal approach: get a general injection point for ART
    IF (l_cosmo_art) THEN
      CALL organize_cosmo_art ('finalize_relaxation', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_cosmo_art: finalize_relaxation')
      ENDIF
    ENDIF
#endif

#ifdef POLLEN
    ! Clean up after Pollen
    IF (l_pollen) THEN
      CALL organize_pollen ('finalize_relaxation', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_pollen: finalize_relaxation')
      ENDIF
    ENDIF
#endif

    ! Final update of temperature and humidity variables due to
    ! cloud microphysics in case of the cloud ice scheme
    IF (lphys) CALL organize_physics ('finish_compute', izerror, yzerrmsg)

    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                     'organize_physics: finish_compute')
    ENDIF

    CALL nullify_tracers

#ifdef COSMOART
    ! CK 20101204 setting minima for tracers
    ! more universal approach: get a general injection point for ART
    IF (l_cosmo_art) THEN
      CALL organize_cosmo_art ('finalize', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_cosmo_art: finalize')
      ENDIF
    ENDIF
#endif

    IF (ltime) CALL get_timings (i_phy_computations, ntstep, dt, izerror)

    !--------------------------------------------------------------------------
    !- Section 6.7: Exchange of boundary data
    !--------------------------------------------------------------------------
 
      ! The calls to exchg_boundaries have to be here in any case, even for a
      ! sequential version, because of possible periodic boundary conditions 
      ! the check, which kind of communication is necessary, is done within 
      ! the subroutine now.

      ! Check, whether additional communication for the convection is
      ! necessary
      lzconv = lconv .AND. lconf_avg .AND.                         &
                ((ntstep < 1) .OR. (MOD(ntstep+2,nincconv)==0))

      IF (ltime_barrier) THEN
        CALL comm_barrier (icomm_cart, ierrstat, yzerrmsg)
        IF (ltime) CALL get_timings (i_barrier_waiting_dyn, ntstep, dt, izerror)
      ENDIF

      IF     ( l2tls .AND. irunge_kutta /= 0 ) THEN
        CALL exchange_runge_kutta
      ELSEIF ( l2tls .AND. irunge_kutta == 0 ) THEN
        CALL exchange_2timelevel
      ELSE ! Leapfrog:
        CALL exchange_leapfrog
      ENDIF

      IF (ltime) CALL get_timings (i_communications_dyn, ntstep, dt, izerror)

    !--------------------------------------------------------------------------
    !- Section 6.8: diagnostics
    !--------------------------------------------------------------------------

    CALL near_surface (nnow)

    !   Analysis of near surface parameters
    !   -----------------------------------

#ifdef NUDGING
    IF (luseobs) CALL organize_assimilation ('surface', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                    'organize_assimilation: surface')
    ENDIF
#endif

    IF (ldiagnos) THEN
      CALL organize_diagnosis ('compute', izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,             &
                                       'organize_diagnosis: compute')
      ENDIF
    ENDIF

#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
    IF (luse_rttov) THEN
      CALL organize_satellites ('compute', izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,             &
                                       'organize_satellites: compute')
      ENDIF
    ENDIF
#endif

    IF (ltime) CALL get_timings (i_add_computations, ntstep, dt, izerror)

    !--------------------------------------------------------------------------
    !- Section 6.9: output of results
    !--------------------------------------------------------------------------

    zgrids_dt(1) = dt
    CALL organize_data ('result', ntstep, 1, 1, zgrids_dt,  &
                         izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                     'result: input-init')
    ENDIF

    IF (ltime) CALL get_timings (i_output, ntstep, dt, izerror)

    !--------------------------------------------------------------------------
    !- Section 6.10: Finalization of this time step
    !--------------------------------------------------------------------------
 
#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
    ! deallocate the satellite variables
    IF (luse_rttov) THEN
      CALL organize_satellites ('dealloc', izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,             &
                                       'organize_satellites: dealloc')
      ENDIF
    ENDIF
#endif

    IF (lbdclim) THEN
      ! record a forecast day
      IF (my_cart_id ==  0) THEN
        zforecasttime = zforecasttime + dt
        IF (zforecasttime >= 86400.0) THEN
          nzdays  = NINT ((ntstep+1)*dt) / 86400
          WRITE (ynote(29:34),'(I6.6)') nzdays
          PRINT *, ynote
          zforecasttime = zforecasttime - 86400.0_ireals
        ENDIF
      ENDIF
    ELSE
      ! record a forecast hour
      IF (my_cart_id ==  0) THEN
        zforecasttime = zforecasttime + dt
        IF (zforecasttime >= 3600.0) THEN
          nzhours = NINT ((ntstep+1)*dt) / 3600
          WRITE (ynote(29:31),'(I3.3)') nzhours
          PRINT *, ynote
          zforecasttime = zforecasttime - 3600.0_ireals
        ENDIF
      ENDIF
    ENDIF

    ! Reset the time step for leapfrog integration
    IF ( ntstep == 0 .AND. (.NOT.l2tls) ) THEN
      dt = 2.0 * dt
    ENDIF

#if defined COUP_OAS_COS
! OASIS4 only
!    CALL oas_cos_update_time(ntstep+1)
#endif

  ENDDO timeloop

  IF (izdebug > 0) THEN
    !PRINT *, 'END OF TIME STEPPING'
    write(*,*) 'advancing cosmo finished'
  ENDIF


  cos_start = cos_start + cos_dt
end subroutine cosmo_advance

subroutine cosmo_finalize() bind(C,name="cosmo_finalize")

use iso_C_binding
use enkf_cosmo_mod

!------------------------------------------------------------------------------
!- Section 7: Final clean up
!------------------------------------------------------------------------------

  IF (izdebug > 0) THEN
    PRINT *, 'CLEAN UP'
  ENDIF

  CALL organize_allocation ('dealloc', ierrstat)
  IF (ldiagnos) THEN
    CALL organize_diagnosis ('dealloc', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,             &
                                     'organize_diagnosis: dealloc')
    ENDIF
  ENDIF

#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
  IF (luse_rttov) THEN
    CALL organize_satellites ('cleanup', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
       CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,             &
          'organize_satellites: cleanup')
    ENDIF
  ENDIF
#endif

#ifdef COSMOART
  IF (l_cosmo_art) THEN
    CALL organize_cosmo_art ('deallocate', ydate_ini, izerror, yzerrmsg)

    IF (izerror /= 0) THEN
      ierrstat = 1005
      yzerrmsg  = ' ERROR    *** Deallocation of space for COSMO_ART failed ***'
      CALL model_abort (my_cart_id, ierrstat, yzerrmsg,             &
                                    'organize_cosmo_art: deallocate')
    ENDIF
  ENDIF
#endif

#ifdef POLLEN
  IF (l_pollen) THEN
    CALL organize_pollen ('deallocate', ydate_ini, izerror, yzerrmsg)

    IF (izerror /= 0) THEN
      ierrstat = 1005
      yzerrmsg  = ' ERROR    *** Deallocation of space for Pollen failed ***'
      CALL model_abort (my_cart_id, ierrstat, yzerrmsg,             &
                                    'organize_pollen: deallocate')
    ENDIF
  ENDIF
#endif

  IF (ltime) THEN
    CALL get_timings (i_cleanup, nstop, dt, izerror)
    CALL collect_timings
  ENDIF

  IF (lasync_io .OR. (num_compute > 1) ) THEN
    CALL mpe_io_shutdown()
  ENDIF

!------------------------------------------------------------------------------
!- Section 8: Part of the IO-PEs
!------------------------------------------------------------------------------

!ELSE comp_pe
!
!   CALL mpe_io_node()
!
!ENDIF comp_pe

!------------------------------------------------------------------------------
!- Section 9: Final MPI-cleanup
!------------------------------------------------------------------------------

  CALL final_environment (ierrstat, yzerrmsg)

!------------------------------------------------------------------------------
!- End of the main program
!------------------------------------------------------------------------------

end subroutine cosmo_finalize
