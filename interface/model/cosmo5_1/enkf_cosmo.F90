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
  yroutine = 'lmorg'

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
! AF Save the value of l_cosmo_art in l_cosmo_art_nl
  IF (l_cosmo_art) THEN
    l_cosmo_art_nl=l_cosmo_art
  ENDIF

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

  IF (ltraj) THEN
    CALL organize_traj('input', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort(my_world_id, 100+izerror, yzerrmsg,'organize_traj: input')
    ENDIF
  ENDIF

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
  IF (leps .OR. lsppt) THEN
    CALL organize_eps ('input', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_world_id, 100+izerror, yzerrmsg,             &
                                      'organize_eps: input')
    ENDIF
  ENDIF

#ifdef RADARFWO
  ! Initialize parameters in data_radar.f90, which have to do with the model
  !  configuration (cartesian grid, PEs, MPI), with the respective parameters of the COSMO-model.
  ! This is done in any case, because it also concerns the "normal" DBZ gridpoint output
  !  defined in the GRIBOUT namelists, which is possible also in case of luse_radarfwo = .FALSE.
  ! This has to be done BEFORE the GRIBOUT-namelist(s) is/are read.
  CALL get_model_config_for_radar ()
#endif

  ! Input of the namelists for the I/O-package
  CALL organize_data ('input', 0, izerror, yzerrmsg)
  IF (izerror /= 0_iintegers) THEN
    CALL model_abort (my_world_id, 100+izerror, yzerrmsg,               &
                                   'organize_data: input')
  ENDIF

#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
  IF (luse_rttov) THEN
    ! Input of the namelists for the RTTOV-package
    CALL organize_satellites ('input', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_world_id, 100+izerror, yzerrmsg,             &
                                     'organize_satellites: input')
    ENDIF
  ENDIF
#endif

  ! Initialize I/O (must be called by all PEs)
  CALL organize_data ('init', 0, izerror, yzerrmsg)
  IF (izerror /= 0_iintegers) THEN
    CALL model_abort (my_world_id, 100+izerror, yzerrmsg,               &
                                   'organize_data: init')
  ENDIF

  ! Initialize tracer module
  CALL trcr_init( izerror )
  IF ( izerror /= 0_iintegers ) THEN
    yzerrmsg = trcr_errorstr( izerror )
    CALL model_abort( my_world_id, izerror, yzerrmsg, 'trcr_init' )
  ENDIF

  ! Initialize the timings
  CALL get_free_unit (nutiming)
  CALL init_timings (nstart, nstop, dt, itype_timing, ldfi,             &
                     lphys, luseobs, l2tls, lsemi_imp, l_cosmo_art,     &
                     l_2mom, luse_radarfwo, izerror)
  IF (izerror /= 0) THEN
    ! no system clock present
    ltime = .FALSE.
  ENDIF

#ifdef MESSY
  CALL messy_initialize
  CALL messy_new_tracer
#endif

!------------------------------------------------------------------------------
! Section 2: Allocation of space and computation of constant fields
!------------------------------------------------------------------------------

! Now comes the part for the compute PEs. The part for the IO-PEs is in
! the ELSE-part at the end of the program.

! comp_pe: IF (lcompute_pe) THEN

  !------------------------------------------------------------------------------
  ! 2.1:  Tracer definition
  !------------------------------------------------------------------------------

  ! Setup tracers for the physics (e.g. microphysics)
  CALL organize_physics( 'tracer', izerror, yzerrmsg )
  IF ( izerror /= 0_iintegers ) THEN
    CALL model_abort( my_world_id, 100+izerror, yzerrmsg,                     &
                                   'organize_physics: tracer' )
  ENDIF

  ! Setup other artificial tracer substances
  CALL gen_trcr_data( 'define', izerror, yzerrmsg )
  IF ( izerror /= 0_iintegers ) THEN
    CALL model_abort( my_world_id, 101+izerror, yzerrmsg,                     &
                                   'gen_trcr_data: define' )
  ENDIF

  !------------------------------------------------------------------------------
  ! 2.2:  Space allocation
  !------------------------------------------------------------------------------

  ! allocate space
  IF (izdebug > 0) THEN
    PRINT *,'    ALLOCATE SPACE'
  ENDIF

#ifdef NETCDF
  IF( lasync_io .AND. nc_asyn_io>0 ) THEN 
    CALL allocate_io_sendbuffer(yzerrmsg, izerror)
    IF (izerror /= 0) THEN
      ierrstat = 3321
      yzerrmsg = ' ERROR    *** Allocation of space for isend_buffer failed *** ' // yzerrmsg(1:100)
      CALL model_abort(my_cart_id, ierrstat, yzerrmsg,'allocate_io_sendbuffer' )
    ENDIF
  ENDIF
#endif

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

  ! tracer fields
  CALL trcr_alloc( izerror )
  IF ( izerror /= 0 ) THEN
    ierrstat = izerror
    yzerrmsg = trcr_errorstr( izerror )
    CALL model_abort( my_cart_id, ierrstat, yzerrmsg, 'allocation: tracers' )
  ENDIF

  ! block fields allocation
  IF (lphys) THEN
     CALL block_fields_allocate( izerror)
     IF ( izerror /= 0 ) THEN
        ierrstat = 1006
        yzerrmsg = 'block field allocation failed'
        CALL model_abort( my_cart_id, ierrstat, yzerrmsg, &
             'src_block_fields_org: block_fields_allocate' )
     ENDIF
  END IF

  !------------------------------------------------------------------------------
  ! 2.3:  Computation of constant fields
  !------------------------------------------------------------------------------

  CALL constant_fields

#ifdef RADARFWO

  IF (ltime) CALL get_timings (i_initializations, ntstep, dt, izerror)

  !------------------------------------------------------------------------------
  ! 2.4:  Initialization of radar forward operator on all compute PEs:
  !         - reading the namelist
  !         - reading radar meta informations from Radar data files
  !         - initializing lookup tables for Mie scattering if necessary
  !         - setting up auxiliary grids ("azimutal slices" for online propag.)
  !         - in case of not using the full 3D polar operator (luse_radarfwo=.false.),
  !           the COSMO model uses still the gridpoint reflectivity calculation
  !           of EMRADSCOPE. This needs an extra initialization "init_only_radar_gridpoint_calc"
  !------------------------------------------------------------------------------

  IF (luse_radarfwo) THEN
    CALL organize_radar ('init', nnew)
  END IF

  IF (ltime) CALL get_timings (i_radarsim, ntstep, dt, izerror)

#endif

!------------------------------------------------------------------------------
!- Section 3: Input of first data sets
!------------------------------------------------------------------------------

  ! Read or generate initial data and the first boundary data sets
  CALL organize_data ('start', 0, izerror, yzerrmsg)
  IF (izerror /= 0_iintegers) THEN
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,                &
                                   'start: input-init')
  ENDIF

#ifdef MESSY
  ! messy_init_memory is called inside organize_data
#endif

  !------------------------------------------------------------------------------
  ! 3.1:  Tracer summary and retrieval of metadata
  !------------------------------------------------------------------------------

  IF (izdebug > 2) THEN
    ! Print tracer list
    CALL trcr_print( izerror )
    IF ( izerror /= 0_iintegers ) THEN
      yzerrmsg = trcr_errorstr( izerror )
      CALL model_abort( my_world_id, izerror, yzerrmsg, 'trcr_print' )
    ENDIF
  ENDIF

  ! Retrieve the required metadata
  ALLOCATE (izlbc(trcr_get_ntrcr()), STAT=izerror)
  CALL trcr_meta_get(izerror, T_LBC_ID, izlbc)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yroutine)
  ENDIF

  ALLOCATE (izclp(trcr_get_ntrcr()), STAT=izerror)
  CALL trcr_meta_get(izerror, T_CLP_ID, izclp)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yroutine)
  ENDIF

  ALLOCATE (izbd_forced(trcr_get_ntrcr()), STAT=izerror)
  CALL trcr_meta_get(izerror, "BD_SET_FORCED", izbd_forced)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yroutine)
  ENDIF

!------------------------------------------------------------------------------
!- Section 4: Initializations and allocation of extra space
!------------------------------------------------------------------------------

  IF (izdebug > 0) THEN
    PRINT *, '  INITIALIZATIONS'
  ENDIF

  !------------------------------------------------------------------------------
  ! 4.1:  Initialization of different packages
  !------------------------------------------------------------------------------

  CALL organize_dynamics ('init', izerror, yzerrmsg, dt, .FALSE.)
  IF (izerror /= 0_iintegers) THEN
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,                &
                                   'organize_dynamics: init')
  ENDIF

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

    ! Initialization of the stochastics physics
    IF (lsppt) THEN
      CALL organize_eps ('init', izerror, yzerrmsg)
      IF ( izerror /= 0 ) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_eps: init')
      ENDIF
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

  IF (ltraj) THEN
    CALL organize_traj('init',izerror,yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort(my_world_id, 100+izerror, yzerrmsg,'organize_traj: init')
    ENDIF
  ENDIF

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

  IF (luse_rttov) THEN
#ifdef RTTOV10
    IF (lobsrad) THEN
      ! Read satpp files
      CALL organize_satellites('input-satpp',izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_world_id, 100+izerror, yzerrmsg,             &
                          'organize_satellites: input-satpp')
      ENDIF
    ENDIF
#endif

#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
    IF (lsynsat) THEN
      ! initialization of variables for the synthetic satellite computations
      CALL organize_satellites ('init', izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,              &
                                       'organize_satellites: init')
      ENDIF
    ENDIF
#endif
  ENDIF

  ! Initialization of blocks fields 
  ! This needs to be done after all arrays have been allocated
  IF (lphys) THEN
     ! Register block fields
     CALL block_fields_register_all( izerror )
     IF (izerror /= 0_iintegers) THEN
        yzerrmsg = 'Block field registration failed'
        CALL model_abort (my_cart_id, izerror, yzerrmsg,           &
             'src_block_fields_org: block_fields_register_all')
     ENDIF

     ! Initialize copy for each physics scheme
     CALL organize_physics ('init_copy', izerror, yzerrmsg)
     IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,              &
             'organize_physics: init_copy')
     ENDIF
  END IF

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
    !   switch off SPPT during DFI (as DFI calls organize_dynamics, "_physics)
    lzspptd = lsppt
    lsppt   = .FALSE.
    CALL dfi_initialization (lbd_frame, undef, izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,            &
                                     'dfi_initialization')
    ENDIF
    lsppt   = lzspptd
  ENDIF
  
  ! Close file for control output
  IF (my_cart_id == 0) THEN
    CLOSE (nuspecif, STATUS='KEEP')
  ENDIF

  ! halve the time step, if ntstep = 0
  ! (in src_setup: nstep = nstart)
  IF (.NOT. l2tls) THEN
    IF (ntstep == 0) THEN
      zforecasttime = 0.0_wp
      dt = 0.5_wp * dt
    ELSE
      zforecasttime = 0.5_wp*dt + (ntstep-1)*dt
      nzdiv         = INT (zforecasttime / 3600.0_wp, iintegers)
      zforecasttime = zforecasttime - nzdiv * 3600.0_wp
    ENDIF
  ELSE
    IF (ntstep == 0) THEN
      zforecasttime = 0.0_wp
    ELSE
      zforecasttime = ntstep*dt
      nzdiv         = INT (zforecasttime / 3600.0_wp, iintegers)
      zforecasttime = zforecasttime - nzdiv * 3600.0_wp
    ENDIF
  ENDIF

  IF (lbdclim) THEN
    ynote         = '...... FORECAST TIME IS NOW xxxxxx DAYS ......'
  ELSE
    ynote         = '...... FORECAST TIME IS NOW xxx HOURS   ......'
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

    ! Copy the modified fields T and QV into timelevel nnow for leapfrog integration:
    IF (.NOT. l2tls) THEN
      CALL trcr_get(izerror, 'QV', ptr_tlev = nnew, ptr = qv_new)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, 'lmorg: T/QV disturbance(s) at model start')
      ENDIF
      IF (izerror /= 0) THEN
      CALL trcr_get(izerror, 'QV', ptr_tlev = nnow, ptr = qv_now)
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, 'lmorg: T/QV disturbance(s) at model start')
      ENDIF
      t(:,:,:,nnow) = t(:,:,:,nnew)
      qv_now(:,:,:) = qv_new(:,:,:)
    END IF

    ! Initial condition on t_so (takes only effect if lsoil=.true.)
    CALL set_tempdist_tso(nnew)

    ! Copy soil temperatures into timelevel nnow for leapfrog integration:
    IF (.NOT. l2tls) THEN
      t_s   (:,:,nnow) = t_s (:,: ,nnew)
      t_g   (:,:,nnow) = t_g (:,: ,nnew)
      IF (lmulti_layer) THEN
        t_so(:,:,:,nnow) = t_so(:,:,:,nnew)
      ELSE
        t_m(:,:,nnow)  = t_m(:,:,nnew)
      ENDIF
    END IF
  ENDIF

#if defined COUP_OAS
  CALL oas_cos_define
!OASIS4 only
!  CALL oas_cos_update_time(0)
#endif

#ifdef MESSY
 CALL messy_init_coupling
 IF (nstart /= 0) CALL messy_channel_read_restart
 CALL messy_init_tracer
#endif

  IF (ltime) CALL get_timings (i_initializations, ntstep, dt, izerror)

  !------------------------------------------------------------------------------
  ! 5b. For periodic BCs we may need an exchange here, otherwise some fields
  !     might not be periodic in the first timestep:
  !------------------------------------------------------------------------------

  IF (lperi_x .OR. lperi_y) THEN

    ! Check, whether additional communication for the convection is
    ! necessary
    lzconv = lconv .AND. lconf_avg .AND.                         &
         ((ntstep < 1) .OR. (MOD(ntstep+2,nincconv)==0))

    IF (ltime_barrier) THEN
      CALL comm_barrier (icomm_cart, ierrstat, yzerrmsg)
      IF (ltime) CALL get_timings (i_barrier_waiting_dyn, ntstep, dt, izerror)
    ENDIF

    IF     ( l2tls ) THEN
      CALL exchange_runge_kutta
    ELSE ! Leapfrog:
      CALL exchange_leapfrog
    ENDIF

    IF (ltime) CALL get_timings (i_communications_dyn, ntstep, dt, izerror)

  ENDIF

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

  ! timeloop: DO ntstep = nstart , nstop
  write(*,*)'advancing cosmo from ',cos_start,' to ',(cos_start+cos_dt-1)

  timeloop: DO ntstep = cos_start,(cos_start+cos_dt-1)
    IF ( (izdebug > 1) .AND. (.NOT. lbdclim)) THEN
      PRINT *, '    STEP ',ntstep
    ENDIF

#ifdef COSMOART
! AF Put l_cosmo_art to false if you don't want to calculate the ART part
    IF (l_cosmo_art_nl) THEN
      l_cosmo_art=(l_cosmo_art_nl .AND. ((ntstep*dt/3600.0_wp) >= artstart ))
    ENDIF
#endif

    !--------------------------------------------------------------------------
    !- Section 6.1: Initialization of this time step
    !--------------------------------------------------------------------------

    ! Set nexch_tag dependend on the time step
    nexch_tag = MOD (ntstep, INT(24.0_wp*3600.0_wp/dt))

    IF (l2tls) THEN
      nnow = 3 - nnow
      nnew = 3 - nnew
    ELSE
      nsp    = nold
      nold   = nnow
      nnow   = nnew
      nnew   = nsp
    ENDIF

    CALL initialize_loop (ntstep, nbd1, nbd2, nold, nnow, nnew)

    IF (ltime) CALL get_timings (i_add_computations, ntstep, dt, izerror)

    !--------------------------------------------------------------------------
    !- Section 6.2.1: physics 
    !--------------------------------------------------------------------------

#ifdef COSMOART
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
      CALL organize_pollen ('prepare_transport', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_pollen: prepare_transport')
      ENDIF
    ENDIF
#endif

#ifdef MESSY
    CALL messy_global_start
    CALL messy_local_start
    CALL messy_vdiff
#endif

    IF (lsppt .AND. lphys) THEN
      CALL organize_eps ('compute', izerror, yzerrmsg)
      IF ( izerror /= 0 ) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                         'organize_eps: compute')
      END IF
    END IF

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
    IF (l_cosmo_art) THEN
      CALL organize_cosmo_art ('finalize_physics', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_cosmo_art: finalize_physics')
      ENDIF
    ENDIF
#endif

!#ifdef POLLEN
!    ! Clean up after Pollen
!    IF (l_pollen) THEN
!      CALL organize_pollen ('finalize_physics', ydate_ini, izerror, yzerrmsg)
!      IF (izerror /= 0_iintegers) THEN
!        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
!                                       'organize_pollen: finalize_physics')
!      ENDIF
!    ENDIF
!#endif

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

!#ifdef POLLEN
!    ! Preparations for Pollen
!    IF (l_pollen) THEN
!      CALL organize_pollen ('prepare_dynamics', ydate_ini, izerror, yzerrmsg)
!      IF (izerror /= 0_iintegers) THEN
!        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
!                                       'organize_pollen: prepare_dynamics')
!      ENDIF
!    ENDIF
!#endif

#ifdef MESSY
    ! must be called here and not later than organize_dynamics
    ! (otherwise tendencies must be applied "by-hand"
    CALL messy_physc
    CALL messy_local_end
    CALL messy_global_end(1)

! moved into organize dynamics
!    CALL main_tracer_beforeadv
! NOTE: IN ORGANIZE DYNAMICS THE "FINAL INTEGRATION" takes place !!!!
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
      CALL organize_pollen ('washout', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_pollen: washout')
      ENDIF
    ENDIF
#endif

    CALL set_trcr_special_bc

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

    !--------------------------------------------------------------------------
    !- Section 6.4: nudging
    !--------------------------------------------------------------------------

#ifdef NUDGING
    IF (luseobs) CALL organize_assimilation ('nudge', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                    'organize_assimilation: nudge')
    ENDIF

    IF (ltime) CALL get_timings (i_nud_computations, ntstep, dt, izerror)
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

      IF (ltime) CALL get_timings (i_add_computations, ntstep, dt, izerror)
    ENDIF

    !--------------------------------------------------------------------------
    !- Section 6.6: spectral nudging and relaxation
    !--------------------------------------------------------------------------

    IF (lspecnudge .AND. ((ntstep < 2) .OR. (MOD(ntstep+1,nincsn) == 0))) THEN
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

!#ifdef POLLEN
!    ! Preparations for Pollen
!    IF (l_pollen) THEN
!      CALL organize_pollen ('prepare_relaxation', ydate_ini, izerror, yzerrmsg)
!      IF (izerror /= 0_iintegers) THEN
!        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
!                                       'organize_pollen: prepare_relaxation')
!      ENDIF
!    ENDIF
!#endif

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
      CALL organize_pollen ('finalize_transport', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_pollen: finalize_transport')
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

      IF     ( l2tls ) THEN
        CALL exchange_runge_kutta
      ELSE ! Leapfrog:
        CALL exchange_leapfrog
      ENDIF

      IF (ltime) CALL get_timings (i_communications_dyn, ntstep, dt, izerror)

    !--------------------------------------------------------------------------
    !- Section 6.8: diagnostics
    !--------------------------------------------------------------------------

    CALL near_surface (nnow)
    IF (ltime) CALL get_timings (i_add_computations, ntstep, dt, izerror)

    !   Analysis of near surface parameters
    !   -----------------------------------

#ifdef NUDGING
    IF (luseobs) THEN
      CALL organize_assimilation ('surface', izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                      'organize_assimilation: surface')
      ENDIF
      IF (ltime) CALL get_timings (i_nud_computations, ntstep, dt, izerror)
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

    !   Trajectories computation
    !   ------------------------

    IF (ltraj) THEN
      CALL organize_traj('compute',izerror,yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort(my_world_id, 100+izerror, yzerrmsg,'organize_traj: compute')
      ENDIF
    ENDIF

    ! put a barrier here to have a clean separation for the timing
    CALL comm_barrier (icomm_cart, ierrstat, yzerrmsg)
    IF (ltime) CALL get_timings (i_add_computations, ntstep, dt, izerror)

#ifdef RADARFWO
    !--------------------------------------------------------------------------
    !- Section 6.8b: radar forward operator and optionally preparing 
    !                of radar feedback files for radar data assimilation
    !                (the latter needs radar observation files)
    !--------------------------------------------------------------------------

    IF (luse_radarfwo) THEN

      CALL organize_radar ('compute', nnow)

      IF (ltime) CALL get_timings (i_radarsim, ntstep, dt, izerror)

    END IF
#endif

    !--------------------------------------------------------------------------
    !- Section 6.9: output of results
    !--------------------------------------------------------------------------

#ifdef MESSY
    CALL messy_global_end(2)
    CALL messy_write_output(1)
#endif

#ifdef COSMOART
! AF Put l_cosmo_art back to initial value for output
    IF (l_cosmo_art_nl) THEN
      l_cosmo_art=l_cosmo_art_nl
    ENDIF
#endif 

    CALL organize_data ('result', ntstep, izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                     'result: input-init')
    ENDIF

#ifdef MESSY
    IF ( (l2tls) .OR.  (.NOT. lstop .AND. .NOT. lbreak) ) &
         CALL messy_write_output(2)
#endif

#ifdef COSMOART
! AF Put l_cosmo_art to false if you don't want to calculate the ART part
    IF (l_cosmo_art_nl) THEN
      l_cosmo_art=(l_cosmo_art_nl .AND. ((ntstep*dt/3600.0_wp) >= artstart ))
    ENDIF

    ! more universal approach: get a general injection point for ART
    IF (l_cosmo_art) THEN
      CALL organize_cosmo_art ('endoftimestep', ydate_ini, izerror, yzerrmsg)
      IF (izerror /= 0_iintegers) THEN
        CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,               &
                                       'organize_cosmo_art: endoftimestep')
      ENDIF
    ENDIF
#endif

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
        IF (zforecasttime >= 86400.0_wp) THEN
          nzdays  = NINT ((ntstep+1)*dt) / 86400
          WRITE (ynote(29:34),'(I6.6)') nzdays
          PRINT *, ynote
          zforecasttime = zforecasttime - 86400.0_wp
        ENDIF
      ENDIF
    ELSE
      ! record a forecast hour
      IF (my_cart_id ==  0) THEN
        zforecasttime = zforecasttime + dt
        IF (zforecasttime >= 3600.0_wp) THEN
          nzhours = NINT ((ntstep+1)*dt) / 3600
          WRITE (ynote(29:31),'(I3.3)') nzhours
          PRINT *, ynote
          zforecasttime = zforecasttime - 3600.0_wp
        ENDIF
      ENDIF
    ENDIF

    ! Reset the time step for leapfrog integration
    IF ( ntstep == 0 .AND. (.NOT.l2tls) ) THEN
      dt = 2.0_wp * dt
    ENDIF

#if defined COUP_OAS
! OASIS4 only
!    CALL oas_cos_update_time(ntstep+1)
#endif

#ifdef MESSY
    IF (lbreak .OR. lstop) EXIT
    CALL messy_timer_reset_time
#endif

  ENDDO timeloop

#ifdef COSMOART
! AF Put l_cosmo_art back to initial value for final clean-up
  IF (l_cosmo_art_nl) THEN
    l_cosmo_art=l_cosmo_art_nl
  ENDIF
#endif

  IF (izdebug > 0) THEN
    ! PRINT *, 'END OF TIME STEPPING'
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

#ifdef MESSY
  CALL messy_free_memory
#endif

  CALL organize_allocation ('dealloc', ierrstat)

  CALL organize_dynamics   ('cleanup', izerror, yzerrmsg, dt, .FALSE.)

  IF (ldiagnos) THEN
    CALL organize_diagnosis ('dealloc', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,             &
                                     'organize_diagnosis: dealloc')
    ENDIF
  ENDIF

  ! Deallocate tracers and metadata
  DEALLOCATE ( izlbc    )
  DEALLOCATE ( izclp    )
  DEALLOCATE ( izbd_forced )

  CALL trcr_cleanup(izerror)
  IF ( izerror /= 0_iintegers ) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, 'trcr_cleanup:')
  ENDIF

  ! Block physics cleanup
  CALL block_fields_cleanup(izerror)

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

  IF (ltraj) THEN
    CALL organize_traj('finalize', izerror, yzerrmsg)
    IF (izerror /= 0_iintegers) THEN
      CALL model_abort(my_world_id, 100+izerror, yzerrmsg,'organize_traj: finalize')
    ENDIF
  ENDIF

  IF (ltime) THEN
    CALL get_timings (i_cleanup, nstop, dt, izerror)
  ENDIF

#ifdef NETCDF
  IF ( nc_asyn_io > 0 .AND. lasync_io ) THEN
    CALL shutdown_netcdfio_sendbuffers()
    IF (ltime) THEN
      CALL get_timings (i_asynio_wait, ntstep, dt, izerror)
    ENDIF
  ENDIF
#endif

  IF (ltime) THEN
    CALL collect_timings
  ENDIF

  IF ( (lasync_io .OR. (num_compute > 1)) .AND. (nc_asyn_io < 1) ) THEN
    CALL mpe_io_shutdown(izerror)
  ENDIF

!------------------------------------------------------------------------------
!- Section 8: Part of the IO-PEs
!------------------------------------------------------------------------------

! ELSE comp_pe

!   IF( nc_asyn_io > 0 ) THEN
! #ifdef NETCDF
!     CALL start_ionode( yzerrmsg, izerror)
!     IF( izerror /= 0 ) THEN
!       CALL model_abort(my_cart_id, izerror, yzerrmsg,  &
!               'start_ionode')
!     ENDIF
!     CALL shutdown_io()
! #endif
!   ELSE
!    CALL mpe_io_node(izerror)
!    IF (izerror /= 0) THEN
!      ierrstat = 1015
!      yzerrmsg  = ' ERROR    *** Running the asynchronous I/O failed ***'
!      CALL model_abort (my_cart_id, ierrstat, yzerrmsg, 'start mpe_io_node')
!    ENDIF
!  ENDIF

! ENDIF comp_pe

!------------------------------------------------------------------------------
!- Section 9: Final MPI-cleanup
!------------------------------------------------------------------------------

#ifdef MESSY
   IF (lstop) THEN
    ! WRITE file 'END' to break rerun chain
    CALL messy_blather_endfile_bi('Simulation finished.', ' ')
  ELSE
    CALL info_bi('Simulation stopped.', ' ')
    ! Notes:
    !    - simulation is stopped (lbreak) and a rerun is started
    !      (continue rerun chain)
  END IF
#endif

  CALL final_environment (ierrstat, yzerrmsg)

!------------------------------------------------------------------------------
!- End of the main program
!------------------------------------------------------------------------------

end subroutine cosmo_finalize
