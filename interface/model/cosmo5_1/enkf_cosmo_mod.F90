module enkf_cosmo_mod
USE info_lm_f90,         ONLY: info_define, info_readnl, info_print

USE data_parameters,     ONLY:   wp, iintegers
USE data_constants,      ONLY:   b1, b2w, b3, b4w, b2i, b4i, rdv, o_m_rdv,   &
                                 rvd_m_o, r_d
USE data_soil,           ONLY:   cf_snow
USE data_fields,         ONLY:                                               &
       dp0, p0, rho, rho0, qrs, llandmask, t_g, vmax_10m, dqvdt, qvsflx,     &
       u, v, w, t, tke, pp, ps, u_bd, v_bd, w_bd, t_bd, pp_bd,               &
       t_snow,    t_s,    qv_s,    t_m,    w_snow,    w_g1,    w_g2,         &
       t_snow_bd, t_s_bd, qv_s_bd, t_m_bd, w_snow_bd, w_g1_bd, w_g2_bd,      &
       plcov_bd, lai_bd, rootdp_bd, vio3_bd, hmo3_bd, t_cl_bd, w_cl_bd,      &
       plcov,    lai,    rootdp,    vio3,    hmo3,    t_cl,    w_cl,         &
       utens, vtens, wtens, ttens, pptens, w_g3, w_g3_bd,                    &
       fr_lake, depth_lk, h_ice, vgust_con, vgust_dyn, vabsmx_10m,           &
       prr_gsp, prr_con, prne_con, bas_con, t_snow_mult, t_ice,              &
       aer_su   , aer_du   , aer_bc   , aer_or   , aer_ss   ,                &
       aer_su_bd, aer_du_bd, aer_bc_bd, aer_or_bd, aer_ss_bd, t_so


USE data_modelconfig,    ONLY:   ie, je, ke, ke1, jstartpar, jendpar,        &
                                 istartpar, iendpar,                         &
                                 dt, dt2, ed2dt, dtdeh, nehddt,              &
                                 jstartu, jendu, jstartv, jendv,             &
                                 istart, iend, jstart, jend, ie_tot, je_tot, &
                                 idt_qv, idt_qc, idt_qi, idt_qni,            &
                                 lalloc_t_cl, idt_qr, idt_qs, idt_qg, idt_qh

USE data_runcontrol,     ONLY:                                               &
       nstart, nstop, ntstep, nlastbound, nincbound, lmulti_layer,           &
       ltime, itype_timing, newbcdt, newbc, nincboufac, nlgw, itype_gscp,    &
       nbd1, nbd2, nold, nnow, nnew, ntke, lartif_data, ldiagnos, luse_rttov,&
       lphys, lconv, luseobs, l2tls, lsemi_imp, lgsp, lsppt,                 &
       yakdat1, yakdat2, nuspecif, ldfi, nincconv, lconf_avg, l2dim,         &
       nbl_exchg, lprog_tke, lspecnudge, lseaice, llake,                     &
       lw_freeslip, leps, idbg_level, lprintdeb_all, itype_calendar,         &
       hlastmxu, hnextmxu, hincmxu, nlastmxu, nnextmxu, l_cosmo_art,         &
       l_pollen, hstart, itype_lbc_qrsg, lmulti_snow, lperi_x, lperi_y,      &
       l2dim, nincsn, itype_aerosol, l_2mom, itype_turb, luse_radarfwo,      &
       ltraj

USE data_parallel,       ONLY:                                               &
       num_compute, nc_asyn_io, icomm_cart, my_cart_id, sendbuf, isendbuflen,&
       iexch_req, imp_reals, nboundlines, my_cart_neigh, my_world_id,        &
       lcompute_pe, lasync_io, ncomm_type, ldatatypes, ltime_barrier,        &
       nexch_tag, cosmo_input_suffix

USE data_io,             ONLY:   ydate_ini, lbdclim, lbdsst, lbd_frame, undef

USE data_flake,          ONLY:   h_Ice_min_flk

#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
USE data_satellites,     ONLY:   lsynsat, lobsrad
#endif

USE mpe_io2,             ONLY:   mpe_io_node, mpe_io_shutdown

#ifdef NETCDF
USE netcdf_io,       ONLY:   start_ionode, shutdown_io, &
    shutdown_netcdfio_sendbuffers, allocate_io_sendbuffer
#endif

USE environment,         ONLY:   exchg_boundaries, comm_barrier,             &
                                 final_environment, model_abort, get_free_unit
USE meteo_utilities,     ONLY:   calrho, calps, tgcom
USE time_utilities,      ONLY:   nutiming, init_timings, get_timings,        &
                  i_initializations, i_add_computations, i_phy_computations, &
                  i_dyn_computations,i_lhn_computations, i_spectr_nudging,   &
                  i_relaxation, i_input, i_output, i_barrier_waiting_dyn,    &
                  i_communications_dyn, i_cleanup, i_asynio_wait,            &
                  i_radarsim, i_nud_computations,                            &
                  collect_timings
USE utilities,           ONLY:   get_utc_date, check_field_NaNs

USE src_setup,           ONLY:   organize_setup, constant_fields

USE src_allocation,      ONLY:   organize_allocation

USE src_artifdata,       ONLY:   artif_heatrate_dist, artif_heatrate_dist_tso, set_tempdist, &
                                 set_tempdist_tso, set_tempdist_bbc_ts, add_noise_tw,        &
                                 gen_trcr_data

USE data_tracer,         ONLY:   T_CLP_ID, T_CLP_ON, T_INI_ID, T_INI_FILE,    &
                                 T_LBC_ID, T_LBC_FILE, T_LBC_CST, T_LBC_ZERO, &
                                 T_LBC_ZEROGRAD, T_LBC_USER, T_ERR_NOTFOUND

USE src_tracer,          ONLY:   trcr_init, trcr_errorstr, trcr_alloc,        &
                                 trcr_print, trcr_meta_get, trcr_cleanup,     &
                                 trcr_get, trcr_get_ntrcr, trcr_get_block,    &
                                 trcr_swap

USE src_traj,            ONLY:   organize_traj

USE src_block_fields_org, ONLY: block_fields_allocate,                       &
                                block_fields_deallocate,                     &
                                block_fields_register_all,                   &
                                block_fields_cleanup

#ifdef COSMOART
USE data_cosmo_art,      ONLY:   lgas, laero, ldust, lseas,                  &
                                 lrad_dust, lrad_aero,                       &
                                 lgasbd, laerobd, lgasini, laeroini,         &
                                 isp_gas, isp_aero, isp_aerotrans,           &
                                 nlastbound_art, nincbound_art,              &
                                 artstart, l_cosmo_art_nl,                   &
                                 nbd1_art, nbd2_art, trcr_idx_gas,           &
                                 trcr_idx_aero
USE art_species_data,    ONLY:   lho, lho2
USE art_aerosol_const,   ONLY:   aerostart
#endif

#ifdef POLLEN
USE data_pollen,         ONLY:   isp_pollen, dtpollen, trcr_idx_pollen
#endif

#ifdef MESSY
! MESSy/BMIL
USE messy_main_channel_bi, ONLY: messy_channel_read_restart
USE messy_main_timer_bi,   ONLY: messy_timer_reset_time
USE messy_main_blather_bi, ONLY: error_bi, info_bi, messy_blather_endfile_bi
USE messy_main_tracer_bi,  ONLY: main_tracer_beforeadv, main_tracer_afteradv
USE messy_main_data_bi,    ONLY: L_IS_CLIENT

! MESSy/SMCL
USE messy_main_timer,      ONLY: lstop, lbreak
#endif

#ifdef TWOMOM_SB
USE data_fields,             ONLY:   reffc_out, reffi_out,                    &
                                     odepthw_so, odepthi_so,                  &
                                     odepthw_th, odepthi_th
USE src_twomom_sb_interface, ONLY:   set_qni_from_qi_sb
#endif

#ifdef RADARFWO
USE src_radar,           ONLY:   organize_radar
USE radar_cosmo,         ONLY:   get_model_config_for_radar
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Local Scalars:

INTEGER  (KIND=iintegers) ::                          &
  ierrstat,      & ! error status variable
  izerror,       & !
  izlbcqi,       & ! type of lateral BC for QI
  izanaqi,       & ! type of IC for QI
  nzhours,       & ! for recording a forecast hour
  nzdays,        & ! for recording a forecast day
#ifdef COSMOART
  isp,           & !
#endif
  i,j,k, kzdims(24), izdebug, nzdiv, iztrcr, nsp

REAL (KIND=wp)            ::                          &
  zforecasttime    ! for recording a forecast hour

LOGICAL                   ::                          &
  lzconv,        & ! to determine whether extra communication for convection
                   ! is needed
  lzspptd          ! dummy switch, to exclude SPPT during digital filter init.

CHARACTER (LEN= 46)       :: ynote
CHARACTER (LEN=255)       :: yzerrmsg
CHARACTER (LEN= 25)       :: yroutine

! Tracer pointers:
REAL (KIND=wp),     POINTER ::                        &
  ztrcr      (:,:,:)  => NULL(), & ! tracer tendency variable
  ztrcr_now  (:,:,:)  => NULL(), & ! tracer variable at nnow
  ztrcr_bd   (:,:,:,:)=> NULL(), & ! tracer boundaries variable
  qv_new     (:,:,:)  => NULL(), & ! QV at nnew
  qv_now     (:,:,:)  => NULL(), & ! QV at nnow
  qc_new     (:,:,:)  => NULL(), & ! QC at nnew
  qc_now     (:,:,:)  => NULL(), & ! QC at nnow
  qc_nx0     (:,:,:)  => NULL(), & ! QC at nx0
  qi_new     (:,:,:)  => NULL(), & ! QI at nnew
  qi_now     (:,:,:)  => NULL(), & ! QI at nnow
  qi_nx0     (:,:,:)  => NULL(), & ! QI at nx0
  qr_now     (:,:,:)  => NULL(), & ! QR at nnow
  qs_now     (:,:,:)  => NULL(), & ! QS at nnow
  qg_now     (:,:,:)  => NULL()    ! QG at nnow

#ifdef TWOMOM_SB
REAL (KIND=wp),     POINTER ::                        &
  qni_now    (:,:,:)  => NULL(), & ! NCICE at nnow
  qni_nx0    (:,:,:)  => NULL()    ! NCICE at nx0
#endif

INTEGER (KIND=iintegers), ALLOCATABLE:: &
  izclp   (:) , & ! clipping type for all tracers
  izlbc   (:) , & ! array containing the lateral BC
                  ! type for all tracers
  izbd_forced(:)  ! BD_SET_FORCED for all tracers

!kuw: time step control for cosmo
integer :: cos_start

!==============================================================================

!==============================================================================
! Internal procedures in lmorg
!==============================================================================

CONTAINS

!==============================================================================
!+ Internal procedure in "lmorg" for initializing each time step
!------------------------------------------------------------------------------

!option! -pvctl ifopt
SUBROUTINE initialize_loop (ntstep, nbd1, nbd2, nold, nnow, nnew)

!------------------------------------------------------------------------------
!
! Description:
!   This routine initializes each time step. It checks whether certain 
!   actions have to be performed and sets the logical variables from
!   the parameterlist. Organizational variables are updated.
!
! Method:
!   - input of new boundary values, if necessary
!   - check whether diagnostic computations have to be performed
!   - check whether outputs have to be done
!   - update organiziational variables
!   - initialize the new time level with boundary values
!   - calculate the density of moist air (rho) for this time level
!
!==============================================================================
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER (KIND=iintegers), INTENT (IN)  ::       &
  ntstep             ! actual time step

! Scalar arguments with intent(inout):
INTEGER (KIND=iintegers), INTENT (IN)  ::       &
  nbd1, nbd2,      & ! indices for the boundary levels
  nold, nnow, nnew   ! indices for the time levels

!------------------------------------------------------------------------------
! Local variables

INTEGER (KIND=iintegers)   ::  &
  nx0, i, j, k,                & ! support variable
  nzjulianday                    ! day of the year

REAL (KIND=wp)             ::  &
  z1, z2,                      & ! factors for time interpolation
  fpvsw, fpvsi, fqvs, zst,     & ! for the statement functions
  zsge, zsp, zacthour            ! actual hour of the forecast

#ifdef COSMOART
! CK 20101204
! ART bd update frequency does not need to be the same as meteo.
! therefore, the weights calculated for meteo might be wrong.

REAL (KIND=wp)             ::        &
  z1_art, z2_art                       ! factors for time interpolation
#endif

REAL (KIND=wp)             ::        &
  zt_s(ie,je)                    ! = t_s   on land and sea
                                 ! = t_ice on sea ice (if present)

#ifdef COSMOART
REAL (KIND=wp),     POINTER ::                        &
  cgas_bd  (:,:,:,:,:)  => NULL(),                    &! cgas_bd
  cgas_now   (:,:,:,:)  => NULL(),                    &! cgas at nnow
  cgas_new   (:,:,:,:)  => NULL(),                    &! cgas at nnew
  caero_bd  (:,:,:,:,:)  => NULL(),                   &! caero_bd
  caero_now   (:,:,:,:)  => NULL(),                   &! caero at nnow
  caero_new   (:,:,:,:)  => NULL()                     ! caero at nnew
#endif

#ifdef POLLEN
REAL (KIND=wp),     POINTER ::                        &
  cpollen_new   (:,:,:,:)  => NULL(),                 &  ! cpollen at nnew
  cpollen_now   (:,:,:,:)  => NULL(),                 &  ! cpollen at nnow
  cpollen_bd    (:,:,:,:,:)  => NULL()                     ! cpollen_bd
#endif

CHARACTER (LEN=25)         :: yzroutine

LOGICAL :: lfound

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE initialize_loop
!------------------------------------------------------------------------------

! Statement functions
!--------------------
  fpvsw(zst)       = b1 * EXP( b2w * (zst-b3)/(zst-b4w) )
  fpvsi(zst)       = b1 * EXP( b2i * (zst-b3)/(zst-b4i) )
  fqvs (zsge, zsp) = rdv * zsge / ( zsp - o_m_rdv * zsge )

! define routine name
  yzroutine = 'initialize_loop'

! get new actual utc date
  CALL get_utc_date(ntstep, ydate_ini, dt, itype_calendar, yakdat1, yakdat2, &
                    nzjulianday, zacthour)

!------------------------------------------------------------------------------
!- Section 1: input of new boundary values, if necessary
!------------------------------------------------------------------------------

  IF (ltime) CALL get_timings (i_add_computations, ntstep, dt, izerror)

  ! get new boundary data file or generate artificial data
  CALL organize_data ('boundary', ntstep, izerror, yzerrmsg)
  IF (izerror /= 0_iintegers) THEN
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,                &
                                   'boundary: input-init')
  ENDIF

  IF (ltime) CALL get_timings (i_input, ntstep, dt, izerror)

#ifdef MESSY
  ! get new boundary data from MMDCLNT of required
  CALL messy_init_loop
#endif

!------------------------------------------------------------------------------
!- Section 2: update organizational variables
!------------------------------------------------------------------------------

  ! cyclic changing of the indices for the time levels
  IF (l2tls) THEN
    nx0  = nnow
  ELSE
    nx0    = nold
  ENDIF

  ! variables concerned with the time step
  dt2    = 2.0_wp * dt
  ed2dt  = 1.0_wp / dt2
  dtdeh  = dt / 3600.0_wp
  nehddt = NINT ( 3600.0_wp / dt )

  ! swap timelevels in case tracers use static fields
  CALL trcr_swap(izerror)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  ! retrieve the required microphysics tracers
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnew, ptr = qv_new)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnow, ptr = qv_now)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnew, ptr = qc_new)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnow, ptr = qc_now)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nx0, ptr = qc_nx0)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qi, ptr_tlev = nnew, ptr = qi_new)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qi, ptr_tlev = nnow, ptr = qi_now)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qi, ptr_tlev = nx0, ptr = qi_nx0)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

#ifdef TWOMOM_SB
  CALL trcr_get(izerror, idt_qni, ptr_tlev = nnow, ptr = qni_now)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qni, ptr_tlev = nx0,  ptr = qni_nx0)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
#endif 

!------------------------------------------------------------------------------
!- Section 3: Initialize new time level with boundary values
!------------------------------------------------------------------------------

  ! factors for linear time interpolation
  z2 = REAL (ntstep+1-nlastbound, wp) / REAL (nincbound, wp)
  z2 = MIN ( 1.0_wp , z2 )
  z1 = 1.0_wp - z2

  ! fields of the atmosphere
  IF (lbd_frame) THEN
!CDIR COLLAPSE
    WHERE (t_bd (:,:,:,nbd2) /= undef)
!CDIR COLLAPSE
      u (:,:,:,nnew) = z1 * u_bd (:,:,:,nbd1) + z2 * u_bd (:,:,:,nbd2)
!CDIR COLLAPSE
      v (:,:,:,nnew) = z1 * v_bd (:,:,:,nbd1) + z2 * v_bd (:,:,:,nbd2)
!CDIR COLLAPSE
      t (:,:,:,nnew) = z1 * t_bd (:,:,:,nbd1) + z2 * t_bd (:,:,:,nbd2)
!CDIR COLLAPSE
      pp(:,:,:,nnew) = z1 * pp_bd(:,:,:,nbd1) + z2 * pp_bd(:,:,:,nbd2)
    ELSEWHERE
!CDIR COLLAPSE
      u (:,:,:,nnew) = u (:,:,:,nnow)
!CDIR COLLAPSE
      v (:,:,:,nnew) = v (:,:,:,nnow)
!CDIR COLLAPSE
      t (:,:,:,nnew) = t (:,:,:,nnow)
!CDIR COLLAPSE
      pp(:,:,:,nnew) = pp(:,:,:,nnow)
    ENDWHERE
    IF (.NOT. lw_freeslip) THEN
!CDIR COLLAPSE
      WHERE (t_bd (:,:,:,nbd2) /= undef)
!CDIR COLLAPSE
        w(:,:,:,nnew) = z1 *  w_bd(:,:,:,nbd1) + z2 *  w_bd(:,:,:,nbd2)
      ELSEWHERE
!CDIR COLLAPSE
        w(:,:,:,nnew) =  w(:,:,:,nnow)
      ENDWHERE
    ENDIF
  ELSE
!CDIR COLLAPSE
    u (:,:,:,nnew) = z1 * u_bd (:,:,:,nbd1) + z2 * u_bd (:,:,:,nbd2)
!CDIR COLLAPSE
    v (:,:,:,nnew) = z1 * v_bd (:,:,:,nbd1) + z2 * v_bd (:,:,:,nbd2)
!CDIR COLLAPSE
    t (:,:,:,nnew) = z1 * t_bd (:,:,:,nbd1) + z2 * t_bd (:,:,:,nbd2)
!CDIR COLLAPSE
    pp(:,:,:,nnew) = z1 * pp_bd(:,:,:,nbd1) + z2 * pp_bd(:,:,:,nbd2)
    IF (.NOT. lw_freeslip) THEN
!CDIR COLLAPSE
      w (:,:,:,nnew) = z1 * w_bd (:,:,:,nbd1) + z2 * w_bd (:,:,:,nbd2)
    END IF
  ENDIF

  ! Tracers
  !---------
  ! loop over tracers
  DO iztrcr = 1, trcr_get_ntrcr()

    ! get pointer to tracer (at nnew)
    CALL trcr_get(izerror, iztrcr, ptr_tlev=nnew, ptr=ztrcr,              &
                  ptr_bd=ztrcr_bd)
    IF (izerror /= 0_iintegers) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF

    ! get pointer to tracer (at nnow)
    CALL trcr_get(izerror, iztrcr, ptr_tlev=nnow, ptr=ztrcr_now)
    IF (izerror /= 0_iintegers) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF

    IF ( ANY(izlbc(iztrcr) == (/T_LBC_FILE, T_LBC_USER/) ) ) THEN

      ! if boundary from file, time interpolation between 2 boundary values
      IF ( lbd_frame ) THEN
        WHERE (ztrcr_bd (:,:,:,nbd2) /= undef)
          ztrcr(:,:,:) = z1 * ztrcr_bd(:,:,:,nbd1) +                    &
                         z2 * ztrcr_bd(:,:,:,nbd2)    
        ELSEWHERE
          ztrcr(:,:,:) = ztrcr_now(:,:,:)
        ENDWHERE
      ELSE
        ztrcr(:,:,:) = z1 * ztrcr_bd(:,:,:,nbd1) +                      &
                       z2 * ztrcr_bd(:,:,:,nbd2) 
      ENDIF

    ELSEIF ( izlbc(iztrcr) == T_LBC_CST ) THEN

      ! if constant boundary, simply copy values at tlev=nnow into values at nnew
      ztrcr(:,:,:) = ztrcr_now(:,:,:)   

    ELSEIF ( izlbc(iztrcr) == T_LBC_ZERO ) THEN

      ztrcr(:,:,:) = 0.0_wp

    ELSEIF ( izlbc(iztrcr) == T_LBC_ZEROGRAD ) THEN

      ! nothing to do since the values have to be computed first

    ELSE
      !error
      yzerrmsg = 'this type of LBC does not exist'
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF

  ENDDO ! loop over tracers

  ! Special case of QI
  !--------------------

  ! field for the cloud ice scheme
  CALL trcr_meta_get(izerror, idt_qi, T_LBC_ID, izlbcqi)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  ! if qi is a prognostic variables (i.e. the above call to meta_get did
  ! not return with a T_ERR_NOTFOUND but we have no boundary values from
  ! file, we need a special treatment...
  IF (ASSOCIATED(qi_new) .AND. (izerror==0) .AND. izlbcqi/=T_LBC_FILE) THEN

    ! Boundary values of cloud ice are interpreted from qc and
    ! qv is recalculated from relative humidity over ice below
    ! a threshold temperature. 
    DO k = 1, ke
!CDIR COLLAPSE
      qi_new (:,:,k) = 0.0_wp
!CDIR COLLAPSE
      IF ( MINVAL(t(:,:,k,nnew)) < 248.15_wp ) THEN
        DO j = 1, je
!CDIR COLLAPSE
          DO i = 1, ie
            IF ( t(i,j,k,nnew) < 248.15_wp ) THEN
              qi_new(i,j,k)  = qc_new(i,j,k)
              qc_new(i,j,k)  = 0.0_wp
              ! this reduction is only useful for GME data without cloud ice
              ! put all these computations into the INT2LM later
              !               qv(i,j,k,nnew) = qv(i,j,k,nnew) &
              !                      *fqvs( fpvsi(t(i,j,k,nnew)), p0(i,j,k)+pp(i,j,k,nnew) ) &
              !                      /fqvs( fpvsw(t(i,j,k,nnew)), p0(i,j,k)+pp(i,j,k,nnew) )
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  ! the initial values are reformulated similarily
  CALL trcr_meta_get(izerror, idt_qi, T_INI_ID, izanaqi)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  IF (ASSOCIATED(qi_now) .AND. (ntstep == 0) .AND. (izerror==0) .AND.         &
      (izanaqi /= T_INI_FILE) ) THEN
    DO k = 1, ke
!CDIR COLLAPSE
      IF ( MINVAL(t(:,:,k,nnow)) < 248.15_wp ) THEN
        DO j = 1, je
!CDIR COLLAPSE
          DO i = 1, ie
            IF ( t(i,j,k,nnow) < 248.15_wp ) THEN
              qi_now(i,j,k)  = qc_now(i,j,k)
              qi_nx0(i,j,k)  = qi_now(i,j,k)
              qc_now(i,j,k)  = 0.0_wp
              qc_nx0(i,j,k)  = 0.0_wp
              ! this reduction is only useful for GME data without cloud ice
              ! put all these computations into the INT2LM later
              !               qv(i,j,k,nnow) = qv(i,j,k,nnow) &
              !                      *fqvs( fpvsi(t(i,j,k,nnow)), p0(i,j,k)+pp(i,j,k,nnow) ) &
              !                      /fqvs( fpvsw(t(i,j,k,nnow)), p0(i,j,k)+pp(i,j,k,nnow) )
              !               qv(i,j,k,nx0)  = qv(i,j,k,nnow)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDDO
#ifdef TWOMOM_SB
    IF (itype_gscp >= 100) THEN
      qni_now(:,:,:) = set_qni_from_qi_sb(qi_now(:,:,:))
      IF (.NOT. l2tls) THEN
        qni_nx0(:,:,:) = set_qni_from_qi_sb(qi_nx0(:,:,:))
      END IF
    END IF
  ELSEIF ( (ntstep == 0) .and. (lartif_data) ) THEN
    IF (my_cart_id == 0) &
         WRITE (*,*) 'Warning initialize_loop: No initial adjustment of qv, qc, and qi'
#endif 
  ENDIF

#ifdef COSMOART
  IF (l_cosmo_art)  THEN
    ! CK 20101204 ART calculates its own interpolation weights
    IF(lgasbd .OR. laerobd) THEN
      ! factors for linear time interpolation
      z2_art = REAL (ntstep+1-nlastbound_art, wp) / &
               REAL (nincbound_art, wp)
      z2_art = MIN ( 1.0_wp , z2_art )
      z1_art = 1.0_wp - z2_art
    ENDIF
    IF (lgas) THEN
      CALL trcr_get_block(izerror, idx_start=trcr_idx_gas(1), idx_end=trcr_idx_gas(isp_gas), &
                 ptr_tlev = nnew, ptr = cgas_new)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
      CALL trcr_get_block(izerror, idx_start=trcr_idx_gas(1), idx_end=trcr_idx_gas(isp_gas), &
                 ptr_tlev = nnow, ptr = cgas_now, ptr_bd = cgas_bd)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
      IF (lgasbd)  THEN
!CDIR COLLAPSE
        ! CK 20101204 ART uses its own interpolation weights
        cgas_new(:,:,:,:) = z1_art*cgas_bd (:,:,:,:,nbd1_art) + &
                            z2_art*cgas_bd (:,:,:,:,nbd2_art)
      ELSE
!CDIR COLLAPSE
        cgas_new(:,:,:,:) = cgas_now(:,:,:,:)
      ENDIF
!CDIR COLLAPSE 
      ! CK 20101204 no more hardcoded references to certain species
      cgas_new(:,:,:,lho)   = cgas_now(:,:,:,lho)
!CDIR COLLAPSE
      cgas_new(:,:,:,lho2)   = cgas_now(:,:,:,lho2)
    ENDIF

    IF (laero)  THEN
      CALL trcr_get_block(izerror, idx_start=trcr_idx_aero(1), idx_end=trcr_idx_aero(isp_aero), &
                 ptr_tlev = nnew, ptr = caero_new)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
      CALL trcr_get_block(izerror, idx_start=trcr_idx_aero(1), idx_end=trcr_idx_aero(isp_aero), &
                 ptr_tlev = nnow, ptr = caero_now, ptr_bd = caero_bd)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
      IF (laerobd)  THEN
!CDIR COLLAPSE
        ! CK 20101204 ART uses its own interpolation weights
        DO isp=1,isp_aerotrans
          caero_new(:,:,:,isp) = z1_art*caero_bd (:,:,:,isp,nbd1_art) + &
                                 z2_art*caero_bd (:,:,:,isp,nbd2_art)
        ENDDO
        DO isp=1,isp_aerotrans+1,isp_aero
          caero_new(:,:,:,isp)  = caero_now(:,:,:,isp)
        ENDDO
      ELSE
!CDIR COLLAPSE
        caero_new(:,:,:,:)  = caero_now(:,:,:,:)
      ENDIF
    ENDIF
  ENDIF
#endif

#ifdef POLLEN
  IF (l_pollen)  THEN
    CALL trcr_get_block(izerror, idx_start=trcr_idx_pollen(1), idx_end=trcr_idx_pollen(isp_pollen), &
               ptr_tlev = nnow, ptr = cpollen_now, ptr_bd = cpollen_bd)
    IF (izerror /= 0) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get_block(izerror, idx_start=trcr_idx_pollen(1), idx_end=trcr_idx_pollen(isp_pollen), &
               ptr_tlev = nnew, ptr = cpollen_new)
    IF (izerror /= 0) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF
!CDIR COLLAPSE
    cpollen_new(:,:,:,:) = cpollen_now(:,:,:,:)
  ENDIF
#endif

  !--------------------------------------------------------------
  ! 3.1 initialize tendency fields
  !--------------------------------------------------------------

  ! initialize tendency fields with 0.0
  utens  (:,:,:) = 0.0_wp
  vtens  (:,:,:) = 0.0_wp
  wtens  (:,:,:) = 0.0_wp
  ttens  (:,:,:) = 0.0_wp
  pptens (:,:,:) = 0.0_wp

  ! loop over tracers
  DO iztrcr = 1, trcr_get_ntrcr()

    ! get pointer to tendency field
    CALL trcr_get( izerror, iztrcr, ptr_tens=ztrcr )

    ! set tendency to zero
    ztrcr(:,:,:) = 0.0_wp

  ENDDO

  !--------------------------------------------------------------
  ! 3.2 surface fields and its boundary settings
  !--------------------------------------------------------------

  ! Calculate the surface pressure ps for the new time level nnew
  CALL calps ( ps(:,:   ,nnew), pp(:,:,ke,nnew), t(:,:,ke,nnew),     &
               qv_new(:,:,ke ), qc_new(:,:,ke ), qrs(:,:,ke)   ,     &
               rho0(:,:,ke), p0(:,:,ke), dp0(:,:,ke),                &
               ie, je, rvd_m_o, r_d,                                 &
               istartpar, iendpar, jstartpar, jendpar )

  ! surface fields
  IF (lbdclim) THEN
    DO j = 1,je
!CDIR COLLAPSE
      DO i = 1,ie
        vio3(i,j) = z1 * vio3_bd(i,j,nbd1) + z2 * vio3_bd(i,j,nbd2)
        hmo3(i,j) = z1 * hmo3_bd(i,j,nbd1) + z2 * hmo3_bd(i,j,nbd2)

        IF ( itype_aerosol == 2 ) THEN
          aer_su(i,j) = z1 * aer_su_bd(i,j,nbd1) + z2 * aer_su_bd(i,j,nbd2)
          aer_du(i,j) = z1 * aer_du_bd(i,j,nbd1) + z2 * aer_du_bd(i,j,nbd2)
          aer_bc(i,j) = z1 * aer_bc_bd(i,j,nbd1) + z2 * aer_bc_bd(i,j,nbd2)
          aer_or(i,j) = z1 * aer_or_bd(i,j,nbd1) + z2 * aer_or_bd(i,j,nbd2)
          aer_ss(i,j) = z1 * aer_ss_bd(i,j,nbd1) + z2 * aer_ss_bd(i,j,nbd2)
        ENDIF

        IF ( llandmask(i,j) .EQV. .TRUE. ) THEN
          ! in this case no distinction between undef/defined (for frames)
          ! is done, because we assume that for climate simulations always
          ! the full field is given
          lai   (i,j) = z1 * lai_bd   (i,j,nbd1) + z2 * lai_bd   (i,j,nbd2)
          rootdp(i,j) = z1 * rootdp_bd(i,j,nbd1) + z2 * rootdp_bd(i,j,nbd2)
          plcov (i,j) = z1 * plcov_bd (i,j,nbd1) + z2 * plcov_bd (i,j,nbd2)
          IF (.NOT. lmulti_layer) THEN
            t_cl  (i,j) = z1 * t_cl_bd  (i,j,nbd1) + z2 * t_cl_bd  (i,j,nbd2)
            w_cl  (i,j) = z1 * w_cl_bd  (i,j,nbd1) + z2 * w_cl_bd  (i,j,nbd2)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  DO j = 1,je
!CDIR COLLAPSE
    DO i = 1,ie
      IF ( llandmask(i,j) .EQV. .TRUE. ) THEN

        IF (t_snow_bd   (i,j,nbd2) /= undef) THEN

          IF(lmulti_layer .AND. lmulti_snow) THEN
            t_snow_mult(i,j,1,nnew) = z1*t_snow_bd(i,j,nbd1) + z2*t_snow_bd(i,j,nbd2)
          ELSE
            t_snow(i,j,nnew) = z1*t_snow_bd(i,j,nbd1) + z2*t_snow_bd(i,j,nbd2)
          ENDIF

          qv_s  (i,j,nnew) = z1*qv_s_bd  (i,j,nbd1) + z2*qv_s_bd  (i,j,nbd2)
          w_snow(i,j,nnew) = z1*w_snow_bd(i,j,nbd1) + z2*w_snow_bd(i,j,nbd2)
          IF(.NOT.lmulti_layer) THEN
            t_s   (i,j,nnew) = z1*t_s_bd   (i,j,nbd1) + z2*t_s_bd (i,j,nbd2)
            t_m   (i,j,nnew) = z1*t_m_bd   (i,j,nbd1) + z2*t_m_bd (i,j,nbd2)
            w_g1  (i,j,nnew) = z1*w_g1_bd  (i,j,nbd1) + z2*w_g1_bd(i,j,nbd2)
            w_g2  (i,j,nnew) = z1*w_g2_bd  (i,j,nbd1) + z2*w_g2_bd(i,j,nbd2)
            IF ( nlgw == 3 ) THEN
              w_g3(i,j,nnew) = z1*w_g3_bd  (i,j,nbd1) + z2*w_g3_bd(i,j,nbd2)
            ENDIF
          ELSE
            t_s   (i,j,nnew) = t_s   (i,j,nx0) ! should only be used, if lphys=.F.
          ENDIF
        ELSE

          IF(lmulti_layer .AND. lmulti_snow) THEN
            t_snow_mult(i,j,1,nnew) = t_snow_mult(i,j,1,nnow)
          ELSE
            t_snow(i,j,nnew) = t_snow(i,j,nnow)
          ENDIF

          qv_s  (i,j,nnew) = qv_s  (i,j,nnow)
          w_snow(i,j,nnew) = w_snow(i,j,nnow)
          IF(.NOT.lmulti_layer) THEN
            t_s   (i,j,nnew) = t_s (i,j,nnow)
            t_m   (i,j,nnew) = t_m (i,j,nnow)
            w_g1  (i,j,nnew) = w_g1(i,j,nnow)
            w_g2  (i,j,nnew) = w_g2(i,j,nnow)
            IF ( nlgw == 3 ) THEN
              w_g3(i,j,nnew) = w_g3(i,j,nnow)
            ENDIF
          ELSE
            t_s   (i,j,nnew) = t_s   (i,j,nx0) ! should only be used, if lphys=.F.
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDDO

  DO j = 1,je
!CDIR COLLAPSE
    DO i = 1,ie
      IF ( llandmask(i,j) .EQV. .FALSE. ) THEN

        ! For lakes: compute saturation specific humidity over the lake
        ! surface at "nnow".
        ! This is required to compute surface fluxes in "flake_interface".

        IF (llake) THEN
          IF (depth_lk(i,j) > 0.0_wp) THEN
            ! Lake model is used and this is a lake point
            t_s   (i,j,nnew) = t_s   (i,j,nnow)
            IF(lmulti_layer .AND. lmulti_snow) THEN
              t_snow_mult(i,j,1,nnew) = t_snow_mult(i,j,1,nnow)
            ELSE
              t_snow(i,j,nnew) = t_snow(i,j,nnow)
            ENDIF
            IF (h_ice(i,j,nnow) < h_Ice_min_flk) THEN
              ! Water surface
              qv_s(i,j,nnow) = fqvs ( fpvsw ( t_s(i,j,nnow) ) , ps (i,j,nnow) )
            ELSE
              ! Ice surface
              qv_s(i,j,nnow) = fqvs ( fpvsi ( t_s(i,j,nnow) ) , ps (i,j,nnow) )
            END IF
            qv_s  (i,j,nnew) = qv_s  (i,j,nnow)
          ELSE
            ! This is a sea point
            IF (lbdclim) THEN
              ! in this case the sea surface temperature must not be constant!!
              IF(lmulti_layer .AND. lmulti_snow) THEN
                t_snow_mult(i,j,1,nnew) = z1*t_snow_bd(i,j,nbd1) + z2*t_snow_bd(i,j,nbd2)
              ELSE
                t_snow(i,j,nnew) = z1*t_snow_bd(i,j,nbd1) + z2*t_snow_bd(i,j,nbd2)
              ENDIF
              t_s   (i,j,nnew) = z1*t_s_bd   (i,j,nbd1) + z2*t_s_bd   (i,j,nbd2)
            ELSE
              IF(lmulti_layer .AND. lmulti_snow) THEN
                t_snow_mult(i,j,1,nnew) = t_snow_mult(i,j,1,nx0)
              ELSE
                t_snow(i,j,nnew) = t_snow(i,j,nx0)
              ENDIF
              IF (lbdsst) THEN
                t_s   (i,j,nnew) = z1*t_s_bd   (i,j,nbd1) + z2*t_s_bd   (i,j,nbd2)
              ELSE
                t_s   (i,j,nnew) = t_s   (i,j,nx0)
              ENDIF
            ENDIF
            IF (lseaice) THEN
              IF (h_ice(i,j,nnow) > 0.0_wp) THEN
                ! Ice surface
                qv_s(i,j,nnow) = fqvs ( fpvsi ( t_ice(i,j,nnow) ) , ps (i,j,nnow) )
              ELSE
                ! Water surface
                qv_s(i,j,nnow) = fqvs ( fpvsw ( t_s(i,j,nnow) ) , ps (i,j,nnow) )
              ENDIF
            ELSE
              qv_s  (i,j,nx0 ) = fqvs ( fpvsw ( t_s(i,j,nx0) ) , ps (i,j,nx0) )
              qv_s  (i,j,nnow) = qv_s  (i,j,nx0)
            ENDIF
            qv_s  (i,j,nnew) = qv_s  (i,j,nnow)
          ENDIF
        ELSE
          ! Lake model is not used
          IF (lbdclim) THEN
            ! in this case the sea surface temperature must not be constant!!
            IF(lmulti_layer .AND. lmulti_snow) THEN
              t_snow_mult(i,j,1,nnew) = z1*t_snow_bd(i,j,nbd1) + z2*t_snow_bd(i,j,nbd2)
            ELSE
              t_snow(i,j,nnew) = z1*t_snow_bd(i,j,nbd1) + z2*t_snow_bd(i,j,nbd2)
            ENDIF
            t_s   (i,j,nnew) = z1*t_s_bd   (i,j,nbd1) + z2*t_s_bd   (i,j,nbd2)
          ELSE
            IF(lmulti_layer .AND. lmulti_snow) THEN
              t_snow_mult(i,j,1,nnew) = t_snow_mult(i,j,1,nx0)
            ELSE
              t_snow(i,j,nnew) = t_snow(i,j,nx0)
            ENDIF
            IF (lbdsst) THEN
              t_s   (i,j,nnew) = z1*t_s_bd   (i,j,nbd1) + z2*t_s_bd   (i,j,nbd2)
            ELSE
              t_s   (i,j,nnew) = t_s   (i,j,nx0)
            ENDIF
          ENDIF
          IF (lseaice) THEN
            IF (h_ice(i,j,nnow) > 0.0_wp) THEN
              ! Ice surface
              qv_s(i,j,nnow) = fqvs ( fpvsi ( t_ice(i,j,nnow) ) , ps (i,j,nnow) )
            ELSE
              ! Water surface
              qv_s(i,j,nnow) = fqvs ( fpvsw ( t_s(i,j,nnow) ) , ps (i,j,nnow) )
            ENDIF
          ELSE
            qv_s  (i,j,nx0 ) = fqvs ( fpvsw ( t_s(i,j,nx0) ) , ps (i,j,nx0) )
            qv_s  (i,j,nnow) = qv_s  (i,j,nx0)
          ENDIF
          qv_s  (i,j,nnew) = qv_s  (i,j,nnow)
        ENDIF
      ENDIF
    ENDDO
  ENDDO

  ! compute the temperature at the boundary soil-atmosphere
  DO j = 1,je
!CDIR COLLAPSE
    DO i = 1,ie
      zt_s(i,j) = t_s(i,j,nnew)
      IF (lseaice) THEN
        IF (.NOT. llandmask(i,j) .AND. h_ice(i,j,nx0) > 0.0_wp) THEN
          zt_s(i,j) = t_ice(i,j,nx0)
        ENDIF
      ENDIF
    ENDDO
  ENDDO

  ! compute t_g on the full domain including boundaries
  ! (all fields are given on full domain)
  IF(lmulti_layer .AND. lmulti_snow) THEN
    CALL tgcom ( t_g (:,:,nnew), t_snow_mult(:,:,1,nnew), &
                 zt_s(:,:)     , w_snow(:,:,nnew), &
                 llandmask(:,:) , ie, je, cf_snow, &
                 1, ie, 1, je)
!US              istartpar, iendpar, jstartpar, jendpar )
  ELSE
    CALL tgcom ( t_g (:,:,nnew), t_snow(:,:,nnew), &
                 zt_s(:,:)     , w_snow(:,:,nnew), &
                 llandmask(:,:) , ie, je, cf_snow, &
                 1, ie, 1, je)
!US              istartpar, iendpar, jstartpar, jendpar )
  ENDIF

  ! compute density of moist air for time-level nnow        
  CALL calrho ( t(:,:,:,nnow), pp(:,:,:,nnow), qv_now(:,:,:), qc_now(:,:,:),  &
                qrs, p0, rho, ie, je, ke, r_d, rvd_m_o)

! Set one or more heating rate disturbances (up to 50). Time(s), location(s) and intensity(ies) are
! determined by namelist parameters of namelist IDEAL
  IF (lartif_data) THEN
    ! Set possible heating rate disturbance(s) in the atmosphere (affects ttens and qtens):
    CALL artif_heatrate_dist(nnow)
    ! Set possible disturbance(s) in the bottom boundary condition for t_s (takes effect if lsoil=.false.):
    ! in case the soil model is turned off:
    CALL set_tempdist_bbc_ts( )
    ! Copy soil temperature into timelevel nx0 for leapfrog integration:
    IF (.NOT. l2tls .AND. ntstep == 0) THEN
      t_s   (:,:,nx0) = t_s (:,: ,nnow)
      t_g   (:,:,nx0) = t_g (:,: ,nnow)
      t_snow(:,:,nx0) = t_snow (:,: ,nnow)
    END IF
    ! Add white noise to the T and W - fields:
    CALL  add_noise_tw(nnow)
  END IF

!-------------------------------------------------------------------------------
!  Section 4: Reinitialize vmax_10m
!-------------------------------------------------------------------------------

  IF (ntstep-1 == nnextmxu) THEN
    vmax_10m  (:,:) =   0.0_wp
    vabsmx_10m(:,:) =   0.0_wp
    vgust_dyn (:,:) =   0.0_wp
    vgust_con (:,:) =   0.0_wp

    ! Determine next step for re-initializing
    hlastmxu = hnextmxu
    hnextmxu = hlastmxu + hincmxu
    nlastmxu = NINT (hlastmxu * 3600.0_wp / dt)
    nnextmxu = NINT (hnextmxu * 3600.0_wp / dt)
  ENDIF

!-------------------------------------------------------------------------------
!  Section 5: Check for NaN's
!-------------------------------------------------------------------------------o

  IF (izdebug > 2) THEN
    lfound = .false.
    call check_field_NaNs(u(:,:,:,nnow),  'u::nnow', lfound, my_cart_id)
    call check_field_NaNs(v(:,:,:,nnow),  'v::nnow', lfound, my_cart_id)
    call check_field_NaNs(w(:,:,:,nnow),  'w::nnow', lfound, my_cart_id)
    call check_field_NaNs(t(:,:,:,nnow),  't::nnow', lfound, my_cart_id)
    call check_field_NaNs(pp(:,:,:,nnow),'pp::nnow', lfound, my_cart_id)
    IF (ASSOCIATED(qv_now)) THEN
      call check_field_NaNs(qv_now(:,:,:), 'qv::nnow', lfound, my_cart_id)
    END IF
    IF (ASSOCIATED(qc_now)) THEN
      call check_field_NaNs(qc_now(:,:,:), 'qc::nnow', lfound, my_cart_id)
    END IF
    IF (ASSOCIATED(qi_now)) THEN
      call check_field_NaNs(qi_now(:,:,:), 'qi::nnow', lfound, my_cart_id)
    END IF
    IF (ASSOCIATED(qr_now)) THEN
      call check_field_NaNs(qr_now(:,:,:), 'qr::nnow', lfound, my_cart_id)
    END IF
    IF (ASSOCIATED(qs_now)) THEN
      call check_field_NaNs(qs_now(:,:,:), 'qs::nnow', lfound, my_cart_id)
    END IF
    IF (ASSOCIATED(qg_now)) THEN
      call check_field_NaNs(qg_now(:,:,:), 'qg::nnow', lfound, my_cart_id)
    END IF
    IF (lfound) THEN
      izerror = 4242
      yzerrmsg = 'NaN encountered in prognostic field'
      CALL model_abort (my_world_id, 100+izerror, yzerrmsg, 'lmorg: initialize_loop')
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE initialize_loop

!==============================================================================
!==============================================================================

SUBROUTINE exchange_leapfrog

CHARACTER(LEN=25) :: yzroutine = 'exchange_leapfrog'

  IF (lzconv) THEN
    kzdims(1:24) =                                                          &
       (/ke,ke,ke,ke,ke1,ke1,ke,ke,ke,ke,ke,ke,                             &
         1,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                   &
       (nnew+39, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
        ie, je, kzdims, jstartpar, jendpar,                                 &
        nbl_exchg, nboundlines, my_cart_neigh,                              &
        lperi_x, lperi_y, l2dim,                                            &
        17000+nexch_tag, .FALSE.,    ncomm_type, izerror, yzerrmsg,         &
        u (:,:,:,nnow), u (:,:,:,nnew), v (:,:,:,nnow), v (:,:,:,nnew),     &
        w (:,:,:,nnow), w (:,:,:,nnew), t (:,:,:,nnow), t (:,:,:,nnew),     &
        pp(:,:,:,nnow), pp(:,:,:,nnew), qrs(:,:,:)    , dqvdt(:,:,:)  ,     &
        qvsflx(:,:) )
  ELSEIF (.NOT. lzconv) THEN
    kzdims(1:24) =                                                          &
       (/ke,ke,ke,ke,ke1,ke1,ke,ke,ke,ke,ke,0,                              &
         0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                   &
       (nnew+36, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
        ie, je, kzdims, jstartpar, jendpar,                                 &
        nbl_exchg, nboundlines, my_cart_neigh,                              &
        lperi_x, lperi_y, l2dim,                                            &
        17000+nexch_tag, .FALSE.,    ncomm_type, izerror, yzerrmsg,         &
        u (:,:,:,nnow), u (:,:,:,nnew), v (:,:,:,nnow), v (:,:,:,nnew),     &
        w (:,:,:,nnow), w (:,:,:,nnew), t (:,:,:,nnow), t (:,:,:,nnew),     &
        pp(:,:,:,nnow), pp(:,:,:,nnew), qrs(:,:,:)    )
  ENDIF

  ! loop over tracers
  DO iztrcr = 1, trcr_get_ntrcr()

    ! get pointer to tracer (at nnow)
    CALL trcr_get(izerror, iztrcr, ptr_tlev=nnow, ptr=ztrcr_now)
    IF (izerror /= 0_iintegers) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF

    ! get pointer to tracer (at nnew)
    CALL trcr_get(izerror, iztrcr, ptr_tlev=nnew, ptr=ztrcr)
    IF (izerror /= 0_iintegers) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF

    kzdims(1:24) =                                                     &
          (/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    ! Final boundary exchange of tracers at tlev=nnow and nnew (leapfrog)
    CALL exchg_boundaries                                              &
         (nnew+55, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
          ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,         &
          my_cart_neigh, lperi_x, lperi_y, l2dim, 18000+nexch_tag+iztrcr,     &
          ldatatypes, ncomm_type, izerror, yzerrmsg,                          &
          ztrcr_now(:,:,:), ztrcr(:,:,:))

  ENDDO

END SUBROUTINE exchange_leapfrog

!==============================================================================
!==============================================================================

SUBROUTINE exchange_runge_kutta
  
CHARACTER (LEN=25) :: yzroutine='exchange_runge_kutta'
INTEGER(KIND=iintegers) :: zntke

  kzdims(1:24)=(/ke,ke,ke1,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
  CALL exchg_boundaries                                                  &
   (50+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
    ie, je, kzdims, jstartpar, jendpar,                                  &
    nbl_exchg, nboundlines, my_cart_neigh,                               &
    lperi_x, lperi_y, l2dim,                                             &
    17000+nexch_tag, ldatatypes, ncomm_type, izerror, yzerrmsg,          &
    u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),      &
    pp(:,:,:,nnew), qrs(:,:,:) )          

  IF ( lzconv ) THEN
    IF ( lprog_tke ) THEN
      IF (itype_turb /= 3 .OR. ntke == 0) THEN
        zntke = nnew
      ELSE
        zntke = ntke
      ENDIF
      kzdims(1:24)=(/ke,1,ke1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
       (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,           &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        18000+nexch_tag, .FALSE., ncomm_type, izerror, yzerrmsg,               &
        dqvdt(:,:,:), qvsflx(:,:), tke(:,:,:,zntke) )
    ELSE
      kzdims(1:24)=(/ke,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
       (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,           &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        18000+nexch_tag, .FALSE., ncomm_type, izerror, yzerrmsg,               &
        dqvdt(:,:,:), qvsflx(:,:) )
    END IF
  ELSE
    IF ( lprog_tke ) THEN
      IF (itype_turb /= 3 .OR. ntke == 0) THEN
        zntke = nnew
      ELSE
        zntke = ntke
      ENDIF
      kzdims(1:24)=(/ke1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
       (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,           &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        18000+nexch_tag, .FALSE., ncomm_type, izerror, yzerrmsg,               &
        tke(:,:,:,zntke) )
    END IF
  END IF

  ! loop over tracers
  DO iztrcr = 1, trcr_get_ntrcr()

    ! get pointer to tracer (at nnew)
    CALL trcr_get(izerror, iztrcr, ptr_tlev=nnew, ptr=ztrcr)
    IF (izerror /= 0_iintegers) THEN
      yzerrmsg = trcr_errorstr(izerror)
      CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
    ENDIF

    ! halo-update
    kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                     &
     (55+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,      &
      ie, je, kzdims, jstartpar, jendpar,                                     &
      nbl_exchg, nboundlines, my_cart_neigh,                                  &
      lperi_x, lperi_y, l2dim,                                                &
      19000+nexch_tag+iztrcr, ldatatypes, ncomm_type, izerror, yzerrmsg,      &
      ztrcr(:,:,:) )

  ENDDO

END SUBROUTINE exchange_runge_kutta

!==============================================================================
!==============================================================================

SUBROUTINE set_trcr_special_bc

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure defines the lateral boundary conditions
!   in case of zero gradient or zero value.
!
! Method:
!
!   The values at the outermost points of the forecast domain (inner domain) 
!   are assigned to the halo points around the forecast domain in case of 
!   zero gradient and 0 is simply assigned to the halo points in case of zero
!   boundaries.
!------------------------------------------------------------------------------

CHARACTER(LEN=25) :: yzroutine = 'set_trcr_special_bc'

  ! loop over tracers 
  DO iztrcr = 1, trcr_get_ntrcr()

  ! check if zero-gradient boundaries are required
    IF (izlbc(iztrcr) == T_LBC_ZEROGRAD                                       &
        .OR. izbd_forced(iztrcr) == 1_iintegers) THEN


      ! get pointer to tracer (at nnew)
      CALL trcr_get(izerror, iztrcr, ptr_tlev=nnew, ptr=ztrcr)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF

      ! western boundary
      IF (my_cart_neigh(1) == -1) THEN
        DO i = 1, nboundlines
          DO  k = 1, ke
            DO  j = jstart, jend
              ztrcr(i,j,k) = ztrcr(istart,j,k)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      ! eastern boundary
      IF (my_cart_neigh(3) == -1) THEN
        DO i = ie-nboundlines+1, ie
          DO  k = 1, ke
            DO  j = jstart, jend
              ztrcr(i,j,k) = ztrcr(iend  ,j,k)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      ! southern boundary
      IF (my_cart_neigh(4) == -1) THEN
        DO  k = 1, ke
          DO  j = 1, nboundlines
            ztrcr(:,j,k) = ztrcr(:,jstart,k)
          ENDDO
        ENDDO
      ENDIF

      ! northern boundary
      IF (my_cart_neigh(2) == -1) THEN
        DO  k = 1, ke
          DO  j = je-nboundlines+1, je
            ztrcr(:,j,k) = ztrcr(:,jend  ,k)
          ENDDO
        ENDDO
      ENDIF

    ! check if zero-value boundary conditions are required
    ELSEIF ( izlbc(iztrcr) == T_LBC_ZERO .AND.                                &
             izbd_forced(iztrcr) == 2_iintegers ) THEN


      ! get pointer to tracer (at nnew)
      CALL trcr_get(izerror, iztrcr, ptr_tlev=nnew, ptr=ztrcr)
      IF (izerror /= 0) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF

      ! western boundary
      IF (my_cart_neigh(1) == -1) THEN
        DO i = 1, nboundlines
          DO  k = 1, ke
            DO  j = jstart, jend
              ztrcr(i,j,k) = 0.0_wp
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      ! eastern boundary
      IF (my_cart_neigh(3) == -1) THEN
        DO i = ie-nboundlines+1, ie
          DO  k = 1, ke
            DO  j = jstart, jend
              ztrcr(i,j,k) = 0.0_wp

            ENDDO
          ENDDO
        ENDDO
      ENDIF

      ! southern boundary
      IF (my_cart_neigh(4) == -1) THEN
        DO  k = 1, ke
          DO  j = 1, nboundlines
            ztrcr(:,j,k) = 0.0_wp
          ENDDO
        ENDDO
      ENDIF

      ! northern boundary
      IF (my_cart_neigh(2) == -1) THEN
        DO  k = 1, ke
          DO  j = je-nboundlines+1, je
            ztrcr(:,j,k) = 0.0_wp
          ENDDO
        ENDDO
      ENDIF

    ENDIF

  ENDDO ! loop over tracers

END SUBROUTINE set_trcr_special_bc

!==============================================================================
!==============================================================================

SUBROUTINE nullify_tracers

CHARACTER(LEN=25) :: yzroutine='nullify_tracers'

!UB This is not consistent with the below clipping of the single hydrometeors
!UB  and can cause small inconsistencies in restart. Therefore we
!UB  recompute qrs after hydrometeor clipping:
!UB ! If the humidity variables are rather small, set them to 0
!UB WHERE (qrs(:,:,:) < 1.0E-15_wp)
!UB   qrs(:,:,:) = 0.0_wp
!UB ENDWHERE

#ifndef MESSY
    ! no clipping for MESSy Tracers this is done and budgeted in
    ! in the MESSy subsubmodel tracer_pdef
    ! loop over tracers
    DO i = 1,trcr_get_ntrcr()

      ! check if clipping is required
      IF ( izclp(i) == T_CLP_ON ) THEN
        ! get pointer to tracer (at nnew)
        CALL trcr_get(izerror, i, ptr_tlev=nnew, ptr=ztrcr)
        IF (izerror /= 0_iintegers) THEN
          yzerrmsg = trcr_errorstr(izerror)
          CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
        ENDIF

        ! clip to zero
        WHERE (ztrcr(:,:,:) < 1.0E-15_wp)
          ztrcr(:,:,:) = 0.0_wp
        ENDWHERE

      ENDIF

    ENDDO ! loop over tracers

    ! Recompute qrs from clipped variables:
    qrs = 0.0_wp
    IF (idt_qr /= -99_iintegers) THEN
      CALL trcr_get(izerror, idt_qr, ptr_tlev=nnew, ptr=ztrcr)
      IF (izerror /= 0_iintegers) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
      qrs(:,:,:) = qrs(:,:,:) + ztrcr(:,:,:)
    END IF
    IF (idt_qi /= -99_iintegers) THEN
      CALL trcr_get(izerror, idt_qi, ptr_tlev=nnew, ptr=ztrcr)
      IF (izerror /= 0_iintegers) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
      qrs(:,:,:) = qrs(:,:,:) + ztrcr(:,:,:)
    END IF
    IF (idt_qs /= -99_iintegers) THEN
      CALL trcr_get(izerror, idt_qs, ptr_tlev=nnew, ptr=ztrcr)
      IF (izerror /= 0_iintegers) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
      qrs(:,:,:) = qrs(:,:,:) + ztrcr(:,:,:)
    END IF
    IF (idt_qg /= -99_iintegers) THEN
      CALL trcr_get(izerror, idt_qg, ptr_tlev=nnew, ptr=ztrcr)
      IF (izerror /= 0_iintegers) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
      qrs(:,:,:) = qrs(:,:,:) + ztrcr(:,:,:)
    END IF
#ifdef TWOMOM_SB
    IF (idt_qh /= -99_iintegers) THEN
      CALL trcr_get(izerror, idt_qh, ptr_tlev=nnew, ptr=ztrcr)
      IF (izerror /= 0_iintegers) THEN
        yzerrmsg = trcr_errorstr(izerror)
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      ENDIF
      qrs(:,:,:) = qrs(:,:,:) + ztrcr(:,:,:)
    END IF
#endif
#endif

END SUBROUTINE nullify_tracers


end module enkf_cosmo_mod
