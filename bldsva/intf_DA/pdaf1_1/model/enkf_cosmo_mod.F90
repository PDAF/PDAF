module enkf_cosmo_mod
USE info_lm_f90,         ONLY: info_define, info_readnl, info_print

USE data_parameters,     ONLY:   ireals, iintegers
USE data_constants,      ONLY:   b1, b2w, b3, b4w, b2i, b4i, rdv, o_m_rdv,   &
                                 rvd_m_o, r_d
USE data_soil,           ONLY:   cf_snow
USE data_fields,         ONLY:                                               &
       dp0, p0, rho, rho0, qrs, llandmask, t_g, vmax_10m, dqvdt, qvsflx,     &
       u, v, w, t, qv, qc, qi, qr, qs, qg, tke, pp, ps, u_bd, v_bd, w_bd,    &
       t_bd, qv_bd, qc_bd, qi_bd, qr_bd, qs_bd, qg_bd, pp_bd,                &
       t_snow,    t_s,    qv_s,    t_m,    w_snow,    w_g1,    w_g2,         &
       t_snow_bd, t_s_bd, qv_s_bd, t_m_bd, w_snow_bd, w_g1_bd, w_g2_bd,      &
       plcov_bd, lai_bd, rootdp_bd, vio3_bd, hmo3_bd, t_cl_bd, w_cl_bd,      &
       plcov,    lai,    rootdp,    vio3,    hmo3,    t_cl,    w_cl,         &
       utens, vtens, wtens, ttens, qvtens, qctens, pptens, w_g3, w_g3_bd,    &
       fr_lake, depth_lk, h_ice, qvt_diff, qitens, vgust_con, vgust_dyn,     &
       prr_gsp, prr_con, prne_con, bas_con, t_snow_mult, t_ice

!AK (20.03.12)
USE data_tracer,         ONLY:   tracer, ntracer, ltracer, tracer_bd,        &
                                 tracertens
!AK (20.03.12)

USE data_modelconfig,    ONLY:   ie, je, ke, ke1, jstartpar, jendpar,        &
                                 istartpar, iendpar, ivctype, vcoord,        &
                                 dt, dt2, ed2dt, dtdeh, nehddt, hhlr, sigmr, &
                                 vcflat, kflat, p0sl, t0sl, dt0lp,           &
                                 jstartu, jendu, jstartv, jendv,             &
                                 istart, iend, jstart, jend, ie_tot, je_tot

USE data_runcontrol,     ONLY:                                               &
       nstart, nstop, ntstep, nlastbound, nincbound, lmulti_layer,           &
       ltime, itype_timing, newbcdt, newbc, nincboufac, nlgw, itype_gscp,    &
       nbd1, nbd2, nold, nnow, nnew, lartif_data, ldiagnos, luse_rttov,      &
       lphys, lconv, luseobs, l2tls, lsemi_imp, lgsp, lprog_qi, lprogprec,   &
       yakdat1, yakdat2, nuspecif, ldfi, nincconv, lconf_avg, l2dim,         &
       irunge_kutta, nbl_exchg, lprog_tke, lspecnudge, lseaice, llake,       &
       lw_freeslip, leps, idbg_level, lprintdeb_all, itype_calendar,         &
       hlastmxu, hnextmxu, hincmxu, nlastmxu, nnextmxu, l_cosmo_art,         &
       l_pollen, hstart, itype_lbc_qrsg, lmulti_snow, lperi_x, lperi_y, l2dim

USE data_parallel,       ONLY:                                               &
       num_compute, icomm_cart, my_cart_id, sendbuf, isendbuflen, iexch_req, &
       imp_reals, nboundlines, my_cart_neigh, nprocx, nprocy, my_world_id,   &
       lcompute_pe, lasync_io, ncomm_type, ldatatypes, ltime_barrier

USE data_io,             ONLY:   ydate_ini, pp_nl, root, lbdclim,            &
                                 lana_qi, llb_qi, llb_qr_qs, llb_qg,         &
                                 lbd_frame, undef

USE data_flake,          ONLY:   h_Ice_min_flk

USE mpe_io,              ONLY:   mpe_io_node, mpe_io_shutdown

USE environment,         ONLY:   exchg_boundaries, comm_barrier,             &
                                 final_environment, model_abort, get_free_unit
USE meteo_utilities,     ONLY:   calrho, calps, tgcom
USE time_utilities,      ONLY:   nutiming, init_timings, get_timings,        &
                  i_initializations, i_add_computations, i_phy_computations, &
                  i_dyn_computations,i_lhn_computations, i_spectr_nudging,   &
                  i_relaxation, i_input, i_output, i_barrier_waiting_dyn,    &
                  i_communications_dyn, i_cleanup, collect_timings
USE utilities,           ONLY:   get_utc_date

USE src_setup,           ONLY:   organize_setup, constant_fields

USE src_allocation,      ONLY:   organize_allocation

USE src_artifdata,       ONLY:   artif_heatrate_dist, artif_heatrate_dist_tso, set_tempdist, &
                                 set_tempdist_tso, set_tempdist_bbc_ts, add_noise_tw

#ifdef COSMOART
USE data_cosmo_art,      ONLY:   lgas, laero, ldust, lseas,                  &
                                 lrad_dust, lrad_aero,                       &
                                 lgasbd, laerobd, lgasini, laeroini,         &
                                 cgas, caero, cgas_bd, caero_bd, cgastens,   &
                                 caerotens, isp_gas, isp_aero,               &
                                 nlastbound_art, nincbound_art,              &
                                 nbd1_art, nbd2_art
USE art_utilities,       ONLY:   gm3ppm_a, ppmgm3_a
USE art_gas_const,       ONLY:   lho, lho2
USE art_aerosol_const,   ONLY:   aerostart
#endif

#ifdef POLLEN
USE data_pollen,         ONLY:   cpollen, cpollentens, cpollen_bd,           &
                                 isp_pollen, dtpollen
#endif

!AK (20.03.12)
USE src_tracer_supply,   ONLY:   organize_tracer_init, organize_tracer,       &
                                 organize_tracer_bound, organize_tracer_source
!AK (20.03.12) 
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
  nzhours,       & ! for recording a forecast hour
  nzdays,        & ! for recording a forecast day
  i,j,k, kzdims(24), izdebug, nzdiv, isp, kstart,     &
!AK (20.03.2012)
  iig,           & ! loop index for ltracer
  iprog,         & ! loop index for prognostic processes of tracers
  nprog            ! number of prognostic processes of tracer
!AK (20.03.2012)

REAL (KIND=ireals)        ::                          &
  zgrids_dt(1),  & ! structure for interface to organize_data
                   ! (would be necessary only for the nesting version)
  zforecasttime, & ! for recording a forecast hour
  zrealdiff        ! for time measurement

LOGICAL                   ::                          &
  lopen,         & ! check whether YUCHKDAT is opened
  lzconv           ! to determine whether extra communication for convection
                   ! is needed

CHARACTER (LEN=46)        :: ynote
CHARACTER (LEN=80)        :: yzerrmsg

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
INTEGER (KIND=iintegers), INTENT (IN)     ::       &
  ntstep             ! actual time step

! Scalar arguments with intent(inout):
INTEGER (KIND=iintegers), INTENT (INOUT)  ::       &
  nbd1, nbd2,      & ! indices for the boundary levels
  nold, nnow, nnew   ! indices for the time levels

!------------------------------------------------------------------------------
! Local variables

INTEGER (KIND=iintegers)   ::  &
  nsp, nx0, i, j, k, i1,       & ! support variable
  nzjulianday                    ! day of the year

REAL (KIND=ireals)         ::        &
  tact, tlbc, tivl,            & ! temporal quantities for control output
  z1, z2,                      & ! factors for time interpolation
  fpvsw, fpvsi, fqvs, zst, zsge, zsp,  & ! for the statement functions
  zacthour                               ! actual hour of the forecast

#ifdef COSMOART
! CK 20101204
! ART bd update frequency does not need to be the same as meteo.
! therefore, the weights calculated for meteo might be wrong.

REAL (KIND=ireals)         ::        &
  z1_art, z2_art                       ! factors for time interpolation
#endif

REAL (KIND=ireals)         ::        &
  zt_s(ie,je)                    ! = t_s   on land and sea
                                 ! = t_ice on sea ice (if present)

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

! get new actual utc date
  CALL get_utc_date(ntstep, ydate_ini, dt, itype_calendar, yakdat1, yakdat2, &
                    nzjulianday, zacthour)

!------------------------------------------------------------------------------
!- Section 1: input of new boundary values, if necessary
!------------------------------------------------------------------------------

  IF (ltime) CALL get_timings (i_add_computations, ntstep, dt, izerror)

  ! get new boundary data file or generate artificial data
  zgrids_dt(1) = dt
  CALL organize_data ('boundary', ntstep, 1, 1, zgrids_dt, &
                       izerror, yzerrmsg)
  IF (izerror /= 0_iintegers) THEN
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg,                &
                                   'boundary: input-init')
  ENDIF

  IF (ltime) CALL get_timings (i_input, ntstep, dt, izerror)

!------------------------------------------------------------------------------
!- Section 2: update organizational variables
!------------------------------------------------------------------------------

  ! cyclic changing of the indices for the time levels
  IF (l2tls) THEN
    nnow = 3 - nnow
    nnew = 3 - nnew
    nx0  = nnow
  ELSE
    nsp    = nold
    nold   = nnow
    nnow   = nnew
    nnew   = nsp
    nx0    = nold
  ENDIF

  ! variables concerned with the time step
  dt2    = 2.0 * dt
  ed2dt  = 1.0 / dt2
  dtdeh  = dt / 3600.0
  nehddt = NINT ( 3600.0_ireals / dt )

!------------------------------------------------------------------------------
!- Section 3: Initialize new time level with boundary values
!------------------------------------------------------------------------------

  ! factors for linear time interpolation
  z2 = REAL (ntstep+1-nlastbound, ireals) / REAL (nincbound, ireals)
  z2 = MIN ( 1.0_ireals , z2 )
  z1 = 1.0 - z2

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
!CDIR COLLAPSE
      qv(:,:,:,nnew) = z1 * qv_bd(:,:,:,nbd1) + z2 * qv_bd(:,:,:,nbd2)
!CDIR COLLAPSE
      qc(:,:,:,nnew) = z1 * qc_bd(:,:,:,nbd1) + z2 * qc_bd(:,:,:,nbd2)
    ELSEWHERE
!CDIR COLLAPSE
      u (:,:,:,nnew) = u (:,:,:,nnow)
!CDIR COLLAPSE
      v (:,:,:,nnew) = v (:,:,:,nnow)
!CDIR COLLAPSE
      t (:,:,:,nnew) = t (:,:,:,nnow)
!CDIR COLLAPSE
      pp(:,:,:,nnew) = pp(:,:,:,nnow)
!CDIR COLLAPSE
      qv(:,:,:,nnew) = qv(:,:,:,nnow)
!CDIR COLLAPSE
      qc(:,:,:,nnew) = qc(:,:,:,nnow)
!CDIR COLLAPSE
    ENDWHERE

!AK (20.03.12)
    DO iig=1, ntracer
      WHERE (t_bd (:,:,:,nbd2) /= undef)
        tracer(:,:,:,nnew,iig) = z1 * tracer_bd(:,:,:,nbd1,iig) + z2 * tracer_bd(:,:,:,nbd2,iig)
      ELSEWHERE
        tracer(:,:,:,nnew,iig) = tracer(:,:,:,nnow,iig)
      ENDWHERE
    ENDDO
!AK (20.03.12)

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
!CDIR COLLAPSE
    qv(:,:,:,nnew) = z1 * qv_bd(:,:,:,nbd1) + z2 * qv_bd(:,:,:,nbd2)
!CDIR COLLAPSE
    qc(:,:,:,nnew) = z1 * qc_bd(:,:,:,nbd1) + z2 * qc_bd(:,:,:,nbd2)

!AK (20.03.12)
    DO iig=1, ntracer
      tracer(:,:,:,nnew,iig) = z1 * tracer_bd(:,:,:,nbd1,iig) + z2 * tracer_bd(:,:,:,nbd2,iig)
    ENDDO
!AK (20.03.12)

 IF (.NOT. lw_freeslip) THEN
!CDIR COLLAPSE
      w (:,:,:,nnew) = z1 * w_bd (:,:,:,nbd1) + z2 * w_bd (:,:,:,nbd2)
    END IF
  ENDIF

#ifdef COSMOART
  IF (l_cosmo_art)  THEN
    ! CK 20101204 ART calculates its own interpolation weights
    IF(lgasbd .OR. laerobd) THEN
      ! factors for linear time interpolation
      z2_art = REAL (ntstep+1-nlastbound_art, ireals) / &
               REAL (nincbound_art, ireals)
      z2_art = MIN ( 1.0_ireals , z2_art )
      z1_art = 1.0 - z2_art
    ENDIF
    IF (lgas) THEN
      IF (lgasbd)  THEN
!CDIR COLLAPSE
        ! CK 20101204 ART uses its own interpolation weights
        cgas(:,:,:,:,nnew) = z1_art*cgas_bd (:,:,:,:,nbd1_art) + &
                             z2_art*cgas_bd (:,:,:,:,nbd2_art)
      ELSE
!CDIR COLLAPSE
        cgas(:,:,:,:,nnew)   = cgas(:,:,:,:,nnow)
      ENDIF
!CDIR COLLAPSE
      ! CK 20101204 no more hardcoded references to certain species
      cgas(:,:,:,lho,nnew)   = cgas(:,:,:,lho,nnow)
!CDIR COLLAPSE
      cgas(:,:,:,lho2,nnew)   = cgas(:,:,:,lho2,nnow)
    ENDIF

    IF (laero)  THEN
      IF (laerobd)  THEN
!CDIR COLLAPSE
        ! CK 20101204 ART uses its own interpolation weights
        caero(:,:,:,:,nnew) = z1_art*caero_bd (:,:,:,:,nbd1_art) + &
                              z2_art*caero_bd (:,:,:,:,nbd2_art)
      ELSE
!CDIR COLLAPSE
        caero(:,:,:,:,nnew)  = caero(:,:,:,:,nnow)
      ENDIF
    ENDIF
  ENDIF
#endif

#ifdef POLLEN
!CDIR COLLAPSE
  IF (l_pollen) cpollen(:,:,:,:,nnew) = cpollen(:,:,:,:,nnow)
#endif

  ! field for the cloud ice scheme
  IF (lprog_qi) THEN
    IF (llb_qi) THEN
      IF (lbd_frame) THEN
!CDIR COLLAPSE
        WHERE (qi_bd (:,:,:,nbd2) /= undef)
!CDIR COLLAPSE
          qi(:,:,:,nnew) = z1 * qi_bd(:,:,:,nbd1) + z2 * qi_bd(:,:,:,nbd2)
        ELSEWHERE
!CDIR COLLAPSE
          qi(:,:,:,nnew) = qi(:,:,:,nnow)
        END WHERE
      ELSE
!CDIR COLLAPSE
        qi(:,:,:,nnew) = z1 * qi_bd(:,:,:,nbd1) + z2 * qi_bd(:,:,:,nbd2)
      ENDIF
    ELSE
      ! Boundary values of cloud ice are interpreted from qc and
      ! qv is recalculated from relative humidity over ice below
      ! a threshold temperature. 
      DO k = 1, ke
!CDIR COLLAPSE
        qi(:,:,k,nnew) = 0.0_ireals
!CDIR COLLAPSE
        IF ( MINVAL(t(:,:,k,nnew)) < 248.15_ireals ) THEN
          DO j = 1, je
!CDIR COLLAPSE
            DO i = 1, ie
              IF ( t(i,j,k,nnew) < 248.15_ireals ) THEN
                qi(i,j,k,nnew) = qc(i,j,k,nnew)
                qc(i,j,k,nnew) = 0.0_ireals
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
  ENDIF

  ! the initial values are reformulated similarily
  IF (lprog_qi) THEN
    IF ( (ntstep == 0) .AND. (.NOT. lana_qi) ) THEN
      DO k = 1, ke
!CDIR COLLAPSE
        IF ( MINVAL(t(:,:,k,nnow)) < 248.15_ireals ) THEN
          DO j = 1, je
!CDIR COLLAPSE
            DO i = 1, ie
              IF ( t(i,j,k,nnow) < 248.15_ireals ) THEN
                qi(i,j,k,nnow) = qc(i,j,k,nnow)
                qi(i,j,k,nx0)  = qi(i,j,k,nnow)
                qc(i,j,k,nnow) = 0.0_ireals
                qc(i,j,k,nx0)  = 0.0_ireals
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
    ENDIF
  ENDIF

  ! prognostic precipitation fields
  IF (llb_qr_qs) THEN
    IF (lbd_frame) THEN
!CDIR COLLAPSE
      WHERE (qr_bd (:,:,:,nbd2) /= undef)
!CDIR COLLAPSE
        qr(:,:,:,nnew) = z1 * qr_bd(:,:,:,nbd1) + z2 * qr_bd(:,:,:,nbd2)
!CDIR COLLAPSE
        qs(:,:,:,nnew) = z1 * qs_bd(:,:,:,nbd1) + z2 * qs_bd(:,:,:,nbd2)
      ELSEWHERE
!CDIR COLLAPSE
        qr(:,:,:,nnew) = qr(:,:,:,nnow)
!CDIR COLLAPSE
        qs(:,:,:,nnew) = qs(:,:,:,nnow)
      END WHERE
    ELSE
!CDIR COLLAPSE
      qr(:,:,:,nnew) = z1 * qr_bd(:,:,:,nbd1) + z2 * qr_bd(:,:,:,nbd2)
!CDIR COLLAPSE
      qs(:,:,:,nnew) = z1 * qs_bd(:,:,:,nbd1) + z2 * qs_bd(:,:,:,nbd2)
    ENDIF
  ENDIF
  IF (llb_qg) THEN
    IF (lbd_frame) THEN
!CDIR COLLAPSE
      WHERE (qg_bd (:,:,:,nbd2) /= undef)
!CDIR COLLAPSE
        qg(:,:,:,nnew) = z1 * qg_bd(:,:,:,nbd1) + z2 * qg_bd(:,:,:,nbd2)
      ELSEWHERE
!CDIR COLLAPSE
        qg(:,:,:,nnew) = qg(:,:,:,nnow)
      END WHERE
    ELSE
!CDIR COLLAPSE
      qg(:,:,:,nnew) = z1 * qg_bd(:,:,:,nbd1) + z2 * qg_bd(:,:,:,nbd2)
    ENDIF
  ENDIF

  ! initialize tendency fields with 0.0
  utens  (:,:,:) = 0.0_ireals
  vtens  (:,:,:) = 0.0_ireals
  wtens  (:,:,:) = 0.0_ireals
  ttens  (:,:,:) = 0.0_ireals
  qvtens (:,:,:) = 0.0_ireals
  qctens (:,:,:) = 0.0_ireals
!WS (01.06.2012) tracertens is already set to zero before call of organize_tracer_source
!  DO iig = 1, ntracer
!    tracertens(:,:,:,iig) = 0.0_ireals
!  ENDDO
!WS (01.06.2012)

  ! Initialise qi tendency so that cloud ice can be treated
  ! in the same way as cloud water.
  qitens (:,:,:) = 0.0_ireals
  pptens (:,:,:) = 0.0_ireals
  ! ... also for humidity tendency due to diffusion
  qvt_diff (:,:,:) = 0.0_ireals

#ifdef COSMOART
  IF (l_cosmo_art) THEN
    IF (lgas)   cgastens (:,:,:,:) = 0.0_ireals
    IF (laero)  caerotens(:,:,:,:) = 0.0_ireals
  ENDIF
#endif

#ifdef POLLEN
  IF (l_pollen) cpollentens(:,:,:,:) = 0.0_ireals
#endif

  ! Calculate the surface pressure ps for the new time level nnew
  CALL calps ( ps(:,:   ,nnew), pp(:,:,ke,nnew), t(:,:,ke,nnew),     &
               qv(:,:,ke,nnew), qc(:,:,ke,nnew), qrs(:,:,ke)   ,     &
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
          IF (depth_lk(i,j) > 0.0_ireals) THEN
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
              t_s   (i,j,nnew) = t_s   (i,j,nx0)
            ENDIF
            IF (lseaice) THEN
              IF (h_ice(i,j,nnow) > 0.0_ireals) THEN
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
            t_s   (i,j,nnew) = t_s   (i,j,nx0)
          ENDIF
          IF (lseaice) THEN
            IF (h_ice(i,j,nnow) > 0.0_ireals) THEN
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
        IF (.NOT. llandmask(i,j) .AND. h_ice(i,j,nx0) > 0.0_ireals) THEN
          zt_s(i,j) = t_ice(i,j,nx0)
        ENDIF
      ENDIF
    ENDDO
  ENDDO

#if !defined COUP_OAS_COS
  IF(lmulti_layer .AND. lmulti_snow) THEN
    CALL tgcom ( t_g (:,:,nnew), t_snow_mult(:,:,1,nnew), &
                 zt_s(:,:)     , w_snow(:,:,nnew), &
                 llandmask(:,:) , ie, je, cf_snow, &
                 istartpar, iendpar, jstartpar, jendpar )
  ELSE
    CALL tgcom ( t_g (:,:,nnew), t_snow(:,:,nnew), &
                 zt_s(:,:)     , w_snow(:,:,nnew), &
                 llandmask(:,:) , ie, je, cf_snow, &
                 istartpar, iendpar, jstartpar, jendpar )
  ENDIF
#endif
  ! compute density of moist air for time-level nnow        
  CALL calrho ( t(:,:,:,nnow), pp(:,:,:,nnow), qv(:,:,:,nnow), qc(:,:,:,nnow),&
                qrs, p0, rho, ie, je, ke, r_d, rvd_m_o)

! Set one or more heating rate disturbances (up to 50). Time(s), location(s) and intensity(ies) are
! determined by namelist parameters of namelist IDEAL
  IF (lartif_data) THEN
    ! Set possible heating rate disturbance(s) in the atmosphere (affects ttens and qtens):
    CALL artif_heatrate_dist(nnow)
    ! Set possible disturbance(s) in the bottom boundary condition for t_s (takes effect if lsoil=.false.):
    ! in case the soil model is turned off:
    CALL set_tempdist_bbc_ts( )
    ! Add white noise to the T and W - fields:
    CALL  add_noise_tw(nnow)
  END IF

!-------------------------------------------------------------------------------
!  Section 4: Reinitialize vmax_10m
!-------------------------------------------------------------------------------

  IF (ntstep-1 == nnextmxu) THEN
    vmax_10m (:,:) =   0.0_ireals
    vgust_dyn(:,:) =   0.0_ireals
    vgust_con(:,:) =   0.0_ireals

    ! Determine next step for re-initializing
    hlastmxu = hnextmxu
    hnextmxu = hlastmxu + hincmxu
    nlastmxu = NINT (hlastmxu * 3600.0_ireals / dt)
    nnextmxu = NINT (hnextmxu * 3600.0_ireals / dt)
  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE initialize_loop

!==============================================================================
!==============================================================================

SUBROUTINE exchange_leapfrog

  IF (lprog_qi .AND. lzconv .AND. .NOT. lprogprec) THEN
    kzdims(1:24) =                                                          &
       (/ke,ke,ke,ke,ke1,ke1,ke,ke,ke,ke,ke,ke,                             &
         ke,ke,ke,ke,ke,ke,1,0,0,0,0,0/)
    CALL exchg_boundaries                                                   &
       (nnew+39, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
        ie, je, kzdims, jstartpar, jendpar,                                 &
        nbl_exchg, nboundlines, my_cart_neigh, lperi_x, lperi_y, l2dim,     &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,            &
        u (:,:,:,nnow), u (:,:,:,nnew), v (:,:,:,nnow), v (:,:,:,nnew),     &
        w (:,:,:,nnow), w (:,:,:,nnew), t (:,:,:,nnow), t (:,:,:,nnew),     &
        qv(:,:,:,nnow), qv(:,:,:,nnew), qc(:,:,:,nnow), qc(:,:,:,nnew),     &
        qi(:,:,:,nnow), qi(:,:,:,nnew), pp(:,:,:,nnow), pp(:,:,:,nnew),     &
        qrs(:,:,:)    , dqvdt(:,:,:)  , qvsflx(:,:) )
  ELSEIF (lprog_qi .AND. lzconv .AND. lprogprec) THEN
    kzdims(1:24) =                                                          &
       (/ke,ke,ke,ke,ke1,ke1,ke,ke,ke,ke,ke,ke,                             &
         ke,ke,ke,ke,ke,ke,ke,ke,ke,ke,1,0/)
    CALL exchg_boundaries                                                   &
       (nnew+39, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
        ie, je, kzdims, jstartpar, jendpar,                                 &
        nbl_exchg, nboundlines, my_cart_neigh,                              &
        lperi_x, lperi_y, l2dim,                                            &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,            &
        u (:,:,:,nnow), u (:,:,:,nnew), v (:,:,:,nnow), v (:,:,:,nnew),     &
        w (:,:,:,nnow), w (:,:,:,nnew), t (:,:,:,nnow), t (:,:,:,nnew),     &
        qv(:,:,:,nnow), qv(:,:,:,nnew), qc(:,:,:,nnow), qc(:,:,:,nnew),     &
        qi(:,:,:,nnow), qi(:,:,:,nnew), qr(:,:,:,nnow), qr(:,:,:,nnew),     &
        qs(:,:,:,nnow), qs(:,:,:,nnew), pp(:,:,:,nnow), pp(:,:,:,nnew),     &
        qrs(:,:,:)    , dqvdt(:,:,:)  , qvsflx(:,:) )
    IF (itype_gscp==4) THEN
      kzdims(1:24) =(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                 &
       ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,       &
        ie, je, kzdims, jstartpar, jendpar,                                 &
        nbl_exchg, nboundlines, my_cart_neigh,                              &
        lperi_x, lperi_y, l2dim,                                            &
        20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,               &
        qg (:,:,:,nnow), qg (:,:,:,nnew) )
    ENDIF
  ELSEIF (lprog_qi .AND. .NOT. lzconv .AND. .NOT. lprogprec) THEN
    kzdims(1:24) =                                                          &
       (/ke,ke,ke,ke,ke1,ke1,ke,ke,ke,ke,ke,ke,                             &
         ke,ke,ke,ke,ke,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                   &
       (nnew+36, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
        ie, je, kzdims, jstartpar, jendpar,                                 &
        nbl_exchg, nboundlines, my_cart_neigh,                              &
        lperi_x, lperi_y, l2dim,                                            &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,            &
        u (:,:,:,nnow), u (:,:,:,nnew), v (:,:,:,nnow), v (:,:,:,nnew),     &
        w (:,:,:,nnow), w (:,:,:,nnew), t (:,:,:,nnow), t (:,:,:,nnew),     &
        qv(:,:,:,nnow), qv(:,:,:,nnew), qc(:,:,:,nnow), qc(:,:,:,nnew),     &
        qi(:,:,:,nnow), qi(:,:,:,nnew), pp(:,:,:,nnow), pp(:,:,:,nnew),     &
        qrs(:,:,:)    )
  ELSEIF (lprog_qi .AND. .NOT. lzconv) THEN
    kzdims(1:24) =                                                          &
       (/ke,ke,ke,ke,ke1,ke1,ke,ke,ke,ke,ke,ke,                             &
         ke,ke,ke,ke,ke,ke,ke,ke,ke,0,0,0/)
    CALL exchg_boundaries                                                   &
       (nnew+36, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
        ie, je, kzdims, jstartpar, jendpar,                                 &
        nbl_exchg, nboundlines, my_cart_neigh,                              &
        lperi_x, lperi_y, l2dim,                                            &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,            &
        u (:,:,:,nnow), u (:,:,:,nnew), v (:,:,:,nnow), v (:,:,:,nnew),     &
        w (:,:,:,nnow), w (:,:,:,nnew), t (:,:,:,nnow), t (:,:,:,nnew),     &
        qv(:,:,:,nnow), qv(:,:,:,nnew), qc(:,:,:,nnow), qc(:,:,:,nnew),     &
        qi(:,:,:,nnow), qi(:,:,:,nnew), qr(:,:,:,nnow), qr(:,:,:,nnew),     &
        qs(:,:,:,nnow), qs(:,:,:,nnew), pp(:,:,:,nnow), pp(:,:,:,nnew),     &
        qrs(:,:,:)    )
    IF (itype_gscp==4) THEN
      kzdims(1:24) =(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                 &
       ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,       &
        ie, je, kzdims, jstartpar, jendpar,                                 &
        nbl_exchg, nboundlines, my_cart_neigh,                              &
        lperi_x, lperi_y, l2dim,                                            &
        20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,               &
        qg (:,:,:,nnow), qg (:,:,:,nnew) )
    ENDIF
  ELSEIF (.NOT. lprog_qi .AND. lzconv) THEN
    kzdims(1:24) =                                                          &
       (/ke,ke,ke,ke,ke1,ke1,ke,ke,ke,ke,ke,ke,                             &
         ke,ke,ke,1,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                   &
       (nnew+33, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
        ie, je, kzdims, jstartpar, jendpar,                                 &
        nbl_exchg, nboundlines, my_cart_neigh,                              &
        lperi_x, lperi_y, l2dim,                                            &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,            &
        u (:,:,:,nnow), u (:,:,:,nnew), v (:,:,:,nnow), v (:,:,:,nnew),     &
        w (:,:,:,nnow), w (:,:,:,nnew), t (:,:,:,nnow), t (:,:,:,nnew),     &
        qv(:,:,:,nnow), qv(:,:,:,nnew), qc(:,:,:,nnow), qc(:,:,:,nnew),     &
        pp(:,:,:,nnow), pp(:,:,:,nnew), dqvdt(:,:,:)  , qvsflx(:,:) )
    IF (lprogprec) THEN
      IF (itype_gscp > 1) THEN
        kzdims(1:24) =                                                      &
           (/ke,ke,ke,ke,ke,0,0,0,0,0,0,0,                                  &
             0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                               &
         ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
          ie, je, kzdims, jstartpar, jendpar,                               &
          nbl_exchg, nboundlines, my_cart_neigh,                            &
          lperi_x, lperi_y, l2dim,                                          &
          20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,             &
          qr(:,:,:,nnow), qr(:,:,:,nnew),                                   &
          qs(:,:,:,nnow), qs(:,:,:,nnew), qrs(:,:,:) )
      ELSE
        kzdims(1:24) =                                                      &
           (/ke,ke,ke,0,0,0,0,0,0,0,0,0,                                    &
             0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                               &
         ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
          ie, je, kzdims, jstartpar, jendpar,                               &
          nbl_exchg, nboundlines, my_cart_neigh,                            &
          lperi_x, lperi_y, l2dim,                                          &
          20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,             &
          qr(:,:,:,nnow), qr(:,:,:,nnew), qrs(:,:,:) )
      ENDIF
    ENDIF
  ELSE
    kzdims(1:24) =                                                          &
       (/ke,ke,ke,ke,ke1,ke1,ke,ke,ke,ke,ke,ke,                             &
         ke,ke,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                   &
       (nnew+30, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,  &
        ie, je, kzdims, jstartpar, jendpar,                                 &
        nbl_exchg, nboundlines, my_cart_neigh,                              &
        lperi_x, lperi_y, l2dim,                                            &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,            &
        u (:,:,:,nnow), u (:,:,:,nnew), v (:,:,:,nnow), v (:,:,:,nnew),     &
        w (:,:,:,nnow), w (:,:,:,nnew), t (:,:,:,nnow), t (:,:,:,nnew),     &
        qv(:,:,:,nnow), qv(:,:,:,nnew), qc(:,:,:,nnow), qc(:,:,:,nnew),     &
        pp(:,:,:,nnow), pp(:,:,:,nnew))
    IF (lprogprec) THEN
      IF (itype_gscp > 1) THEN
        kzdims(1:24) =                                                      &
           (/ke,ke,ke,ke,ke,0,0,0,0,0,0,0,                                  &
             0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                               &
         ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
          ie, je, kzdims, jstartpar, jendpar,                               &
          nbl_exchg, nboundlines, my_cart_neigh,                            &
          lperi_x, lperi_y, l2dim,                                          &
          20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,             &
          qr(:,:,:,nnow), qr(:,:,:,nnew),                                   &
          qs(:,:,:,nnow), qs(:,:,:,nnew), qrs(:,:,:) )
      ELSE
        kzdims(1:24) =                                                      &
           (/ke,ke,ke,0,0,0,0,0,0,0,0,0,                                    &
             0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                               &
         ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
          ie, je, kzdims, jstartpar, jendpar,                               &
          nbl_exchg, nboundlines, my_cart_neigh,                            &
          lperi_x, lperi_y, l2dim,                                          &
          20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,             &
          qr (:,:,:,nnow), qr (:,:,:,nnew), qrs(:,:,:) )
      ENDIF
    ENDIF
  ENDIF

#ifdef COSMOART
  IF (l_cosmo_art) THEN
    IF (lgas) THEN
      DO isp = 1,isp_gas
        kzdims(1:24) =                                                      &
           (/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                               &
           ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
            ie, je, kzdims, jstartpar, jendpar,                             &
            nbl_exchg, nboundlines, my_cart_neigh,                          &
            lperi_x, lperi_y, l2dim,                                        &
            20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,           &
            cgas(:,:,:,isp,nnow), cgas(:,:,:,isp,nnew))
      ENDDO
    ENDIF
    IF (laero) THEN
      DO isp = 1,isp_aero
        kzdims(1:24) =                                                      &
           (/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                               &
           ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
            ie, je, kzdims, jstartpar, jendpar,                             &
            nbl_exchg, nboundlines, my_cart_neigh,                          &
            lperi_x, lperi_y, l2dim,                                        &
            20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,           &
            caero(:,:,:,isp,nnow), caero(:,:,:,isp,nnew))
      ENDDO
    ENDIF
  ENDIF
#endif

#ifdef POLLEN
  IF (l_pollen) THEN
    DO isp = 1,isp_pollen
      kzdims(1:24) =                                                        &
         (/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                 &
         ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
          ie, je, kzdims, jstartpar, jendpar,                               &
          nbl_exchg, nboundlines, my_cart_neigh,                            &
          lperi_x, lperi_y, l2dim,                                          &
          20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,             &
          cpollen(:,:,:,isp,nnow), cpollen(:,:,:,isp,nnew))
    ENDDO
  ENDIF
#endif

!AK (20.03.2012)
  DO iig=1, ntracer
    nprog = 0_iintegers
    DO iprog=1, 7
      nprog = nprog + ltracer(iprog,iig)
    ENDDO
    IF (nprog .GE. 1) THEN
      kzdims(1:24)=(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                &
        (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,      &
                                                               ie, je,     &
        kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
        lperi_x, lperi_y, l2dim,                                           &
        20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,              &
        tracer (:,:,:,nnow,iig), tracer (:,:,:,nnew,iig) )
    ENDIF
  ENDDO
!AK (20.03.2012)

END SUBROUTINE exchange_leapfrog

!==============================================================================
!==============================================================================

SUBROUTINE exchange_runge_kutta
  
  IF (lprog_qi) THEN
    IF (lprogprec) THEN
      ! this is former itype_gscp = 5
      IF (itype_gscp == 3) THEN
        kzdims(1:24)=(/ke,ke,ke1,ke,ke,ke,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                  &
         (50+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
          ie, je, kzdims, jstartpar, jendpar,                                  &
          nbl_exchg, nboundlines, my_cart_neigh,                               &
          lperi_x, lperi_y, l2dim,                                             &
          20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,             &
          u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),      &
          qv(:,:,:,nnew), qc(:,:,:,nnew), qi(:,:,:,nnew), qr(:,:,:,nnew),      &
          qs(:,:,:,nnew), pp(:,:,:,nnew), qrs(:,:,:) )
      END IF
      IF (itype_gscp == 4) THEN
        kzdims(1:24)=(/ke,ke,ke1,ke,ke,ke,ke,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                  &
         (50+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
          ie, je, kzdims, jstartpar, jendpar,                                  &
          nbl_exchg, nboundlines, my_cart_neigh,                               &
          lperi_x, lperi_y, l2dim,                                             &
          20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,             &
          u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),      &
          qv(:,:,:,nnew), qc(:,:,:,nnew), qi(:,:,:,nnew), qr(:,:,:,nnew),      &
          qs(:,:,:,nnew), qg(:,:,:,nnew), pp(:,:,:,nnew), qrs(:,:,:) )
      ENDIF
    ELSE ! .NOT. lprogprec:
      ! this is former itype_gscp = 3
      kzdims(1:24)=(/ke,ke,ke1,ke,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
       (50+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,               &
        u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),        &
        qv(:,:,:,nnew), qc(:,:,:,nnew), qi(:,:,:,nnew), pp(:,:,:,nnew),        &
        qrs(:,:,:) )
    ENDIF
  ELSE ! .NOT. lprog_qi:
    IF (lprogprec) THEN
      IF (itype_gscp > 1) THEN        
        ! this is former itype_gscp = 4
        kzdims(1:24)=(/ke,ke,ke1,ke,ke,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                  &
         (50+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
          ie, je, kzdims, jstartpar, jendpar,                                  &
          nbl_exchg, nboundlines, my_cart_neigh,                               &
          lperi_x, lperi_y, l2dim,                                             &
          20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,             &
          u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),      &
          qv(:,:,:,nnew), qc(:,:,:,nnew), qr(:,:,:,nnew), qs(:,:,:,nnew),      &
          pp(:,:,:,nnew), qrs(:,:,:) )
      ELSE ! kessler_pp:
        kzdims(1:24)=(/ke,ke,ke1,ke,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                  &
         (50+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,   &
          ie, je, kzdims, jstartpar, jendpar,                                  &
          nbl_exchg, nboundlines, my_cart_neigh,                               &
          lperi_x, lperi_y, l2dim,                                             &
          20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,             &
          u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),      &
          qv(:,:,:,nnew), qc(:,:,:,nnew), qr(:,:,:,nnew), pp(:,:,:,nnew),      &
          qrs(:,:,:) )          
      ENDIF          
    ELSE ! .NOT. lprogprec:
      kzdims(1:24)=(/ke,ke,ke1,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
       (50+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,               &
        u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),        &
        qv(:,:,:,nnew), qc(:,:,:,nnew), pp(:,:,:,nnew), qrs(:,:,:) )
    ENDIF
  END IF
  
  IF ( lzconv ) THEN
    IF ( lprog_tke ) THEN
      kzdims(1:24)=(/ke,1,ke1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
       (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,           &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,                  &
        dqvdt(:,:,:), qvsflx(:,:), tke(:,:,:,nnew) )
    ELSE
      kzdims(1:24)=(/ke,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
       (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,           &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,                  &
        dqvdt(:,:,:), qvsflx(:,:) )
    END IF
  ELSE
    IF ( lprog_tke ) THEN
      kzdims(1:24)=(/ke1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
       (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,           &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,                  &
        tke(:,:,:,nnew) )
    END IF
  END IF

#ifdef COSMOART
! CK 20101204 boundary exchange also for ART
  IF (l_cosmo_art) THEN
    IF (lgas) THEN
      DO isp = 1,isp_gas
        kzdims(1:24) =                                                         &
           (/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                  &
           (2, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,       &
            ie, je, kzdims, jstartpar, jendpar,                                &
            nbl_exchg, nboundlines, my_cart_neigh,                             &
            lperi_x, lperi_y, l2dim,                                           &
            20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,              &
            cgas(:,:,:,isp,nnew))
      ENDDO
    ENDIF
    IF (laero) THEN
      DO isp = 1,isp_aero
        kzdims(1:24) =                                                         &
           (/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                  &
           (2, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,       &
            ie, je, kzdims, jstartpar, jendpar,                                &
            nbl_exchg, nboundlines, my_cart_neigh,                             &
            lperi_x, lperi_y, l2dim,                                           &
            20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,              &
            caero(:,:,:,isp,nnew))
      ENDDO
    ENDIF
  ENDIF
#endif

#ifdef POLLEN
  IF (l_pollen) THEN
    DO isp = 1,isp_pollen
      kzdims(1:24) =                                                           &
         (/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
         (2, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,         &
          ie, je, kzdims, jstartpar, jendpar,                                  &
          nbl_exchg, nboundlines, my_cart_neigh,                               &
          lperi_x, lperi_y, l2dim,                                             &
          20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,                &
          cpollen(:,:,:,isp,nnew))
    ENDDO
  ENDIF
#endif

!AK (20.03.2012)
  DO iig=1, ntracer
    nprog = 0_iintegers
    DO iprog=1, 7
      nprog = nprog + ltracer(iprog,iig)
    ENDDO
    IF (nprog .GE. 1) THEN
      kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                   &
        (2, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,ie, je,  &
        kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh,    &
        lperi_x, lperi_y, l2dim,                                              &
        20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,                 &
        tracer (:,:,:,nnew,iig) )
    ENDIF
  ENDDO
!AK (20.03.2012)


END SUBROUTINE exchange_runge_kutta

!==============================================================================
!==============================================================================

SUBROUTINE exchange_2timelevel

  IF (lprog_qi .AND. lzconv) THEN
    kzdims(1:24) =(/ke,ke,ke1,ke,ke,ke,ke,ke,ke,ke,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                      &
       (nnew+39, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,               &
        u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),        &
        qv(:,:,:,nnew), qc(:,:,:,nnew), qi(:,:,:,nnew), pp(:,:,:,nnew),        &
        qrs(:,:,:)    , dqvdt(:,:,:)  , qvsflx(:,:) )
  ELSEIF (lprog_qi .AND. .NOT. lzconv) THEN
    kzdims(1:24) =(/ke,ke,ke1,ke,ke,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                      &
       (nnew+36, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,               &
        u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),        &
        qv(:,:,:,nnew), qc(:,:,:,nnew), qi(:,:,:,nnew), pp(:,:,:,nnew),        &
        qrs(:,:,:)    )
  ELSEIF (.NOT. lprog_qi .AND. lzconv) THEN
    kzdims(1:24) =(/ke,ke,ke1,ke,ke,ke,ke,ke,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                      &
       (nnew+33, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,               &
        u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),        &
        qv(:,:,:,nnew), qc(:,:,:,nnew), pp(:,:,:,nnew),                        &
        dqvdt(:,:,:)  , qvsflx(:,:) )
  ELSE
    kzdims(1:24) =(/ke,ke,ke1,ke,ke,ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                      &
       (nnew+30, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,     &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, ldatatypes, ncomm_type, izerror, yzerrmsg,               &
        u (:,:,:,nnew), v (:,:,:,nnew), w (:,:,:,nnew), t (:,:,:,nnew),        &
        qv(:,:,:,nnew), qc(:,:,:,nnew), pp(:,:,:,nnew))
  ENDIF

!AK (20.03.2012)
  DO iig=1, ntracer
    nprog = 0_iintegers
    DO iprog=1, 7
      nprog = nprog + ltracer(iprog,iig)
    ENDDO
    IF (nprog .GE. 1) THEN
      kzdims(1:24)=(/ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                &
        (nnew+30, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,&
        kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh,   &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,           &
        tracer (:,:,:,nnew,iig) )
    ENDIF
  ENDDO
!AK (20.03.2012)

  IF (lprogprec) THEN
    IF (itype_gscp == 4) THEN
      kzdims(1:24) =                                                           &
         (/ke,ke,ke,ke,0,0,0,0,0,0,0,0,                                        &
           0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
       ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,          &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
       20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,                   &
        qr(:,:,:,nnew), qs(:,:,:,nnew), qg(:,:,:,nnew), qrs(:,:,:) )
    ELSEIF (itype_gscp > 1) THEN
      kzdims(1:24) =                                                           &
         (/ke,ke,ke,0,0,0,0,0,0,0,0,0,                                         &
           0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
       ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,          &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
       20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,                   &
        qr(:,:,:,nnew), qs(:,:,:,nnew), qrs(:,:,:) )
    ELSE
      kzdims(1:24) =                                                           &
         (/ke,ke,0,0,0,0,0,0,0,0,0,0,                                          &
           0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                    &
       ( 0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,          &
        ie, je, kzdims, jstartpar, jendpar,                                    &
        nbl_exchg, nboundlines, my_cart_neigh,                                 &
        lperi_x, lperi_y, l2dim,                                               &
        20000+ntstep, .FALSE., ncomm_type, izerror, yzerrmsg,                  &
        qr(:,:,:,nnew), qrs(:,:,:) )
    ENDIF
  ENDIF

END SUBROUTINE exchange_2timelevel

!==============================================================================
!==============================================================================

SUBROUTINE exchange_l2dim

  DO k = 1, ke
    DO j = 1,nboundlines
          
      t (:,jstart-j,k,nnew) = t (:,jend  +1-j,k,nnew)
      pp(:,jstart-j,k,nnew) = pp(:,jend  +1-j,k,nnew)
      qv(:,jstart-j,k,nnew) = qv(:,jend  +1-j,k,nnew)
      qc(:,jstart-j,k,nnew) = qc(:,jend  +1-j,k,nnew)
      qrs(:,jstart-j,k)     = qrs(:,jend  +1-j,k)
      t (:,jend  +j,k,nnew) = t (:,jstart-1+j,k,nnew)
      pp(:,jend  +j,k,nnew) = pp(:,jstart-1+j,k,nnew)
      qv(:,jend  +j,k,nnew) = qv(:,jstart-1+j,k,nnew)
      qc(:,jend  +j,k,nnew) = qc(:,jstart-1+j,k,nnew)
      qrs (:,jend  +j,k)    = qrs (:,jstart-1+j,k)

!AK (20.03.2012)
      DO iig=1, ntracer
        tracer(:,jstart-j,k,nnew,iig) = tracer(:,jend  +1-j,k,nnew,iig)
        tracer(:,jend  +j,k,nnew,iig) = tracer(:,jstart-1+j,k,nnew,iig)
      ENDDO
!AK (20.03.2012)

      IF (lprogprec) THEN
        qr(:,jstart-j,k,nnew) = qr(:,jend  +1-j,k,nnew)
        qr(:,jend  +j,k,nnew) = qr(:,jstart-1+j,k,nnew)
        IF (itype_gscp > 1) THEN
          qs(:,jstart-j,k,nnew) = qs(:,jend  +1-j,k,nnew)
          qs(:,jend  +j,k,nnew) = qs(:,jstart-1+j,k,nnew)
          IF (itype_gscp == 4) THEN
            qg(:,jstart-j,k,nnew) = qg(:,jend  +1-j,k,nnew)
            qg(:,jend  +j,k,nnew) = qg(:,jstart-1+j,k,nnew)
          ENDIF
        ENDIF
      ENDIF
      IF (lprog_qi) THEN
        qi(:,jstart-j,k,nnew) = qi(:,jend  +1-j,k,nnew)
        qi(:,jend  +j,k,nnew) = qi(:,jstart-1+j,k,nnew)
      ENDIF

      IF (lzconv) THEN
        dqvdt(:,jstart-j,k)   = dqvdt(:,jend  +1-j,k)
        dqvdt(:,jend  +j,k)   = dqvdt(:,jstart-1+j,k)
      ENDIF

      u(:,jstartu-j,k,nnew) = u(:,jendu  +1-j,k,nnew)
      v(:,jstartv-j,k,nnew) = v(:,jendv  +1-j,k,nnew)
      u(:,jendu  +j,k,nnew) = u(:,jstartu-1+j,k,nnew)
      v(:,jendv  +j,k,nnew) = v(:,jstartv-1+j,k,nnew)
      
    ENDDO
  ENDDO

  DO k = 1, ke1
    DO j = 1,nboundlines
      w (:,jstart-j,k,nnew) = w (:,jend  +1-j,k,nnew)
      w (:,jend  +j,k,nnew) = w (:,jstart-1+j,k,nnew)
    ENDDO
  ENDDO

  IF (lzconv) THEN
    DO j = 1,nboundlines
      qvsflx (:,jstart-j) = qvsflx (:,jend  +1-j)
      qvsflx (:,jend  +j) = qvsflx (:,jstart-1+j)
    ENDDO
  ENDIF
     
  IF ( .NOT.l2tls ) THEN

    DO k = 1, ke
      DO j = 1,nboundlines
        
        t (:,jstart-j,k,nnow) = t (:,jend  +1-j,k,nnow)
        pp(:,jstart-j,k,nnow) = pp(:,jend  +1-j,k,nnow)
        qv(:,jstart-j,k,nnow) = qv(:,jend  +1-j,k,nnow)
        qc(:,jstart-j,k,nnow) = qc(:,jend  +1-j,k,nnow)
        t (:,jend  +j,k,nnow) = t (:,jstart-1+j,k,nnow)
        pp(:,jend  +j,k,nnow) = pp(:,jstart-1+j,k,nnow)
        qv(:,jend  +j,k,nnow) = qv(:,jstart-1+j,k,nnow)
        qc(:,jend  +j,k,nnow) = qc(:,jstart-1+j,k,nnow)

        IF (lprogprec) THEN
          qr(:,jstart-j,k,nnow) = qr(:,jend  +1-j,k,nnow)
          qr(:,jend  +j,k,nnow) = qr(:,jstart-1+j,k,nnow)
          IF (itype_gscp > 1) THEN
            qs(:,jstart-j,k,nnow) = qs(:,jend  +1-j,k,nnow)
            qs(:,jend  +j,k,nnow) = qs(:,jstart-1+j,k,nnow)
            IF (itype_gscp == 4) THEN
              qg(:,jstart-j,k,nnow) = qg(:,jend  +1-j,k,nnow)
              qg(:,jend  +j,k,nnow) = qg(:,jstart-1+j,k,nnow)
            ENDIF
          ENDIF
        ENDIF
        IF (lprog_qi) THEN
          qi(:,jstart-j,k,nnow) = qi(:,jend  +1-j,k,nnow)
          qi(:,jend  +j,k,nnow) = qi(:,jstart-1+j,k,nnow)
        ENDIF

        u(:,jstartu-j,k,nnow) = u(:,jendu  +1-j,k,nnow)
        v(:,jstartv-j,k,nnow) = v(:,jendv  +1-j,k,nnow)
        u(:,jendu  +j,k,nnow) = u(:,jstartu-1+j,k,nnow)
        v(:,jendv  +j,k,nnow) = v(:,jstartv-1+j,k,nnow)
        
      ENDDO
    ENDDO

    DO k = 1, ke1
      DO j = 1,nboundlines
        w (:,jstart-j,k,nnow) = w (:,jend  +1-j,k,nnow)
        w (:,jend  +j,k,nnow) = w (:,jstart-1+j,k,nnow)
      ENDDO
    ENDDO

  ENDIF
      
END SUBROUTINE exchange_l2dim

!==============================================================================

SUBROUTINE set_qrqsqg_boundaries

  ! Now we have to set the nnew values for qr and qs in a consistent way:
  ! this is an intermediate solution, as long as no better treatment of 
  ! the boundary values is found

  ! Treatment of rain and snow
  ! --------------------------

!US-JF:  IF (.NOT.llb_qr_qs) THEN
    IF     (itype_lbc_qrsg == 1) THEN
      ! set boundary value to first interior row

      ! qr
      IF (ALLOCATED(qr)) THEN
        ! western boundary
        IF (my_cart_neigh(1) == -1) THEN
          DO i = 1, nboundlines
            DO  k = 1, ke
!CDIR NOLOOPCHG
              DO  j = jstart, jend
                qr(i,j,k,nnew) = qr(istart,j,k,nnew)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! eastern boundary
        IF (my_cart_neigh(3) == -1) THEN
          DO i = ie-nboundlines+1, ie
            DO  k = 1, ke
!CDIR NOLOOPCHG
              DO  j = jstart, jend
                qr(i,j,k,nnew) = qr(iend  ,j,k,nnew)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! southern boundary
        IF (my_cart_neigh(4) == -1) THEN
          DO  k = 1, ke
            DO  j = 1, nboundlines
              qr(:,j,k,nnew) = qr(:,jstart,k,nnew)
            ENDDO
          ENDDO
        ENDIF
        ! northern boundary
        IF (my_cart_neigh(2) == -1) THEN
          DO  k = 1, ke
            DO  j = je-nboundlines+1, je
              qr(:,j,k,nnew) = qr(:,jend  ,k,nnew)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      ! qs
      IF (ALLOCATED(qs)) THEN
        ! western boundary
        IF (my_cart_neigh(1) == -1) THEN
          DO  k = 1, ke
            DO i = 1, nboundlines
!CDIR NOLOOPCHG
              DO  j = jstart, jend
                qs(i,j,k,nnew) = qs(istart,j,k,nnew)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! eastern boundary
        IF (my_cart_neigh(3) == -1) THEN
          DO  k = 1, ke
            DO i = ie-nboundlines+1, ie
!CDIR NOLOOPCHG
              DO  j = jstart, jend
                qs(i,j,k,nnew) = qs(iend  ,j,k,nnew)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! southern boundary
        IF (my_cart_neigh(4) == -1) THEN
          DO  k = 1, ke
            DO  j = 1, nboundlines
              qs(:,j,k,nnew) = qs(:,jstart,k,nnew)
            ENDDO
          ENDDO
        ENDIF
        ! northern boundary
        IF (my_cart_neigh(2) == -1) THEN
          DO  k = 1, ke
            DO  j = je-nboundlines+1, je
              qs(:,j,k,nnew) = qs(:,jend  ,k,nnew)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

    ELSEIF (itype_lbc_qrsg == 2) THEN

      ! set all values to 0

      ! qr
      IF (ALLOCATED(qr)) THEN
        ! western boundary
        IF (my_cart_neigh(1) == -1) THEN
          DO i = 1, nboundlines
            DO  k = 1, ke
!CDIR NOLOOPCHG
              DO  j = jstart, jend
                qr(i,j,k,nnew) = 0.0_ireals
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! eastern boundary
        IF (my_cart_neigh(3) == -1) THEN
          DO i = ie-nboundlines+1, ie
            DO  k = 1, ke
!CDIR NOLOOPCHG
              DO  j = jstart, jend
                qr(i,j,k,nnew) = 0.0_ireals
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! southern boundary
        IF (my_cart_neigh(4) == -1) THEN
          DO  k = 1, ke
            DO  j = 1, nboundlines
              qr(:,j,k,nnew) = 0.0_ireals
            ENDDO
          ENDDO
        ENDIF
        ! northern boundary
        IF (my_cart_neigh(2) == -1) THEN
          DO  k = 1, ke
            DO  j = je-nboundlines+1, je
              qr(:,j,k,nnew) = 0.0_ireals
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      ! qs
      IF (ALLOCATED(qs)) THEN
        ! western boundary
        IF (my_cart_neigh(1) == -1) THEN
          DO  k = 1, ke
            DO i = 1, nboundlines
!CDIR NOLOOPCHG
              DO  j = jstart, jend
                qs(i,j,k,nnew) = 0.0_ireals
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! eastern boundary
        IF (my_cart_neigh(3) == -1) THEN
          DO  k = 1, ke
            DO i = ie-nboundlines+1, ie
!CDIR NOLOOPCHG
              DO  j = jstart, jend
                qs(i,j,k,nnew) = 0.0_ireals
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! southern boundary
        IF (my_cart_neigh(4) == -1) THEN
          DO  k = 1, ke
            DO  j = 1, nboundlines
              qs(:,j,k,nnew) = 0.0_ireals
            ENDDO
          ENDDO
        ENDIF
        ! northern boundary
        IF (my_cart_neigh(2) == -1) THEN
          DO  k = 1, ke
            DO  j = je-nboundlines+1, je
              qs(:,j,k,nnew) = 0.0_ireals
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!   ELSEIF (itype_lbc_qrsg == 3) THEN
!     ! nothing is done for this option
    ENDIF
!US-JF:  ENDIF      ! .NOT.llb_qr_qs

  ! Treatment of graupel
  ! --------------------

!US-JF:  IF (.NOT.llb_qg) THEN
    IF (ALLOCATED(qg)) THEN
      IF     (itype_lbc_qrsg == 1) THEN

        ! western boundary
        IF (my_cart_neigh(1) == -1) THEN
          DO  k = 1, ke
            DO i = 1, nboundlines
!CDIR NOLOOPCHG 
              DO  j = jstart, jend
                qg(i,j,k,nnew) = qg(istart,j,k,nnew)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! eastern boundary 
        IF (my_cart_neigh(3) == -1) THEN
          DO  k = 1, ke
            DO i = ie-nboundlines+1, ie
!CDIR NOLOOPCHG
              DO  j = jstart, jend
                qg(i,j,k,nnew) = qg(iend  ,j,k,nnew)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! southern boundary
        IF (my_cart_neigh(4) == -1) THEN
          DO  k = 1, ke
            DO  j = 1, nboundlines
              qg(:,j,k,nnew) = qg(:,jstart,k,nnew)
            ENDDO
          ENDDO
        ENDIF
          ! northern boundary
        IF (my_cart_neigh(2) == -1) THEN
          DO  k = 1, ke
            DO  j = je-nboundlines+1, je
              qg(:,j,k,nnew) = qg(:,jend  ,k,nnew)
            ENDDO
          ENDDO
        ENDIF

      ELSEIF (itype_lbc_qrsg == 2) THEN

        ! set all values to 0
        ! western boundary
        IF (my_cart_neigh(1) == -1) THEN
          DO  k = 1, ke
            DO i = 1, nboundlines
!CDIR NOLOOPCHG 
              DO  j = jstart, jend
                qg(i,j,k,nnew) = 0.0_ireals
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! eastern boundary 
        IF (my_cart_neigh(3) == -1) THEN
          DO  k = 1, ke
            DO i = ie-nboundlines+1, ie
!CDIR NOLOOPCHG
              DO  j = jstart, jend
                qg(i,j,k,nnew) = 0.0_ireals
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! southern boundary
        IF (my_cart_neigh(4) == -1) THEN
          DO  k = 1, ke
            DO  j = 1, nboundlines
              qg(:,j,k,nnew) = 0.0_ireals
            ENDDO
          ENDDO
        ENDIF
        ! northern boundary
        IF (my_cart_neigh(2) == -1) THEN
          DO  k = 1, ke
            DO  j = je-nboundlines+1, je
              qg(:,j,k,nnew) = 0.0_ireals
            ENDDO
          ENDDO
        ENDIF

!     ELSEIF (itype_lbc_qrsg == 3) THEN
!       ! nothing has to be done for this option
      ENDIF
    ENDIF
!US-JF:  ENDIF

END SUBROUTINE set_qrqsqg_boundaries

!==============================================================================

SUBROUTINE nullify_tracers

    ! If the humidity variables are rather small, set them to 0
    IF (ALLOCATED(qr)) THEN
      WHERE (qr (:,:,:,nnew) < 1.0E-15_ireals)
        qr (:,:,:,nnew) = 0.0_ireals
      ENDWHERE
    ENDIF
    IF (ALLOCATED(qs)) THEN
      WHERE (qs (:,:,:,nnew) < 1.0E-15_ireals)
        qs (:,:,:,nnew) = 0.0_ireals
      ENDWHERE
    ENDIF
    IF (ALLOCATED(qg)) THEN
      WHERE (qg (:,:,:,nnew) < 1.0E-15_ireals)
        qg (:,:,:,nnew) = 0.0_ireals
      ENDWHERE
    ENDIF
    IF (ALLOCATED(qi)) THEN
      WHERE (qi (:,:,:,nnew) < 1.0E-15_ireals)
        qi (:,:,:,nnew) = 0.0_ireals
      ENDWHERE
    ENDIF
    WHERE (qc (:,:,:,nnew) < 1.0E-15_ireals)
      qc (:,:,:,nnew) = 0.0_ireals
    ENDWHERE
    WHERE (qv (:,:,:,nnew) < 1.0E-15_ireals)
      qv (:,:,:,nnew) = 0.0_ireals
    ENDWHERE
    WHERE (qrs(:,:,:) < 1.0E-15_ireals)
      qrs(:,:,:) = 0.0_ireals
    ENDWHERE

#ifdef POLLEN
    IF (l_pollen) THEN
      WHERE (cpollen(:,:,:,:,nnew) < 1.0E-15_ireals)
        cpollen(:,:,:,:,nnew) = 1.0E-14_ireals
      ENDWHERE
    ENDIF
#endif

END SUBROUTINE nullify_tracers


end module enkf_cosmo_mod
