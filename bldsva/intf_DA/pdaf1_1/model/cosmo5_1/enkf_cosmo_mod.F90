!+ Main program for the LM
!------------------------------------------------------------------------------

PROGRAM lmorg

!------------------------------------------------------------------------------
!
! Description:
!  Lmorg is the main program that organizes the forecast done by the LM.
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Ulrich Schaettler
!  Initial release
! 1.2        1998/03/30 Ulrich Schaettler
!  Introduction of digital filtering for data initialization
! 1.4        1998/05/22 Guenther Doms
!  Adaptions for the two time-level time integration scheme
! 1.5        1998/06/29 Guenther Doms
!  Call of new routine 'print_vertcoord' for printing the vertical coordinates
! 1.8        1998/08/03 Ulrich Schaettler
!  Added parameterlist to init_meanvalues
! 1.9        1998/09/16 Guenther Doms
!  Call of routine 'vertical_diffusion' in case of two_time_level integration.
! 1.10       1998/09/29 Ulrich Schaettler
!  Change parameterlist of routine collect_timings.
! 1.12       1998/10/19 Ulrich Schaettler
!  Introduced possibility of boundary values from different runs.
! 1.13       1998/10/22 Christoph Schraff
!  Introduction of data assimilation by the nudging method.
! 1.17       1998/11/17 Ulrich Schaettler
!  Writing ready files after output steps, if required.
! 1.19       1998/12/11 Christoph Schraff
!  Bug fix for writing ready files for analyses.
! 1.20       1999/01/07 Guenhter Doms
!  Renaming of some global variables.
! 1.22       1999/02/08 Michael Buchhold
!  Introduction of 2D analyses of near surface parameters (e.g. t2m, rh2m)
! 1.24       1999/03/01 Guenther Doms
!  Call of the cloud ice parameterization schema (optionally).
! 1.29       1999/05/11 Ulrich Schaettler
!  Adapted interfaces to utility-modules and excluded several internal 
!  procedures
! 1.33       1999/10/14 Reinhold Hess
!  Call to budget diagnosis 'organize_diabudget', additional use of ldiag.
! 1.34       1999/12/10 Ulrich Schaettler
!  Introduced organization routines to most components which changes most
!  interfaces.
! 1.39       2000/05/03 Ulrich Schaettler
!  Adaptations to include asynchronous IO and the interactive nesting.
! 1.40       2000/05/23 Ulrich Schaettler
!  Correction for using MPE_IO on sequential platforms.
! 2.6        2001/06/12 Guenther Doms
!  Correction in an data exchange related to the cloud-ice scheme
! 2.8        2001/07/06 Ulrich Schaettler
!  Eliminated non-necessary variables from the USE-lists
!  Inquire the ASCII files now by UNIT-number instead of filename
! 2.11       2001/09/28 Ulrich Schaettler
!  Changed treatment of cloud ice, if initial and boundary fields are provided
! 2.17       2002/05/08 Ulrich Schaettler
!  Initializations of new time level for multi-layer soil model variables
! 2.18       2002/07/16 Ulrich Schaettler
!  Eliminated variable rhde, use cf_snow instead (by Reinhold Schrodin)
!  Exchange of variables qr, qs for 2 time level scheme (by Almut Gassmann)
!  Changed initialization of new time level, if working with boundary data
!  defined on frames (by Lucio Torrisi).
!  Use lmulti_layer to define boundary values for old soil model correctly
! 2.19       2002/10/24 Ulrich Schaettler
!  Moved re-initialization of vmax_10m from near_surface to initialize_loop.
! 3.2        2003/02/07 Ulrich Schaettler
!  Moved some communications from the Convection scheme to lmorg.
! 3.5        2003/09/02 Ulrich Schaettler
!  Adaptations for the interface of exchg_boundaries
!  Corrections for the digital filtering
! 3.6        2003/12/11 Ulrich Schaettler
!  Changed interface for calling the digital filtering dfi_initialization
! 3.7        2004/02/18 Ulrich Schaettler
!  Introduced CALLs to new component for computing synthetic satellite images.
!  Adaptations in the boundary exchange for prognostic precipitation.
!  Treatment of humidity-variables (are set to 0.0, if they are too small)
! 3.8        2004/03/23 Ulrich Schaettler
!  Correction in initialize_loop for the 2 timelevel scheme.
! 3.11       2004/07/21 Michael Gertz
!  FPE trapping for IBM machines introduced
! 3.13       2004/12/03 Ulrich Schaettler
!  Introduced CALL for Latent Heat Nudging (LHN) within organize_assimilation
!                                                      (Klaus Stephan)
!  Adaptations in boundary exchange for new graupel scheme; 
!  Modification for 2D runs (l2dim) (Thorsten Reinhardt)
! 3.14       2005/01/25 Ulrich Schaettler
!  Adapted treatment of t_s for new multi-layer soil model 
! 3.15       2005/03/03 Ulrich Schaettler
!  Replaced FLOAT by REAL
! 3.16       2005/07/22 Ulrich Schaettler
!  Adaptations for the 2 time level Runge-Kutta scheme
! 3.18       2006/03/03 Ulrich Schaettler
!  Changed treatment of ASCII files to introduce restart possibility
!  Additional boundary variables in case of lbdclim=.TRUE.
!  Added get_timing routine for LHN   (Klaus Stephan)
!  Changed computation of surface specific humidity for lakes (D. Mironov)
!  Added treatment of vertical velocity w for lw_freeslip=.FALSE.
! 3.21       2006/12/04 Ulrich Schaettler, Jochen Foerstner
!  Read Namelist Input for physics in any case (also for "dry" runs)
!  Read Namelist Input for EPS mode, if leps=.TRUE. 
!  Call near_surface with time level nnow
!  Set new time level for qr, qs, qg 
! V3_23        2007/03/30 Ulrich Schaettler
!  Adapting print outs to idbg_level and flexible dt
! V3_24        2007/04/26 Ulrich Schaettler
!  Bug correction in call to dfi_initialization
!  Eliminated nincmxu and introduced control as for other increments
! V4_3         2008/02/25 Ulrich Schaettler 
!  Correction for computing zforecasttime in case of restart runs
!  Adapted call to init_timings
! V4_4         2008/07/16 Ulrich Schaettler, Dmitrii Mironov
!  Use qitens so that cloud ice can be treated in the same way as cloud water (DM)
!  Eliminated timing variables which are unused
!  Adapted interface of get_timings
!  Changed NL parameter lyear_360 to itype_calendar, to have several options
! V4_8         2009/02/16 Ulrich Schaettler
!  Renamed switch for artificial data from lgen to lartif_data
!  Reset also variables for dynamical and convective gusts
! V4_9         2009/07/16 Ulrich Schaettler, Heike Vogel
!  Implement a first version of COSMO_ART
!  Introduced internal subroutine set_qrqsqg_boundaries. Set these boundaries
!   for the total field to 0, if no boundary fields are given.
!  Replace hardcoded 0 in call to get_timings for cleanup by nstop
!  Read diagnosis namelist even if ldiagnos=.FALSE.; (some settings are needed)
!  Inserted NEC Compiler directives
!  Eliminate reduction of qv in case lana_qi, llb_qi = .TRUE.
! V4_11        2009/11/30 Ekaterina Machulskaya
!  Adaptations to run with the multi-layer snow model
! V4_12        2010/05/11 Jan-Peter Schulz
!  Adapted interface parameter t_s for SR tgcom to account for seaice (JPS)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
!  New information module info_lm_f90 added
! V4_14        2010/06/14 Ulrich Schaettler
!  Putting the call for the info-module to the right place
! V4_17        2011/02/24 Ulrich Blahak
!  Adapted interface of exchg_boundaries; corrected kzdims(1:20) -> kzdims(1:24);
!  eliminated my_peri_neigh; added calls for artificial disturbances on
!  initial temperature (atmosphere and/or soil) or artificial heating rates
!  (atmosphere and/or soil) or artificial t_s disturbances in case of lsoil=.false.;
! V4_18        2011/05/26 Ulrich Schaettler
!  Introduced conditional compilation for Nudging and synthetic satellite images
!  from Christoph Knote (for COSMO-ART)
!  lphot removed, not necessary anymore
!  lho, lho2 used instead of hardcoded cgas species index
!  ART boundaries use and calculate their own relaxation weights
! V4_20        2011/08/31 Ulrich Schaettler
!  Implemented interface to OASIS coupler using conditional compilation with -DCOUP_OAS
!   (by CLM Community)
! V4_21        2011/12/06 Oliver Fuhrer
!  Bug Fix in call to exchg_boundaries for POLLEN runs
! V4_23        2012/05/10 Ulrich Schaettler, Oli Fuhrer, Burkhardt Rockel
!                         Lucio Torrisi
!  Removed src_2timelevel and related exchange routine
!  Removed switch lprogprec
!  Removed SR exchange_l2dim (this is now done within exchg_boundaries)
!  Removed qvt_diff (Oli Fuhrer)
!  Splitted call to organize_data 'input-init' into two calls for preparation
!    of the new tracer module (Oli Fuhrer)
!  Introduced time step increment nincsn for calling spectral nudging (BR)
!  Initialize t_s (new) with boundary values depending on lbdsst (LT)
!  Implementation of time dependent boundary values for aerosol optical depths
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Ulrich Blahak
!                         Florian Prill, Hans-Juergen Panitz, Carlos Osuna
!  Introduced the tracer module (AR, OF)
!  Introduced 2-moment microphysics: define and set pointer for additional
!    tracer in the 2-moment scheme (in SR initialize_loop) (UB)
!  Modified interface to init_timings (UB)
!  Added exchange in case of periodic BCs right before
!   time stepping. This is done for safety, otherwise some fields
!   may not be periodic in the first timestep, if the user does specialized
!   things in src_artifdata.f90.  (UB)
!  Adapted interfaces to SR mpe_io_node, mpe_io_shutdown from mpe_io2 (FP)
!  Introduced nexch_tag for MPI boundary exchange tag to replace ntstep (HJP)
! Carlos Osuna:
!  Introduce fork in those PE that will act as I/O PE for asynchronous netcdf.
!  In compute PE prepare buffers used to send data to I/O PE.
!  In I/O PE control the start and shutdown of I/O PE.
! V4_26        2012/12/06 Andreas Messer, Anne Roches
!  Extensions to read satpp-files if model is compiled with RTTOV10 (AM)
!  Changes and technical adaptations to the tracer handling (by Anne Roches)
! V4_27        2013/03/19 Michael Baldauf, Astrid Kerkweg, Ulrich Schaettler
!  Introduced another call to organize_dynamics for cleanup and deallocation (MB)
!  MESSy interface introduced (AK)
! V4_28        2013/07/12 Ulrich Schaettler, KIT, Hans-Juergen Panitz (CLM)
!  Adapted COSMO-ART to new tracer version
!  Eliminated reference to ivctype, vcoord (not used) (US)
!  Reinitialization of vabsmx_10m (which has been forgotten before, HJP)
! V4_29        2013-10-02 Astrid Kerkweg, Ulrich Blahak
!  Unification of MESSy interfaces and COSMO Tracer structure
!  Moved call to trcr_alloc at the end of the allocation block
!  Added some missing initializations for leapfrog integration for lartif_data and
!   2-moment-scheme (UB)
! V4_30        2013/11/08 Ulrich Schaettler
!  Changed use of datatypes to FALSE in exchange_leapfrog (had problems on the IBM)
!  Some technical changes for COSMOART and POLLEN
! V5_1         2014-11-28 Ulrich Schaettler, Damian Wojcik, Ulrich Blahak, Oliver Fuhrer
!                         Michael Baldauf, Anne Roches, Xavier Lapillonne, Lucio Torrisi
!  Print information for tracers only if idbg_level > 2
!  Embraced USE data_satellites with if defined RTTOVx, because of unresolved
!  references when compiling purpar, nudpar (US)
!  Changed computation of nexch_tag by using INT(24*3600/dt) as second argument
!   for the MOD function: then also dt < 1 will work (DW)
!  Use the correct timelevel to exchange TKE in case it is advected (UB)
!  Added interface to radar forward operator and respective timer calls.
!   (luse_radarfwo, organize_radar(), -DRADARFWO).
!  Replaced ireals by wp (working precision) (OF)
!  Added check for NaNs in every timestep, but only for idbg_level > 2 (OF)
!  Reshuffle time indices now in lmorg, not in initialize_loop (MB)
!  Removed BD_SET_FORCED (forced boundary conditions for all tracers) in 
!    SR set_trcr_special_bc (AR)
!  Implemented calls to Online Trajectory Module (AR)
!  Implemented block data structure: registration, initialization and clean-up
!    for blocked data fields (XL)
!  Calls added for stochastic perturbation of physics tendencies (SPPT) (LT)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Declarations:
! Modules used:

! Load the library information data:
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
       nexch_tag

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

!==============================================================================

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

comp_pe: IF (lcompute_pe) THEN

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

!------------------------------------------------------------------------------
!- Section 6: Time stepping
!------------------------------------------------------------------------------

  IF (izdebug > 0) THEN
    PRINT *, '  TIME STEPPING'
  ENDIF

  timeloop: DO ntstep = nstart , nstop

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
    PRINT *, 'END OF TIME STEPPING'
  ENDIF

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

ELSE comp_pe

  IF( nc_asyn_io > 0 ) THEN
#ifdef NETCDF
    CALL start_ionode( yzerrmsg, izerror)
    IF( izerror /= 0 ) THEN
      CALL model_abort(my_cart_id, izerror, yzerrmsg,  &
              'start_ionode')
    ENDIF
    CALL shutdown_io()
#endif
  ELSE
   CALL mpe_io_node(izerror)
   IF (izerror /= 0) THEN
     ierrstat = 1015
     yzerrmsg  = ' ERROR    *** Running the asynchronous I/O failed ***'
     CALL model_abort (my_cart_id, ierrstat, yzerrmsg, 'start mpe_io_node')
   ENDIF
 ENDIF

ENDIF comp_pe

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

!==============================================================================
! End of the program
!==============================================================================

END PROGRAM lmorg

!------------------------------------------------------------------------------
