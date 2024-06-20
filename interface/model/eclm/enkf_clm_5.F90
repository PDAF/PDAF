!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz and Guowei He (Forschungszentrum Juelich GmbH)
!
!This file is part of TSMP-PDAF
!
!TSMP-PDAF is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TSMP-PDAF is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU LesserGeneral Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public License
!along with TSMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------
!enkf_clm_5.F90: Wrapper functions for CLM 5 - based on
! cime/src/mct/main/cime_driver.F90
!-------------------------------------------------------------------------------------------
module enkf_clm_5

  use shr_kind_mod,      only: r8 => SHR_KIND_R8
  use shr_kind_mod,      only: i8 => SHR_KIND_I8
  use shr_kind_mod,      only: cs => SHR_KIND_CS
  use shr_kind_mod,      only: cl => SHR_KIND_CL
  use shr_sys_mod,       only: shr_sys_abort, shr_sys_flush
  use shr_const_mod,     only: shr_const_cday
  use shr_file_mod,      only: shr_file_setLogLevel, shr_file_setLogUnit
  use shr_file_mod,      only: shr_file_setIO, shr_file_getUnit, shr_file_freeUnit
  use shr_scam_mod,      only: shr_scam_checkSurface
  use shr_map_mod,       only: shr_map_setDopole
  use shr_mpi_mod,       only: shr_mpi_min, shr_mpi_max
  use shr_mpi_mod,       only: shr_mpi_bcast, shr_mpi_commrank, shr_mpi_commsize
  use shr_mem_mod,       only: shr_mem_init, shr_mem_getusage
  use shr_cal_mod,       only: shr_cal_date2ymd, shr_cal_ymd2date, shr_cal_advdateInt
  use shr_cal_mod,       only: shr_cal_ymds2rday_offset
  use shr_orb_mod,       only: shr_orb_params
  use shr_frz_mod,       only: shr_frz_freezetemp_init
  use shr_reprosum_mod,  only: shr_reprosum_setopts
  use mct_mod            ! mct_ wrappers for mct lib
  use perf_mod
  use ESMF
  use cime_comp_mod 
  !----------------------------------------------------------------------------
  ! component model interfaces (init, run, final methods)
  !----------------------------------------------------------------------------
  
  use atm_comp_mct  , only: atm_init=>atm_init_mct, atm_run=>atm_run_mct, atm_final=>atm_final_mct
  use lnd_comp_mct  , only: lnd_init=>lnd_init_mct, lnd_run=>lnd_run_mct, lnd_final=>lnd_final_mct
  use ocn_comp_mct  , only: ocn_init=>ocn_init_mct, ocn_run=>ocn_run_mct, ocn_final=>ocn_final_mct
  use ice_comp_mct  , only: ice_init=>ice_init_mct, ice_run=>ice_run_mct, ice_final=>ice_final_mct
  use glc_comp_mct  , only: glc_init=>glc_init_mct, glc_run=>glc_run_mct, glc_final=>glc_final_mct
  use wav_comp_mct  , only: wav_init=>wav_init_mct, wav_run=>wav_run_mct, wav_final=>wav_final_mct
  use rof_comp_mct  , only: rof_init=>rof_init_mct, rof_run=>rof_run_mct, rof_final=>rof_final_mct
  use esp_comp_mct  , only: esp_init=>esp_init_mct, esp_run=>esp_run_mct, esp_final=>esp_final_mct
  
  !----------------------------------------------------------------------------
  ! cpl7 modules
  !----------------------------------------------------------------------------
  
  ! mpi comm data & routines, plus logunit and loglevel
  use seq_comm_mct, only: CPLID, GLOID, logunit, loglevel
  use seq_comm_mct, only: ATMID, LNDID, OCNID, ICEID, GLCID, ROFID, WAVID, ESPID
  use seq_comm_mct, only: ALLATMID,ALLLNDID,ALLOCNID,ALLICEID,ALLGLCID,ALLROFID,ALLWAVID,ALLESPID
  use seq_comm_mct, only: CPLALLATMID,CPLALLLNDID,CPLALLOCNID,CPLALLICEID
  use seq_comm_mct, only: CPLALLGLCID,CPLALLROFID,CPLALLWAVID,CPLALLESPID
  use seq_comm_mct, only: CPLATMID,CPLLNDID,CPLOCNID,CPLICEID,CPLGLCID,CPLROFID,CPLWAVID,CPLESPID
  use seq_comm_mct, only: num_inst_atm, num_inst_lnd, num_inst_rof
  use seq_comm_mct, only: num_inst_ocn, num_inst_ice, num_inst_glc
  use seq_comm_mct, only: num_inst_wav, num_inst_esp
  use seq_comm_mct, only: num_inst_xao, num_inst_frc, num_inst_phys
  use seq_comm_mct, only: num_inst_total, num_inst_max
  use seq_comm_mct, only: seq_comm_iamin, seq_comm_name, seq_comm_namelen
  use seq_comm_mct, only: seq_comm_init, seq_comm_setnthreads, seq_comm_getnthreads
  use seq_comm_mct, only: seq_comm_getinfo => seq_comm_setptrs
  use seq_comm_mct, only: cpl_inst_tag
  
  ! clock & alarm routines and variables
  use seq_timemgr_mod, only: seq_timemgr_type
  use seq_timemgr_mod, only: seq_timemgr_clockInit
  use seq_timemgr_mod, only: seq_timemgr_clockAdvance
  use seq_timemgr_mod, only: seq_timemgr_clockPrint
  use seq_timemgr_mod, only: seq_timemgr_EClockGetData
  use seq_timemgr_mod, only: seq_timemgr_alarmIsOn
  use seq_timemgr_mod, only: seq_timemgr_histavg_type
  use seq_timemgr_mod, only: seq_timemgr_type_never
  use seq_timemgr_mod, only: seq_timemgr_alarm_restart
  use seq_timemgr_mod, only: seq_timemgr_alarm_stop
  use seq_timemgr_mod, only: seq_timemgr_alarm_datestop
  use seq_timemgr_mod, only: seq_timemgr_alarm_history
  use seq_timemgr_mod, only: seq_timemgr_alarm_atmrun
  use seq_timemgr_mod, only: seq_timemgr_alarm_lndrun
  use seq_timemgr_mod, only: seq_timemgr_alarm_ocnrun
  use seq_timemgr_mod, only: seq_timemgr_alarm_icerun
  use seq_timemgr_mod, only: seq_timemgr_alarm_glcrun
  use seq_timemgr_mod, only: seq_timemgr_alarm_glcrun_avg
  use seq_timemgr_mod, only: seq_timemgr_alarm_ocnnext
  use seq_timemgr_mod, only: seq_timemgr_alarm_tprof
  use seq_timemgr_mod, only: seq_timemgr_alarm_histavg
  use seq_timemgr_mod, only: seq_timemgr_alarm_rofrun
  use seq_timemgr_mod, only: seq_timemgr_alarm_wavrun
  use seq_timemgr_mod, only: seq_timemgr_alarm_esprun
  use seq_timemgr_mod, only: seq_timemgr_alarm_barrier
  use seq_timemgr_mod, only: seq_timemgr_alarm_pause
  use seq_timemgr_mod, only: seq_timemgr_pause_active
  use seq_timemgr_mod, only: seq_timemgr_pause_component_active
  use seq_timemgr_mod, only: seq_timemgr_pause_component_index
  
  ! "infodata" gathers various control flags into one datatype
  use seq_infodata_mod, only: seq_infodata_putData, seq_infodata_GetData
  use seq_infodata_mod, only: seq_infodata_init, seq_infodata_exchange
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_orb_variable_year
  use seq_infodata_mod, only: seq_infodata_print, seq_infodata_init2
  
  ! domain related routines
  use seq_domain_mct, only : seq_domain_check
  
  ! history file routines
  use seq_hist_mod, only : seq_hist_write, seq_hist_writeavg, seq_hist_writeaux
  
  ! restart file routines
  use seq_rest_mod, only : seq_rest_read, seq_rest_write
  
  ! flux calc routines
  use seq_flux_mct, only: seq_flux_init_mct, seq_flux_initexch_mct, seq_flux_ocnalb_mct
  use seq_flux_mct, only: seq_flux_atmocn_mct, seq_flux_atmocnexch_mct, seq_flux_readnl_mct
  
  ! domain fraction routines
  use seq_frac_mct, only : seq_frac_init, seq_frac_set
  
  ! i/o subroutines
  use seq_io_mod, only : seq_io_cpl_init
  
  ! rearrange type routines
  use cplcomp_exchange_mod, only: seq_mctext_decomp
  
  ! diagnostic routines
  use seq_diag_mct, only : seq_diag_zero_mct, seq_diag_avect_mct, seq_diag_lnd_mct
  use seq_diag_mct, only : seq_diag_rof_mct, seq_diag_ocn_mct, seq_diag_atm_mct
  use seq_diag_mct, only : seq_diag_ice_mct, seq_diag_accum_mct, seq_diag_print_mct
  
  ! list of fields transferred between components
  use seq_flds_mod, only : seq_flds_a2x_fluxes, seq_flds_x2a_fluxes
  use seq_flds_mod, only : seq_flds_i2x_fluxes, seq_flds_x2i_fluxes
  use seq_flds_mod, only : seq_flds_l2x_fluxes, seq_flds_x2l_fluxes
  use seq_flds_mod, only : seq_flds_o2x_fluxes, seq_flds_x2o_fluxes
  use seq_flds_mod, only : seq_flds_g2x_fluxes, seq_flds_x2g_fluxes
  use seq_flds_mod, only : seq_flds_w2x_fluxes, seq_flds_x2w_fluxes
  use seq_flds_mod, only : seq_flds_r2x_fluxes, seq_flds_x2r_fluxes
  use seq_flds_mod, only : seq_flds_set
  
  ! component type and accessor functions
  use component_type_mod , only: component_get_iamin_compid, component_get_suffix
  use component_type_mod , only: component_get_name, component_get_c2x_cx
  use component_type_mod , only: atm, lnd, ice, ocn, rof, glc, wav, esp
  use component_mod      , only: component_init_pre
  use component_mod      , only: component_init_cc, component_init_cx, component_run, component_final
  use component_mod      , only: component_init_areacor, component_init_aream
  use component_mod      , only: component_exch, component_diag
  
  ! prep routines (includes mapping routines between components and merging
  ! routines)
  use prep_lnd_mod
  use prep_ice_mod
  use prep_wav_mod
  use prep_rof_mod
  use prep_glc_mod
  use prep_ocn_mod
  use prep_atm_mod
  use prep_aoflux_mod
  
  !--- mapping routines ---
  use seq_map_type_mod
  use seq_map_mod      ! generic mapping
  
  ! --- timing routines ---
  use t_drv_timers_mod

!#include <mpif.h>  
!  !----------------------------------------------------------------------------
!  ! temporary variables
!  !----------------------------------------------------------------------------
!
!  !- from prep routines (arrays of instances)
!  type(mct_aVect) , pointer :: a2x_ox(:) => null()
!  type(mct_aVect) , pointer :: o2x_ax(:) => null()
!  type(mct_aVect) , pointer :: xao_ox(:) => null()
!  type(mct_aVect) , pointer :: xao_ax(:) => null()
!
!  !- from component type (single instance inside array of components)
!  type(mct_aVect) , pointer :: o2x_ox => null()
!  type(mct_aVect) , pointer :: a2x_ax => null()
!
!  character(len=CL) :: inst_suffix
!  logical           :: iamin_id
!  character(len=seq_comm_namelen) :: compname
!
!  !----------------------------------------------------------------------------
!  ! domains & related
!  !----------------------------------------------------------------------------
!
!  !--- domain fractions (only defined on cpl pes) ---
!  type(mct_aVect) , pointer :: fractions_ax(:)   ! Fractions on atm grid, cpl processes
!  type(mct_aVect) , pointer :: fractions_lx(:)   ! Fractions on lnd grid, cpl processes
!  type(mct_aVect) , pointer :: fractions_ix(:)   ! Fractions on ice grid, cpl processes
!  type(mct_aVect) , pointer :: fractions_ox(:)   ! Fractions on ocn grid, cpl processes
!  type(mct_aVect) , pointer :: fractions_gx(:)   ! Fractions on glc grid, cpl processes
!  type(mct_aVect) , pointer :: fractions_rx(:)   ! Fractions on rof grid, cpl processes
!  type(mct_aVect) , pointer :: fractions_wx(:)   ! Fractions on wav grid, cpl processes
!
!  !--- domain equivalent 2d grid size ---
!  integer  :: atm_nx, atm_ny  ! nx, ny of 2d grid, if known
!  integer  :: lnd_nx, lnd_ny
!  integer  :: ice_nx, ice_ny
!  integer  :: ocn_nx, ocn_ny
!  integer  :: rof_nx, rof_ny
!  integer  :: glc_nx, glc_ny
!  integer  :: wav_nx, wav_ny
!
!  !----------------------------------------------------------------------------
!  ! Infodata: inter-model control flags, domain info
!  !----------------------------------------------------------------------------
!
!  type (seq_infodata_type), target :: infodata ! single instance for cpl and all comps
!
!  !----------------------------------------------------------------------------
!  ! time management
!  !----------------------------------------------------------------------------
!
!  type (seq_timemgr_type), SAVE :: seq_SyncClock ! array of all clocks & alarm
!  type (ESMF_Clock), target :: EClock_d      ! driver clock
!  type (ESMF_Clock), target :: EClock_a      ! atmosphere clock
!  type (ESMF_Clock), target :: EClock_l      ! land clock
!  type (ESMF_Clock), target :: EClock_o      ! ocean clock
!  type (ESMF_Clock), target :: EClock_i      ! ice clock
!  type (ESMF_Clock), target :: EClock_g      ! glc clock
!  type (ESMF_Clock), target :: EClock_r      ! rof clock
!  type (ESMF_Clock), target :: EClock_w      ! wav clock
!  type (ESMF_Clock), target :: EClock_e      ! esp clock
!
!  logical  :: restart_alarm          ! restart alarm
!  logical  :: history_alarm          ! history alarm
!  logical  :: histavg_alarm          ! history alarm
!  logical  :: stop_alarm             ! stop alarm
!  logical  :: atmrun_alarm           ! atm run alarm
!  logical  :: lndrun_alarm           ! lnd run alarm
!  logical  :: icerun_alarm           ! ice run alarm
!  logical  :: ocnrun_alarm           ! ocn run alarm
!  logical  :: ocnnext_alarm          ! ocn run alarm on next timestep
!  logical  :: glcrun_alarm           ! glc run alarm
!  logical  :: glcrun_avg_alarm       ! glc run averaging alarm
!  logical  :: rofrun_alarm           ! rof run alarm
!  logical  :: wavrun_alarm           ! wav run alarm
!  logical  :: esprun_alarm           ! esp run alarm
!  logical  :: tprof_alarm            ! timing profile alarm
!  logical  :: barrier_alarm          ! barrier alarm
!  logical  :: t1hr_alarm             ! alarm every hour
!  logical  :: t2hr_alarm             ! alarm every two hours
!  logical  :: t3hr_alarm             ! alarm every three hours
!  logical  :: t6hr_alarm             ! alarm every six hours
!  logical  :: t12hr_alarm            ! alarm every twelve hours
!  logical  :: t24hr_alarm            ! alarm every twentyfour hours
!  logical  :: t1yr_alarm             ! alarm every year, at start of year
!  logical  :: pause_alarm            ! pause alarm
!  logical  :: write_hist_alarm       ! alarm to write a history file under multiple conditions
!  integer  :: drv_index              ! seq_timemgr index for driver
!
!  real(r8) :: days_per_year = 365.0  ! days per year
!
!  integer  :: dtime                  ! dt of one coupling interval
!  integer  :: ncpl                   ! number of coupling intervals per day
!  integer  :: ymd                    ! Current date (YYYYMMDD)
!  integer  :: year                   ! Current date (YYYY)
!  integer  :: month                  ! Current date (MM)
!  integer  :: day                    ! Current date (DD)
!  integer  :: tod                    ! Current time of day (seconds)
!  integer  :: ymdtmp                 ! temporary date (YYYYMMDD)
!  integer  :: todtmp                 ! temporary time of day (seconds)
!  character(CL) :: orb_mode          ! orbital mode
!  character(CS) :: tfreeze_option    ! Freezing point calculation
!  integer  :: orb_iyear              ! orbital year
!  integer  :: orb_iyear_align        ! associated with model year
!  integer  :: orb_cyear              ! orbital year for current orbital computation
!  integer  :: orb_nyear              ! orbital year associated with currrent model year
!  real(r8) :: orb_eccen              ! orbital eccentricity
!  real(r8) :: orb_obliq              ! obliquity in degrees
!  real(r8) :: orb_mvelp              ! moving vernal equinox long
!  real(r8) :: orb_obliqr             ! Earths obliquity in rad
!  real(r8) :: orb_lambm0             ! Mean long of perihelion at vernal equinox (radians)
!  real(r8) :: orb_mvelpp             ! moving vernal equinox long
!  real(r8) :: wall_time_limit        ! wall time limit in hours
!  real(r8) :: wall_time              ! current wall time used
!  character(CS) :: force_stop_at     ! force stop at next (month, day, etc)
!  logical  :: force_stop             ! force the model to stop
!  integer  :: force_stop_ymd         ! force stop ymd
!  integer  :: force_stop_tod         ! force stop tod
!
!  !--- for documenting speed of the model ---
!  character(8) :: dstr              ! date string
!  character(10) :: tstr              ! time string
!  integer       :: begStep, endStep  ! Begining and ending step number
!  character(CL) :: calendar          ! calendar name
!  real(r8)      :: simDays           ! Number of simulated days
!  real(r8)      :: SYPD              ! Simulated years per day
!  real(r8)      :: Time_begin        ! Start time
!  real(r8)      :: Time_end          ! Ending time
!  real(r8)      :: Time_bstep        ! Start time
!  real(r8)      :: Time_estep        ! Ending time
!  real(r8)      :: time_brun         ! Start time
!  real(r8)      :: time_erun         ! Ending time
!  real(r8)      :: cktime            ! delta time
!  real(r8)      :: cktime_acc(10)    ! cktime accumulator array 1 = all, 2 = atm, etc
!  integer       :: cktime_cnt(10)    ! cktime counter array
!  real(r8)      :: max_cplstep_time
!  character(CL) :: timing_file       ! Local path to tprof filename
!  character(CL) :: timing_dir        ! timing directory
!  character(CL) :: tchkpt_dir        ! timing checkpoint directory
!
!  !----------------------------------------------------------------------------
!  ! control flags
!  !----------------------------------------------------------------------------
!
!  logical  :: atm_present            ! .true.  => atm is present
!  logical  :: lnd_present            ! .true.  => land is present
!  logical  :: ice_present            ! .true.  => ice is present
!  logical  :: ocn_present            ! .true.  => ocn is present
!  logical  :: glc_present            ! .true.  => glc is present
!  logical  :: glclnd_present         ! .true.  => glc is computing land coupling
!  logical  :: glcocn_present         ! .true.  => glc is computing ocean runoff
!  logical  :: glcice_present         ! .true.  => glc is computing icebergs
!  logical  :: rofice_present         ! .true.  => rof is computing icebergs
!  logical  :: rof_present            ! .true.  => rof is present
!  logical  :: flood_present          ! .true.  => rof is computing flood
!  logical  :: wav_present            ! .true.  => wav is present
!  logical  :: esp_present            ! .true.  => esp is present
!
!  logical  :: atm_prognostic         ! .true.  => atm comp expects input
!  logical  :: lnd_prognostic         ! .true.  => lnd comp expects input
!  logical  :: ice_prognostic         ! .true.  => ice comp expects input
!  logical  :: iceberg_prognostic     ! .true.  => ice comp can handle iceberg input
!  logical  :: ocn_prognostic         ! .true.  => ocn comp expects input
!  logical  :: ocnrof_prognostic      ! .true.  => ocn comp expects runoff input
!  logical  :: glc_prognostic         ! .true.  => glc comp expects input
!  logical  :: rof_prognostic         ! .true.  => rof comp expects input
!  logical  :: wav_prognostic         ! .true.  => wav comp expects input
!  logical  :: esp_prognostic         ! .true.  => esp comp expects input
!
!  logical  :: atm_c2_lnd             ! .true.  => atm to lnd coupling on
!  logical  :: atm_c2_ocn             ! .true.  => atm to ocn coupling on
!  logical  :: atm_c2_ice             ! .true.  => atm to ice coupling on
!  logical  :: atm_c2_wav             ! .true.  => atm to wav coupling on
!  logical  :: lnd_c2_atm             ! .true.  => lnd to atm coupling on
!  logical  :: lnd_c2_rof             ! .true.  => lnd to rof coupling on
!  logical  :: lnd_c2_glc             ! .true.  => lnd to glc coupling on
!  logical  :: ocn_c2_atm             ! .true.  => ocn to atm coupling on
!  logical  :: ocn_c2_ice             ! .true.  => ocn to ice coupling on
!  logical  :: ocn_c2_wav             ! .true.  => ocn to wav coupling on
!  logical  :: ice_c2_atm             ! .true.  => ice to atm coupling on
!  logical  :: ice_c2_ocn             ! .true.  => ice to ocn coupling on
!  logical  :: ice_c2_wav             ! .true.  => ice to wav coupling on
!  logical  :: rof_c2_lnd             ! .true.  => rof to lnd coupling on
!  logical  :: rof_c2_ocn             ! .true.  => rof to ocn coupling on
!  logical  :: rof_c2_ice             ! .true.  => rof to ice coupling on
!  logical  :: glc_c2_lnd             ! .true.  => glc to lnd coupling on
!  logical  :: glc_c2_ocn             ! .true.  => glc to ocn coupling on
!  logical  :: glc_c2_ice             ! .true.  => glc to ice coupling on
!  logical  :: wav_c2_ocn             ! .true.  => wav to ocn coupling on
!
!  logical  :: dead_comps             ! .true.  => dead components
!  logical  :: esmf_map_flag          ! .true.  => use esmf for mapping
!
!  logical  :: areafact_samegrid      ! areafact samegrid flag
!  logical  :: single_column          ! scm mode logical
!  real(r8) :: scmlon                 ! single column lon
!  real(r8) :: scmlat                 ! single column lat
!  logical  :: aqua_planet            ! aqua planet mode
!  real(r8) :: nextsw_cday            ! radiation control
!  logical  :: atm_aero               ! atm provides aerosol data
!
!  character(CL) :: cpl_seq_option    ! coupler sequencing option
!  logical  :: skip_ocean_run         ! skip the ocean model first pass
!  logical  :: cpl2ocn_first          ! use to call initial cpl2ocn timer
!  logical  :: run_barriers           ! barrier the component run calls
!
!  character(CS) :: aoflux_grid       ! grid for a/o flux calc: atm xor ocn
!  character(CS) :: vect_map          ! vector mapping type
!
!  character(CL) :: atm_gnam          ! atm grid
!  character(CL) :: lnd_gnam          ! lnd grid
!  character(CL) :: ocn_gnam          ! ocn grid
!  character(CL) :: ice_gnam          ! ice grid
!  character(CL) :: rof_gnam          ! rof grid
!  character(CL) :: glc_gnam          ! glc grid
!  character(CL) :: wav_gnam          ! wav grid
!
!  logical  :: samegrid_ao            ! samegrid atm and ocean
!  logical  :: samegrid_al            ! samegrid atm and land
!  logical  :: samegrid_lr            ! samegrid land and rof
!  logical  :: samegrid_oi            ! samegrid ocean and ice
!  logical  :: samegrid_ro            ! samegrid runoff and ocean
!  logical  :: samegrid_aw            ! samegrid atm and wave
!  logical  :: samegrid_ow            ! samegrid ocean and wave
!  logical  :: samegrid_lg            ! samegrid glc and land
!  logical  :: samegrid_og            ! samegrid glc and ocean
!  logical  :: samegrid_ig            ! samegrid glc and ice
!  logical  :: samegrid_alo           ! samegrid atm, lnd, ocean
!
!  logical       :: read_restart      ! local read restart flag
!  character(CL) :: rest_file         ! restart file path + filename
!
!  logical  :: shr_map_dopole         ! logical for dopole in shr_map_mod
!  logical  :: domain_check           ! .true.  => check consistency of domains
!  logical  :: reprosum_use_ddpdd     ! setup reprosum, use ddpdd
!  real(r8) :: reprosum_diffmax       ! setup reprosum, set rel_diff_max
!  logical  :: reprosum_recompute     ! setup reprosum, recompute if tolerance exceeded
!
!  logical  :: output_perf = .false.  ! require timing data output for this pe
!  logical  :: in_first_day = .true.  ! currently simulating first day
!
!  !--- history & budgets ---
!  logical :: do_budgets              ! heat/water budgets on
!  logical :: do_histinit             ! initial hist file
!  logical :: do_histavg              ! histavg on or off
!  logical :: do_hist_r2x             ! create aux files: r2x
!  logical :: do_hist_l2x             ! create aux files: l2x
!  logical :: do_hist_a2x24hr         ! create aux files: a2x
!  logical :: do_hist_l2x1yrg         ! create aux files: l2x 1yr glc forcings
!  logical :: do_hist_a2x             ! create aux files: a2x
!  logical :: do_hist_a2x3hrp         ! create aux files: a2x 3hr precip
!  logical :: do_hist_a2x3hr          ! create aux files: a2x 3hr states
!  logical :: do_hist_a2x1hri         ! create aux files: a2x 1hr instantaneous
!  logical :: do_hist_a2x1hr          ! create aux files: a2x 1hr
!  integer :: budget_inst             ! instantaneous budget flag
!  integer :: budget_daily            ! daily budget flag
!  integer :: budget_month            ! monthly budget flag
!  integer :: budget_ann              ! annual budget flag
!  integer :: budget_ltann            ! long term budget flag for end of year writing
!  integer :: budget_ltend            ! long term budget flag for end of run writing
!
!  character(CL) :: hist_a2x_flds     = &
!       'Faxa_swndr:Faxa_swvdr:Faxa_swndf:Faxa_swvdf'
!
!  character(CL) :: hist_a2x3hrp_flds = &
!       'Faxa_rainc:Faxa_rainl:Faxa_snowc:Faxa_snowl'
!
!  character(CL) :: hist_a2x24hr_flds = &
!       'Faxa_bcphiwet:Faxa_bcphodry:Faxa_bcphidry:Faxa_ocphiwet:Faxa_ocphidry:&
!       &Faxa_ocphodry:Faxa_dstwet1:Faxa_dstdry1:Faxa_dstwet2:Faxa_dstdry2:Faxa_dstwet3:&
!       &Faxa_dstdry3:Faxa_dstwet4:Faxa_dstdry4:Sa_co2prog:Sa_co2diag'
!
!  character(CL) :: hist_a2x1hri_flds = &
!       'Faxa_swndr:Faxa_swvdr:Faxa_swndf:Faxa_swvdf'
!
!  character(CL) :: hist_a2x1hr_flds  = &
!       'Sa_u:Sa_v'
!
!  character(CL) :: hist_a2x3hr_flds  = &
!       'Sa_z:Sa_topo:Sa_u:Sa_v:Sa_tbot:Sa_ptem:Sa_shum:Sa_dens:Sa_pbot:Sa_pslv:Faxa_lwdn:&
!       &Faxa_rainc:Faxa_rainl:Faxa_snowc:Faxa_snowl:&
!       &Faxa_swndr:Faxa_swvdr:Faxa_swndf:Faxa_swvdf:&
!       &Sa_co2diag:Sa_co2prog'
!
!  ! --- other ---
!
!  integer  :: ocnrun_count           ! number of times ocn run alarm went on
!  logical  :: exists                 ! true if file exists
!  integer  :: ierr                   ! MPI error return
!
!  character(*), parameter :: NLFileName = "drv_in"  ! input namelist filename
!
!  integer  :: info_debug = 0         ! local info_debug level
!
!  !----------------------------------------------------------------------------
!  ! memory monitoring
!  !----------------------------------------------------------------------------
!  real(r8) :: msize,msize0,msize1     ! memory size (high water)
!  real(r8) :: mrss ,mrss0 ,mrss1      ! resident size (current memory use)
!
!  !----------------------------------------------------------------------------
!  ! threading control
!  !----------------------------------------------------------------------------
!  integer  :: nthreads_GLOID         ! OMP global number of threads
!  integer  :: nthreads_CPLID         ! OMP cpl number of threads
!  integer  :: nthreads_ATMID         ! OMP atm number of threads
!  integer  :: nthreads_LNDID         ! OMP lnd number of threads
!  integer  :: nthreads_ICEID         ! OMP ice number of threads
!  integer  :: nthreads_OCNID         ! OMP ocn number of threads
!  integer  :: nthreads_GLCID         ! OMP glc number of threads
!  integer  :: nthreads_ROFID         ! OMP glc number of threads
!  integer  :: nthreads_WAVID         ! OMP wav number of threads
!  integer  :: nthreads_ESPID         ! OMP esp number of threads
!
!  integer  :: pethreads_GLOID        ! OMP number of threads per task
!
!  logical  :: drv_threading          ! driver threading control
!
!  !----------------------------------------------------------------------------
!  ! communicator groups and related
!  !----------------------------------------------------------------------------
!  integer  :: global_comm
!  integer  :: mpicom_GLOID          ! MPI global communicator
!  integer  :: mpicom_CPLID          ! MPI cpl communicator
!  integer  :: mpicom_OCNID          ! MPI ocn communicator for ensemble member 1
!
!  integer  :: mpicom_CPLALLATMID    ! MPI comm for CPLALLATMID
!  integer  :: mpicom_CPLALLLNDID    ! MPI comm for CPLALLLNDID
!  integer  :: mpicom_CPLALLICEID    ! MPI comm for CPLALLICEID
!  integer  :: mpicom_CPLALLOCNID    ! MPI comm for CPLALLOCNID
!  integer  :: mpicom_CPLALLGLCID    ! MPI comm for CPLALLGLCID
!  integer  :: mpicom_CPLALLROFID    ! MPI comm for CPLALLROFID
!  integer  :: mpicom_CPLALLWAVID    ! MPI comm for CPLALLWAVID
!
!  integer  :: iam_GLOID             ! pe number in global id
!  logical  :: iamin_CPLID           ! pe associated with CPLID
!  logical  :: iamroot_GLOID         ! GLOID masterproc
!  logical  :: iamroot_CPLID         ! CPLID masterproc
!
!  logical  :: iamin_CPLALLATMID     ! pe associated with CPLALLATMID
!  logical  :: iamin_CPLALLLNDID     ! pe associated with CPLALLLNDID
!  logical  :: iamin_CPLALLICEID     ! pe associated with CPLALLICEID
!  logical  :: iamin_CPLALLOCNID     ! pe associated with CPLALLOCNID
!  logical  :: iamin_CPLALLGLCID     ! pe associated with CPLALLGLCID
!  logical  :: iamin_CPLALLROFID     ! pe associated with CPLALLROFID
!  logical  :: iamin_CPLALLWAVID     ! pe associated with CPLALLWAVID
!
!
!  !----------------------------------------------------------------------------
!  ! complist: list of comps on this pe
!  !----------------------------------------------------------------------------
!
!  ! allow enough room for names of all physical components + coupler,
!  ! where each string can be up to (max_inst_name_len+1) characters
!  ! long (+1 allows for a space before each name)
!  character(len=(seq_comm_namelen+1)*(num_inst_phys+1)) :: complist
!
!  !----------------------------------------------------------------------------
!  ! comp_num_<comp>: unique component number for each component type
!  !----------------------------------------------------------------------------
!  integer, parameter :: comp_num_atm = 1
!  integer, parameter :: comp_num_lnd = 2
!  integer, parameter :: comp_num_ice = 3
!  integer, parameter :: comp_num_ocn = 4
!  integer, parameter :: comp_num_glc = 5
!  integer, parameter :: comp_num_rof = 6
!  integer, parameter :: comp_num_wav = 7
!  integer, parameter :: comp_num_esp = 8
!
!  !----------------------------------------------------------------------------
!  ! misc
!  !----------------------------------------------------------------------------
!
!  integer, parameter :: ens1=1         ! use first instance of ensemble only
!  integer, parameter :: fix1=1         ! temporary hard-coding to first ensemble, needs to be fixed
!  integer :: eai, eli, eoi, eii, egi, eri, ewi, eei, exi, efi  ! component instance counters
!
!  !----------------------------------------------------------------------------
!  ! formats
!  !----------------------------------------------------------------------------
!  character(*), parameter :: subname = '(seq_mct_drv)'
!  character(*), parameter :: F00 = "('"//subname//" : ', 4A )"
!  character(*), parameter :: F0L = "('"//subname//" : ', A, L6 )"
!  character(*), parameter :: F01 = "('"//subname//" : ', A, 2i8, 3x, A )"
!  character(*), parameter :: F0R = "('"//subname//" : ', A, 2g23.15 )"
!  character(*), parameter :: FormatA = '(A,": =============== ", A44,          " ===============")'
!  character(*), parameter :: FormatD = '(A,": =============== ", A20,I10.8,I8,8x,   " ===============")'
!  character(*), parameter :: FormatR = '(A,": =============== ", A31,F12.3,1x, " ===============")'
!  character(*), parameter :: FormatQ = '(A,": =============== ", A20,2F10.2,4x," ===============")'

  contains

subroutine clm_init(finname, pdaf_id, pdaf_max, mype) bind(C,name="clm_init")
  use cime_comp_mod, only : cime_pre_init1
  use cime_comp_mod, only : cime_pre_init2
  use cime_comp_mod, only : cime_init
  use spmdMod, only: mpicom, comp_id, spmd_init
  use iso_C_binding
  use enkf_clm_mod

  implicit none

  character(kind=c_char,len=1),dimension(100),intent(in) :: finname 
  integer(c_int), intent(in) :: pdaf_id
  integer(c_int), intent(in) :: pdaf_max
  integer(c_int), intent(in) :: mype
  integer(c_int) :: counter

  !--------------------------------------------------------------------------
  ! timing variables
  !--------------------------------------------------------------------------
  integer(i8) :: beg_count, end_count, irtc_rate
  real(r8)    :: cime_pre_init1_time, ESMF_Initialize_time, &
       cime_pre_init2_time, cime_init_time_adjustment

  !--------------------------------------------------------------------------
  ! For ESMF logging
  !--------------------------------------------------------------------------
  character(len=CS)       :: esmf_logfile_option
  type(ESMF_LogKind_Flag) :: esmf_logfile_kind

  ! finname input is already input file name with ensemble id attached
  ! this is done in wrapper_tsmp.c this loop just converts c into f string

  loop_string: do counter=1,clmprefixlen
        nlfilename(counter:counter) = finname(counter)
  end do loop_string

! spmd_init and mct_world init are called in cime_pre_init1 now.

!#if (defined COUP_OAS_COS || defined COUP_OAS_PFL)
!  write(6,*) 'Not yet supported model combination.'
!#elif (defined CLMSA)
!  !call spmd_init(COMM_model_oas)
!  !call mct_world_init(1,COMM_model_oas,mpicom,comp_id)
!#else
!  write(6,*) 'Not yet supported model combination.'
!#endif
  !write(*,*) "before cime_pre_init1", pdaf_id, pdaf_max
  !--------------------------------------------------------------------------
  ! Setup and initialize the communications and logging.
  !--------------------------------------------------------------------------
  beg_count = shr_sys_irtc(irtc_rate)
  call cime_pre_init1(esmf_logfile_option, COMM_model_clm, &
                      pdaf_id=pdaf_id, pdaf_max=pdaf_max)

  end_count = shr_sys_irtc(irtc_rate)

  cime_pre_init1_time = real( (end_count-beg_count), r8)/real(irtc_rate, r8)

  !--------------------------------------------------------------------------
  ! Initialize ESMF.  This is done outside of the ESMF_INTERFACE ifdef
  ! because it is needed for the time manager, even if the ESMF_INTERFACE
  ! is not used.
  !--------------------------------------------------------------------------
  beg_count = shr_sys_irtc(irtc_rate)


  select case(esmf_logfile_option)
  case('ESMF_LOGKIND_SINGLE')
     esmf_logfile_kind = ESMF_LOGKIND_SINGLE
  case('ESMF_LOGKIND_MULTI')
     esmf_logfile_kind = ESMF_LOGKIND_MULTI
  case('ESMF_LOGKIND_MULTI_ON_ERROR')
     esmf_logfile_kind = ESMF_LOGKIND_MULTI_ON_ERROR
  case('ESMF_LOGKIND_NONE')
     esmf_logfile_kind = ESMF_LOGKIND_NONE
  case default
     call shr_sys_abort('CIME ERROR: invalid ESMF logfile kind'//trim(esmf_logfile_option))
  end select
  write(6,*) "esmf_initialize"
  call ESMF_Initialize(logkindflag=esmf_logfile_kind)

  end_count = shr_sys_irtc(irtc_rate)
  ESMF_Initialize_time = real( (end_count-beg_count), r8)/real(irtc_rate, r8)

  !--------------------------------------------------------------------------
  ! Read in the configuration information and initialize the time manager.
  !--------------------------------------------------------------------------
  ! Timer initialization has to be after determination of the maximum number
  ! of threads used across all components, so called inside of
  ! cime_pre_init2, as are t_startf and t_stopf for CPL:INIT and
  ! cime_pre_init2.
  !--------------------------------------------------------------------------
  beg_count = shr_sys_irtc(irtc_rate)
  write(6,*) "cime-pre-init2"
  call cime_pre_init2()

  end_count = shr_sys_irtc(irtc_rate)
  cime_pre_init2_time = real( (end_count-beg_count), r8)/real(irtc_rate, r8)

  !--------------------------------------------------------------------------
  ! Call the initialize, run and finalize routines.
  !--------------------------------------------------------------------------

  call t_startf('CPL:INIT')
  call t_adj_detailf(+1)

  call t_startstop_valsf('CPL:cime_pre_init1',  walltime=cime_pre_init1_time)
  call t_startstop_valsf('CPL:ESMF_Initialize', walltime=ESMF_Initialize_time)
  call t_startstop_valsf('CPL:cime_pre_init2',  walltime=cime_pre_init2_time)

  call cime_init()

  call t_adj_detailf(-1)
  call t_stopf('CPL:INIT')

  cime_init_time_adjustment = cime_pre_init1_time  &
       + ESMF_Initialize_time &
       + cime_pre_init2_time
  call t_startstop_valsf('CPL:INIT',  walltime=cime_init_time_adjustment, &
       callcount=0)

#if defined CLMSA
  call define_clm_statevec(mype)
#endif 


end subroutine clm_init

!===============================================================================

!subroutine cime_comp_barriers(mpicom, timer)
!  integer         , intent(in) :: mpicom
!  character(len=*), intent(in) :: timer
!  integer :: ierr
!
!  if (run_barriers) then
!     call t_drvstartf (trim(timer))
!     call mpi_barrier(mpicom,ierr)
!     call t_drvstopf (trim(timer))
!  endif
!end subroutine cime_comp_barriers


subroutine clm_advance(ntstep, tstartcycle, mype) bind(C,name="clm_advance")
  use seq_comm_mct,   only: atm_layout, lnd_layout, ice_layout, glc_layout,  &
      rof_layout, ocn_layout, wav_layout, esp_layout
  use shr_string_mod, only: shr_string_listGetIndexF
  use seq_comm_mct, only: num_inst_driver
  use ESMF_TimeMod            , only : ESMF_Time
  use clm_time_manager        , only : get_nstep, is_last_step, advance_timestep

  use enkf_clm_mod, only : set_clm_statevec 

  use iso_C_binding

  implicit none

  integer(c_int),intent(in) :: ntstep
  integer(c_int),intent(in) :: tstartcycle
  integer(c_int),intent(in) :: mype
  integer :: counter=0
  integer :: nstep

  ! gptl timer lookup variables
  integer, parameter :: hashcnt=7
  integer      :: hashint(hashcnt)
  ! Driver pause/resume
  logical      :: drv_pause  ! Driver writes pause restart file
  character(len=CL)  :: drv_resume ! Driver resets state from restart file

  type(ESMF_Time)  :: etime_curr      ! Current model time
  real(r8)       :: tbnds1_offset     ! Time offset for call to seq_hist_writeaux
  logical      :: lnd2glc_averaged_now  ! Whether lnd2glc averages were taken this timestep

101 format( A, i10.8, i8, 12A, A, F8.2, A, F8.2 )
102 format( A, i10.8, i8, A, 8L3 )
103 format( 5A )
104 format( A, i10.8, i8)
105 format( A, i10.8, i8, A, f10.2, A, f10.2, A, A, i5, A, A)
108 format( A, f10.2, A, i8.8)
109 format( A, 2f10.3)


  hashint = 0

  call seq_infodata_putData(infodata,atm_phase=1,lnd_phase=1,ocn_phase=1,ice_phase=1)
  call seq_timemgr_EClockGetData( EClock_d, stepno=begstep)
  call seq_timemgr_EClockGetData( EClock_d, dtime=dtime)
  call seq_timemgr_EClockGetData( EClock_d, calendar=calendar)
  ncpl = 86400/dtime
  cktime_acc = 0._r8
  cktime_cnt = 0
  stop_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_stop)
  if (seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_datestop)) then
     if (iamroot_CPLID) then
      write(logunit,*) ' '
      write(logunit,103) subname,' NOTE: Stopping from alarm STOP DATE'
      write(logunit,*) ' '
     endif
     stop_alarm = .true.
  endif
  force_stop = .false.
  force_stop_ymd = -1
  force_stop_tod = -1

  !|----------------------------------------------------------
  !| Beginning of driver time step loop
  !|----------------------------------------------------------

  call t_startf ('CPL:RUN_LOOP_BSTART')
  call mpi_barrier(mpicom_GLOID,ierr)
  call t_stopf ('CPL:RUN_LOOP_BSTART')
  Time_begin = mpi_wtime()
  Time_bstep = mpi_wtime()

  do counter=1,ntstep
    nstep = get_nstep()

    call t_startf('CPL:RUN_LOOP', hashint(1))
    call t_startf('CPL:CLOCK_ADVANCE')

    !----------------------------------------------------------
    !| Advance Clock
    !  (this is time that models should have before they return
    !  to the driver).  Write timestamp and run alarm status
    !----------------------------------------------------------

    call seq_timemgr_clockAdvance(seq_SyncClock, force_stop, &
                                  force_stop_ymd, force_stop_tod)
    call seq_timemgr_EClockGetData( EClock_d, curr_ymd=ymd, curr_tod=tod)
    call shr_cal_date2ymd(ymd,year,month,day)
    stop_alarm    = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_stop)
    atmrun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_atmrun)
    lndrun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_lndrun)
    rofrun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_rofrun)
    icerun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_icerun)
    glcrun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_glcrun)
    wavrun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_wavrun)
    esprun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_esprun)
    ocnrun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_ocnrun)
    ocnnext_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_ocnnext)
    restart_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_restart)
    history_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_history)
    histavg_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_histavg)
    tprof_alarm   = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_tprof)
    barrier_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_barrier)
    pause_alarm   = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_pause)

    ! Does the driver need to pause?
    drv_pause = pause_alarm .and. seq_timemgr_pause_component_active(drv_index)

    if (glc_prognostic) then
    ! Is it time to average fields to pass to glc?
    !
    ! Note that the glcrun_avg_alarm just controls what is passed to glc
    ! in terms
    ! of averaged fields - it does NOT control when glc is called
    ! currently -
    ! glc will be called on the glcrun_alarm setting - but it might not be
    ! passed relevant
    ! info if the time averaging period to accumulate information passed
    ! to glc is greater
    ! than the glcrun interval
      glcrun_avg_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_glcrun_avg)
      if (glcrun_avg_alarm .and. .not. glcrun_alarm) then
        write(logunit,*) 'ERROR: glcrun_avg_alarm is true, but glcrun_alarm is false'
        write(logunit,*) 'Make sure that NCPL_BASE_PERIOD, GLC_NCPL and GLC_AVG_PERIOD'
        write(logunit,*) 'are set so that glc averaging only happens at glc coupling times'
        write(logunit,*) '(It is allowable for glc coupling to be more frequent '
        write(logunit,*) 'than glc averaging,'
        write(logunit,*) 'but not for glc averaging to be more frequent '
        write(logunit,*) 'than glc coupling.)'
        call shr_sys_abort(subname//' glcrun_avg_alarm is true, but glcrun_alarm is false')
      end if
    else
      ! glcrun_avg_alarm shouldn't matter in this case
      glcrun_avg_alarm = .false.
    end if

    ! this probably belongs in seq_timemgr somewhere using proper clocks
    t1hr_alarm = .false.
    t2hr_alarm = .false.
    t3hr_alarm = .false.
    t6hr_alarm = .false.
    t12hr_alarm = .false.
    t24hr_alarm = .false.
    t1yr_alarm = .false.
    if (mod(tod, 3600) == 0) t1hr_alarm = .true.
    if (mod(tod, 7200) == 0) t2hr_alarm = .true.
    if (mod(tod,10800) == 0) t3hr_alarm = .true.
    if (mod(tod,21600) == 0) t6hr_alarm = .true.
    if (mod(tod,43200) == 0) t12hr_alarm = .true.
    if (tod            == 0) t24hr_alarm = .true.
    if (month==1 .and. day==1 .and. tod==0) t1yr_alarm = .true.

    lnd2glc_averaged_now = .false.

    if (seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_datestop)) then
      if (iamroot_CPLID) then
        write(logunit,*) ' '
        write(logunit,103) subname,' NOTE: Stopping from alarm STOP DATE'
        write(logunit,*) ' '
      endif
      stop_alarm = .true.
    endif

    ! update the orbital data as needed
    if (trim(orb_mode) == trim(seq_infodata_orb_variable_year)) then
      orb_nyear =  orb_iyear + (year - orb_iyear_align)
      if (orb_nyear /= orb_cyear) then
        orb_cyear = orb_nyear
        call shr_orb_params(orb_cyear, orb_eccen, orb_obliq, orb_mvelp, &
                            orb_obliqr, orb_lambm0, orb_mvelpp, iamroot_CPLID)
        call seq_infodata_putData(infodata,orb_eccen=orb_eccen,orb_obliqr=orb_obliqr, &
                                  orb_lambm0=orb_lambm0,orb_mvelpp=orb_mvelpp)
      endif
    endif

    ! override ocnrun_alarm and ocnnext_alarm for first ocn run
    ! skip_ocean_run is initialized above to true if it's a startup
    ! if it's not a startup, ignore all of this
    ! stop the overide on the second ocnrun_alarm

    if (ocnrun_alarm) ocnrun_count = ocnrun_count + 1
    if (ocnrun_count > 1) skip_ocean_run = .false.
    if (skip_ocean_run) then
      ocnrun_alarm = .false.
      ocnnext_alarm = .false.
    endif

    if (iamroot_CPLID) then
      if (loglevel > 1) then
        write(logunit,102) ' Alarm_state: model date = ',ymd,tod, &
              ' aliogrw run alarms = ',  atmrun_alarm, lndrun_alarm, &
              icerun_alarm, ocnrun_alarm, glcrun_alarm, &
              rofrun_alarm, wavrun_alarm, esprun_alarm
        write(logunit,102) ' Alarm_state: model date = ',ymd,tod, &
              ' 1.2.3.6.12.24 run alarms = ',  t1hr_alarm, t2hr_alarm, &
              t3hr_alarm, t6hr_alarm, t12hr_alarm, t24hr_alarm
        call shr_sys_flush(logunit)
      endif
    endif

    call t_stopf ('CPL:CLOCK_ADVANCE')

    !----------------------------------------------------------
    !| MAP ATM to OCN
    !  Set a2x_ox as a module variable in prep_ocn_mod
    !  This will be used later in the ice prep and in the
    !  atm/ocn flux calculation
    !----------------------------------------------------------

    if (iamin_CPLID .and. (atm_c2_ocn .or. atm_c2_ice)) then
      call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:OCNPRE1_BARRIER')
      call t_drvstartf('CPL:OCNPRE1', cplrun=.true., &
                       barrier=mpicom_CPLID,hashint=hashint(3))
      if (drv_threading) then
        call seq_comm_setnthreads(nthreads_CPLID)
      end if

      call prep_ocn_calc_a2x_ox(timer='CPL:ocnpre1_atm2ocn')
      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
        call t_drvstopf  ('CPL:OCNPRE1',cplrun=.true.,hashint=hashint(3))
      endif

      !----------------------------------------------------------
      !| ATM/OCN SETUP (rasm_option1)
      !----------------------------------------------------------

      if ((trim(cpl_seq_option) == 'RASM_OPTION1') .and. &
           iamin_CPLID .and. ocn_present) then

        call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:ATMOCN1_BARRIER')
        call t_drvstartf('CPL:ATMOCN1',cplrun=.true., barrier=mpicom_CPLID, &
                         hashint=hashint(4))
        if (drv_threading) then 
          call seq_comm_setnthreads(nthreads_CPLID)
        end if

        if (ocn_prognostic) then
          ! Map ice to ocn
          if (ice_c2_ocn) then
            call prep_ocn_calc_i2x_ox(timer='CPL:atmocnp_ice2ocn')
          end if
          ! Map wav to ocn
          if (wav_c2_ocn) then 
            call prep_ocn_calc_w2x_ox(timer='CPL:atmocnp_wav2ocn')
          end if
        endif

        !----------------------------------------------------------
        !| atm/ocn flux on atm grid (rasm_option1 and aoflux='atm')
        !----------------------------------------------------------

        if (trim(aoflux_grid) == 'atm') then
          ! compute o2x_ax for flux_atmocn, will be updated before atm merge
          ! do not use fractions because fractions here are NOT consistent
          ! with fractions in atm_mrg
          if (ocn_c2_atm) then
            call prep_atm_calc_o2x_ax(timer='CPL:atmoca_ocn2atm')
          end if
          call t_drvstartf ('CPL:atmocna_fluxa',barrier=mpicom_CPLID)
          do exi = 1,num_inst_xao
            eai = mod((exi-1),num_inst_atm) + 1
            eoi = mod((exi-1),num_inst_ocn) + 1
            efi = mod((exi-1),num_inst_frc) + 1
            a2x_ax => component_get_c2x_cx(atm(eai))
            o2x_ax => prep_atm_get_o2x_ax()    ! array over all instances
            xao_ax => prep_aoflux_get_xao_ax() ! array over all instances
            call seq_flux_atmocn_mct(infodata, tod, dtime, a2x_ax, &
                                     o2x_ax(eoi), xao_ax(exi))
          end do
          call t_drvstopf  ('CPL:atmocna_fluxa')

          if (atm_c2_ocn) then 
            call prep_aoflux_calc_xao_ox(timer='CPL:atmocna_atm2ocn')
          end if
        endif  ! aoflux_grid

        !----------------------------------------------------------
        !| atm/ocn flux on ocn grid (rasm_option1 and aoflux='ocn')
        !----------------------------------------------------------

        if (trim(aoflux_grid) == 'ocn') then
          call t_drvstartf('CPL:atmocnp_fluxo',barrier=mpicom_CPLID,hashint=hashint(6))
          do exi = 1,num_inst_xao
            eai = mod((exi-1),num_inst_atm) + 1
            eoi = mod((exi-1),num_inst_ocn) + 1
            efi = mod((exi-1),num_inst_frc) + 1
            a2x_ox => prep_ocn_get_a2x_ox()
            o2x_ox => component_get_c2x_cx(ocn(eoi))
            xao_ox => prep_aoflux_get_xao_ox()
            call seq_flux_atmocn_mct(infodata, tod, dtime,  &
                                     a2x_ox(eai),o2x_ox, xao_ox(exi))
          end do
          call t_drvstopf  ('CPL:atmocnp_fluxo',hashint=hashint(6))
        end if

        !----------------------------------------------------------
        !| ocn prep-merge (rasm_option1)
        !----------------------------------------------------------

        xao_ox => prep_aoflux_get_xao_ox()
        call prep_ocn_mrg(infodata, fractions_ox, &
                          xao_ox=xao_ox, timer_mrg='CPL:atmocnp_mrgx2o')

        ! Accumulate ocn inputs - form partial sum of tavg ocn inputs 
        ! (virtual "send" to ocn)
        call prep_ocn_accum(timer='CPL:atmocnp_accum')

        !----------------------------------------------------------
        !| ocn albedos (rasm_option1)
        ! (MUST BE AFTER prep_ocn_mrg for swnet to ocn to be computed properly)
        !----------------------------------------------------------

        call t_drvstartf ('CPL:atmocnp_ocnalb', barrier=mpicom_CPLID,hashint=hashint(5))
        do exi = 1,num_inst_xao
          efi = mod((exi-1),num_inst_frc) + 1
          eai = mod((exi-1),num_inst_atm) + 1
          xao_ox => prep_aoflux_get_xao_ox()  ! array over all instances
          a2x_ox => prep_ocn_get_a2x_ox()
          call seq_flux_ocnalb_mct(infodata, ocn(1), a2x_ox(eai), &
                                   fractions_ox(efi), xao_ox(exi))
        end do
        call t_drvstopf  ('CPL:atmocnp_ocnalb',hashint=hashint(5))

        !----------------------------------------------------------
        !| ocn budget (rasm_option1)
        !----------------------------------------------------------

        if (do_budgets) then
          call cime_comp_barriers(mpicom=mpicom_CPLID,timer='CPL:BUDGET0_BARRIER')
          call t_drvstartf ('CPL:BUDGET0',budget=.true.,barrier=mpicom_CPLID)
          xao_ox => prep_aoflux_get_xao_ox() ! array over all instances
          call seq_diag_ocn_mct(ocn(ens1), xao_ox(1), fractions_ox(ens1), infodata, &
               do_o2x=.true., do_x2o=.true., do_xao=.true.)
          call t_drvstopf ('CPL:BUDGET0',budget=.true.)
        end if

        if (drv_threading) then 
           call seq_comm_setnthreads(nthreads_GLOID)
        end if
        call t_drvstopf('CPL:ATMOCN1',cplrun=.true.,hashint=hashint(4))
      end if

      !----------------------------------------------------------
      !| ATM/OCN SETUP-SEND (cesm1_orig, cesm1_orig_tight, cesm1_mod,
      !cesm1_mod_tight, or rasm_option1)
      !----------------------------------------------------------

      if ((trim(cpl_seq_option) == 'CESM1_ORIG' .or. &
           trim(cpl_seq_option) == 'CESM1_ORIG_TIGHT' .or. &
           trim(cpl_seq_option) == 'CESM1_MOD'  .or. &
           trim(cpl_seq_option) == 'CESM1_MOD_TIGHT'  .or. &
           trim(cpl_seq_option) == 'RASM_OPTION1'  ) .and. &
           ocn_present .and. ocnrun_alarm) then

        !----------------------------------------------------
        ! "startup" wait (cesm1_orig, cesm1_mod, or rasm_option1)
        !----------------------------------------------------

        if (iamin_CPLALLOCNID) then
          ! want to know the time the ocean pes waited for the cpl pes
          ! at the first ocnrun_alarm, min ocean wait is wait time
          ! do not use t_barrierf here since it can be "off", use mpi_barrier
          do eoi = 1,num_inst_ocn
            if (ocn(eoi)%iamin_compid) then
              call t_drvstartf ('CPL:C2O_INITWAIT')
            end if
          end do
          call mpi_barrier(mpicom_CPLALLOCNID,ierr)
          do eoi = 1,num_inst_ocn
            if (ocn(eoi)%iamin_compid) then 
              call t_drvstopf  ('CPL:C2O_INITWAIT')
            end if
          end do
          cpl2ocn_first = .false.
        end if

        !----------------------------------------------------
        !| ocn average (cesm1_orig, cesm1_orig_tight, cesm1_mod,
        !cesm1_mod_tight, or rasm_option1)
        !----------------------------------------------------

        if (iamin_CPLID .and. ocn_prognostic) then
          call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:OCNPREP_BARRIER')
          call t_drvstartf ('CPL:OCNPREP',cplrun=.true.,barrier=mpicom_CPLID)
          if (drv_threading) then
            call seq_comm_setnthreads(nthreads_CPLID)
          end if 
          ! finish accumulating ocean inputs
          ! reset the value of x2o_ox with the value in x2oacc_ox
          ! (module variable in prep_ocn_mod)
          call prep_ocn_accum_avg(timer_accum='CPL:ocnprep_avg')
          call component_diag(infodata, ocn, flow='x2c', comment= 'send ocn', &
                              info_debug=info_debug, timer_diag='CPL:ocnprep_diagav')

          if (drv_threading) then 
            call seq_comm_setnthreads(nthreads_GLOID)
          end if 
          call t_drvstopf('CPL:OCNPREP',cplrun=.true.)
        end if

        !----------------------------------------------------
        !| cpl -> ocn (cesm1_orig, cesm1_orig_tight, cesm1_mod,
        !cesm1_mod_tight, or rasm_option1)
        !----------------------------------------------------

        if (iamin_CPLALLOCNID .and. ocn_prognostic) then
          call component_exch(ocn, flow='x2c', &
               infodata=infodata, infodata_string='cpl2ocn_run', &
               mpicom_barrier=mpicom_CPLALLOCNID, run_barriers=run_barriers, &
               timer_barrier='CPL:C2O_BARRIER', timer_comp_exch='CPL:C2O', &
               timer_map_exch='CPL:c2o_ocnx2ocno', &
               timer_infodata_exch='CPL:c2o_infoexch')
        end if
      end if ! end of OCN SETUP

      !----------------------------------------------------------
      !| LND SETUP-SEND
      !----------------------------------------------------------

      if (lnd_present .and. lndrun_alarm) then
        !----------------------------------------------------
        !| lnd prep-merge
        !----------------------------------------------------
        if (iamin_CPLID) then
          call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:LNDPREP_BARRIER')
          call t_drvstartf ('CPL:LNDPREP',cplrun=.true.,barrier=mpicom_CPLID)
          if (drv_threading) then
            call seq_comm_setnthreads(nthreads_CPLID)
          end if
          if (atm_c2_lnd) then
            call prep_lnd_calc_a2x_lx(timer='CPL:lndprep_atm2lnd')
          endif

          if (lnd_prognostic) then
            call prep_lnd_mrg(infodata, timer_mrg='CPL:lndprep_mrgx2l')
            call component_diag(infodata, lnd, flow='x2c', comment= 'send lnd', &
                                info_debug=info_debug, timer_diag='CPL:lndprep_diagav')
          endif

          if (drv_threading) then
            call seq_comm_setnthreads(nthreads_GLOID)
          end if
          call t_drvstopf  ('CPL:LNDPREP',cplrun=.true.)
        endif

        !----------------------------------------------------
        !| cpl -> lnd
        !----------------------------------------------------
        if (iamin_CPLALLLNDID) then
          call component_exch(lnd, flow='x2c', &
               infodata=infodata, infodata_string='cpl2lnd_run', &
               mpicom_barrier=mpicom_CPLALLLNDID, run_barriers=run_barriers, &
               timer_barrier='CPL:C2L_BARRIER', timer_comp_exch='CPL:C2L', &
               timer_map_exch='CPL:c2l_lndx2lndl', timer_infodata_exch='CPL:c2l_infoexch')
        endif
      endif

      !----------------------------------------------------------
      !| ICE SETUP-SEND
      !  Note that for atm->ice mapping below will leverage the assumption that
      !  the
      !  ice and ocn are on the same grid and that mapping of atm to ocean is
      !  done already for use by atmocn flux and ice model prep
      !----------------------------------------------------------
      if (ice_present .and. icerun_alarm) then
        !----------------------------------------------------
        !| ice prep-merge
        !----------------------------------------------------
        if (iamin_CPLID .and. ice_prognostic) then
          call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:ICEPREP_BARRIER')
          call t_drvstartf ('CPL:ICEPREP',cplrun=.true.,barrier=mpicom_CPLID)
          if (drv_threading) then 
            call seq_comm_setnthreads(nthreads_CPLID)
          end if

          if (ocn_c2_ice) then
            call prep_ice_calc_o2x_ix(timer='CPL:iceprep_ocn2ice')
          end if

          if (atm_c2_ice) then
            ! This is special to avoid remapping atm to ocn
            ! Note it is constrained that different prep modules cannot
            ! use or call each other
            a2x_ox => prep_ocn_get_a2x_ox() ! array
            call prep_ice_calc_a2x_ix(a2x_ox, timer='CPL:iceprep_atm2ice')
          end if

          call prep_ice_mrg(infodata, timer_mrg='CPL:iceprep_mrgx2i')
          call component_diag(infodata, ice, flow='x2c', comment= 'send ice', &
                              info_debug=info_debug, timer_diag='CPL:iceprep_diagav')

          if (drv_threading) then 
            call seq_comm_setnthreads(nthreads_GLOID)
          end if
          call t_drvstopf  ('CPL:ICEPREP',cplrun=.true.)
        end if

        !----------------------------------------------------
        !| cpl -> ice
        !----------------------------------------------------

        if (iamin_CPLALLICEID .and. ice_prognostic) then
          call component_exch(ice, flow='x2c', &
               infodata=infodata, infodata_string='cpl2ice_run', &
               mpicom_barrier=mpicom_CPLALLICEID, run_barriers=run_barriers, &
               timer_barrier='CPL:C2I_BARRIER', timer_comp_exch='CPL:C2I', &
               timer_map_exch='CPL:c2i_icex2icei', timer_infodata_exch='CPL:ice_infoexch')
        endif
      endif

      !----------------------------------------------------------
      !| WAV SETUP-SEND
      !----------------------------------------------------------
      if (wav_present .and. wavrun_alarm) then
        !----------------------------------------------------------
        !| wav prep-merge
        !----------------------------------------------------------
        if (iamin_CPLID .and. wav_prognostic) then
          call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:WAVPREP_BARRIER')
          call t_drvstartf ('CPL:WAVPREP',cplrun=.true.,barrier=mpicom_CPLID)
          if (drv_threading) then 
            call seq_comm_setnthreads(nthreads_CPLID)
          end if

          if (atm_c2_wav) then
            call prep_wav_calc_a2x_wx(timer='CPL:wavprep_atm2wav')
          endif

          if (ocn_c2_wav) then
            call prep_wav_calc_o2x_wx(timer='CPL:wavprep_ocn2wav')
          endif

          if (ice_c2_wav) then
            call prep_wav_calc_i2x_wx(timer='CPL:wavprep_ice2wav')
          endif

          call prep_wav_mrg(infodata, fractions_wx, timer_mrg='CPL:wavprep_mrgx2w')

          call component_diag(infodata, wav, flow='x2c', comment= 'send wav', &
                              info_debug=info_debug, timer_diag='CPL:wavprep_diagav')

          if (drv_threading) then
            call seq_comm_setnthreads(nthreads_GLOID)
          end if   
          call t_drvstopf  ('CPL:WAVPREP',cplrun=.true.)
        endif

        !----------------------------------------------------------
        !| cpl -> wav
        !----------------------------------------------------------

        if (iamin_CPLALLWAVID .and. wav_prognostic) then
          call component_exch(wav, flow='x2c', &
               infodata=infodata, infodata_string='cpl2wav_run', &
               mpicom_barrier=mpicom_CPLALLWAVID, run_barriers=run_barriers, &
               timer_barrier='CPL:C2W_BARRIER', timer_comp_exch='CPL:C2W', &
               timer_map_exch='CPL:c2w_wavx2wavw', &
               timer_infodata_exch='CPL:c2w_infoexch')
        end if
      end if

      !----------------------------------------------------------
      !| ROF SETUP-SEND
      !----------------------------------------------------------

      if (rof_present .and. rofrun_alarm) then
        !----------------------------------------------------
        !| rof prep-merge
        !----------------------------------------------------
        if (iamin_CPLID .and. rof_prognostic) then
          call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:ROFPREP_BARRIER')
          call t_drvstartf ('CPL:ROFPREP', cplrun=.true., barrier=mpicom_CPLID)
          if (drv_threading) then 
            call seq_comm_setnthreads(nthreads_CPLID)
          end if
          call prep_rof_accum_avg(timer='CPL:rofprep_l2xavg')

          if (lnd_c2_rof) then
            call prep_rof_calc_l2r_rx(fractions_lx, timer='CPL:rofprep_lnd2rof')
          end if

          call prep_rof_mrg(infodata, fractions_rx, timer_mrg='CPL:rofprep_mrgx2r')
          call component_diag(infodata, rof, flow='x2c', comment= 'send rof', &
                              info_debug=info_debug, timer_diag='CPL:rofprep_diagav')

          if (drv_threading) then 
            call seq_comm_setnthreads(nthreads_GLOID)
          end if   
          call t_drvstopf  ('CPL:ROFPREP',cplrun=.true.)
        end if

        !----------------------------------------------------
        !| cpl -> rof
        !----------------------------------------------------

        if (iamin_CPLALLROFID .and. rof_prognostic) then
          call component_exch(rof, flow='x2c', &
               infodata=infodata, infodata_string='cpl2rof_run', &
               mpicom_barrier=mpicom_CPLALLLNDID, run_barriers=run_barriers, &
               timer_barrier='CPL:C2R_BARRIER', timer_comp_exch='CPL:C2R', &
               timer_map_exch='CPL:c2r_rofx2rofr', timer_infodata_exch='CPL:c2r_infoexch')
        end if
      end if

      !----------------------------------------------------------
      !| RUN ICE MODEL
      !----------------------------------------------------------

      if (ice_present .and. icerun_alarm) then
        call component_run(Eclock_i, ice, ice_run, infodata, &
             seq_flds_x2c_fluxes=seq_flds_x2i_fluxes, &
             seq_flds_c2x_fluxes=seq_flds_i2x_fluxes, &
             comp_prognostic=ice_prognostic, comp_num=comp_num_ice, &
             timer_barrier= 'CPL:ICE_RUN_BARRIER', timer_comp_run='CPL:ICE_RUN', &
             run_barriers=run_barriers, ymd=ymd, tod=tod,comp_layout=ice_layout)
      endif

      !----------------------------------------------------------
      !| RUN LND MODEL
      !----------------------------------------------------------

      if (lnd_present .and. lndrun_alarm) then
        call component_run(Eclock_l, lnd, lnd_run, infodata, &
             seq_flds_x2c_fluxes=seq_flds_x2l_fluxes, &
             seq_flds_c2x_fluxes=seq_flds_l2x_fluxes, &
             comp_prognostic=lnd_prognostic, comp_num=comp_num_lnd, &
             timer_barrier= 'CPL:LND_RUN_BARRIER', timer_comp_run='CPL:LND_RUN', &
             run_barriers=run_barriers, ymd=ymd, tod=tod,comp_layout=lnd_layout)
      endif

      !----------------------------------------------------------
      !| RUN ROF MODEL
      !----------------------------------------------------------

      if (rof_present .and. rofrun_alarm) then
        call component_run(Eclock_r, rof, rof_run, infodata, &
             seq_flds_x2c_fluxes=seq_flds_x2r_fluxes, &
             seq_flds_c2x_fluxes=seq_flds_r2x_fluxes, &
             comp_prognostic=rof_prognostic, comp_num=comp_num_rof, &
             timer_barrier= 'CPL:ROF_RUN_BARRIER', timer_comp_run='CPL:ROF_RUN', &
             run_barriers=run_barriers, ymd=ymd, tod=tod,comp_layout=rof_layout)
      endif

      !----------------------------------------------------------
      !| RUN WAV MODEL
      !----------------------------------------------------------

      if (wav_present .and. wavrun_alarm) then
        call component_run(Eclock_w, wav, wav_run, infodata, &
             seq_flds_x2c_fluxes=seq_flds_x2w_fluxes, &
             seq_flds_c2x_fluxes=seq_flds_w2x_fluxes, &
             comp_prognostic=wav_prognostic, comp_num=comp_num_wav, &
             timer_barrier= 'CPL:WAV_RUN_BARRIER', timer_comp_run='CPL:WAV_RUN', &
             run_barriers=run_barriers, ymd=ymd, tod=tod,comp_layout=wav_layout)
      endif

      !----------------------------------------------------------
      !| RUN OCN MODEL (cesm1_orig_tight or cesm1_mod_tight)
      !----------------------------------------------------------

      if ((trim(cpl_seq_option) == 'CESM1_ORIG_TIGHT' .or. &
           trim(cpl_seq_option) == 'CESM1_MOD_TIGHT'   ) .and. &
           ocn_present .and. ocnrun_alarm) then
        call component_run(Eclock_o, ocn, ocn_run, infodata, &
             seq_flds_x2c_fluxes=seq_flds_x2o_fluxes, &
             seq_flds_c2x_fluxes=seq_flds_o2x_fluxes, &
             comp_prognostic=ocn_prognostic, comp_num=comp_num_ocn, &
             timer_barrier= 'CPL:OCNT_RUN_BARRIER', timer_comp_run='CPL:OCNT_RUN', &
               run_barriers=run_barriers, ymd=ymd, tod=tod,comp_layout=ocn_layout)
      endif

      !----------------------------------------------------------
      !| OCN RECV-POST (cesm1_orig_tight or cesm1_mod_tight)
      !----------------------------------------------------------

      if ((trim(cpl_seq_option) == 'CESM1_ORIG_TIGHT' .or. &
           trim(cpl_seq_option) == 'CESM1_MOD_TIGHT'   ) .and. &
           ocn_present .and. ocnnext_alarm) then

        !----------------------------------------------------------
        !| ocn -> cpl (cesm1_orig_tight or cesm1_mod_tight)
        !----------------------------------------------------------

        if (iamin_CPLALLOCNID) then
          call component_exch(ocn, flow='c2x', &
               infodata=infodata, infodata_string='ocn2cpl_run', &
               mpicom_barrier=mpicom_CPLALLOCNID, run_barriers=run_barriers, &
               timer_barrier='CPL:O2CT_BARRIER', timer_comp_exch='CPL:O2CT', &
               timer_map_exch='CPL:o2c_ocno2ocnx', &
               timer_infodata_exch='CPL:o2c_infoexch')
        endif

        !----------------------------------------------------------
        !| ocn post (cesm1_orig_tight or cesm1_mod_tight)
        !----------------------------------------------------------

        if (iamin_CPLID) then
          call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:OCNPOSTT_BARRIER')
          call t_drvstartf('CPL:OCNPOSTT',cplrun=.true.,barrier=mpicom_CPLID)
          if (drv_threading) then
            call seq_comm_setnthreads(nthreads_CPLID)
          end if
          call component_diag(infodata, ocn, flow='c2x', comment= 'recv ocn', &
                              info_debug=info_debug, timer_diag='CPL:ocnpost_diagav')

          if (drv_threading) then
            call seq_comm_setnthreads(nthreads_GLOID)
          end if
          call t_drvstopf('CPL:OCNPOSTT',cplrun=.true.)
        endif
      endif

      !----------------------------------------------------------
      !| ATM/OCN SETUP (cesm1_orig, cesm1_orig_tight, cesm1_mod or
      !cesm1_mod_tight)
      !----------------------------------------------------------
      if ((trim(cpl_seq_option) == 'CESM1_ORIG'       .or. &
           trim(cpl_seq_option) == 'CESM1_ORIG_TIGHT' .or. &
           trim(cpl_seq_option) == 'CESM1_MOD'        .or. &
           trim(cpl_seq_option) == 'CESM1_MOD_TIGHT' ) .and. &
           iamin_CPLID .and. ocn_present) then
        call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:ATMOCNP_BARRIER')
        call t_drvstartf('CPL:ATMOCNP', cplrun=.true., &
             barrier=mpicom_CPLID, hashint=hashint(7))
        if (drv_threading) then
          call seq_comm_setnthreads(nthreads_CPLID)
        end if
        !----------------------------------------------------------
        !| ocn prep-merge (cesm1_orig or cesm1_orig_tight)
        !----------------------------------------------------------

        if (ocn_prognostic) then
          ! Map ice to ocn
          if (ice_c2_ocn) then
            call prep_ocn_calc_i2x_ox(timer='CPL:atmocnp_ice2ocn')
          end if
          ! Map wav to ocn
          if (wav_c2_ocn) then
            call prep_ocn_calc_w2x_ox(timer='CPL:atmocnp_wav2ocn')
          end if
          
          if (cpl_seq_option == 'CESM1_ORIG' .or. &
              cpl_seq_option == 'CESM1_ORIG_TIGHT') then
              xao_ox => prep_aoflux_get_xao_ox()
            call prep_ocn_mrg(infodata, fractions_ox, xao_ox=xao_ox, &
                              timer_mrg='CPL:atmocnp_mrgx2o')

            ! Accumulate ocn inputs - form partial sum of tavg ocn inputs
            ! (virtual "send" to ocn)
            call prep_ocn_accum(timer='CPL:atmocnp_accum')
          endif
        endif

        !----------------------------------------------------------
        !| atm/ocn flux on atm grid ((cesm1_orig, cesm1_orig_tight, cesm1_mod
        !or cesm1_mod_tight) and aoflux='atm')
        !----------------------------------------------------------

        if (trim(aoflux_grid) == 'atm') then
          ! compute o2x_ax for flux_atmocn, will be updated before atm merge
          ! do not use fractions because fractions here are NOT consistent
          ! with fractions in atm_mrg
          if (ocn_c2_atm) then 
            call prep_atm_calc_o2x_ax(timer='CPL:atmoca_ocn2atm')
          end if
          call t_drvstartf ('CPL:atmocna_fluxa',barrier=mpicom_CPLID)
          do exi = 1,num_inst_xao
            eai = mod((exi-1),num_inst_atm) + 1
            eoi = mod((exi-1),num_inst_ocn) + 1
            efi = mod((exi-1),num_inst_frc) + 1
            a2x_ax => component_get_c2x_cx(atm(eai))
            o2x_ax => prep_atm_get_o2x_ax()    ! array over all instances
            xao_ax => prep_aoflux_get_xao_ax() ! array over all instances
            call seq_flux_atmocn_mct(infodata, tod, dtime,  & 
                                     a2x_ax, o2x_ax(eoi), xao_ax(exi))
          end do
          call t_drvstopf  ('CPL:atmocna_fluxa')

          if (atm_c2_ocn) then
            call prep_aoflux_calc_xao_ox(timer='CPL:atmocna_atm2ocn')
          end if
        endif  ! aoflux_grid

          !----------------------------------------------------------
          !| atm/ocn flux on ocn grid ((cesm1_orig, cesm1_orig_tight, cesm1_mod
          !or cesm1_mod_tight) and aoflux='ocn')
          !----------------------------------------------------------

          if (trim(aoflux_grid) == 'ocn') then
            call t_drvstartf ('CPL:atmocnp_fluxo',barrier=mpicom_CPLID)
            do exi = 1,num_inst_xao
               eai = mod((exi-1),num_inst_atm) + 1
               eoi = mod((exi-1),num_inst_ocn) + 1
               efi = mod((exi-1),num_inst_frc) + 1
               a2x_ox => prep_ocn_get_a2x_ox()
               o2x_ox => component_get_c2x_cx(ocn(eoi))
               xao_ox => prep_aoflux_get_xao_ox()
               call seq_flux_atmocn_mct(infodata, tod, dtime, & 
                                        a2x_ox(eai), o2x_ox, xao_ox(exi))
            end do
            call t_drvstopf  ('CPL:atmocnp_fluxo')
            !         else if (trim(aoflux_grid) == 'atm') then
            !            !--- compute later ---
            !
            !         else if (trim(aoflux_grid) == 'exch') then
            !            xao_ax   => prep_aoflux_get_xao_ax()
            !            xao_ox   => prep_aoflux_get_xao_ox()
            !
            !            call t_drvstartf
            !            ('CPL:atmocnp_fluxe',barrier=mpicom_CPLID)
            !            call seq_flux_atmocnexch_mct( infodata, atm(eai),
            !            ocn(eoi), &
            !                 fractions_ax(efi), fractions_ox(efi),
            !                 xao_ax(exi), xao_ox(exi) )
            !            call t_drvstopf  ('CPL:atmocnp_fluxe')
          endif  ! aoflux_grid

          !----------------------------------------------------------
          !| ocn prep-merge (cesm1_mod or cesm1_mod_tight)
          !----------------------------------------------------------

          if (ocn_prognostic) then
            if (cpl_seq_option == 'CESM1_MOD' .or. &
                cpl_seq_option == 'CESM1_MOD_TIGHT') then

              xao_ox => prep_aoflux_get_xao_ox()
              call prep_ocn_mrg(infodata, fractions_ox, xao_ox=xao_ox, & 
                                timer_mrg='CPL:atmocnp_mrgx2o')

              ! Accumulate ocn inputs - form partial sum of tavg ocn inputs
              ! (virtual "send" to ocn)
              call prep_ocn_accum(timer='CPL:atmocnp_accum')
            endif
          endif

          !----------------------------------------------------------
          !| ocn albedos (cesm1_orig, cesm1_orig_tight, cesm1_mod or
          !cesm1_mod_tight)
          !  (MUST BE AFTER prep_ocn_mrg for swnet to ocn to be computed
          !  properly
          !----------------------------------------------------------

          call t_drvstartf ('CPL:atmocnp_ocnalb', barrier=mpicom_CPLID)
          do exi = 1,num_inst_xao
            efi = mod((exi-1),num_inst_frc) + 1
            eai = mod((exi-1),num_inst_atm) + 1
            xao_ox => prep_aoflux_get_xao_ox()   ! array over all instances
            a2x_ox => prep_ocn_get_a2x_ox()
            call seq_flux_ocnalb_mct(infodata, ocn(1), a2x_ox(eai),  &
                                     fractions_ox(efi), xao_ox(exi))
          end do
          call t_drvstopf  ('CPL:atmocnp_ocnalb')

          !----------------------------------------------------------
          !| ocn budget (cesm1_orig, cesm1_orig_tight, cesm1_mod or
          !cesm1_mod_tight)
          !----------------------------------------------------------

          if (do_budgets) then
            call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:BUDGET0_BARRIER')
            call t_drvstartf ('CPL:BUDGET0',budget=.true.,barrier=mpicom_CPLID)
            xao_ox => prep_aoflux_get_xao_ox() ! array over all instances
            call seq_diag_ocn_mct(ocn(ens1), xao_ox(1), fractions_ox(ens1), infodata, &
                 do_o2x=.true., do_x2o=.true., do_xao=.true.)
            call t_drvstopf ('CPL:BUDGET0',budget=.true.)
          end if

          if (drv_threading) then 
            call seq_comm_setnthreads(nthreads_GLOID)
          end if
          call t_drvstopf  ('CPL:ATMOCNP',cplrun=.true.,hashint=hashint(7))
      end if

      !----------------------------------------------------------
      !| LND RECV-POST
      !----------------------------------------------------------

      if (lnd_present .and. lndrun_alarm) then
        !----------------------------------------------------------
        !| lnd -> cpl
        !----------------------------------------------------------
        if (iamin_CPLALLLNDID) then
          call component_exch(lnd, flow='c2x', infodata=infodata, &
                              infodata_string='lnd2cpl_run', &
                              mpicom_barrier=mpicom_CPLALLLNDID, &
                              run_barriers=run_barriers, &
                              timer_barrier='CPL:L2C_BARRIER', & 
                              timer_comp_exch='CPL:L2C', &
                              timer_map_exch='CPL:l2c_lndl2lndx', &
                              timer_infodata_exch='lnd2cpl_run')
        end if

        !----------------------------------------------------------
        !| lnd post
        !----------------------------------------------------------

        if (iamin_CPLID) then
          call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:LNDPOST_BARRIER')
          call t_drvstartf('CPL:LNDPOST',cplrun=.true.,barrier=mpicom_CPLID)
          if (drv_threading) then 
            call seq_comm_setnthreads(nthreads_CPLID)
          end if
          call component_diag(infodata, lnd, flow='c2x', comment='recv lnd', &
               info_debug=info_debug, timer_diag='CPL:lndpost_diagav')

          ! Accumulate rof and glc inputs (module variables in prep_rof_mod
          ! and prep_glc_mod)
          if (lnd_c2_rof) then
            call prep_rof_accum(timer='CPL:lndpost_accl2r')
          end if
          if (lnd_c2_glc) then
            call prep_glc_accum(timer='CPL:lndpost_accl2g' )
          end if

          if (drv_threading) then 
            call seq_comm_setnthreads(nthreads_GLOID)
          end if
          call t_drvstopf  ('CPL:LNDPOST',cplrun=.true.)
        endif
      endif
      !----------------------------------------------------------
      !| GLC SETUP-SEND
      !----------------------------------------------------------
      if (glc_present .and. glcrun_alarm) then
        !----------------------------------------------------
        !| glc prep-merge
        !----------------------------------------------------
        if (iamin_CPLID .and. glc_prognostic) then
          call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:GLCPREP_BARRIER')
          call t_drvstartf ('CPL:GLCPREP',cplrun=.true.,barrier=mpicom_CPLID)
          if (drv_threading) then
            call seq_comm_setnthreads(nthreads_CPLID)
          end if
          if (lnd_c2_glc) then
            ! NOTE - only create appropriate input to glc if the avg_alarm
            ! is on
            if (glcrun_avg_alarm) then
              call prep_glc_accum_avg(timer='CPL:glcprep_avg')
              lnd2glc_averaged_now = .true.
              ! Note that l2x_gx is obtained from mapping the module
              ! variable l2gacc_lx
              call prep_glc_calc_l2x_gx(fractions_lx, timer='CPL:glcprep_lnd2glc')
              call prep_glc_mrg(infodata, fractions_gx, timer_mrg='CPL:glcprep_mrgx2g')
              call component_diag(infodata, glc, flow='x2c', comment='send glc', &
                                  info_debug=info_debug, timer_diag='CPL:glcprep_diagav')
            else
              call prep_glc_zero_fields()
            end if  ! glcrun_avg_alarm
          end if  ! lnd_c2_glc
          if (drv_threading) then
            call seq_comm_setnthreads(nthreads_GLOID)
          end if
          call t_drvstopf  ('CPL:GLCPREP',cplrun=.true.)
        end if  ! iamin_CPLID .and. glc_prognostic
        ! Set the infodata field on all tasks (not just those with
        ! iamin_CPLID).
        if (glc_prognostic) then
          if (glcrun_avg_alarm) then
            call seq_infodata_PutData(infodata, glc_valid_input=.true.)
          else
            call seq_infodata_PutData(infodata, glc_valid_input=.false.)
          end if
        end if
        !----------------------------------------------------
        !| cpl -> glc
        !----------------------------------------------------
        if (iamin_CPLALLGLCID .and. glc_prognostic) then
          call component_exch(glc, flow='x2c', &
               infodata=infodata, infodata_string='cpl2glc_run', &
               mpicom_barrier=mpicom_CPLALLGLCID, run_barriers=run_barriers, &
               timer_barrier='CPL:C2G_BARRIER', timer_comp_exch='CPL:C2G', &
               timer_map_exch='CPL:c2g_glcx2glcg', &
               timer_infodata_exch='CPL:c2g_infoexch')
        end if
      end if
      !----------------------------------------------------------
      !| ROF RECV-POST
      !----------------------------------------------------------
      if (rof_present .and. rofrun_alarm) then
        !----------------------------------------------------------
        !| rof -> cpl
        !----------------------------------------------------------
        if (iamin_CPLALLROFID) then
          call component_exch(rof, flow='c2x', &
               infodata=infodata, infodata_string='rof2cpl_run', &
               mpicom_barrier=mpicom_CPLALLROFID, run_barriers=run_barriers, &
               timer_barrier='CPL:R2C_BARRIER', timer_comp_exch='CPL:R2C', &
               timer_map_exch='CPL:r2c_rofr2rofx', &
               timer_infodata_exch='CPL:r2c_infoexch')
        end if
        !----------------------------------------------------------
        !| rof post
        !----------------------------------------------------------
        if (iamin_CPLID) then
          call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:ROFPOST_BARRIER')
          call t_drvstartf('CPL:ROFPOST',cplrun=.true.,barrier=mpicom_CPLID)
          if (drv_threading) then 
            call seq_comm_setnthreads(nthreads_CPLID)
          end if

          call component_diag(infodata, rof, flow='c2x', comment= 'recv rof', &
               info_debug=info_debug, timer_diag='CPL:rofpost_diagav')

          if (rof_c2_lnd) then
            call prep_lnd_calc_r2x_lx(timer='CPL:rofpost_rof2lnd')
          endif

          if (rof_c2_ice) then
            call prep_ice_calc_r2x_ix(timer='CPL:rofpost_rof2ice')
          endif

          if (rof_c2_ocn) then
            call prep_ocn_calc_r2x_ox(timer='CPL:rofpost_rof2ocn')
          endif

          call t_drvstopf  ('CPL:ROFPOST', cplrun=.true.)
        end if
      end if

      if (rof_present) then
        if (iamin_CPLID) then
          call cime_comp_barriers(mpicom=mpicom_CPLID, timer='DRIVER_ROFPOST_BARRIER')
          call t_drvstartf('DRIVER_ROFPOST',cplrun=.true.,barrier=mpicom_CPLID)
          if (drv_threading) then 
            call seq_comm_setnthreads(nthreads_CPLID)
          end if
          if (do_hist_r2x) then
            call t_drvstartf ('driver_rofpost_histaux', barrier=mpicom_CPLID)
            ! Write coupler's hr2x file at 24 hour marks, 
            ! and at the end of the run interval, even if that's not at a 24
            ! hour mark.
            write_hist_alarm = t24hr_alarm .or. stop_alarm
            do eri = 1,num_inst_rof
              inst_suffix =  component_get_suffix(rof(eri))
              call seq_hist_writeaux(infodata, EClock_d, rof(eri), flow='c2x', &
                   aname='r2x',dname='domrb',inst_suffix=trim(inst_suffix), &
                   nx=rof_nx, ny=rof_ny, nt=1, write_now=write_hist_alarm)
            end do
            call t_drvstopf ('driver_rofpost_histaux')
          end if
          call t_drvstopf  ('DRIVER_ROFPOST', cplrun=.true.)
        end if
      end if

      !----------------------------------------------------------
      !| Budget with old fractions
      !----------------------------------------------------------

      ! WJS (2-17-11): I am just using the first instance for the budgets
      ! because we
      ! don't expect budgets to be conserved for our case (I case). Also note
      ! that we
      ! don't expect budgets to be conserved for the interactive ensemble use
      ! case either.
      ! tcraig (aug 2012): put this after rof->cpl so the budget sees the new
      ! r2x_rx.
      ! it will also use the current r2x_ox here which is the value from the
      ! last timestep
      ! consistent with the ocean coupling

      if (iamin_CPLID .and. do_budgets) then
        call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:BUDGET1_BARRIER')
        call t_drvstartf('CPL:BUDGET1',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
        if (lnd_present) then
          call seq_diag_lnd_mct(lnd(ens1), fractions_lx(ens1), infodata, &
               do_l2x=.true., do_x2l=.true.)
         end if
         if (rof_present) then
           call seq_diag_rof_mct(rof(ens1), fractions_rx(ens1), infodata)
         end if
         if (ice_present) then
           call seq_diag_ice_mct(ice(ens1), fractions_ix(ens1), infodata, &
                do_x2i=.true.)
         end if
         call t_drvstopf  ('CPL:BUDGET1',cplrun=.true.,budget=.true.)
      end if
      !----------------------------------------------------------
      !| ICE RECV-POST
      !----------------------------------------------------------
      if (ice_present .and. icerun_alarm) then
        !----------------------------------------------------------
        !| ice -> cpl
        !----------------------------------------------------------
        if (iamin_CPLALLICEID) then
          call component_exch(ice, flow='c2x', &
               infodata=infodata, infodata_string='ice2cpl_run', &
               mpicom_barrier=mpicom_CPLALLICEID, run_barriers=run_barriers, &
               timer_barrier='CPL:I2C_BARRIER', timer_comp_exch='CPL:I2C', &
               timer_map_exch='CPL:i2c_icei2icex', &
               timer_infodata_exch='CPL:i2c_infoexch')
          end if
          !----------------------------------------------------------
          !| ice post
          !----------------------------------------------------------
          if (iamin_CPLID) then
            call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:ICEPOST_BARRIER')
            call t_drvstartf('CPL:ICEPOST',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) then 
              call seq_comm_setnthreads(nthreads_CPLID)
            end if
            call component_diag(infodata, ice, flow='c2x', comment= 'recv ice', &
                 info_debug=info_debug, timer_diag='CPL:icepost_diagav')

            if (drv_threading) then 
              call seq_comm_setnthreads(nthreads_GLOID)
            end if
            call t_drvstopf('CPL:ICEPOST',cplrun=.true.)
          end if
      end if
      !----------------------------------------------------------
      !| Update fractions based on new ice fractions
      !----------------------------------------------------------
      if (iamin_CPLID) then
        call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:FRACSET_BARRIER')
        call t_drvstartf ('CPL:FRACSET',cplrun=.true., barrier=mpicom_CPLID)
        if (drv_threading) then 
          call seq_comm_setnthreads(nthreads_CPLID)
        end if 
        call t_drvstartf ('CPL:fracset_fracset',barrier=mpicom_CPLID)

        do efi = 1,num_inst_frc
          eii = mod((efi-1),num_inst_ice) + 1
          call seq_frac_set(infodata, ice(eii), &
               fractions_ax(efi), fractions_ix(efi), fractions_ox(efi))
        end do
        call t_drvstopf('CPL:fracset_fracset')
        if (drv_threading) then 
          call seq_comm_setnthreads(nthreads_GLOID)
        end if
        call t_drvstopf('CPL:FRACSET',cplrun=.true.)
      end if
      !----------------------------------------------------------
      !| ATM/OCN SETUP (rasm_option2)
      !----------------------------------------------------------
      if ((trim(cpl_seq_option) == 'RASM_OPTION2') .and. &
           iamin_CPLID .and. ocn_present) then
        call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:ATMOCN2_BARRIER')
        call t_drvstartf ('CPL:ATMOCN2',cplrun=.true.,barrier=mpicom_CPLID)
        if (drv_threading) then 
          call seq_comm_setnthreads(nthreads_CPLID)
        end if
        if (ocn_prognostic) then
          ! Map ice to ocn
          if (ice_c2_ocn) then 
            call prep_ocn_calc_i2x_ox(timer='CPL:atmocnp_ice2ocn')
          end if
          ! Map wav to ocn
          if (wav_c2_ocn) call prep_ocn_calc_w2x_ox(timer='CPL:atmocnp_wav2ocn')
        end if
        !----------------------------------------------------------
        !| atm/ocn flux on atm grid (rasm_option2 and aoflux_grid='atm')
        !----------------------------------------------------------
        if (trim(aoflux_grid) == 'atm') then
          ! compute o2x_ax for flux_atmocn, will be updated before atm merge
          ! can use fractions because fractions here are consistent with
          ! fractions in atm_mrg
          if (ocn_c2_atm) then 
            call prep_atm_calc_o2x_ax(fractions_ox,timer='CPL:atmoca_ocn2atm')
          end if
          call t_drvstartf ('CPL:atmocna_fluxa',barrier=mpicom_CPLID)
          do exi = 1,num_inst_xao
            eai = mod((exi-1),num_inst_atm) + 1
            eoi = mod((exi-1),num_inst_ocn) + 1
            efi = mod((exi-1),num_inst_frc) + 1
            a2x_ax => component_get_c2x_cx(atm(eai))
            o2x_ax => prep_atm_get_o2x_ax()    ! array over all instances
            xao_ax => prep_aoflux_get_xao_ax() ! array over all instances
            call seq_flux_atmocn_mct(infodata, tod, dtime, & 
                                     a2x_ax, o2x_ax(eoi), xao_ax(exi))
          end do
          call t_drvstopf('CPL:atmocna_fluxa')
          if (atm_c2_ocn) then 
            call prep_aoflux_calc_xao_ox(timer='CPL:atmocna_atm2ocn')
          end if
        end if  ! aoflux_grid

        !----------------------------------------------------------
        !| atm/ocn flux on ocn grid (rasm_option2 and aoflux_grid='ocn')
        !----------------------------------------------------------

        if (trim(aoflux_grid) == 'ocn') then
          call t_drvstartf ('CPL:atmocnp_fluxo',barrier=mpicom_CPLID)
          do exi = 1,num_inst_xao
            eai = mod((exi-1),num_inst_atm) + 1
            eoi = mod((exi-1),num_inst_ocn) + 1
            efi = mod((exi-1),num_inst_frc) + 1
            a2x_ox => prep_ocn_get_a2x_ox()
            o2x_ox => component_get_c2x_cx(ocn(eoi))
            xao_ox => prep_aoflux_get_xao_ox()
            call seq_flux_atmocn_mct(infodata, tod, dtime, & 
                                     a2x_ox(eai), o2x_ox, xao_ox(exi))
          end do
          call t_drvstopf  ('CPL:atmocnp_fluxo')
        endif  ! aoflux_grid
        !----------------------------------------------------------
        !| ocn prep-merge (rasm_option2)
        !----------------------------------------------------------
        xao_ox => prep_aoflux_get_xao_ox()
        call prep_ocn_mrg(infodata, fractions_ox, xao_ox=xao_ox, &
                          timer_mrg='CPL:atmocnp_mrgx2o')

        ! Accumulate ocn inputs - form partial sum of tavg ocn inputs (virtual
        ! "send" to ocn)
        call prep_ocn_accum(timer='CPL:atmocnp_accum')
        !----------------------------------------------------------
        !| ocn albedos (rasm_option2)
        !  (MUST BE AFTER prep_ocn_mrg for swnet to ocn to be computed
        !  properly
        !----------------------------------------------------------
        call t_drvstartf ('CPL:atmocnp_ocnalb', barrier=mpicom_CPLID)
        do exi = 1,num_inst_xao
          efi = mod((exi-1),num_inst_frc) + 1
          eai = mod((exi-1),num_inst_atm) + 1
          xao_ox => prep_aoflux_get_xao_ox()   ! array over all instances
          a2x_ox => prep_ocn_get_a2x_ox()
          call seq_flux_ocnalb_mct(infodata, ocn(1), & 
                                   a2x_ox(eai), fractions_ox(efi), xao_ox(exi))
        end do
        call t_drvstopf  ('CPL:atmocnp_ocnalb')
        !----------------------------------------------------------
        !| ocn budget (rasm_option2)
        !----------------------------------------------------------
        if (do_budgets) then
          call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:BUDGET0_BARRIER')
          call t_drvstartf ('CPL:BUDGET0',budget=.true.,barrier=mpicom_CPLID)
          xao_ox => prep_aoflux_get_xao_ox() ! array over all instances
          call seq_diag_ocn_mct(ocn(ens1), xao_ox(1), fractions_ox(ens1), infodata, &
                                do_o2x=.true., do_x2o=.true., do_xao=.true.)
          call t_drvstopf ('CPL:BUDGET0',budget=.true.)
        end if
        if (drv_threading) then
          call seq_comm_setnthreads(nthreads_GLOID)
        end if
        call t_drvstopf  ('CPL:ATMOCN2',cplrun=.true.)
      endif
      !----------------------------------------------------------
      !| OCN SETUP-SEND (rasm_option2)
      !----------------------------------------------------------
      if ((trim(cpl_seq_option) == 'RASM_OPTION2'  ) .and. &
           ocn_present .and. ocnrun_alarm) then
        !----------------------------------------------------
        ! "startup" wait (rasm_option2)
        !----------------------------------------------------
        if (iamin_CPLALLOCNID) then
          ! want to know the time the ocean pes waited for the cpl pes
          ! at the first ocnrun_alarm, min ocean wait is wait time
          ! do not use t_barrierf here since it can be "off", use mpi_barrier
          do eoi = 1,num_inst_ocn
            if (ocn(eoi)%iamin_compid) then 
              call t_drvstartf ('CPL:C2O_INITWAIT')
            end if
          end do
          call mpi_barrier(mpicom_CPLALLOCNID,ierr)
          do eoi = 1,num_inst_ocn
            if (ocn(eoi)%iamin_compid) then 
              call t_drvstopf  ('CPL:C2O_INITWAIT')
            end if
          end do
          cpl2ocn_first = .false.
        end if
        !----------------------------------------------------
        !| ocn average (rasm_option2)
        !----------------------------------------------------
        if (iamin_CPLID .and. ocn_prognostic) then
          call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:OCNPRE2_BARRIER')
          call t_drvstartf ('CPL:OCNPRE2',cplrun=.true.,barrier=mpicom_CPLID)
          if (drv_threading) then 
            call seq_comm_setnthreads(nthreads_CPLID)
          end if
          ! finish accumulating ocean inputs
          ! reset the value of x2o_ox with the value in x2oacc_ox
          ! (module variable in prep_ocn_mod)
          call prep_ocn_accum_avg(timer_accum='CPL:ocnprep_avg')
          call component_diag(infodata, ocn, flow='x2c', comment= 'send ocn', &
                              info_debug=info_debug, timer_diag='CPL:ocnprep_diagav')
          if (drv_threading) then 
            call seq_comm_setnthreads(nthreads_GLOID)
          end if
          call t_drvstopf  ('CPL:OCNPRE2',cplrun=.true.)
        endif
        !----------------------------------------------------
        !| cpl -> ocn (rasm_option2)
        !----------------------------------------------------
        if (iamin_CPLALLOCNID .and. ocn_prognostic) then
          call component_exch(ocn, flow='x2c', &
               infodata=infodata, infodata_string='cpl2ocn_run', &
               mpicom_barrier=mpicom_CPLALLOCNID, run_barriers=run_barriers, &
               timer_barrier='CPL:C2O2_BARRIER', timer_comp_exch='CPL:C2O2', &
               timer_map_exch='CPL:c2o2_ocnx2ocno', &
               timer_infodata_exch='CPL:c2o2_infoexch')
          end if
      end if
      !----------------------------------------------------------
      !| ATM SETUP-SEND
      !----------------------------------------------------------
      if (atm_present .and. atmrun_alarm) then
        !----------------------------------------------------------
        !| atm prep-merge
        !----------------------------------------------------------
        if (iamin_CPLID .and. atm_prognostic) then
          call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:ATMPREP_BARRIER')
          call t_drvstartf ('CPL:ATMPREP',cplrun=.true.,barrier=mpicom_CPLID)
          if (drv_threading) then 
            call seq_comm_setnthreads(nthreads_CPLID)
          end if
          if (ocn_c2_atm) then
            if (trim(aoflux_grid) == 'ocn') then
              ! map xao_ox states and fluxes to xao_ax if fluxes were
              ! computed on ocn grid
              call prep_aoflux_calc_xao_ax(fractions_ox, flds='states_and_fluxes', &
                                           timer='CPL:atmprep_xao2atm')
            end if
            ! recompute o2x_ax now for the merge with fractions associated
            ! with merge
            call prep_atm_calc_o2x_ax(fractions_ox, timer='CPL:atmprep_ocn2atm')
            ! map xao_ox albedos to the atm grid, these are always computed
            ! on the ocean grid
            call prep_aoflux_calc_xao_ax(fractions_ox, flds='albedos', & 
                                         timer='CPL:atmprep_alb2atm')
          end if

          if (ice_c2_atm) then
            call prep_atm_calc_i2x_ax(fractions_ix, timer='CPL:atmprep_ice2atm')
          end if

          if (lnd_c2_atm) then
            call prep_atm_calc_l2x_ax(fractions_lx, timer='CPL:atmprep_lnd2atm')
          end if

          if (associated(xao_ax)) then
            call prep_atm_mrg(infodata, fractions_ax, xao_ax=xao_ax, &
                              timer_mrg='CPL:atmprep_mrgx2a')
          end if
          call component_diag(infodata, atm, flow='x2c', & 
                              comment= 'send atm', &
                              info_debug=info_debug, &
                              timer_diag='CPL:atmprep_diagav')
          call t_drvstopf('CPL:ATMPREP',cplrun=.true.)
          if (drv_threading) then  
            call seq_comm_setnthreads(nthreads_GLOID)
           end if
        end if
        !----------------------------------------------------------
        !| cpl -> atm
        !----------------------------------------------------------
        if (iamin_CPLALLATMID .and. atm_prognostic) then
          call component_exch(atm, flow='x2c', infodata=infodata, &
                              infodata_string='cpl2atm_run', &
                              mpicom_barrier=mpicom_CPLALLATMID, & 
                              run_barriers=run_barriers, &
                              timer_barrier='CPL:C2A_BARRIER', &
                              timer_comp_exch='CPL:C2A', &
                              timer_map_exch='CPL:c2a_atmx2atmg', &
                              timer_infodata_exch='CPL:c2a_infoexch')
        end if
      end if
      !----------------------------------------------------------
      !| RUN OCN MODEL (NOT cesm1_orig_tight or cesm1_mod_tight)
      !----------------------------------------------------------
      if ((trim(cpl_seq_option) /= 'CESM1_ORIG_TIGHT' .and. &
            trim(cpl_seq_option) /= 'CESM1_MOD_TIGHT'   ) .and. &
            ocn_present .and. ocnrun_alarm) then
        call component_run(Eclock_o, ocn, ocn_run, infodata, &
              seq_flds_x2c_fluxes=seq_flds_x2o_fluxes, &
              seq_flds_c2x_fluxes=seq_flds_o2x_fluxes, &
              comp_prognostic=ocn_prognostic, comp_num=comp_num_ocn, &
              timer_barrier= 'CPL:OCN_RUN_BARRIER', timer_comp_run='CPL:OCN_RUN', &
              run_barriers=run_barriers, ymd=ymd, tod=tod,comp_layout=ocn_layout)
      end if
      !----------------------------------------------------------
      !| RUN ATM MODEL
      !----------------------------------------------------------
      if (atm_present .and. atmrun_alarm) then
        call component_run(Eclock_a, atm, atm_run, infodata, &
              seq_flds_x2c_fluxes=seq_flds_x2a_fluxes, &
              seq_flds_c2x_fluxes=seq_flds_a2x_fluxes, &
              comp_prognostic=atm_prognostic, comp_num=comp_num_atm, &
              timer_barrier= 'CPL:ATM_RUN_BARRIER', timer_comp_run='CPL:ATM_RUN', &
               run_barriers=run_barriers, ymd=ymd, tod=tod, comp_layout=atm_layout)
      end if
      !----------------------------------------------------------
      !| RUN GLC MODEL
      !----------------------------------------------------------
      if (glc_present .and. glcrun_alarm) then
        call component_run(Eclock_g, glc, glc_run, infodata, &
             seq_flds_x2c_fluxes=seq_flds_x2g_fluxes, &
             seq_flds_c2x_fluxes=seq_flds_g2x_fluxes, &
             comp_prognostic=glc_prognostic, comp_num=comp_num_glc, &
             timer_barrier= 'CPL:GLC_RUN_BARRIER', timer_comp_run='CPL:GLC_RUN', &
             run_barriers=run_barriers, ymd=ymd, tod=tod,comp_layout=glc_layout)
      end if
      !----------------------------------------------------------
      !| WAV RECV-POST
      !----------------------------------------------------------
      if (wav_present .and. wavrun_alarm) then
        !----------------------------------------------------------
        !| wav -> cpl
        !----------------------------------------------------------
        if (iamin_CPLALLWAVID) then
          call component_exch(wav, flow='c2x', infodata=infodata, &
                              infodata_string='wav2cpl_run', &
                              mpicom_barrier=mpicom_CPLALLWAVID, &
                              run_barriers=run_barriers, &
                              timer_barrier='CPL:W2C_BARRIER', & 
                              timer_comp_exch='CPL:W2C', &
                              timer_map_exch='CPL:w2c_wavw2wavx', &
                              timer_infodata_exch='CPL:w2c_infoexch')
        end if
        !----------------------------------------------------------
        !| wav post
        !----------------------------------------------------------
        if (iamin_CPLID) then
          call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:WAVPOST_BARRIER')
          call t_drvstartf('CPL:WAVPOST',cplrun=.true.,barrier=mpicom_CPLID)
          if (drv_threading) then
            call seq_comm_setnthreads(nthreads_CPLID)
          end if
          call component_diag(infodata, wav, flow='c2x', comment= 'recv wav', &
               info_debug=info_debug, timer_diag='CPL:wavpost_diagav')
          if (drv_threading) then 
            call seq_comm_setnthreads(nthreads_GLOID)
          end if
          call t_drvstopf  ('CPL:WAVPOST',cplrun=.true.)
        end if
      end if
      !----------------------------------------------------------
      !| GLC RECV-POST
      !----------------------------------------------------------
      if (glc_present .and. glcrun_alarm) then
        !----------------------------------------------------------
        !| glc -> cpl
        !----------------------------------------------------------
        if (iamin_CPLALLGLCID) then
          call component_exch(glc, flow='c2x', infodata=infodata, &
                              infodata_string='glc2cpl_run', &
                              mpicom_barrier=mpicom_CPLALLGLCID, & 
                              run_barriers=run_barriers, &
                              timer_barrier='CPL:G2C_BARRIER', & 
                              timer_comp_exch='CPL:G2C', &
                              timer_map_exch='CPL:g2c_glcg2glcx', &
                              timer_infodata_exch='CPL:g2c_infoexch')
        end if
        !----------------------------------------------------------
        !| glc post
        !----------------------------------------------------------
        if (iamin_CPLID) then
          call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:GLCPOST_BARRIER')
          call t_drvstartf('CPL:GLCPOST',cplrun=.true.,barrier=mpicom_CPLID)
          if (drv_threading) then 
            call seq_comm_setnthreads(nthreads_CPLID)
          end if
          call component_diag(infodata, glc, flow='c2x', comment= 'recv glc', &
               info_debug=info_debug, timer_diag='CPL:glcpost_diagav')
          if (glc_c2_lnd) then
            call prep_lnd_calc_g2x_lx(timer='CPL:glcpost_glc2lnd')
          endif
          if (glc_c2_ice) then
            call prep_ice_calc_g2x_ix(timer='CPL:glcpost_glc2ice')
          endif
          if (glc_c2_ocn) then
            call prep_ocn_calc_g2x_ox(timer='CPL:glcpost_glc2ocn')
          endif
          if (drv_threading) then 
            call seq_comm_setnthreads(nthreads_GLOID)
          end if
          call t_drvstopf  ('CPL:GLCPOST',cplrun=.true.)
        end if
      end if
      !----------------------------------------------------------
      !| ATM RECV-POST
      !----------------------------------------------------------
      if (atm_present .and. atmrun_alarm) then
      !----------------------------------------------------------
      !| atm -> cpl
      !----------------------------------------------------------
        if (iamin_CPLALLATMID) then
          call component_exch(atm, flow='c2x', &
              infodata=infodata, infodata_string='atm2cpl_run', &
              mpicom_barrier=mpicom_CPLALLATMID, run_barriers=run_barriers, &
              timer_barrier='CPL:A2C_BARRIER', timer_comp_exch='CPL:A2C', &
              timer_map_exch='CPL:a2c_atma2atmx', timer_infodata_exch='CPL:a2c_infoexch')
        end if
        !----------------------------------------------------------
        !| atm post
        !----------------------------------------------------------
        if (iamin_CPLID) then
          call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:ATMPOST_BARRIER')
          call t_drvstartf ('CPL:ATMPOST',cplrun=.true.,barrier=mpicom_CPLID)
          if (drv_threading) then
            call seq_comm_setnthreads(nthreads_CPLID)
          end if
          call component_diag(infodata, atm, flow='c2x', comment= 'recv atm', &
               info_debug=info_debug, timer_diag='CPL:atmpost_diagav')
          if (drv_threading) then
            call seq_comm_setnthreads(nthreads_GLOID)
          end if
          call t_drvstopf  ('CPL:ATMPOST',cplrun=.true.)
        end if
      end if
      !----------------------------------------------------------
      !| Budget with new fractions
      !----------------------------------------------------------
      if (iamin_CPLID .and. do_budgets) then
        call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:BUDGET2_BARRIER')
        call t_drvstartf('CPL:BUDGET2',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
        if (atm_present) then
          call seq_diag_atm_mct(atm(ens1), fractions_ax(ens1), infodata, &
                                 do_a2x=.true., do_x2a=.true.)
        end if
        if (ice_present) then
          call seq_diag_ice_mct(ice(ens1), fractions_ix(ens1), infodata, &
                                do_i2x=.true.)
        end if
        call t_drvstopf  ('CPL:BUDGET2',cplrun=.true.,budget=.true.)
        call t_drvstartf('CPL:BUDGET3',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
        call seq_diag_accum_mct()
        call t_drvstopf  ('CPL:BUDGET3',cplrun=.true.,budget=.true.)
        call t_drvstartf ('CPL:BUDGETF',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
        if (.not. dead_comps) then
          call seq_diag_print_mct(EClock_d,stop_alarm,budget_inst, &
               budget_daily, budget_month, budget_ann, budget_ltann, budget_ltend)
        end if
        call seq_diag_zero_mct(EClock=EClock_d)
        call t_drvstopf  ('CPL:BUDGETF',cplrun=.true.,budget=.true.)
      end if
      !----------------------------------------------------------
      !| OCN RECV-POST (NOT cesm1_orig_tight and cesm1_mod_tight)
      !----------------------------------------------------------
      if ((trim(cpl_seq_option) /= 'CESM1_ORIG_TIGHT' .and. &
           trim(cpl_seq_option) /= 'CESM1_MOD_TIGHT'   ) .and. &
           ocn_present .and. ocnnext_alarm) then
        !----------------------------------------------------------
        !| ocn -> cpl (NOT cesm1_orig_tight and cesm1_mod_tight)
        !----------------------------------------------------------
        if (iamin_CPLALLOCNID) then
          call component_exch(ocn, flow='c2x', &
               infodata=infodata, infodata_string='ocn2cpl_run', &
               mpicom_barrier=mpicom_CPLALLOCNID, run_barriers=run_barriers, &
               timer_barrier='CPL:O2C_BARRIER', timer_comp_exch='CPL:O2C', &
               timer_map_exch='CPL:o2c_ocno2ocnx', timer_infodata_exch='CPL:o2c_infoexch')
        end if
        !----------------------------------------------------------
        !| ocn post (NOT cesm1_orig_tight and cesm1_mod_tight)
        !----------------------------------------------------------
        if (iamin_CPLID) then
          call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:OCNPOST_BARRIER')
          call t_drvstartf ('CPL:OCNPOST',cplrun=.true.,barrier=mpicom_CPLID)
          if (drv_threading) then
            call seq_comm_setnthreads(nthreads_CPLID)
          end if
          call component_diag(infodata, ocn, flow='c2x', comment= 'recv ocn', &
                              info_debug=info_debug, timer_diag='CPL:ocnpost_diagav')
          if (drv_threading) then 
            call seq_comm_setnthreads(nthreads_GLOID)
          end if
          call t_drvstopf  ('CPL:OCNPOST',cplrun=.true.)
        end if
      end if
      !----------------------------------------------------------
      !| Write driver restart file
      !----------------------------------------------------------
      if ( (restart_alarm .or. drv_pause) .and. iamin_CPLID) then
        call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:RESTART_BARRIER')
        call t_drvstartf ('CPL:RESTART',cplrun=.true.,barrier=mpicom_CPLID)
        if (drv_threading) then
          call seq_comm_setnthreads(nthreads_CPLID)
        end if
        if (iamroot_CPLID) then
          write(logunit,104) ' Write restart file at ',ymd,tod
          call shr_sys_flush(logunit)
        end if
        call seq_rest_write(EClock_d, seq_SyncClock, infodata,       &
             atm, lnd, ice, ocn, rof, glc, wav, esp,                 &
             fractions_ax, fractions_lx, fractions_ix, fractions_ox, &
             fractions_rx, fractions_gx, fractions_wx, trim(cpl_inst_tag))
        if (drv_threading) then 
          call seq_comm_setnthreads(nthreads_GLOID)
        end if
        call t_drvstopf  ('CPL:RESTART',cplrun=.true.)
      end if
      !----------------------------------------------------------
      !| Write history file, only AVs on CPLID
      !----------------------------------------------------------
      if (iamin_CPLID) then
        call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:HISTORY_BARRIER')
        call t_drvstartf ('CPL:HISTORY',cplrun=.true.,barrier=mpicom_CPLID)
        if (history_alarm) then
          if (drv_threading) then
            call seq_comm_setnthreads(nthreads_CPLID)
          end if
          if (iamroot_CPLID) then
            write(logunit,104) ' Write history file at ',ymd,tod
            call shr_sys_flush(logunit)
          end if
          call seq_hist_write(infodata, EClock_d, &
               atm, lnd, ice, ocn, rof, glc, wav, &
               fractions_ax, fractions_lx, fractions_ix, fractions_ox,     &
               fractions_rx, fractions_gx, fractions_wx, trim(cpl_inst_tag))
          if (drv_threading) then 
            call seq_comm_setnthreads(nthreads_GLOID)
          end if
        end if
      if (do_histavg) then
        call seq_hist_writeavg(infodata, EClock_d, &
             atm, lnd, ice, ocn, rof, glc, wav, histavg_alarm, &
             trim(cpl_inst_tag))
      end if
      if (do_hist_a2x) then
        do eai = 1,num_inst_atm
          inst_suffix =  component_get_suffix(atm(eai))
          if (trim(hist_a2x_flds) == 'all') then
             call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                  aname='a2x',dname='doma', inst_suffix=trim(inst_suffix), &
                  nx=atm_nx, ny=atm_ny, nt=ncpl)
          else
            call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                 aname='a2x',dname='doma', inst_suffix=trim(inst_suffix), &
                 nx=atm_nx, ny=atm_ny, nt=ncpl, flds=hist_a2x_flds)
          end if
        end do
      endif

      if (do_hist_a2x1hri .and. t1hr_alarm) then
        do eai = 1,num_inst_atm
          inst_suffix =  component_get_suffix(atm(eai))
          if (trim(hist_a2x1hri_flds) == 'all') then
            call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                 aname='a2x1hi',dname='doma',inst_suffix=trim(inst_suffix), &
                 nx=atm_nx, ny=atm_ny, nt=24)
          else
            call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                 aname='a2x1hi',dname='doma',inst_suffix=trim(inst_suffix), &
                nx=atm_nx, ny=atm_ny, nt=24, flds=hist_a2x1hri_flds)
          end if
        end do
      end if

      if (do_hist_a2x1hr) then
        do eai = 1,num_inst_atm
          inst_suffix =  component_get_suffix(atm(eai))
          if (trim(hist_a2x1hr_flds) == 'all') then
            call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                 aname='a2x1h',dname='doma',inst_suffix=trim(inst_suffix), &
                 nx=atm_nx, ny=atm_ny, nt=24, write_now=t1hr_alarm)
          else
            call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                 aname='a2x1h',dname='doma',inst_suffix=trim(inst_suffix), &
                 nx=atm_nx, ny=atm_ny, nt=24, write_now=t1hr_alarm, flds=hist_a2x1hr_flds)
          end if
        end do
      end if

      if (do_hist_a2x3hr) then
        do eai = 1,num_inst_atm
          inst_suffix =  component_get_suffix(atm(eai))
          if (trim(hist_a2x3hr_flds) == 'all') then
            call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                 aname='a2x3h',dname='doma',inst_suffix=trim(inst_suffix), &
                 nx=atm_nx, ny=atm_ny, nt=8, write_now=t3hr_alarm)
          else
            call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                 aname='a2x3h',dname='doma',inst_suffix=trim(inst_suffix), &
                 nx=atm_nx, ny=atm_ny, nt=8, write_now=t3hr_alarm, flds=hist_a2x3hr_flds)
          end if
        end do
      end if

      if (do_hist_a2x3hrp) then
        do eai = 1,num_inst_atm
          inst_suffix = component_get_suffix(atm(eai))
          if (trim(hist_a2x3hrp_flds) == 'all') then
            call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                 aname='a2x3h_prec',dname='doma',inst_suffix=trim(inst_suffix), &
                 nx=atm_nx, ny=atm_ny, nt=8, write_now=t3hr_alarm)
          else
            call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                 aname='a2x3h_prec',dname='doma',inst_suffix=trim(inst_suffix), &
                 nx=atm_nx, ny=atm_ny, nt=8, write_now=t3hr_alarm, flds=hist_a2x3hrp_flds)
          end if
        end do
      end if

      if (do_hist_a2x24hr) then
        do eai = 1,num_inst_atm
          inst_suffix = component_get_suffix(atm(eai))
          if (trim(hist_a2x24hr_flds) == 'all') then
            call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
            aname='a2x1d',dname='doma',inst_suffix=trim(inst_suffix), &
            nx=atm_nx, ny=atm_ny, nt=1, write_now=t24hr_alarm)
          else
             call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                  aname='a2x1d',dname='doma',inst_suffix=trim(inst_suffix), &
                  nx=atm_nx,ny=atm_ny,nt=1,write_now=t24hr_alarm, flds=hist_a2x24hr_flds)
          end if
        end do
      end if

      if (do_hist_l2x1yrg) then
        ! We use a different approach here than for other aux hist files:
        ! For other files, we let seq_hist_writeaux accumulate fields in time.
        ! However, if we stop in the middle of an accumulation period, these accumulated
        ! fields get reset (because they aren't written to the cpl restart file); this
        ! is potentially a problem for this year-long accumulation. Thus,
        ! here, we use the existing accumulated fields from prep_glc_mod, because those
        ! *do* continue properly through a restart.
        ! The logic here assumes that we average the lnd2glc fields exactly
        ! at the year boundary - no more and no less. If that's not the case,
        ! we're likely to be writing the wrong thing to these aux files, so we check
        ! that assumption here.
        if (t1yr_alarm .and. .not. lnd2glc_averaged_now) then
          write(logunit,*) 'ERROR: histaux_l2x1yrg requested;'
          write(logunit,*) 'it is the year boundary, '
          write(logunit,*) 'but lnd2glc fields were not averaged this time step.'
          write(logunit,*) 'One possible reason is that you are running '
          write(logunit,*) 'with a stub glc model.'
          write(logunit,*) '(It only works to request histaux_l2x1yrg '
          write(logunit,*) 'if running with a prognostic glc model.)'
          call shr_sys_abort(subname// &
          ' do_hist_l2x1yrg and t1yr_alarm are true, but lnd2glc_averaged_now is false')
        end if

        if (lnd2glc_averaged_now .and. .not. t1yr_alarm) then
          ! If we're averaging more frequently than yearly, then just
          ! writing the current values of the averaged fields once per year won't give
          ! the true annual averages.
          write(logunit,*) 'ERROR: histaux_l2x1yrg requested;'
          write(logunit,*) 'lnd2glc fields were averaged this time step, '
          write(logunit,*) 'but it is not the year boundary.'
          write(logunit,*) '(It only works to request histaux_l2x1yrg '
          write(logunit,*) 'if GLC_AVG_PERIOD is yearly.)'
          call shr_sys_abort(subname// &
               ' do_hist_l2x1yrg and lnd2glc_averaged_now are true, &
               & but t1yr_alarm is false')
        end if

        if (t1yr_alarm) then
          call seq_timemgr_EClockGetData( EClock_d, ECurrTime = etime_curr)
          ! We need to pass in tbnds1_offset because (unlike with most
          ! seq_hist_writeaux calls) here we don't call seq_hist_writeaux
          ! every time step, so the automatically determined lower time bound can be
          ! wrong. For typical runs with a noleap calendar, we want tbnds1_offset =
          ! -365. However, to determine this more generally, based on the
          ! calendar we're using, we call this shr_cal routine.
          call shr_cal_ymds2rday_offset(etime=etime_curr, &
                                        rdays_offset = tbnds1_offset, &
                                        years_offset = -1)
          do eli = 1,num_inst_lnd
            inst_suffix = component_get_suffix(lnd(eli))
            ! Use yr_offset=-1 so the file with fields from year 1 has
            ! time stamp 0001-01-01 rather than 0002-01-01, etc.
            call seq_hist_writeaux(infodata, EClock_d, lnd(eli), flow='c2x', &
                 aname='l2x1yr_glc',dname='doml',inst_suffix=trim(inst_suffix), &
                 nx=lnd_nx, ny=lnd_ny, nt=1, write_now=.true., &
                 tbnds1_offset = tbnds1_offset, yr_offset=-1, &
                 av_to_write=prep_glc_get_l2gacc_lx_one_instance(eli))
          end do
        end if
      end if
      if (do_hist_l2x) then
        do eli = 1,num_inst_lnd
          inst_suffix =  component_get_suffix(lnd(eli))
          call seq_hist_writeaux(infodata, EClock_d, lnd(eli), flow='c2x', &
               aname='l2x',dname='doml',inst_suffix=trim(inst_suffix),  &
               nx=lnd_nx, ny=lnd_ny, nt=ncpl)
        end do
      end if
      call t_drvstopf  ('CPL:HISTORY',cplrun=.true.)
    end if
    !----------------------------------------------------------
    !| RUN ESP MODEL
    !----------------------------------------------------------
    if (esp_present .and. esprun_alarm) then
      ! Make sure that all couplers are here in multicoupler mode before
      ! running ESP component
      if (num_inst_driver > 1) then
        call mpi_barrier(global_comm, ierr)
      end if
      call component_run(Eclock_e, esp, esp_run, infodata, &
           comp_prognostic=esp_prognostic, comp_num=comp_num_esp, &
           timer_barrier= 'CPL:ESP_RUN_BARRIER', timer_comp_run='CPL:ESP_RUN', &
           run_barriers=run_barriers, ymd=ymd, tod=tod,comp_layout=esp_layout)
      !---------------------------------------------------------------------
      !| ESP computes resume options for other components -- update everyone
      !---------------------------------------------------------------------
      call seq_infodata_exchange(infodata, CPLALLESPID, 'esp2cpl_run')
    end if
    !----------------------------------------------------------
    !| RESUME (read restart) if signaled
    !----------------------------------------------------------
    call seq_infodata_GetData(infodata, cpl_resume=drv_resume)
    if (len_trim(drv_resume) > 0) then
      if (iamroot_CPLID) then
        write(logunit,103) subname,' Reading restart (resume) file',trim(drv_resume)
        call shr_sys_flush(logunit)
      end if
      if (iamin_CPLID) then
        call seq_rest_read(drv_resume, infodata,                          &
             atm, lnd, ice, ocn, rof, glc, wav, esp,                      &
             fractions_ax, fractions_lx, fractions_ix, fractions_ox,      &
             fractions_rx, fractions_gx, fractions_wx)
      end if
      ! Clear the resume file so we don't try to read it again
      drv_resume = ' '
      call seq_infodata_PutData(infodata, cpl_resume=drv_resume)
    end if
    !----------------------------------------------------------
    !| Timing and memory diagnostics
    !----------------------------------------------------------
    call t_drvstartf ('CPL:TSTAMP_WRITE',cplrun=.true.)
    if (tod == 0 .or. info_debug > 1) then
      if (iamroot_CPLID) then
        call date_and_time(dstr,tstr)
        Time_estep = mpi_wtime()
        cktime = time_estep-time_bstep
        cktime_acc(1) = cktime_acc(1) + cktime
        cktime_cnt(1) = cktime_cnt(1) + 1
#ifndef CPL_BYPASS
        write(logunit,101) ' tStamp_write: model date = ',ymd,tod, &
                           ' wall clock = ',dstr(1:4),'-',dstr(5:6),'-',dstr(7:8),' ',&
                           tstr(1:2),':',tstr(3:4),':',tstr(5:6), &
                           ' avg dt = ',cktime_acc(1)/cktime_cnt(1),' dt = ',cktime
#endif
        Time_bstep = mpi_wtime()
        call shr_sys_flush(logunit)
        if(cktime > max_cplstep_time .and. max_cplstep_time > 0.0) then
          call shr_sys_abort(subname//'Wall clock time exceeds max_cplstep_time')
        else if(max_cplstep_time < -0.05) then
          ! if max_cplstep_time is < 0 we use abs(max_cplstep_time)
          ! times the initial cktime value as a threshhold
          max_cplstep_time = -(max_cplstep_time)*cktime
        end if
      end if
    end if

    if (tod == 0 .and. wall_time_limit > 0.0_r8 .and. .not. force_stop) then
      time_erun = mpi_wtime()
      ! time_*run is seconds, wall_time_limit is hours
      wall_time = (time_erun - time_brun) / 3600._r8   ! convert secs to hrs
      write(logunit,109) subname//' check wall_time_limit: ',wall_time,wall_time_limit
      if (wall_time > wall_time_limit) then
        force_stop = .true.
        force_stop_tod = 0
        if (trim(force_stop_at) == 'month') then
          call shr_cal_date2ymd(ymd,year,month,day)
          month = month + 1
          do while (month > 12)
            month = month - 12
            year = year + 1
          end do
          call shr_cal_ymd2date(year,month,1,force_stop_ymd)
        else if (trim(force_stop_at) == 'year') then  ! next year
          call shr_cal_date2ymd(ymd,year,month,day)
          call shr_cal_ymd2date(year+1,1,1,force_stop_ymd)
        else if (trim(force_stop_at) == 'day') then   ! next day
          ymdtmp = ymd
          call shr_cal_advDateInt(1,'days',ymdtmp,0,force_stop_ymd,todtmp,calendar)
        else    ! day is default
          ymdtmp = ymd
          call shr_cal_advDateInt(1,'days',ymdtmp,0,force_stop_ymd,todtmp,calendar)
        end if
        write(logunit,108) subname//'reached wall_time_limit(hours)=',wall_time_limit, &
              ' :stop at ',force_stop_ymd
      end if
    end if
#ifndef CPL_BYPASS
    if (tod == 0 .or. info_debug > 1) then
      !! Report on memory usage
      !! For now, just look at the first instance of each component
      if ( iamroot_CPLID .or. &
           ocn(ens1)%iamroot_compid .or. &
           atm(ens1)%iamroot_compid .or. &
           lnd(ens1)%iamroot_compid .or. &
           ice(ens1)%iamroot_compid .or. &
           glc(ens1)%iamroot_compid .or. &
           wav(ens1)%iamroot_compid) then
        call shr_mem_getusage(msize,mrss,.true.)
        write(logunit,105) ' memory_write: model date = ',ymd,tod, &
                           ' memory = ',msize,' MB (highwater)    ',mrss,' MB (usage)', &
                           '  (pe=',iam_GLOID,' comps=',trim(complist)//')'
      end if
    end if
#endif
    if (info_debug > 1) then
      if (iamroot_CPLID) then
        call seq_infodata_GetData(infodata,nextsw_cday=nextsw_cday)
        write(logunit,*) '  nextsw_cday = ',nextsw_cday
      end if
    end if
    call t_drvstopf('CPL:TSTAMP_WRITE',cplrun=.true.)
    call t_stopf('CPL:RUN_LOOP', hashint(1))
    ! --- Write out performance data
    call t_startf('CPL:TPROF_WRITE')
    if ((tprof_alarm) .or. ((tod == 0) .and. in_first_day)) then
      if ((tod == 0) .and. in_first_day) then
        in_first_day = .false.
      end if
      call t_adj_detailf(+1)
      call t_startf("CPL:sync1_tprof")
      call mpi_barrier(mpicom_GLOID,ierr)
      call t_stopf("CPL:sync1_tprof")
      write(timing_file,'(a,i8.8,a1,i5.5)') &
            trim(tchkpt_dir)//"/model_timing"//trim(cpl_inst_tag)//"_",ymd,"_",tod
      call t_set_prefixf("CPL:")
      if (output_perf) then
        call t_prf(filename=trim(timing_file), mpicom=mpicom_GLOID, &
                   num_outpe=0, output_thispe=output_perf)
      else
        call t_prf(filename=trim(timing_file), mpicom=mpicom_GLOID, &
                   num_outpe=0)
      end if
      call t_unset_prefixf()
      call t_startf("CPL:sync2_tprof")
      call mpi_barrier(mpicom_GLOID,ierr)
      call t_stopf("CPL:sync2_tprof")
      call t_adj_detailf(-1)
    end if
    call t_stopf  ('CPL:TPROF_WRITE')
    if (barrier_alarm) then
      call t_drvstartf  ('CPL:BARRIERALARM',cplrun=.true.)
      call mpi_barrier(mpicom_GLOID,ierr)
      call t_drvstopf   ('CPL:BARRIERALARM',cplrun=.true.)
    end if
  end do   ! driver run loop
  !|----------------------------------------------------------
  !| End of driver time step loop
  !|---------------------------------------------------------

#if defined CLMSA
  ! Calling PDAF Function to set state vector before assimiliation
  call set_clm_statevec(tstartcycle, mype)
#endif

  call t_startf ('CPL:RUN_LOOP_BSTOP')
  call mpi_barrier(mpicom_GLOID,ierr)
  call t_stopf ('CPL:RUN_LOOP_BSTOP')

  Time_end = mpi_wtime()

  write(6,*)'clm5: completed timestep ',nstep

end subroutine clm_advance

subroutine clm_finalize() bind(C,name="clm_finalize")
  use iso_C_binding

  use ESMF,          only : ESMF_Initialize, ESMF_Finalize
  use cime_comp_mod, only : cime_final

  call cime_final()

  !--------------------------------------------------------------------------
  ! Clean-up
  !--------------------------------------------------------------------------
  ! call ESMF_Finalize( )


end subroutine clm_finalize

end module enkf_clm_5
