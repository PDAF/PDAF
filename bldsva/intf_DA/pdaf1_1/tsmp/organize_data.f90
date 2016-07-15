!+ External procedure for organizing I/O
!------------------------------------------------------------------------------

SUBROUTINE organize_data (yaction, nstep, n_num_now, n_num_tot, &
                          grids_dt, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   This external subroutine organizes the I/O of the model. Before doing this
!   it has to be called with yaction = 'input-init' for performing necessary 
!   initializations:
!     - read NAMELIST-groups for I/O
!     - initialize MPE_IO
!
!   The other actions are:
!     - 'start'  :   read the initial data and the first two boundary data sets
!                    and initialize LM variable table
!     - 'boundary':  read the following boundary data sets
!     - 'result':    write results
!
! Method:
!   IF-construction for checking yaction
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.34       1999/12/10 Ulrich Schaettler
!  Initial release
! 1.37       2000/03/24 Guenther Doms
!  Variables 'setup_vartab' added 
! 1.39       2000/05/03 Ulrich Schaettler
!  Included Subroutines for Namelist Input and initialization of MPE_IO.
! 1.40       2000/05/23 Ulrich Schaettler
!  Correction for using MPE_IO on sequential platforms.
! 1.41       2000/06/02 Ulrich Schaettler
!  Correction in Subroutine input_gribout.
! 2.8        2001/07/06 Ulrich Schaettler
!  Changed the organization of I/O; 
!  Introduced new variable table for meteorol. fields;
!  Introduced new NAMELIST variables yvarini, yvarbd, ymode_read, 
!                                    ymode_write, nvers
! 2.9        2001/07/16 Ulrich Schaettler
!  Corrections for reading variable lists.
! 2.11       2001/09/28 Ulrich Schaettler
!  Added checking of input lists for initial and boundary variables.
!  Corrected treatment of error variable in calc_ngrib
!  Read namelist variables for cloud ice treatment
! 2.12       2001/11/07 Ulrich Schaettler
!  Corrected pointer connection to boundary field for qi (qi_bd)
! 2.14       2002/02/15 Ulrich Schaettler
!  Modifications to allow height of snow fall limit as diagnostic quantity
!  Modifications to allow use of the SLEVE coordinate; bug corrections
! 2.17       2002/05/08 Ulrich Schaettler
!  Added checks for setting of Namelist input; moved namelist variable ncenter
!  to group IOCTL
! 2.18       2002/07/16 Ulrich Schaettler
!  Introduced variables qr, qs, qrs in the grib tables (by Almut Gassmann)
!  Introduced new namelist variables lbd_frame, npstrframe (by Lucio Torrisi)
!  Use lmulti_layer to define boundary values for old soil model correctly
! 2.19       2002/10/24 Ulrich Schaettler
!  Eliminated T_CL, W_CL from output list for lmulti_layer=.TRUE.
! 3.2        2003/02/07 Ulrich Schaettler
!  Read new Namelist parameter ilevbotnoframe
! 3.5        2003/09/02 Ulrich Schaettler
!  Include new variables for integrated cloud water, cloud ice and zenith 
!  delay in the Grib tables
!  Corrected a bug for the digital filtering
! 3.6        2003/12/11 Ulrich Schaettler
!  Changes for the multi-layer soil model;
!  Modifications for checking the IOSTAT-value when reading the NAMELISTs
! 3.7        2004/02/18 Ulrich Schaettler
!  Possibility for specifying additional GRIB tables;
!  Organization for output of synthetic satellite images as GRIB fields;
!  Adaptations in calc_ngrib for writing GRIB files every 15/30 minutes;
!  Adaptations for checking, if cloud-ice is turned on/off
! 3.8        2004/03/23 Ulrich Schaettler
!  Bug-Fix in the treatment of the linked-list for Namelist group /GRIBOUT/
!  Replaced FLOAT by REAL
! 3.10       2004/06/23 Ulrich Schaettler
!  Added check of isynsat_stat, whether synthetic satellite images are written
! 3.12       2004/09/15 Christoph Schraff
!  Adaptation to read prognostic rain / snow from analysis file according
!  to new namelist variable 'lana_qr_qs'.
! 3.13       2004/12/03 Ulrich Schaettler
!  Removed SR setup_vartab and put it to an extra module src_setup_vartab.
!  Adaptations for multi-layer soil model and new graupel scheme.
! 3.14       2005/01/25 Ulrich Schaettler
!  Introduced error message, if variables from I/O-lists cannot be found in
!  the LM Grib tables
! 3.15       2005/03/03 Ulrich Schaettler
!  Introduced new NL Parameter: nunit_of_time; Added check, whether start,
!  end and hincbound can be resolved by dt
! 3.16       2005/07/22 Thorsten Reinhardt, Helmut Frank, Ulrich Schaettler
!  Introduced new NL Parameter lana_qg;                             (T.R.)
!  Implemented test, whether output variables are really allocated; (U.S.)
!  Adapted ynaman to the use of the multi-layer soil model          (H.F.)
! 3.17       2005/12/12 Ulrich Schaettler
!  Corrected check, whether variables appearing in an I/O list are allocated
!  Input of new Namelist parameters lana_rho_snow, lan_rho_snow
! 3.18       2006/03/03 Ulrich Schaettler
!  Additional namelist variables from Climate-LM version
!  Added fields for the lake model FLake in the I/O lists (Dmitrii Mironov)
!  Added structure pp_restart for organization of Restarts
!  Additional components in structure pp_nl for output of subdomain fields
! 3.19       2006/04/25 Ulrich Schaettler
!  Some corrections for writing restart files;
!  In case of llake, read an additional initial field (T_S_LAKE)
! 3.20       2006/06/29 Helmut Frank
!  Limit now%nextstep to values less or equal now%outsteps to avoid that
!  the boundary of array now%ngrib is exceeded.
! 3.21       2006/12/04 Ulrich Schaettler
!  Adaptation of a loop-limit in input_gribout; adaptation of the restart lists
!  In case of restart must not call organize_input with lfirst=.TRUE.
!  Additional modifications for the FLake model
!  Modifications to read boundaries for qr,qs,qg
! 3.22       2007/01/24 Jochen Foerstner
!  Corrected check that restarts are only written after nudging ends
! V3_23        2007/03/30 Ulrich Schaettler
!  Need tinc_lh for restart, if ldiabf_lh is true
!  Adopted I/O to flexible time steps dt, that do not "fit" to a full hour
!  If .NOT. lbdclim, put vio3, hmo3 to the constant fields (M. Raschendorfer)
! V3_24        2007/04/26 Ulrich Schaettler
!  Initialization and distribution of NL variable ydir (for output directory)
! V3_25        2007/05/21 Ulrich Schaettler
!  Eliminated nincbound for Namelist output in YUSPECIF
! V3_26        2007/05/29 Ulrich Schaettler
!  Adaptations for writing analysis files with a flexible dt
!  Computation of IO-control variables (lastbound, nextbound) also for lgen=TRUE
! V3_27        2007/06/29 Ulrich Schaettler
!  Bug Correction for calculating vector ngrib
! V4_1         2007/12/04 Ulrich Schaettler
!  Bug correction for setting up restart files with FLake model
!  If ytunitbd == d, only full hours are allowed for hincbound
!  Check that FI cannot be written on modellevels
! V4_3         2008/02/25 Matthias Raschendorfer, Matteo Buzzi
!  Added additional fields to the restart files ('edr', topo corrections)
! V4_4         2008/07/16 Ulrich Schaettler
!  Read external parameter for sso-scheme, if lsso is set
!  Replaced lshallow by itype_conv
! V4_5         2008/09/10 Guenther Zaengl
!  Adaptations for new reference atmosphere and lateral radiative boundary condition
!  Write external parameters for sso-scheme to constant file, if lsso is set
! V4_7         2008/12/12 Ulrich Schaettler
!  FOR_E and FOR_D are added to default list for constant fields yvarc
! V4_8         2009/02/16 Ulrich Schaettler
!  Moved Namelist variable lgen to src_setup as lartif_data
!  Moved consistency check for artificial data to organize_dynamics
!  Use noutlevels for actual existing output levels
!  Added new fields for restart
!  Write sso- and forest-fields in case of analysis runs
!  Add l_ke_in_gds to partly replace ldwd_grib_use
! V4_9         2009/07/16 Ulrich Schaettler
!  Add COSMO-ART fields to initial, boundary and restart files
!  Use of 3 more GRIB tables
!  Removed unnecessary SSO fields from restarts
!  Included radtopo-fields to initial- and constant- fields list
! V4_10        2009/09/11 Davide Cesari
!  Initialize default values for a subdomain output to avoid compiler problems
!  in case of small numerical discrepancies
! V4_11        2009/11/30 Ekaterina Machulskaya, Jan-Peter Schulz, Guenther Zaengl
!  Put variables for multi-layer snow model to I/O lists (EM)
!  Modifications for the new sea-ice model with respect to T_ICE and H_ICE (JPS)
!  Optional use of reference atmosphere with constant BruntVaisala frequency 
!     for artifdata (GZ)
!  Eliminated option for cold start of FLake model (is put to INT2LM)
! V4_12        2010/05/11 Michael Baldauf, Juergen Helmert
!  Set l_dzeta_d_needed to .TRUE., if vorticity variables are written
!  Include new aerosol fields to initial and restart lists, if itype_aerosol=2 (JH)
!  Eliminated Namelist variables yvarini, yvarbd. These are set by default now (Uli)
!  Include DURSUN_M in restart files (Oli Fuhrer)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Ulrich Schaettler
!  The field t_s_lake is removed again after adaptations in the SST Analysis
! V4_17        2011/02/24 Ulrich Blahak
!  - Eliminated unused parameter lartificial
!  - read namelist for idealized cases;
!  - enabled restart runs for idealized cases
!  - In case of idealized runs, enable arbitrary output time specifications using hcomb 
!    by setting outblock%lhour = .false.
!  - introduced noutst_max (default 1000) instead of hardwired 100
!  - prepare setting of hincbound to 5 min instead of 15, but at the moment commented out.
!  - fixed calculation of ncount in calc_ngrib()
!  - Introduced namelist parameter "itype_vertint" in GRIBOUT namelist group to 
!    specify the type of vertical interpolation to p- and z-levels 
! V4_18        2011/05/26 Ulrich Schaettler
!  Moved NL parameter yform_write to the group(s) /GRIBOUT/ to be able to
!     specify the output format differently for every group. 
!  Introduced conditional compilation for Nudging and SYNSAT
!  Changes to COSMO-ART part: introduced extra variable tables for ART and POLLEN
!   (Christoph Knote)
!  Check that MSG_xx variables are only used in NetCDF output
!  - Fixed generation of BC fields for idealized restart runs. This
!    has to be done in any case, not just for non-periodic runs. (Uli Blahak)
!  - Eliminated lperi_x / lperi_y / l2dim from dependency list. (Uli Blahak)
! V4_19        2011/08/01 Ulrich Schaettler
!  Introduced conditional compilation for NetCDF and GRIBDWD
! V4_20        2011/08/31 Ulrich Schaettler
!  Replaced variable namelist (Fortran Keyword) by outblock
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

USE data_parameters,    ONLY:  ireals, iintegers, iwlength, intgribf

!------------------------------------------------------------------------------

#ifdef POLLEN
USE data_pollen,        ONLY:  lpollenini, lpollenbd
#endif

!------------------------------------------------------------------------------

USE data_parallel,      ONLY:                              &
    lasync_io,       & ! if .TRUE.: the model runs with extra PEs for
                       ! asynchronous IO
!   nprocx,          & ! number of processors in x-direction
!   nprocy,          & ! number of processors in y-direction
    nprocio,         & ! number of extra processors for doing asynchronous IO
    nproc,           & ! total number of processors: nprocx * nprocy
    num_compute,     & ! number of compute PEs
    num_io,          & ! number of IO PEs
    my_world_id,     & ! rank of this subdomain in the global communicator
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    icomm_world,     & ! communicator that belongs to igroup_world, i.e.
                       ! = MPI_COMM_WORLD
    icomm_compute,   & ! communicator for the group of compute PEs
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    imp_integers,    & ! determines the correct INTEGER type used in the
                       ! model for MPI
!   imp_byte,        & ! determines the correct BYTE type used in the model
!                      ! for MPI
    imp_character,   & ! determines the correct CHARACTER type used in the
                       ! model for MPI
    imp_logical,     & ! determines the correct LOGICAL   type used in the
                       ! model for MPI
    lcompute_pe,     & ! indicates whether this is a compute PE or not
    intbuf, realbuf, logbuf, charbuf ! Buffer for distributing the namelists

!------------------------------------------------------------------------------

USE data_io,            ONLY:                              &
    noutst,       & ! maximum number of output timesteps
    ntrip,        & ! maximum number of timing triples
    nzmxin,       & ! maximum number of input model variables
    nzmxml,       & ! maximum number of output model variables
    nzmxpl,       & ! maximum number of pressure-level variables
    nzmxzl,       & ! maximum number of height-level variables
    nzmxc,        & ! maximum number of constant variables
    nlevels,      & ! maximum number of pressure or height levels
    noutlevels,   & ! maximum actual existing number of output levels
    var,          & ! array for LM variable table
    pp_nl,        & ! structure for gribout namelist
    root,         & ! pointer to the root of gribout namelists
    max_gribtabs, & ! maximum number of GRIB tables in LM variable table
    nanfld,       & ! max. number of input fields to be checked for
    lchkini,      & ! checking the initial data
    lchkbd,       & ! checking the boundary data
    ldwd_grib_use,& ! use some DWD specific Grib settings
    l_ke_in_gds,  & ! explicit GDS entry for number of model levels
    lbdana,       & ! boundary data are analysed data
    lbdclim,      & ! boundary data in climate model     ! PIK  (D.Hauffe)
                    ! (in climate mode also some external parameters have
                    !  to be updated, which are held constant in forecast
                    !  mode; e.g. plant cover, root depth)
    lana_qi,      & ! if .TRUE., take qi-values from analysis file
                    ! else, qi is set in the model
    llb_qi,       & ! if .TRUE., take qi_bd-values from lateral boundaries file
                    ! else, qi_bd is set in the model
    lana_qr_qs,   & ! if .TRUE., take qr- and qs-values from analysis file
                    ! else, qr and qs are set in the model
    llb_qr_qs,    & ! if .TRUE., take qr_bd- and qs_bd-values from lateral
                    ! bound. file else, qr_bd and qs_bd are set in the model
    lana_qg,      & ! if .TRUE., take qg-values from analysis file
                    ! else, qg is set in the model
    llb_qg,       & ! if .TRUE., take qg_bd-values from lateral boundaries file
                    ! else, qg_bd is set in the model
    lana_rho_snow,& ! if .TRUE., take rho_snow-values from analysis file
                    ! else, it is set in the model
    lbd_frame,    & ! if .TRUE., boundary data are on a frame
    npstrframe,   & ! width (number of points) of the strip around
                    ! the b.d. frame
    ilevbotnoframe,&! bottom model level with b.d. defined on the whole grid
                    ! (model levels below ilevbotnoframe are defined on a frame)
    lanfld,       & ! contains switches for all fields to be checked for
    ynaman,       & ! name of fields to be checked for time range indicator
    nsma_stat,    & ! status for soil humidity analysis
! ub>>
    noutst_max      ! max number of output timesteps
! ub<<


USE data_io,            ONLY:                              &
    ytrans_in,    & ! directory for reading ready-files
    ytrans_out,   & ! directory for writing ready-files
    nincwait,     & ! if ready-file is not available wait nincwait seconds
                    ! until next attempt
    nmaxwait,     & ! if ready-file is not available after nmaxwait seconds,
                    ! abort the program
    ngribout,     & ! number of GRIBOUT namelist groups
    ncenter,      & ! originating center identification
    num_gribtabs, & ! number of GRIB tables used in LM variable table
    lst_gribtabs, & ! IDs of GRIB tables use
!   ydate_ini,    & ! start of the forecast
!   ydate_bd,     & ! start of the forecast for boundary data
    ydirini,      & ! directory of initial data
    ydirbd,       & ! directory of boundary data
    ytunitbd,     & ! unit of time for the boundary data
    ymode_read,   & ! mode for opening the (read) Grib files
    ymode_write,  & ! mode for opening the (write) Grib files
    yform_read,   & ! format of the (read) files
    nuchkdat,     & ! Unit number for checking the I/O data
    ntrans_out,   & ! Unit Number for writing ready-Files during output
    ar_des,       & ! structure for LM variable table
!   list_description, & ! structure for list descriptions
    list_ini,     & ! variable list for initial data
    list_bd,      & ! variable list for boundary data
    list_res_o,   & ! variable list for restart data (old step)
    list_res_n,   & ! variable list for restart data (now step)
    yvarini,      & ! list of variables for Input
    yvarbd,       & ! list of variables for Output
    nyvar_i,      & ! number of variables for Input
    nyvar_b         ! number of variables for Output

USE data_io,            ONLY:                              &
    nhour_restart,  & ! start-, stop-, inc of writing restart files (tstep)
    pp_restart,     & ! structure of type pp_nl, for restart files
    ydir_restart,   & ! directory of restart file
    ytunit_restart, & ! unit of timescale
    yncglob_institution,   & ! originating center name
    yncglob_title,         & ! title string for the output
    yncglob_source,        & ! program name and version
    yncglob_project_id,    & ! identification of the project of simulation
    yncglob_experiment_id, & ! identification of the experiment of simulation
    yncglob_contact,       & ! contact e.g. email address
    yncglob_references,    & ! URL, report etc.
    ncglob_realization       ! number of the realization of the experiment

!------------------------------------------------------------------------------

USE data_runcontrol,    ONLY:                               &
    nvers,        & ! version number of experiment for documentation
    l2tls,        & ! time integration by two timelevel RK-scheme (.TRUE.)
                    ! else with split-explicit scheme (only for l2tls=FALSE!)
    nuspecif,     & ! output unit for protocolling the task
    lphys,        & ! forecast with physical parametrizations
    lforest,      & ! if .true., run with forest (evergreen and deciduous)
    ltur,         & ! forecast with vertical diffusion
    itype_conv,   & ! type of convection parameterization
    lrad,         & ! Radiation scheme
    lemiss,       & ! surface emissivity map
    lstomata ,    & ! minimum stomata resistance map
    itype_aerosol,& ! type of aerosol map
    lmulti_layer, & ! run multi-layer soil model
    lmulti_snow , & ! run multi-layer snow model
    lseaice,      & ! forecast with sea ice model
    llake,        & ! forecast with lake model FLake
    lsso,         & ! forecast with sub-grid scale orography scheme
    lmelt,        & ! soil model with melting process
    ldiagnos,     & ! perform diagnostic calculations
    newbc,        & ! number of times that boundary update is analysis after 1h
    newbcdt,      & ! time step increment of boundary update being derived from
                    ! the (latest) analysis (rather than forecast) fields
    nincboufac,   & ! factor to 'nincbound' when new boundary update is from
                    ! analysis data
    lprog_qi,     & ! if .TRUE., running with cloud ice
    itype_gscp,   & ! type of microphys. parametrization
    itype_turb,   & ! type of turbulent diffusion parametrization
    lprog_tke,    & ! prognostic treatment of TKE (for itype_turb=5/7)
    lprogprec,    & ! if .TRUE., forecast with prognostic rain and snow (qr, qs)
    ldiniprec,    & ! diagnostic initialisation of prognostic precip (qr, qs)
    lw_freeslip,  & ! if .TRUE.:  with free slip lateral boundary condition and
                    ! if .FALSE.: specified lateral boundary values for w
    lartif_data,  & ! forecast with self-defined artificial data
    luseobs,      & ! on - off switch for using observational data
    luse_rttov,   & ! if rttov-library is used
    l_cosmo_art,  & ! if .TRUE., run the COSMO_ART
    l_pollen,     & ! if .TRUE., run the pollen
    lradtopo,     & ! topographical radiation correction if lradtopo=.TRUE.
    l_dzeta_d_needed  ! metric coeff. dzeta_dlam, dzeta_dphi are needed

USE data_runcontrol,    ONLY:                               &
    idbg_level,   & ! to control the verbosity of debug output
    ldebug_io,    & ! if .TRUE., debug output for I/O
    lprintdeb_all,& ! whether all or only one task prints debug output
    ldiabf_lh,    & ! include diabatic forcing due to latent heat in RK-scheme
    isynsat_stat, & ! status of synthetic satellite images
    nlgw,         & ! number of prognostic soil water levels
    nlgw_bd,      & ! number of prognostic soil water levels in boundary data
    nlgw_ini,     & ! number of prognostic soil water levels in initial data
    nstart,       & ! first time step of the forecast
    nstop,        & ! last time step of the forecast
    nfinalstop,   & ! last time step of the total forecast
                    ! (necessary, if simulation is splitted into periods)
    nbd1,         & ! indices for permutation of the
    nbd2,         & ! two boundary time levels
    hstart,       & ! start of the forecast in full hours
    hstop,        & ! end of the forecast in hours
    hlastbound,   & ! last hour with boundary update
    hincbound,    & ! hour increment for reading boundary update
    hnextbound,   & ! next hour with boundary update
    nlastbound,   & ! time step of the last boundary update
    nnextbound,   & ! time step of the next boundary update
    nincbound,    & ! time step increment of boundary update
    yakdat1         ! actual date (ydate_ini+ntstep/dt) in the form
                    ! ddmmyyyyhh (day, month, year, hour)

!------------------------------------------------------------------------------

USE data_modelconfig,   ONLY:                              &

! 1. vertical coordinate parameters and related variables
! -------------------------------------------------------
    vcoord,       & ! vertical coordinate of LM                    half level
    sigmr,        & ! sigma-coordinate referring to PMSL           half level
    hhlr,         & ! height-coordinate referring to MSL (z=0)     half level
    vcflat,       & ! vertical coordinate where the terrain following system
    ivctype,      & ! type of vertical coordinate
    svc1,         & ! vertical decay rate for large-scale topo part
    svc2,         & ! vertical decay rate for small-scale topo part
    kflat,        & ! level index where (half) levels become flat

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ie_tot,       & ! number of grid points in zonal direction
    je_tot,       & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction

! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------
    dlon,         & ! grid point distance in zonal direction (in degrees)
    dlat,         & ! grid point distance in meridional direction (in degrees)
    startlon_tot, & ! transformed longitude of the lower left grid point
                    ! of the total domain (in degrees, E>0)
    startlat_tot, & ! transformed latitude of the lower left grid point
                    ! of the total domain (in degrees, N>0)

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------
    dt,           & ! long time-step

! 6. variables for the reference atmosphere
! -----------------------------------------
    irefatm,      & ! 1: old reference atm. based on dT/dlnp = const
                    ! 2: new reference atm. with asymptotically isothermal stratosphere
                    ! 3: constant BV-frequency (idealized runs only)
    p0sl,         & ! reference pressure at sea level
    t0sl,         & ! reference temperature at sea level
    dt0lp,        & ! d (t0) / d (ln p0)
    delta_t,      & ! temp. difference between sea level and stratosphere (for irefatm = 2)
    h_scal,       & ! scale height for irefatm = 2
    bvref,        & ! BV-freq. for irefatm=3

! 7. Layer index corresponding to a specified pressure
! ----------------------------------------------------
    klv850,       & ! k index of the LM-mainlevel, on 850 HPa
    klv800,       & ! k index of the LM-mainlevel, on 800 HPa
    klv500,       & ! k index of the LM-mainlevel, on 500 HPa
    klv400,       & ! k index of the LM-mainlevel, on 400 HPa
    klv300          ! k index of the LM-mainlevel, on 300 HPa

!------------------------------------------------------------------------------

#ifdef NUDGING
USE data_nudge_all,     ONLY:  nudgend   ! for checking with restart files
#endif

!------------------------------------------------------------------------------

USE environment,        ONLY:  model_abort, get_free_unit
USE io_utilities,       ONLY:  make_fn
USE parallel_utilities, ONLY:  distribute_values

USE mpe_io,             ONLY:  mpe_io_init, mpe_io_reconfig, mpe_readnldb

USE src_artifdata,      ONLY:  input_artifctl
USE src_input,          ONLY:  organize_input
USE src_output,         ONLY:  init_output, organize_output, makegds
USE src_artifdata,      ONLY:  gen_ini_data, gen_bound_data
USE src_setup_vartab,   ONLY:  setup_vartab

#ifdef COSMOART
USE art_setup_vartab,   ONLY:  setup_vartab_art
#endif

#ifdef POLLEN
USE pol_setup_vartab,   ONLY:  setup_vartab_pol
#endif

!==============================================================================

use data_parallel, only: cosmo_input_suffix

IMPLICIT NONE

!==============================================================================
!
! Parameter list:
CHARACTER (LEN= *),       INTENT(IN)            ::                      &
  yaction      ! action to be performed

INTEGER (KIND=iintegers), INTENT(IN)            ::                      &
  nstep,     & ! actual time-step
  n_num_now, & ! current nest number
  n_num_tot    ! total number of nests

REAL (KIND=ireals), INTENT(IN)                  ::                      &
  grids_dt(n_num_tot) ! array with the time steps on all grids

INTEGER (KIND=iintegers), INTENT(OUT)           ::                      &
  ierror       ! error status

CHARACTER (LEN=  *),      INTENT(OUT)           ::                      &
  yerrmsg      ! error message

!==============================================================================

! Local variables: 
INTEGER (KIND=iintegers)   ::                                           &
  nuin,            & ! Namelist INPUT file
  izerrstat,       & ! error status variable
  izerrstat1,      & ! error status variable
  izmemstat,       & ! error status for memory allocation
  izreadstat,      & ! error status for reading NAMELIST groups
  izmplcode,       & ! error status for MPE_IO or MPI-routines
  iztest_compute,  & ! test communicator for MPE_IO
  ibits,           & ! number of bits in GRIB integer variable
  n, n1, iz1, iz2, iz3, istat, ngout, izdebug

INTEGER (KIND=iintegers), SAVE   ::                                     &
  ibdf,            & ! interval between boundary fields 
                     ! (has to be saved for next call)
  nsetupvar = 0      ! number of nest vartab is actually set for

REAL    (KIND=ireals)      ::                                           &
  tact,tlbc,tivl,  & ! temporal quantities for control output
  fc_hour,         & !
  zendlon,         & ! intermediate storage
  zendlat            !

LOGICAL                   ::                          &
  lwritefc,        & ! check whether Forecasts are written in an output step
  lwriteana,       & ! check whether Analysis are written in an output step
  lwriteready        ! check whether Ready files shall be written

CHARACTER (LEN=260)        ::       &
  yzname             ! full path- and file-name of ready files

CHARACTER (LEN=10)           ::             &
  yzloclist(nzmxml)  ! local copy of the list

CHARACTER (LEN= 3)         ::       &
  yzhead,          & ! Header for creating ready-file name
  c_nnum             ! String for specifying the file name for nest output

CHARACTER (LEN= 14)         ::       &
  yinput             ! Namelist INPUT file

CHARACTER (LEN=90)         ::       &
  yzerrmsg           ! error message variable

CHARACTER (LEN=15)         ::       &
  yzroutine          ! name of this routine for error messages

! Pointer for the linked list
TYPE(pp_nl) , POINTER      ::  now, next

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
!- Begin Subroutine organize_data
!------------------------------------------------------------------------------

ierror     = 0
istat      = 0
yerrmsg    = '   '
izerrstat  = 0
izerrstat1 = 0
izmemstat  = 0
izreadstat = 0
izmplcode  = 0
yzroutine  = 'organize_data'

c_nnum     = '   '
IF ( n_num_now > 1 ) THEN
  IF ( n_num_now < 10 ) THEN
    WRITE ( c_nnum, '( A1, I1, A1 )' ) '.', n_num_now, ' '
  ELSE
    WRITE ( c_nnum, '( A1, I2.2 )'   ) '.', n_num_now
  ENDIF
ENDIF

! Initialize, whether additional debug output shall be done
IF (ldebug_io) THEN
  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF
ELSE
  izdebug = 0
ENDIF

!------------------------------------------------------------------------------
! Section 1: Initialization of I/O
!------------------------------------------------------------------------------

IF (yaction == 'input-init') THEN

  !----------------------------------------------------------------------------
  ! Section 1.1: Read in NAMELIST groups for I/O
  !----------------------------------------------------------------------------

  ! Open NAMELIST-INPUT file for IO
  IF (my_world_id == 0) THEN
    IF (idbg_level > 0) THEN
      PRINT *,'    INPUT OF THE NAMELISTS FOR GRIB-IO'
    ENDIF
    !yinput   = 'INPUT_IO'
write(yinput,'(a,i5.5)') 'INPUT_IO_',cosmo_input_suffix

    nuin     =  1

    OPEN(nuin   , FILE=yinput  , FORM=  'FORMATTED', STATUS='UNKNOWN',  &
         IOSTAT=izerrstat)
    IF(izerrstat /= 0) THEN
      yerrmsg = ' ERROR    *** Error while opening file INPUT_IO *** '
      ierror  = 6001
      RETURN
    ENDIF
  ENDIF

  CALL input_ioctl  (nuspecif, nuin, izerrstat )

  IF (izerrstat < 0) THEN
    yerrmsg = ' ERROR    *** while reading NAMELIST Group /IOCTL/ ***'
    ierror  = 6002
    RETURN
  ENDIF

  ! In mpe_io, izerrstat is set to 0 again, therefore use izerrstat1
  CALL input_dbase  (nuspecif, nuin, izerrstat1)
  izerrstat = izerrstat + izerrstat1

  IF (izerrstat1 < 0) THEN
    yerrmsg = ' ERROR    *** while reading NAMELIST Group /DATABASE/ ***'
    ierror  = 6003
    RETURN
  ENDIF

  CALL input_gribin (nuspecif, nuin, izerrstat)

  IF (izerrstat < 0) THEN
    yerrmsg = ' ERROR    *** while reading NAMELIST Group /GRIBIN/ ***'
    ierror  = 6004
    RETURN
  ENDIF

  ! There can be several groups for Grib output: A linked list is created
  ! with an entry for every group. The number of groups is specified by
  ! ngribout in the NAMELIST group IOCTL.

  ALLOCATE(root, STAT = izmemstat)      ! root is declared in data_iocontrol
  IF (izmemstat /= 0) THEN
    yerrmsg = ' ERROR allocating root pointer'
    ierror   = 6010
    RETURN
  ENDIF

  now => root
  noutst = 0

  ! Initialize noutlevels with maximum of number of modellevels and 
  ! "satellite" levels (which are 32: see organize_satellites, 
  ! variables synme7, synmsg)
  noutlevels = MAX (ke+1, 32)

  DO ngout = 1, ngribout

    CALL input_gribout(now, nuspecif, nuin, izreadstat, izerrstat)

    IF (izerrstat < 0) THEN
      yerrmsg = ' ERROR    *** while reading namelist /GRIBOUT/:    ***'
      WRITE (yerrmsg(49:50),'(I2)') ngout
      ierror   = 6005
      RETURN
    ELSEIF (izerrstat > 0) THEN
      yzerrmsg = ' ERROR reading gribout namelist in loop'
      ierror   = 6006
      RETURN
    ENDIF

    ! Adapt noutlevels
    noutlevels = MAX (noutlevels, now%kezin, now%kepin)

    ! Loop exit or next namelist group
    IF (ngout == ngribout) THEN
      NULLIFY (now%next)
    ELSE
      ALLOCATE(next, STAT=izmemstat)
      IF (izmemstat /= 0) THEN
        yerrmsg = ' ERROR allocating new pointer'
        ierror   = 6011
        RETURN
      ENDIF
      now%next => next
      now => next
    ENDIF
  ENDDO

  IF (my_world_id == 0) THEN
    ! Close file for input of the NAMELISTS
    CLOSE (nuin    , STATUS='KEEP')

    IF (izerrstat /= 0) THEN
      ierror    = 6002
      yerrmsg  = ' ERROR *** Wrong values occured in NAMELIST INPUT_IO ***'
      RETURN
    ENDIF
  ENDIF

! UB>>
! Read in the namelist for idealized runs. This is done in any case to
! enable the possibility to mix idealized features (e.g., warm bubbles)
! into real case runs. If namelist file 'INPUT_IDEAL' does not exist,
! nothing is read and all namelist parameters take on their default
! values, which are neutral to operational runs.
! HOWEVER, WE DISABLE IT IN CASE OF LARTIF_DATA=.FALSE. FOR NOW.
! IF THE USER WOULD LIKE TO MIX ARTIFICIAL ELEMENTS INTO REAL CASE
! SIMULATIONS, COMMENT OUT THE lartif_data-IF-CLAUSE!
! Also, do the same in **organize_physics.f90** and **src_output.f90**!
  IF (lartif_data) THEN
    CALL input_artifctl (nuspecif, izerrstat)
    IF (izerrstat /= 0) THEN
      ierror    = 6007
      RETURN
    ENDIF
  END IF
! UB<<

  !----------------------------------------------------------------------------
  ! Section 1.2: Initialize MPE_IO (has to be done by all PEs
  !----------------------------------------------------------------------------

  ! In LM the communicator for the compute PEs has already been build in
  ! init_procgrid, so here only a test is done whether mpe_io_init builds
  ! the same communicator
  IF (lasync_io .OR. (num_compute > 1) ) THEN
    CALL mpe_io_init (icomm_world, num_compute, num_io, iztest_compute,   &
                      izmplcode)
  ENDIF

  IF (izmplcode /= 0) THEN
    ierror  = 6013
    yerrmsg = '*** ERROR: initializing MPE_IO ***'
    RETURN
  ENDIF

  IF (lcompute_pe .AND. (lasync_io .OR. (num_compute > 1) )) THEN
    CALL mpe_io_reconfig (icomm_compute, izmplcode)
    IF (izmplcode /= 0) THEN
      ierror  = 6014
      yerrmsg = '*** ERROR: reconfiguring communicator for MPE_IO ***'
      RETURN
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 1.3: Initialize organizational variables for GRIB
  !----------------------------------------------------------------------------

  ! initialize variable iwlength
  ibits       = BIT_SIZE (1_intgribf)
  iwlength    = INT (ibits, iintegers) / 8

  !----------------------------------------------------------------------------
  ! Section 1.4: Initialize variables for ASCII-file handling
  !----------------------------------------------------------------------------

  CALL get_free_unit (ntrans_out)
  CALL get_free_unit (nuchkdat)

  !----------------------------------------------------------------------------
  ! Section 1.5: Initialize variables for restart files
  !----------------------------------------------------------------------------

  ALLOCATE(pp_restart, STAT = izmemstat)
  IF (izmemstat /= 0) THEN
    yerrmsg = ' ERROR allocating pp_restart pointer'
    ierror   = 6014
    RETURN
  ENDIF

  zendlon = startlon_tot  + REAL(ie_tot-1, ireals) * dlon
  IF (zendlon > 180.0_ireals) THEN
    zendlon = zendlon - 360.0_ireals
  ENDIF

  zendlat = startlat_tot  + REAL(je_tot-1, ireals) * dlat

  ! Initialize necessary components of pp_restart
  pp_restart%yvarml(:) = '          '    ! will be used for 'o'-list
  pp_restart%yvarpl(:) = '          '    ! will be used for 'n'-list (leapfrog)
  pp_restart%yvarzl(:)        = ''       ! not used
  pp_restart%yvarsl(:)        = ''       ! not used
  pp_restart%yvarc (:)        = ''       ! not used
  pp_restart%ilist_ml(:,:)    =  0       ! will be set later
  pp_restart%ilist_pl(:,:)    =  0       ! will be set later
  pp_restart%ilist_zl(:,:)    = -1       ! not used
  pp_restart%ilist_sl(:,:)    = -1       ! not used
  pp_restart%ilist_c (:,:)    = -1       ! not used
  pp_restart%n_num            =  1       ! nr. of grid in nesting
  pp_restart%nyvar_m          =  0       ! will be set later
  pp_restart%nyvar_p          =  0       ! will be set later
  pp_restart%nyvar_z          =  0       ! not used
  pp_restart%nyvar_s          =  0       ! not used
  pp_restart%nyvar_c          =  0       ! not used

  AllOCATE (pp_restart%ngrib(1), STAT=istat)
  pp_restart%ngrib(:)         =  0       ! not used
  pp_restart%outsteps         =  0       ! not used
  pp_restart%nexthour         =  nhour_restart(1)
  pp_restart%nextstep         =  NINT (3600.0_ireals * nhour_restart(1) / dt)-1
                                            ! first step for restart files
  pp_restart%lhour            = .TRUE.
  pp_restart%nprocess_ini_out = 0           ! not used
  pp_restart%nprocess_bd_out  = 0           ! not used
  pp_restart%nunit_of_time    = 1           ! for hourly output
  pp_restart%slon             = startlon_tot! not used
  pp_restart%slat             = startlat_tot! not used
  pp_restart%elon             = zendlon     ! not used
  pp_restart%elat             = zendlat     ! not used
  pp_restart%i_out_start      = 1           ! not used
  pp_restart%j_out_start      = 1           ! not used
  pp_restart%i_out_end        = ie_tot      ! not used
  pp_restart%j_out_end        = je_tot      ! not used
  pp_restart%ie_out_tot       = ie_tot      ! not used
  pp_restart%je_out_tot       = je_tot      ! not used
  pp_restart%yform_write      = 'bina'      ! not 'grb1','ncdf'
  pp_restart%ysystem          = 'FILE'      ! default
  pp_restart%ydir             = ydir_restart! has to be set by namelist
  pp_restart%ysuffix          = ' '         ! 'o' or 'n'
  pp_restart%ytunit           = ytunit_restart ! has to be set by namelist
  pp_restart%ydomain          = 'f'         ! full domain
  pp_restart%nrbit            = 0           ! not used
  pp_restart%ydbid            = ' '         ! not used
  pp_restart%ydbtype          = ' '         ! not used
  pp_restart%ydbpw            = ' '         ! not used
  pp_restart%plev(:)          = -1.0        ! not used
  pp_restart%zlev(:)          = -1.0        ! not used
  pp_restart%kepin            =    0        ! not used
  pp_restart%kezin            =    0        ! not used
  pp_restart%lcheck           = .TRUE.      ! not used
  pp_restart%lwrite_const     = .FALSE.     ! not used
  pp_restart%luvmasspoint     = .FALSE.     ! not used
  pp_restart%lanalysis        = .FALSE.     ! not used
  pp_restart%l_p_filter       = .FALSE.     ! not used
  pp_restart%l_z_filter       = .FALSE.     ! not used
  pp_restart%l_fi_ps_smooth   = .FALSE.     ! not used

  NULLIFY (pp_restart%next)

!------------------------------------------------------------------------------
! Section 2: Input of initial data and the first two boundary data sets
!------------------------------------------------------------------------------

ELSEIF (yaction == 'start') THEN

  !----------------------------------------------------------------------------
  ! Section 2.1: Initialize LM-variable table
  !----------------------------------------------------------------------------

  CALL setup_vartab

#ifdef COSMOART
  CALL setup_vartab_art
#endif

#ifdef POLLEN
  CALL setup_vartab_pol
#endif

  nsetupvar = n_num_now

  ! Allocate and define the description lists for input, output
 
  ALLOCATE(list_ini    (nyvar_i), STAT=istat)
  ALLOCATE(list_bd     (nyvar_b), STAT=istat)

  IF(istat /= 0) THEN
    CALL model_abort(my_cart_id, 2029,                                    &
                 'Error allocating description lists', 'organize_data')
  ENDIF

  ! list of initial variables
  loop_list_i: DO n = 1, nyvar_i
    list_ini(n)%name = yvarini(n)
    DO iz3 = 1, num_gribtabs
      DO iz2 = 0, 255
        DO iz1 = 1,4
          IF (var(iz1,iz2,iz3)%name(1:LEN_TRIM(var(iz1,iz2,iz3)%name)) ==  &
                   list_ini(n)%name(1:LEN_TRIM(     list_ini(n)%name)) ) THEN
            ! Set location in variable table
            list_ini(n)%iloc1 = iz1
            list_ini(n)%iloc2 = iz2
            list_ini(n)%iloc3 = iz3

            ! Set dimension of this variable (1, ke or ke+1)
            SELECT CASE (var(iz1,iz2,iz3)%rank)
            CASE (4)
              list_ini(n)%idimvert = UBOUND (var(iz1,iz2,iz3)%p4,3)
            CASE (3)
              SELECT CASE (var(iz1,iz2,iz3)%name)
              CASE ('PS        ','T_G       ','QV_S      ','W_SNOW    ',  &
                    'T_S       ','T_M       ','W_G1      ','W_G2      ',  &
                    'W_G3      ','W_I       ','T_SNOW    ','RHO_SNOW  ',  &
                    'H_SNOW    ','H_ICE     ','T_ICE     ','T_MNW_LK  ',  &
                    'T_WML_LK  ','T_BOT_LK  ','T_B1_LK   ','C_T_LK    ',  &
                    'H_ML_LK   ','H_B1_LK   ')
                ! these are 2D variables with rank=3 because of time dependency
                list_ini(n)%idimvert = 1
              CASE DEFAULT
                ! these are real 3D variables
                list_ini(n)%idimvert = UBOUND (var(iz1,iz2,iz3)%p3,3)
              END SELECT
            CASE (2)
              list_ini(n)%idimvert = 1
            END SELECT
            CYCLE loop_list_i
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO loop_list_i

  ! list of boundary variables
  loop_list_b: DO n = 1, nyvar_b
    list_bd(n)%name = yvarbd(n)
    DO iz3 = 1, num_gribtabs
      DO iz2 = 0, 255
        DO iz1 = 1,4
          IF (var(iz1,iz2,iz3)%name(1:LEN_TRIM(var(iz1,iz2,iz3)%name)) ==  &
                   list_bd (n)%name(1:LEN_TRIM(     list_bd (n)%name)) ) THEN
            ! Set location in variable table
            list_bd (n)%iloc1 = iz1
            list_bd (n)%iloc2 = iz2
            list_bd (n)%iloc3 = iz3

            ! Set dimension of this variable (1, ke or ke+1)
            SELECT CASE (var(iz1,iz2,iz3)%rank)
            CASE (4)
              list_bd (n)%idimvert = UBOUND (var(iz1,iz2,iz3)%p4,3)
            CASE (3)
              SELECT CASE (var(iz1,iz2,iz3)%name)
              CASE ('PS        ','T_G       ','QV_S      ','W_SNOW    ',  &
                    'T_S       ','T_M       ','W_G1      ','W_G2      ',  &
                    'W_G3      ','W_I       ','T_SNOW    ','RHO_SNOW  ',  &
                    'H_SNOW    ','H_ICE     ','T_ICE     ','T_MNW_LK  ',  &
                    'T_WML_LK  ','T_BOT_LK  ','T_B1_LK   ','C_T_LK    ',  &
                    'H_ML_LK   ','H_B1_LK   ')
                ! these are 2D variables with rank=3 because of time dependency
                list_bd(n)%idimvert = 1
              CASE DEFAULT
                ! these are real 3D variables
                list_bd(n)%idimvert = UBOUND (var(iz1,iz2,iz3)%p3,3)
              END SELECT
            CASE (2)
              list_bd (n)%idimvert = 1
            END SELECT

            CYCLE loop_list_b
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO loop_list_b

  ! organization indices for timestepping
  nbd1 = 1
  nbd2 = 2

  ! list of restart variables: "old" timestep
  ! constant fields:
  pp_restart%yvarml( 1)  = 'HSURF     '
  pp_restart%yvarml( 2)  = 'FR_LAND   '
  pp_restart%yvarml( 3)  = 'Z0        '
  pp_restart%yvarml( 4)  = 'SOILTYP   '
  pp_restart%yvarml( 5)  = 'PLCOV     '
  pp_restart%yvarml( 6)  = 'LAI       '
  pp_restart%yvarml( 7)  = 'ROOTDP    '
  pp_restart%yvarml( 8)  = 'VIO3      '
  pp_restart%yvarml( 9)  = 'HMO3      '

  ! prognostic variables and dynamics:
  pp_restart%yvarml(10)  = 'U         '
  pp_restart%yvarml(11)  = 'V         '
  pp_restart%yvarml(12)  = 'W         '
  pp_restart%yvarml(13)  = 'T         '
  pp_restart%yvarml(14)  = 'QV        '
  pp_restart%yvarml(15)  = 'QC        '
  pp_restart%yvarml(16)  = 'PP        '
  pp_restart%yvarml(17)  = 'T_SNOW    '
  pp_restart%yvarml(18)  = 'W_I       '
  pp_restart%yvarml(19)  = 'QV_S      '
  pp_restart%yvarml(20)  = 'W_SNOW    '
  pp_restart%yvarml(21)  = 'T_S       '

  ! fields from the parameterizations
  ! microphysics:  nothing special

  ! radiation:
  pp_restart%yvarml(22)  = 'SOHR_RAD  '
  pp_restart%yvarml(23)  = 'THHR_RAD  '
  pp_restart%yvarml(24)  = 'ALB_RAD   '
  pp_restart%yvarml(25)  = 'SOBS_RAD  '
  pp_restart%yvarml(26)  = 'THBS_RAD  '
  pp_restart%yvarml(27)  = 'PABS_RAD  '
  pp_restart%yvarml(28)  = 'SOBT_RAD  '
  pp_restart%yvarml(29)  = 'THBT_RAD  '
  pp_restart%yvarml(30)  = 'CLCH      '
  pp_restart%yvarml(31)  = 'CLCM      '
  pp_restart%yvarml(32)  = 'CLCL      '
  pp_restart%yvarml(33)  = 'CLCT      '

  ! turbulence:
  pp_restart%yvarml(34)  = 'TKVM      '
  pp_restart%yvarml(35)  = 'TKVH      '
  pp_restart%yvarml(36)  = 'TCM       '
  pp_restart%yvarml(37)  = 'TCH       '

  ! convection:
  pp_restart%yvarml(38)  = 'CLC_CON   '
  pp_restart%yvarml(39)  = 'CLW_CON   '
  pp_restart%yvarml(40)  = 'PRR_CON   '
  pp_restart%yvarml(41)  = 'PRS_CON   '
  pp_restart%yvarml(42)  = 'TOP_CON   '
  pp_restart%yvarml(43)  = 'BAS_CON   '
  pp_restart%yvarml(44)  = 'DU_CON    '
  pp_restart%yvarml(45)  = 'DV_CON    '
  pp_restart%yvarml(46)  = 'DT_CON    '
  pp_restart%yvarml(47)  = 'DQV_CON   '
  pp_restart%yvarml(48)  = 'DQC_CON   '
  pp_restart%yvarml(49)  = 'DQI_CON   '
  pp_restart%yvarml(50)  = 'MFLX_CON  '
  pp_restart%yvarml(51)  = 'CAPE_CON  '
  pp_restart%yvarml(52)  = 'QCVG_CON  '
  pp_restart%yvarml(53)  = 'TKE_CON   '

  ! soil model:

  ! and all the rest:
  pp_restart%yvarml(54)  = 'TMAX_2M   '
  pp_restart%yvarml(55)  = 'VMAX_10M  '
  pp_restart%yvarml(56)  = 'TMIN_2M   '
  pp_restart%yvarml(57)  = 'RUNOFF_S  '
  pp_restart%yvarml(58)  = 'RUNOFF_G  '
  pp_restart%yvarml(59)  = 'ASOB_S    '
  pp_restart%yvarml(60)  = 'ATHB_S    '
  pp_restart%yvarml(61)  = 'ASOB_T    '
  pp_restart%yvarml(62)  = 'ATHB_T    '
  pp_restart%yvarml(63)  = 'ALHFL_S   '
  pp_restart%yvarml(64)  = 'ASHFL_S   '
  pp_restart%yvarml(65)  = 'AUMFL_S   '
  pp_restart%yvarml(66)  = 'AVMFL_S   '
  pp_restart%yvarml(67)  = 'APAB_S    '
  pp_restart%yvarml(68)  = 'DQVDT     '
  pp_restart%yvarml(69)  = 'QVSFLX    '
  pp_restart%yvarml(70)  = 'QRS       '
  pp_restart%yvarml(71)  = 'SNOW_GSP  '
  pp_restart%yvarml(72)  = 'RAIN_GSP  '
  pp_restart%yvarml(73)  = 'RAIN_CON  '
  pp_restart%yvarml(74)  = 'SNOW_CON  '
  pp_restart%yvarml(75)  = 'PRR_GSP   '
  pp_restart%yvarml(76)  = 'PRS_GSP   '
  pp_restart%yvarml(77)  = 'T_G       '
  pp_restart%yvarml(78)  = 'PS        '
  pp_restart%yvarml(79)  = 'VGUST_CON '
  pp_restart%yvarml(80)  = 'VGUST_DYN '
  pp_restart%yvarml(81)  = 'SOD_T     '
  pp_restart%yvarml(82)  = 'ASOD_T    '

  pp_restart%nyvar_m     = 82

  IF (ldiabf_lh) THEN
    pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'TINC_LH   '
    pp_restart%nyvar_m     = pp_restart%nyvar_m + 1
  ENDIF

  IF ( lphys .AND. ltur ) THEN
    SELECT CASE (itype_turb)
    CASE (3)
      pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'TKE       '
      pp_restart%yvarml(pp_restart%nyvar_m + 2)  = 'TKETENS   '
      pp_restart%yvarml(pp_restart%nyvar_m + 3)  = 'RCLD      '
      pp_restart%yvarml(pp_restart%nyvar_m + 4)  = 'EDR       '
      pp_restart%nyvar_m     = pp_restart%nyvar_m + 4
    CASE (5:8)
      IF (lprog_tke) THEN
        pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'TKE       '
        pp_restart%yvarml(pp_restart%nyvar_m + 2)  = 'TKETENS   '
        pp_restart%nyvar_m     = pp_restart%nyvar_m + 2
      ELSE
        pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'TKE       '
        pp_restart%nyvar_m     = pp_restart%nyvar_m + 1
      ENDIF
    CASE DEFAULT
      pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'TKE       '
      pp_restart%yvarml(pp_restart%nyvar_m + 2)  = 'EDR       '
      pp_restart%nyvar_m     = pp_restart%nyvar_m + 2
    END SELECT
  ENDIF

  ! variables for the soil model
  IF (lmulti_layer) THEN
    pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'T_SO      '
    pp_restart%yvarml(pp_restart%nyvar_m + 2)  = 'W_SO      '
    pp_restart%yvarml(pp_restart%nyvar_m + 3)  = 'W_SO_ICE  '
    pp_restart%yvarml(pp_restart%nyvar_m + 4)  = 'FRESHSNW  '
    pp_restart%yvarml(pp_restart%nyvar_m + 5)  = 'RHO_SNOW  '
    pp_restart%nyvar_m     = pp_restart%nyvar_m + 5
  ELSE
    pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'T_M       '
    pp_restart%yvarml(pp_restart%nyvar_m + 2)  = 'T_CL      '
    pp_restart%yvarml(pp_restart%nyvar_m + 3)  = 'W_G1      '
    pp_restart%yvarml(pp_restart%nyvar_m + 4)  = 'W_G2      '
    pp_restart%yvarml(pp_restart%nyvar_m + 5)  = 'W_CL      '
    pp_restart%nyvar_m     = pp_restart%nyvar_m + 5
  ENDIF

  ! variables for topographic radiation correction
  IF (lradtopo) THEN
    pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'SKYVIEW   '
    pp_restart%yvarml(pp_restart%nyvar_m + 2)  = 'SLO_ASP   '
    pp_restart%yvarml(pp_restart%nyvar_m + 3)  = 'SLO_ANG   '
    pp_restart%yvarml(pp_restart%nyvar_m + 4)  = 'HORIZON   '
    pp_restart%nyvar_m     = pp_restart%nyvar_m + 4
  ENDIF

  ! sunshine duration and some more fields for radiation
  pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'DURSUN    '
  pp_restart%yvarml(pp_restart%nyvar_m + 2)  = 'DURSUN_M  '
  pp_restart%nyvar_m     = pp_restart%nyvar_m + 2

  ! variables for the sea ice and/or the FLake model
  IF (lseaice .OR. llake) THEN
    pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'T_ICE     '
    pp_restart%yvarml(pp_restart%nyvar_m + 2)  = 'H_ICE     '
    pp_restart%nyvar_m     = pp_restart%nyvar_m + 2
  ENDIF

  ! variables for the FLake model
  IF (llake) THEN
    pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'FR_LAKE   '
    pp_restart%yvarml(pp_restart%nyvar_m + 2)  = 'DEPTH_LK  '
    pp_restart%yvarml(pp_restart%nyvar_m + 3)  = 'T_MNW_LK  '
    pp_restart%yvarml(pp_restart%nyvar_m + 4)  = 'T_WML_LK  '
    pp_restart%yvarml(pp_restart%nyvar_m + 5)  = 'T_BOT_LK  '
    pp_restart%yvarml(pp_restart%nyvar_m + 6)  = 'C_T_LK    '
    pp_restart%yvarml(pp_restart%nyvar_m + 7)  = 'H_ML_LK   '
    pp_restart%nyvar_m     = pp_restart%nyvar_m + 7

    !_cdm The following FLake external parameters are assigned their default 
    !     values that are kept constant in space and time.
    !     These external-parameter fields are handled internally by LM
    !     and are not included into the LM IO lists
    !     (this should be done in the future).
    !       'FETCH_LK  '
    !       'DP_BS_LK  '
    !       'T_BS_LK   '
    !       'GAMSO_LK  '
    !     In the present configuration,
    !     the bottom-sediment module of FLake is switched off.
    !     The respective fields are handled internally by LM
    !     and are not included into the LM IO lists.
    !       'T_B1_LK   '
    !       'H_B1_LK   '
  ENDIF

  ! sso scheme:
  IF (lsso) THEN
    pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'SSO_STDH  '
    pp_restart%yvarml(pp_restart%nyvar_m + 2)  = 'SSO_GAMMA '
    pp_restart%yvarml(pp_restart%nyvar_m + 3)  = 'SSO_THETA '
    pp_restart%yvarml(pp_restart%nyvar_m + 4)  = 'SSO_SIGMA '
    pp_restart%yvarml(pp_restart%nyvar_m + 5)  = 'AUSTRSSO  '
    pp_restart%yvarml(pp_restart%nyvar_m + 6)  = 'AVSTRSSO  '
    pp_restart%yvarml(pp_restart%nyvar_m + 7)  = 'AVDISSSO  '
    pp_restart%nyvar_m     = pp_restart%nyvar_m + 7
  ENDIF

  IF (lemiss) THEN
    pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'EMIS_RAD  '
    pp_restart%nyvar_m     = pp_restart%nyvar_m + 1
  END IF

  IF (lstomata) THEN
    pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'RSMIN     '
    pp_restart%nyvar_m     = pp_restart%nyvar_m + 1
  END IF


  IF ( itype_aerosol == 2 ) THEN
    pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'AER_SO4   '
    pp_restart%yvarml(pp_restart%nyvar_m + 2)  = 'AER_DUST  '
    pp_restart%yvarml(pp_restart%nyvar_m + 3)  = 'AER_ORG   '
    pp_restart%yvarml(pp_restart%nyvar_m + 4)  = 'AER_BC    '
    pp_restart%yvarml(pp_restart%nyvar_m + 5)  = 'AER_SS    '
    pp_restart%nyvar_m     = pp_restart%nyvar_m + 5
  ENDIF

  IF (lforest) THEN
    pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'FOR_E     '
    pp_restart%yvarml(pp_restart%nyvar_m + 2)  = 'FOR_D     '
    pp_restart%nyvar_m     = pp_restart%nyvar_m + 2
  ENDIF

  ! additional variables depending on chosen configuration
  IF (lprog_qi) THEN
    pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'QI        '
    pp_restart%nyvar_m     = pp_restart%nyvar_m + 1
  ENDIF
  IF (lprogprec) THEN
    pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'QR        '
    pp_restart%nyvar_m     = pp_restart%nyvar_m + 1

    IF (itype_gscp > 1) THEN
      pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'QS        '
      pp_restart%nyvar_m     = pp_restart%nyvar_m + 1
    ENDIF
    IF (itype_gscp == 4) THEN
      pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'QG        '
      pp_restart%yvarml(pp_restart%nyvar_m + 2)  = 'PRG_GSP   '
      pp_restart%yvarml(pp_restart%nyvar_m + 3)  = 'GRAU_GSP  '
      pp_restart%nyvar_m     = pp_restart%nyvar_m + 3
    ENDIF
  ENDIF

  ! HJP: accumulated (ntri=4) and averaged (ntri=3) quantities
  ! that might also be an output variable for climate mode

  ! ntri=4 variables
  IF (ldiagnos) THEN
    pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'TDIV_HUM  '
    pp_restart%yvarml(pp_restart%nyvar_m + 2)  = 'AEVAP_S   '
    pp_restart%nyvar_m     = pp_restart%nyvar_m + 2
  ENDIF

  ! ntri=3 variables
  IF (lbdclim) THEN
    pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'T_2M_AV   '
    pp_restart%yvarml(pp_restart%nyvar_m + 2)  = 'TD_2M_AV  '
    pp_restart%yvarml(pp_restart%nyvar_m + 3)  = 'U_10M_AV  '
    pp_restart%yvarml(pp_restart%nyvar_m + 4)  = 'V_10M_AV  '
    pp_restart%nyvar_m     = pp_restart%nyvar_m + 4
  ENDIF
  pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'ASWDIR_S  '
  pp_restart%yvarml(pp_restart%nyvar_m + 2)  = 'ASWDIFD_S '
  pp_restart%yvarml(pp_restart%nyvar_m + 3)  = 'ASWDIFU_S '
  pp_restart%yvarml(pp_restart%nyvar_m + 4)  = 'ALWD_S    '
  pp_restart%yvarml(pp_restart%nyvar_m + 5)  = 'ALWU_S    '
  pp_restart%nyvar_m     = pp_restart%nyvar_m + 5

#ifdef COSMOART
  ! Restart fields for COSMOART
  IF (l_cosmo_art) THEN
    ! CK 20101213 restart fields now set within organize_cosmo_art
    CALL art_restart_fields("old",istat)
  ENDIF
#endif

#ifdef POLLEN
  IF (l_pollen) THEN
    pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'CNC_BETU '
    pp_restart%yvarml(pp_restart%nyvar_m + 2)  = 'FR_BETU  '
    pp_restart%yvarml(pp_restart%nyvar_m + 3)  = 'QP0_BETU '
    pp_restart%yvarml(pp_restart%nyvar_m + 4)  = 'TSS_BETU '
    pp_restart%yvarml(pp_restart%nyvar_m + 5)  = 'TSE_BETU '
    pp_restart%yvarml(pp_restart%nyvar_m + 6)  = 'TTS_BETU '
    pp_restart%yvarml(pp_restart%nyvar_m + 7)  = 'TTE_BETU '
    pp_restart%yvarml(pp_restart%nyvar_m + 8)  = 'HCEM_BETU'
    pp_restart%yvarml(pp_restart%nyvar_m + 9)  = 'SDES_BETU'
    pp_restart%nyvar_m     = pp_restart%nyvar_m + 9
  ENDIF
#endif

  ALLOCATE(list_res_o  (pp_restart%nyvar_m), STAT=istat)

  ! determine the indices in the list of variable table
  ! and initialize list_res_o
  ! list 1 of restart variables
  loop_list_r1: DO n = 1, pp_restart%nyvar_m
    list_res_o(n)%name = pp_restart%yvarml(n)
    DO iz3 = 1, num_gribtabs
      DO iz2 = 0, 255
        DO iz1 = 1,4
          IF (TRIM(var(iz1,iz2,iz3)%name) == TRIM(pp_restart%yvarml(n))) THEN
            ! Set location in variable table
            pp_restart%ilist_ml(1,n) = iz1
            pp_restart%ilist_ml(2,n) = iz2
            pp_restart%ilist_ml(3,n) = iz3

            list_res_o(n)%iloc1      = iz1
            list_res_o(n)%iloc2      = iz2
            list_res_o(n)%iloc3      = iz3

            ! Set dimension of this variable (1, ke or ke+1)
            SELECT CASE (var(iz1,iz2,iz3)%rank)
            CASE (4)
              list_res_o(n)%idimvert = UBOUND (var(iz1,iz2,iz3)%p4,3)
            CASE (3)
              SELECT CASE (var(iz1,iz2,iz3)%name)
              CASE ('PS        ','T_G       ','QV_S      ','W_SNOW    ',  &
                    'T_S       ','T_M       ','W_G1      ','W_G2      ',  &
                    'W_G3      ','W_I       ','T_SNOW    ','RHO_SNOW  ',  &
                    'H_SNOW    ','H_ICE     ','T_ICE     ','T_MNW_LK  ',  &
                    'T_WML_LK  ','T_BOT_LK  ','T_B1_LK   ','C_T_LK    ',  &
                    'H_ML_LK   ','H_B1_LK   ')
                ! these are 2D variables with rank=3 because of time dependency
                list_res_o(n)%idimvert = 1
              CASE DEFAULT
                ! these are real 3D variables
                list_res_o(n)%idimvert = UBOUND (var(iz1,iz2,iz3)%p3,3)
              END SELECT
            CASE (2)
              list_res_o(n)%idimvert = 1
            END SELECT

            CYCLE loop_list_r1
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO loop_list_r1



  ! list of restart variables: "now" timestep
  IF (.NOT. l2tls) THEN
    ! atmospheric variables
    pp_restart%yvarpl( 1)  = 'U         '
    pp_restart%yvarpl( 2)  = 'V         '
    pp_restart%yvarpl( 3)  = 'W         '
    pp_restart%yvarpl( 4)  = 'T         '
    pp_restart%yvarpl( 5)  = 'QV        '
    pp_restart%yvarpl( 6)  = 'QC        '
    pp_restart%yvarpl( 7)  = 'PP        '
    pp_restart%yvarpl( 8)  = 'T_SNOW    '
    pp_restart%yvarpl( 9)  = 'W_I       '
    pp_restart%yvarpl(10)  = 'QV_S      '
    pp_restart%yvarpl(11)  = 'W_SNOW    '
    pp_restart%yvarpl(12)  = 'T_S       '
    pp_restart%yvarpl(13)  = 'T_G       '
    pp_restart%yvarpl(14)  = 'PS        '

    pp_restart%nyvar_p     =  14

    SELECT CASE (itype_turb)
    CASE (3)
      pp_restart%yvarpl(pp_restart%nyvar_p + 1)  = 'TKE       '
      pp_restart%nyvar_p     = pp_restart%nyvar_p + 1
    CASE (5:8)
      IF (lprog_tke) THEN
        pp_restart%yvarpl(pp_restart%nyvar_p + 1)  = 'TKE       '
        pp_restart%nyvar_p     = pp_restart%nyvar_p + 1
      ENDIF
    END SELECT

    IF (lprog_qi) THEN
      pp_restart%yvarpl(pp_restart%nyvar_p + 1)  = 'QI        '
      pp_restart%nyvar_p     = pp_restart%nyvar_p + 1
    ENDIF
    IF (lprogprec) THEN
      pp_restart%yvarpl(pp_restart%nyvar_p + 1)  = 'QR        '
      pp_restart%nyvar_p     = pp_restart%nyvar_p + 1

      IF (itype_gscp > 1) THEN
        pp_restart%yvarpl(pp_restart%nyvar_p + 1)  = 'QS        '
        pp_restart%nyvar_p     = pp_restart%nyvar_p + 1
      ENDIF
      IF (itype_gscp == 4) THEN
        pp_restart%yvarpl(pp_restart%nyvar_p + 1)  = 'QG        '
        pp_restart%nyvar_p     = pp_restart%nyvar_p + 1
      ENDIF
    ENDIF

    ! variables for the soil model
    IF (lmulti_layer) THEN
      pp_restart%yvarpl(pp_restart%nyvar_p + 1)  = 'T_SO      '
      pp_restart%yvarpl(pp_restart%nyvar_p + 2)  = 'W_SO      '
      pp_restart%yvarpl(pp_restart%nyvar_p + 3)  = 'W_SO_ICE  '
      pp_restart%yvarpl(pp_restart%nyvar_p + 4)  = 'RHO_SNOW  '
      pp_restart%nyvar_p     = pp_restart%nyvar_p + 4
    ELSE
      pp_restart%yvarpl(pp_restart%nyvar_p + 1)  = 'T_M       '
      pp_restart%yvarpl(pp_restart%nyvar_p + 2)  = 'W_G1      '
      pp_restart%yvarpl(pp_restart%nyvar_p + 3)  = 'W_G2      '
      pp_restart%nyvar_p     = pp_restart%nyvar_p + 3
    ENDIF

    ! variables for the sea ice and/or the FLake model
    IF (lseaice .OR. llake) THEN
      pp_restart%yvarpl(pp_restart%nyvar_p + 1)  = 'T_ICE     '
      pp_restart%yvarpl(pp_restart%nyvar_p + 2)  = 'H_ICE     '
      pp_restart%nyvar_p     = pp_restart%nyvar_p + 2
    ENDIF

    ! variables for the FLake model
    IF (llake) THEN
      pp_restart%yvarpl(pp_restart%nyvar_p + 1)  = 'T_MNW_LK  '
      pp_restart%yvarpl(pp_restart%nyvar_p + 2)  = 'T_WML_LK  '
      pp_restart%yvarpl(pp_restart%nyvar_p + 3)  = 'T_BOT_LK  '
      pp_restart%yvarpl(pp_restart%nyvar_p + 4)  = 'C_T_LK    '
      pp_restart%yvarpl(pp_restart%nyvar_p + 5)  = 'H_ML_LK   '
      pp_restart%nyvar_p     = pp_restart%nyvar_p + 5

      !_cdm In the present configuration,
      !     the bottom-sediment module of FLake is switched off.
      !     The respective fields are handled internally by LM
      !     and are not included into the LM IO lists.
      !       'T_B1_LK   '
      !       'H_B1_LK   '
    ENDIF

#ifdef COSMOART
    ! Restart fields for COSMO_ART
    IF (l_cosmo_art) THEN
      ! CK 20101213 restart fields now set within organize_cosmo_art
      CALL art_restart_fields("now",istat)
    ENDIF
#endif

#ifdef POLLEN
    IF (l_pollen) THEN
      pp_restart%yvarpl(pp_restart%nyvar_p + 1)  = 'CNC_BETU   '
      pp_restart%nyvar_p     = pp_restart%nyvar_p + 1
    ENDIF
#endif

    ALLOCATE(list_res_n(pp_restart%nyvar_p), STAT=istat)

    ! determine the indices in the list of variable table
    ! and initialize list_res_n
    ! list 2 of restart variables
    loop_list_r2: DO n = 1, pp_restart%nyvar_p
      list_res_n(n)%name = pp_restart%yvarpl(n)
      DO iz3 = 1, num_gribtabs
        DO iz2 = 0, 255
          DO iz1 = 1,4
            IF ( TRIM(var(iz1,iz2,iz3)%name) == TRIM(pp_restart%yvarpl(n))) THEN
              ! Set location in variable table
              pp_restart%ilist_pl(1,n) = iz1
              pp_restart%ilist_pl(2,n) = iz2
              pp_restart%ilist_pl(3,n) = iz3

              list_res_n(n)%iloc1      = iz1
              list_res_n(n)%iloc2      = iz2
              list_res_n(n)%iloc3      = iz3

              ! Set dimension of this variable (1, ke or ke+1)
              SELECT CASE (var(iz1,iz2,iz3)%rank)
              CASE (4)
                list_res_n(n)%idimvert = UBOUND (var(iz1,iz2,iz3)%p4,3)
              CASE (3)
                SELECT CASE (var(iz1,iz2,iz3)%name)
              CASE ('PS        ','T_G       ','QV_S      ','W_SNOW    ',  &
                    'T_S       ','T_M       ','W_G1      ','W_G2      ',  &
                    'W_G3      ','W_I       ','T_SNOW    ','RHO_SNOW  ',  &
                    'H_SNOW    ','H_ICE     ','T_ICE     ','T_MNW_LK  ',  &
                    'T_WML_LK  ','T_BOT_LK  ','T_B1_LK   ','C_T_LK    ',  &
                    'H_ML_LK   ','H_B1_LK   ')
                  ! these are 2D variables with rank=3 because of time dependency
                  list_res_n(n)%idimvert = 1
                CASE DEFAULT
                  ! these are real 3D variables
                  list_res_n(n)%idimvert = UBOUND (var(iz1,iz2,iz3)%p3,3)
                END SELECT
              CASE (2)
                list_res_n(n)%idimvert = 1
              END SELECT

              CYCLE loop_list_r2
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO loop_list_r2
  ENDIF

  IF (n_num_now > 1) RETURN

  !----------------------------------------------------------------------------
  ! Section 2.2: Read or create the data
  !----------------------------------------------------------------------------

  ! Find the hours and steps for the first two boundary data sets
  ! (is 0 and hincbound for non-restarts)
  hlastbound = 0.0_ireals
  endless: DO
    IF ( (hlastbound <= hstart) .AND. (hstart < hlastbound+hincbound) ) THEN
      EXIT endless
    ENDIF
    hlastbound = hlastbound + hincbound
  ENDDO endless

  ! the starting hour is also a boundary update hour
  hnextbound = hlastbound + hincbound
  nlastbound = NINT (3600.0_ireals * hlastbound / dt)
  nnextbound = NINT (3600.0_ireals * hnextbound / dt)
  nincbound  = nnextbound - nlastbound

  IF (izdebug >= 10) THEN
    PRINT *, '       Computation of control variables for boundary update: '
    PRINT *, '        (lastbound, nextbound in h/n:)   ',    &
              hlastbound, hnextbound, nlastbound, nnextbound, nincbound
  ENDIF

  ibdf = 1
  IF (lartif_data .EQV. .FALSE.) THEN
    ! Set correct levtop, levbot for soil moisture content
    IF (.NOT. lmulti_layer) THEN
      IF (nlgw_ini == 2) THEN
        var(1, 86,1)%levtop =    0;    var(1, 86,1)%levbot =  10
        var(2, 86,1)%levtop =   10;    var(2, 86,1)%levbot = 100
        var(3, 86,1)%levtop =    0;    var(3, 86,1)%levbot =   0
      ELSEIF (nlgw_ini == 3) THEN
        var(1, 86,1)%levtop =    0;    var(1, 86,1)%levbot =   2
        var(2, 86,1)%levtop =    2;    var(2, 86,1)%levbot =  10
        var(3, 86,1)%levtop =   10;    var(3, 86,1)%levbot = 100
      ENDIF
    ENDIF

    IF (nstart == 0) THEN
      ! read the initial file
      IF (izdebug >= 10) THEN
        PRINT *, '     Read initial data for hstart = ', hstart
      ENDIF

      CALL organize_input ('initial', nlgw, nlgw_ini, list_ini, nyvar_i,    &
                           yform_read, ' ')
    ELSE
      ! read the restart file(s)
      IF (izdebug >= 10) THEN
        PRINT *, '     Read restart data for hstart = ', hstart
      ENDIF

      CALL organize_input ('restart', nlgw, nlgw_ini, list_res_o,           &
                           pp_restart%nyvar_m, 'bina', 'o')
      IF (.NOT. l2tls) THEN
        CALL organize_input ('restart', nlgw, nlgw_ini, list_res_n,         &
                             pp_restart%nyvar_p, 'bina', 'n')
      ENDIF
    ENDIF

    IF (      (MOD(nstart,newbcdt) >= newbcdt-nincboufac*nincbound)         &
        .AND. (nstart/newbcdt <= newbc-1))                                  &
       ibdf = ((1+ nlastbound/newbcdt) *newbcdt - nlastbound) / nincbound
    nincbound  = nincbound * ibdf
    IF (izdebug >= 1) THEN
      tact = dt / 3600.0_ireals * REAL(nstep, ireals)
      tlbc = dt / 3600.0_ireals * REAL(nlastbound + nincbound, ireals)
      tivl = dt / 3600.0_ireals * REAL(nincbound, ireals)
      WRITE (*,'(A,F9.2,A,F9.2,A,F4.2,I4)')                                 &
             ' LATERAL BC: act. hour ', tact, 'hour of NEW BC field ',tlbc, &
             'interval ', tivl, ibdf
    ENDIF

    ! Set correct levtop, levbot for soil moisture content
    IF (.NOT. lmulti_layer) THEN
      IF (nlgw_bd == 2) THEN
        var(1, 86,1)%levtop =    0;    var(1, 86,1)%levbot =  10
        var(2, 86,1)%levtop =   10;    var(2, 86,1)%levbot = 100
        var(3, 86,1)%levtop =    0;    var(3, 86,1)%levbot =   0
      ELSEIF (nlgw_bd == 3) THEN
        var(1, 86,1)%levtop =    0;    var(1, 86,1)%levbot =   2
        var(2, 86,1)%levtop =    2;    var(2, 86,1)%levbot =  10
        var(3, 86,1)%levtop =   10;    var(3, 86,1)%levbot = 100
      ENDIF
    ENDIF

    IF (nstart == 0) THEN
      CALL organize_input ('boundary', nlgw, nlgw_bd, list_bd, nyvar_b,     &
                           yform_read, ' ', .TRUE.,  1, nlastbound)
    ELSE
      CALL organize_input ('boundary', nlgw, nlgw_bd, list_bd, nyvar_b,     &
                           yform_read, ' ', .FALSE.,  1, nlastbound)
    ENDIF
    CALL organize_input ('boundary', nlgw, nlgw_bd, list_bd, nyvar_b,     &
                         yform_read, ' ', .FALSE., 2, nlastbound+nincbound)
  ELSE
    IF (nstart == 0) THEN
      IF (my_cart_id == 0) THEN
        PRINT *,'    GENERATE ARTIFICIAL DATA'
      ENDIF

      CALL gen_ini_data
      CALL gen_bound_data (1)
    ELSE
! UB>>
! Has been commented to allow for restart runs in idealized studies.
!!$      ierror  = 6101
!!$      yerrmsg = 'No generation of artifical data for nstart > 0'
!!$      RETURN
! UB<<

      !.. We need to generate the boundary fields here,
      !   which are not part of the restart list.
      !   Therefore, we recover them as in the first original
      !   model time step and compute them also for
      !   the restart-timestep, if they are time dependent.
      ! NOTE: we have to do this also for periodic
      !   runs, because bd-fields for timelevel nnew are
      !   needed to initialize the nnew-
      !   values at the beginning of the first model time step.
      !   However, they have no further effect in periodic runs.

      IF (my_cart_id == 0) THEN
        WRITE (*,*) &
        '    RECOVERING THE TIME-CONSTANT BOUNDARY FIELDS FOR THE IDEALIZED RUN'
      ENDIF

      !    1) generating initial data like in the first time step
      !       of the original run:
      CALL gen_ini_data
      !    2) copy them to the boundary fields like in the
      !       first timestep of the original run:
      CALL gen_bound_data (1)
      !    3) Generate boundary fields for later timesteps (if time dependent):
      CALL gen_bound_data (2)


      !.. Read in restart-fields (which overwrite the "wrong" actual model fields
      !      which have been initialized by the above call to gen_ini_data():

      IF (my_cart_id == 0) THEN
        WRITE (*,*) '    READ IN RESTART DATA'
      ENDIF

      ! read the restart file(s)
      CALL organize_input ('restart', nlgw, nlgw_ini, list_res_o,           &
                           pp_restart%nyvar_m, 'bina', 'o')
      IF (.NOT. l2tls) THEN
        CALL organize_input ('restart', nlgw, nlgw_ini, list_res_n,         &
                             pp_restart%nyvar_p, 'bina', 'n')
      ENDIF

    ENDIF
  ENDIF

  ! Print vertical coordinate parameters in control output
  CALL print_vertcoord

  ! Initialize the output
  now => root
  gribout_loop_init: DO

    ! list of constant variables
    yzloclist(:) = '          '
    n1 = 0_iintegers
    loop_list_c: DO n = 1, now%nyvar_c
      DO iz3 = 1, num_gribtabs
        DO iz2 = 0, 255
          DO iz1 = 1,4
            IF (var(iz1,iz2,iz3)%name(1:LEN_TRIM(var(iz1,iz2,iz3)%name)) ==  &
                        now%yvarc (n)(1:LEN_TRIM(now%yvarc (n))) ) THEN
              IF (var(iz1,iz2,iz3)%idef_stat < 0) THEN
                IF (my_cart_id == 0) THEN
                  PRINT *,'variable ', now%yvarc(n)(1:LEN_TRIM(now%yvarc(n))),&
                           ' is not allocated and is removed from the list'
                ENDIF
                CYCLE loop_list_c
              ENDIF

              n1 = n1+1
              now%ilist_c (1,n1) = iz1
              now%ilist_c (2,n1) = iz2
              now%ilist_c (3,n1) = iz3
              yzloclist(n1)(1:LEN_TRIM(now%yvarc (n))) =                &
                                 now%yvarc (n)(1:LEN_TRIM(now%yvarc (n)))
              CYCLE loop_list_c
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      ! If this point is reached, no variable with name now%yvarc(n) was found
      IF (my_cart_id == 0) THEN
        PRINT *, 'Variable ', now%yvarc(n), ' was not found in the GRIB table'
      ENDIF
    ENDDO loop_list_c
  
    ! now adjust the list
    now%nyvar_c = n1
    DO n = 1, now%nyvar_c
      now%yvarc (n)(1:10) = yzloclist(n)(1:10)
    ENDDO

    ! list of model level variables
    ! check the list of model level variables, if there are e.g. variables,
    ! that are not computed during the LM run
    yzloclist(:) = '          '
    n1 = 0_iintegers
    loop_list_m: DO n = 1, now%nyvar_m
      DO iz3 = 1, num_gribtabs
        DO iz2 = 0, 255
          DO iz1 = 1,4
            IF (var(iz1,iz2,iz3)%name(1:LEN_TRIM(var(iz1,iz2,iz3)%name)) ==  &
                        now%yvarml(n)(1:LEN_TRIM(now%yvarml(n))) ) THEN
              IF (var(iz1,iz2,iz3)%idef_stat < 0) THEN
                IF (my_cart_id == 0) THEN
                  PRINT *,'variable ', now%yvarml(n)(1:LEN_TRIM(now%yvarml(n))),&
                           ' is not allocated and is removed from the list'
                ENDIF
                CYCLE loop_list_m
              ENDIF

              ! Check for special variables that cannot be written on model levels
              SELECT CASE (TRIM(now%yvarml(n)))
              CASE ('FI')
                IF (my_cart_id == 0) THEN
                  PRINT *,'variable ', TRIM(now%yvarml(n)),        &
                          ' cannot be written on modellevels'
                ENDIF
                CYCLE loop_list_m
              END SELECT

              ! if this point is reached, the variable can be written
              n1 = n1+1
              now%ilist_ml(1,n1) = iz1
              now%ilist_ml(2,n1) = iz2
              now%ilist_ml(3,n1) = iz3
              yzloclist(n1)(1:LEN_TRIM(now%yvarml(n))) =                &
                                 now%yvarml(n)(1:LEN_TRIM(now%yvarml(n)))
!             IF (lzreset) THEN
!               IF ( (var(iz1,iz2,iz3)%ntri == 3) .OR.         &
!                    (var(iz1,iz2,iz3)%ntri == 4) ) THEN
!                 IF (ASSOCIATED(var(iz1,iz2,iz3)%p2)) THEN
!                   var(iz1,iz2,iz3)%p2(:,:) = 0.0_ireals
!                 ENDIF
!               ENDIF
!             ENDIF

              CYCLE loop_list_m
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      ! If this point is reached, no variable with name now%yvarml(n) was found
      IF (my_cart_id == 0) THEN
        PRINT *, 'Variable ', now%yvarml(n), ' was not found in the GRIB table'
      ENDIF
    ENDDO loop_list_m

    ! now adjust the list
    now%nyvar_m = n1
    DO n = 1, now%nyvar_m
      now%yvarml(n)(1:10) = yzloclist(n)(1:10)
    ENDDO

    ! list of pressure level variables
    ! check the list of model level variables, if there are e.g. variables,
    ! that are not computed during the LM run
    yzloclist(:) = '          '
    n1 = 0_iintegers
    loop_list_p: DO n = 1, now%nyvar_p
      DO iz3 = 1, num_gribtabs
        DO iz2 = 0, 255
          DO iz1 = 1,4
            IF (var(iz1,iz2,iz3)%name(1:LEN_TRIM(var(iz1,iz2,iz3)%name)) ==  &
                        now%yvarpl(n)(1:LEN_TRIM(now%yvarpl(n))) ) THEN
              IF (var(iz1,iz2,iz3)%idef_stat < 0) THEN
                IF (my_cart_id == 0) THEN
                  PRINT *,'variable ', now%yvarpl(n)(1:LEN_TRIM(now%yvarpl(n))),&
                           ' is not allocated and is removed from the list'
                ENDIF
                CYCLE loop_list_p
              ENDIF

              ! if this point is reached, the variable can be written
              n1 = n1+1
              now%ilist_pl(1,n1) = iz1
              now%ilist_pl(2,n1) = iz2
              now%ilist_pl(3,n1) = iz3
              yzloclist(n1)(1:LEN_TRIM(now%yvarpl(n))) =                &
                                 now%yvarpl(n)(1:LEN_TRIM(now%yvarpl(n)))
              CYCLE loop_list_p
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      ! If this point is reached, no variable with name now%yvarpl(n) was found
      IF (my_cart_id == 0) THEN
        PRINT *, 'Variable ', now%yvarpl(n), ' was not found in the GRIB table'
      ENDIF
    ENDDO loop_list_p

    ! now adjust the list
    now%nyvar_p = n1
    DO n = 1, now%nyvar_p
      now%yvarpl(n)(1:10) = yzloclist(n)(1:10)
    ENDDO

    ! list of z-level variables
    yzloclist(:) = '          '
    n1 = 0_iintegers
    loop_list_z: DO n = 1, now%nyvar_z
      DO iz3 = 1, num_gribtabs
        DO iz2 = 0, 255
          DO iz1 = 1,4
            IF (var(iz1,iz2,iz3)%name(1:LEN_TRIM(var(iz1,iz2,iz3)%name)) ==  &
                        now%yvarzl(n)(1:LEN_TRIM(now%yvarzl(n))) ) THEN
              IF (var(iz1,iz2,iz3)%idef_stat < 0) THEN
                IF (my_cart_id == 0) THEN
                  PRINT *,'variable ', now%yvarzl(n)(1:LEN_TRIM(now%yvarzl(n))),&
                           ' is not allocated and is removed from the list'
                ENDIF
                CYCLE loop_list_z
              ENDIF

              ! if this point is reached, the variable can be written
              n1 = n1+1
              now%ilist_zl(1,n1) = iz1
              now%ilist_zl(2,n1) = iz2
              now%ilist_zl(3,n1) = iz3
              yzloclist(n1)(1:LEN_TRIM(now%yvarzl(n))) =                &
                                 now%yvarzl(n)(1:LEN_TRIM(now%yvarzl(n)))
              CYCLE loop_list_z
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      ! If this point is reached, no variable with name now%yvarzl(n) was found
      IF (my_cart_id == 0) THEN
        PRINT *, 'Variable ', now%yvarzl(n), ' was not found in the GRIB table'
      ENDIF
    ENDDO loop_list_z

    ! now adjust the list
    now%nyvar_z = n1
    DO n = 1, now%nyvar_z
      now%yvarzl(n)(1:10) = yzloclist(n)(1:10)
    ENDDO

    ! list of satellite channel variables
    yzloclist(:) = '          '
    n1 = 0_iintegers
    loop_list_s: DO n = 1, now%nyvar_s
      DO iz3 = 1, num_gribtabs
        DO iz2 = 0, 255
          DO iz1 = 1,4
            IF (var(iz1,iz2,iz3)%name(1:LEN_TRIM(var(iz1,iz2,iz3)%name)) ==  &
                        now%yvarsl(n)(1:LEN_TRIM(now%yvarsl(n))) ) THEN

! we cannot check this here, because these variables are allocated later
!             IF (var(iz1,iz2,iz3)%idef_stat < 0) THEN
!               IF (my_cart_id == 0) THEN
!                 PRINT *,'variable ', now%yvarsl(n)(1:LEN_TRIM(now%yvarsl(n))),&
!                          ' is not allocated and is removed from the list'
!               ENDIF
!               CYCLE loop_list_s
!             ENDIF

              ! check that the variables MSG_xx are NOT written for
              ! grib output; they can only be written in NetCDF
              IF (now%yform_write /= 'ncdf') THEN
                IF ( (var(iz1,iz2,iz3)%name == 'MSG_TB    ') .OR.  &
                     (var(iz1,iz2,iz3)%name == 'MSG_TBC   ') .OR.  &
                     (var(iz1,iz2,iz3)%name == 'MSG_RAD   ') .OR.  &
                     (var(iz1,iz2,iz3)%name == 'MSG_RADC  ') ) THEN
                  IF (my_cart_id == 0) THEN
                    PRINT *,'variable ', TRIM(now%yvarsl(n)), &
                            ' can only be written for Grib Output'
                  ENDIF
                  CYCLE loop_list_s
                ENDIF
              ENDIF

              ! if this point is reached, the variable can be written
              n1 = n1+1
              now%ilist_sl(1,n) = iz1
              now%ilist_sl(2,n) = iz2
              now%ilist_sl(3,n) = iz3
              yzloclist(n1)(1:LEN_TRIM(now%yvarsl(n))) =                &
                                 now%yvarsl(n)(1:LEN_TRIM(now%yvarsl(n)))
              CYCLE loop_list_s
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      ! If this point is reached, no variable with name now%yvarsl(n) was found
      IF (my_cart_id == 0) THEN
        PRINT *, 'Variable ', now%yvarsl(n), ' was not found in the GRIB table'
      ENDIF
    ENDDO loop_list_s

    ! now adjust the list
    now%nyvar_s = n1
    DO n = 1, now%nyvar_s
      now%yvarsl(n)(1:10) = yzloclist(n)(1:10)
    ENDDO

    IF (ASSOCIATED(now%next)) THEN
      now => now%next
    ELSE
      EXIT gribout_loop_init
    ENDIF
  ENDDO gribout_loop_init

  CALL init_output (root, n_num_now)

!------------------------------------------------------------------------------
! Section 3: Input of boundary data
!------------------------------------------------------------------------------

ELSEIF (yaction == 'boundary') THEN

  IF (lartif_data) THEN
    !    3) Generate time dependent boundary fields 
    !      (nothing implemented in gen_bound_data so far ...):
    CALL gen_bound_data (2)
  END IF

  IF ( (nstep+1 > nnextbound) .AND. (nstep < nstop) ) THEN

    ! Initialize LM variable table again, if the nest has changed to
    ! get the cross reference for the pointers right
    IF ( nsetupvar /= n_num_now ) THEN
      CALL setup_vartab

#ifdef COSMOART
      CALL setup_vartab_art
#endif

#ifdef POLLEN
      CALL setup_vartab_pol
#endif

      nsetupvar = n_num_now
    ENDIF

    nbd1 = 3 - nbd1
    nbd2 = 3 - nbd2

    hlastbound = hnextbound
    hnextbound = hlastbound + hincbound
    nlastbound = nnextbound
    nnextbound = NINT (3600.0_ireals * hnextbound / dt)
    nincbound  = nnextbound - nlastbound

    IF (izdebug >= 10) THEN
      PRINT *, '       Computation of control variables for boundary update: '
      PRINT *, '        (lastbound, nextbound in h/n):  ',               &
                hlastbound, hnextbound, nlastbound, nnextbound, nincbound
    ENDIF

    IF (ibdf > 1) nincbound = nincbound / ibdf
    ibdf = 1
    IF (      (MOD(nstep,newbcdt) >= newbcdt-nincboufac*nincbound)          &
        .AND. (nstep/newbcdt <= newbc-1))                                   &
      ibdf = ((1+ nlastbound/newbcdt) *newbcdt - nlastbound) / nincbound
    nincbound  = nincbound * ibdf
    IF (izdebug >= 1) THEN
      tact = dt / 3600.0_ireals * REAL(nstep, ireals)
      tlbc = dt / 3600.0_ireals * REAL(nlastbound + nincbound, ireals)
      tivl = dt / 3600.0_ireals * REAL(nincbound, ireals)
      WRITE (*,'(A,F6.2,A,F6.2,A,F6.2,I4)')                                 &
             ' LATERAL BC: act. hour ', tact, 'hour of NEW BC field ',tlbc, &
             'interval ', tivl, ibdf
    ENDIF
    IF (lartif_data .EQV. .FALSE.) THEN
      ! Set correct levtop, levbot for soil moisture content
      IF (.NOT. lmulti_layer) THEN
        IF (nlgw_bd == 2) THEN
          var(1, 86,1)%levtop =    0;    var(1, 86,1)%levbot =  10
          var(2, 86,1)%levtop =   10;    var(2, 86,1)%levbot = 100
          var(3, 86,1)%levtop =    0;    var(3, 86,1)%levbot =   0
        ELSEIF (nlgw_bd == 3) THEN
          var(1, 86,1)%levtop =    0;    var(1, 86,1)%levbot =   2
          var(2, 86,1)%levtop =    2;    var(2, 86,1)%levbot =  10
          var(3, 86,1)%levtop =   10;    var(3, 86,1)%levbot = 100
        ENDIF
      ENDIF

      CALL organize_input ('boundary', nlgw, nlgw_bd, list_bd, nyvar_b,      &
                           yform_read, ' ', .FALSE., nbd2, nnextbound)
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 4: Output of data
!------------------------------------------------------------------------------
  
ELSEIF (yaction == 'result') THEN

  ! define startpoint of the list
  now => root
  lwritefc    = .FALSE.
  lwriteana   = .FALSE.
  lwriteready = .FALSE.

  ! Set correct levtop, levbot for soil moisture content
  IF (.NOT. lmulti_layer) THEN
    IF (nlgw == 2) THEN
      var(1, 86,1)%levtop =    0;    var(1, 86,1)%levbot =  10
      var(2, 86,1)%levtop =   10;    var(2, 86,1)%levbot = 100
      var(3, 86,1)%levtop =    0;    var(3, 86,1)%levbot =   0
    ELSEIF (nlgw == 3) THEN
      var(1, 86,1)%levtop =    0;    var(1, 86,1)%levbot =   2
      var(2, 86,1)%levtop =    2;    var(2, 86,1)%levbot =  10
      var(3, 86,1)%levtop =   10;    var(3, 86,1)%levbot = 100
    ENDIF
  ENDIF

  gribout_loop: DO

    ! this is necessary for running with digital filtering (ndfi=1),
    ! because timestepping then starts with nstart > 0.
    ! (but it would be better to adjust calc_ngrib)
    IF ( (nstart > 0) .AND. (now%ngrib(now%nextstep) == 0) ) THEN
      now%nextstep = MIN( now%nextstep + 1, now%outsteps)
    ENDIF

    ! current timestep is the next output timestep of the namelist
    IF ( now%n_num == n_num_now ) THEN
      IF (   (nstep == now%ngrib(now%nextstep))              .AND.  &
           ( (nstep /= nstop) .OR. (nstep == nfinalstop) ) )        THEN
        ! Do output only, if it is not the last step of a simulation period
        ! But do it, if it is the last step of the whole simulation

        ! Initialize LM variable table again, if the nest has changed to
        ! get the cross reference for the pointers right and
        ! set grid description section for this grid
        IF ( nsetupvar /= n_num_now ) THEN
          CALL setup_vartab

#ifdef COSMOART
          CALL setup_vartab_art
#endif

#ifdef POLLEN
          CALL setup_vartab_pol
#endif

          CALL makegds
          nsetupvar = n_num_now
        ENDIF

#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
        ! output of satellite channels variables
        IF (luse_rttov .AND. now%nyvar_s > 0 .AND. isynsat_stat == 0) THEN
          CALL organize_output (now, 's', now%nyvar_s, now%yvarsl,        &
                                now%ilist_sl)
          lwriteready = .TRUE.
        ENDIF
#endif

        ! interpolation/output of data on constant pressure levels
        IF (now%nyvar_p > 0) THEN
          CALL organize_output (now, 'p', now%nyvar_p, now%yvarpl,        &
                                now%ilist_pl)
          lwriteready = .TRUE.
        ENDIF

        ! interpolation/output of data on constant height levels
        IF (now%nyvar_z > 0) THEN
          CALL organize_output (now, 'z', now%nyvar_z, now%yvarzl,        &
                                now%ilist_zl)
          lwriteready = .TRUE.
        ENDIF

        ! output of all model variables in the output lists
        IF (now%nyvar_m > 0) THEN
          CALL organize_output (now, ' ', now%nyvar_m, now%yvarml,        &
                                now%ilist_ml)
          lwriteready = .TRUE.
        ENDIF

        now%nextstep = MIN( now%nextstep + 1, now%outsteps)
        IF (izdebug > 10) THEN
          PRINT *, '     Next output will be done in step:  ',     &
                     now%ngrib(now%nextstep), now%nextstep
        ENDIF

        ! check, whether forecasts or analyses are written in this output step
        lwritefc  = (lwritefc ) .OR. (.NOT. now%lanalysis)
        lwriteana = (lwriteana) .OR. (      now%lanalysis)

      ENDIF
    ENDIF

    IF (ASSOCIATED(now%next)) THEN
      ! exchange current namelist with the next of the list
      now => now%next
    ELSE
      EXIT gribout_loop
    ENDIF
  ENDDO gribout_loop


  ! Write restart-files, if necessary
  IF ((nstep == pp_restart%nextstep) .AND.                                   &
      (pp_restart%nextstep <= NINT (3600.0_ireals * nhour_restart(2) / dt))) &
                                                                         THEN
    fc_hour = REAL(pp_restart%nextstep, ireals)
    ! write restart files
    CALL organize_output (pp_restart, 'o', pp_restart%nyvar_m,               &
                          pp_restart%yvarml, pp_restart%ilist_ml)

    IF (.NOT. l2tls) THEN
      CALL organize_output (pp_restart, 'n', pp_restart%nyvar_p,               &
                            pp_restart%yvarpl, pp_restart%ilist_pl)
    ENDIF

    ! determine next output hour and step
    pp_restart%nexthour = pp_restart%nexthour + nhour_restart(3)
    pp_restart%nextstep = NINT (3600.0_ireals * pp_restart%nexthour / dt) - 1
  ENDIF

  ! Write ready-files, if required
  IF ((lwriteready) .AND. (ytrans_out /= '   ') .AND. (my_cart_id == 0)) THEN
    IF (lwriteana) THEN
      yzhead = 'LMA'
    ELSEIF (lwritefc ) THEN
      yzhead = 'LMF'
    ENDIF

    ! Create the filename LMF_forecasttime
    CALL make_fn (yzhead, yakdat1, 'f',' ', nstep, dt, .TRUE., ytrans_out,   &
                  yzname, izdebug, izerrstat)

    ! Add the nest specification (except for the mother nest)
    IF (now%n_num > 1) THEN
      yzname = yzname(1:LEN_TRIM(yzname)) // c_nnum
    ENDIF

    ! Write the file
    OPEN  (ntrans_out, FILE=yzname, FORM='FORMATTED')
    WRITE (ntrans_out, '(A)') 'ready'
    CLOSE (ntrans_out)

  ENDIF

!------------------------------------------------------------------------------
! Section 5: All other actions are wrong
!------------------------------------------------------------------------------

ELSE

    ierror  = 1
    yerrmsg = 'ERROR *** No valid action for organize_data ***'

ENDIF
  
!------------------------------------------------------------------------------
! Internal procedures
!------------------------------------------------------------------------------

CONTAINS

!==============================================================================
!+ Print the vertical coordinate parameters on file yuspecif
!------------------------------------------------------------------------------

SUBROUTINE print_vertcoord

!------------------------------------------------------------------------------
!
! Description:
!
! Method:
!
!==============================================================================

INTEGER (KIND=iintegers)   :: k  

!- End of header
!==============================================================================

! Print vertical coordinate parameters on unit 'yuspecif'

  IF(my_cart_id == 0) THEN

    IF (ivctype == 1 ) THEN
      write(nuspecif, '(A)' )                                             &
             '0     Vertical coordinate: type 1 (pressure-based hybrid)' 
    ELSEIF (ivctype == 2 ) THEN
      write(nuspecif, '(A)' )                                             &
             '0     Vertical coordinate: type 2 (Gal-Chen hybrid)' 
    ELSEIF (ivctype == 3 ) THEN
      write(nuspecif, '(A)' )                                             &
             '0     Vertical coordinate: type 3 (SLEVE coordinate)'
    ENDIF

    WRITE(nuspecif,'(A)')         ' '
    WRITE(nuspecif,'(A)')   '    k    VCOORD(k)      Z(k)        P0(k) '
    DO k = 1, ke+1
      write(nuspecif,'(1X,I5,F12.4,F12.1,F12.4)' ) &
                  k, vcoord(k), hhlr(k), sigmr(k)*p0sl*0.01
    ENDDO

    WRITE(nuspecif,'(A)')         ' '
    WRITE(nuspecif,'(A,F10.4)')                                            &
      '   Boundary for change from z to terrain-following: vcflat = ', vcflat
    WRITE(nuspecif,'(A,I3)')                                            &
      '   Half-level index where levels become flat: k = ' , kflat
    WRITE(nuspecif,'(A)')         ' '
    IF (irefatm == 1) THEN
      WRITE(nuspecif,'(A38)')       ' Linear vertical gradient w.r.t. ln p '
      WRITE(nuspecif,'(A)')         ' '
      WRITE(nuspecif,'(A37)')       ' Values of the Reference Atmosphere: '
      WRITE(nuspecif,'(A37,F15.4)')                                            &
       '   Pressure at Sea-Level:            ', p0sl * 0.01_ireals
      WRITE(nuspecif,'(A37,F15.4)')                                            &
       '   Temperature at Sea-Level:         ', t0sl-273.15_ireals
      WRITE(nuspecif,'(A37,F15.4)')                                            &
       '   Temperature-Rate of decrease:     ', dt0lp
    ELSE IF (irefatm == 2) THEN
      WRITE(nuspecif,'(A61)')       ' Exponential profile with asymptotic isothermal stratosphere '
      WRITE(nuspecif,'(A)')         ' '
      WRITE(nuspecif,'(A37)')       ' Values of the Reference Atmosphere: '
      WRITE(nuspecif,'(A37,F15.4)')                                            &
       '   Pressure at Sea-Level:            ', p0sl * 0.01_ireals
      WRITE(nuspecif,'(A37,F15.4)')                                            &
       '   Temperature at Sea-Level:         ', t0sl-273.15_ireals
      WRITE(nuspecif,'(A37,F15.4)')                                            &
       '   Temp. diff. sea-level-stratosph.: ', delta_t
      WRITE(nuspecif,'(A37,F15.4)')                                            &
       '  Scale height for temp. decrease:   ', h_scal
    ELSE IF (irefatm == 3) THEN
      WRITE(nuspecif,'(A61)')       ' Reference atmosphere with constant Brunt-Vaisala frequency'
      WRITE(nuspecif,'(A)')         ' '
      WRITE(nuspecif,'(A37)')       ' Values of the Reference Atmosphere: '
      WRITE(nuspecif,'(A37,F15.4)')                                            &
       '   Pressure at Sea-Level:            ', p0sl * 0.01_ireals
      WRITE(nuspecif,'(A37,F15.4)')                                            &
       '   Temperature at Sea-Level:         ', t0sl-273.15_ireals
      WRITE(nuspecif,'(A37,F15.4)')                                            &
       '   Brunt-Vaisala-frequency           ', bvref
    ENDIF
    WRITE(nuspecif,'(A)')         ' '
    IF (ivctype == 3) THEN
       WRITE(nuspecif,'(A,F10.4)')                                            &
         '   Decay Rate for Large-Scale Topography: svc1 = ', svc1
       WRITE(nuspecif,'(A,F10.4)')                                            &
         '   Decay Rate for Small-Scale Topography: svc2 = ', svc2
    ENDIF

    WRITE(nuspecif,'(A)')                                                    &
           ' Full-level index corresponding to specific pressure levels:'
    WRITE(nuspecif,'(A,I3)') '   300 hPa: k = ', klv300
    WRITE(nuspecif,'(A,I3)') '   400 hPa: k = ', klv400
    WRITE(nuspecif,'(A,I3)') '   500 hPa: k = ', klv500
    WRITE(nuspecif,'(A,I3)') '   800 hPa: k = ', klv800
    WRITE(nuspecif,'(A,I3)') '   850 hPa: k = ', klv850

  ENDIF

END SUBROUTINE print_vertcoord

!==============================================================================
!+ Internal procedure in "organize_data" for the input of NAMELIST ioctl
!------------------------------------------------------------------------------

SUBROUTINE input_ioctl (nuspecif, nuin, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the NAMELIST-group ioctl.
!   The group contains variables for the organization of the input and output.
!
! Method:
!
!------------------------------------------------------------------------------

! Parameters
  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    nuspecif,     & ! Unit number for protocolling the task
    nuin            ! Unit number for Namelist INPUT file

  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable

! Variables for default values

  INTEGER (KIND=iintegers)   ::       &
    nincwait_d,    & ! seconds to wait until next attempt if a ready file is
                     ! not available
    nmaxwait_d,    & ! maximum seconds to wait until abort if a ready file is
                     ! not available
    nsma_stat_d,   & ! status for soil humidity analysis
    nvers_d,       & ! number of experiment for documentation
    ngribout_d,    & ! number of GRIBOUT namelists
    ncenter_d,     & ! originating center identification
    num_gribtabs_d,& ! number of GRIB tables used
    lst_gribtabs_d(max_gribtabs), & ! list of different GRIB tables
    ncglob_realization_d, & ! nr. of realization of the experiment
    nhour_restart_d(3) ! for start/stop/increment of writing restart files

  LOGICAL                    ::       &
    ldwd_grib_use_d,&! use some DWD specific Grib settings
    l_ke_in_gds_d, & ! explicit GDS entry for number of model levels
    lasync_io_d,   & ! if .TRUE., running with extra PEs for asynchronous IO
    lbdclim_d        ! if .TRUE., reading additional boundary data for 
                     ! climate mode

  REAL (KIND=ireals)         ::       &
    hrestart  (3), & ! tripe for start/stop/increment of writing restart files
    hrestart_d(3)    ! tripe for start/stop/increment of writing restart files

  CHARACTER (LEN=  4) ::  &
    yform_read_d     ! format of the (read) files

  CHARACTER (LEN=200)  ::    &
    yncglob_institution_d,   & ! originating center name
    yncglob_title_d,         & ! title string for the output
    yncglob_source_d,        & ! program name and version
    yncglob_project_id_d,    & ! identification of the project of simulation
    yncglob_experiment_id_d, & ! identification of the experiment of simulation
    yncglob_contact_d,       & ! contact e.g. email address
    yncglob_references_d       ! URL, report etc.

  CHARACTER (LEN=250)  ::    &
    ydir_restart_d             ! directory for restart files

  CHARACTER (LEN=  1) :: ytunit_restart_d ! for restart files

  CHARACTER (LEN=100) :: ytrans_in_d  ! directory for reading ready-files
  CHARACTER (LEN=100) :: ytrans_out_d ! directory for writing ready-files
  CHARACTER (LEN=  3) :: ymode_read_d ! mode for opening the (read) Grib files
  CHARACTER (LEN=  3) :: ymode_write_d! mode for opening the (write) Grib files

  INTEGER (KIND=iintegers)   :: ierr, iz_err, nc

! Define the namelist group
  NAMELIST /ioctl/ lasync_io, nsma_stat, ytrans_in, ytrans_out,            &
                   nincwait, nmaxwait, nvers, ymode_read, ymode_write,     &
                   ngribout, ncenter, num_gribtabs, lst_gribtabs,          &
                   nhour_restart, lbdclim, yform_read,                     &
                   yncglob_institution, yncglob_title, yncglob_source,     &
                   yncglob_project_id, yncglob_experiment_id,              &
                   yncglob_contact, yncglob_references, ncglob_realization,&
                   ydir_restart, ytunit_restart, ldwd_grib_use, l_ke_in_gds


!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE input_ioctl
!-------------------------------------------------------------------------------

iz_err = 0_iintegers

IF (my_world_id == 0) THEN

!-------------------------------------------------------------------------------
!- Section 1: Initialize the default variables
!-------------------------------------------------------------------------------

  nsma_stat_d    = 0_iintegers
  nincwait_d     = 0_iintegers
  nmaxwait_d     = 0_iintegers
  nvers_d        = 1_iintegers
  ytrans_in_d    = '   '
  ytrans_out_d   = '   '
  ymode_read_d   = 'r  '
  ymode_write_d  = 'w  '
  lasync_io_d    = .FALSE.
  ldwd_grib_use_d= .TRUE.
  l_ke_in_gds_d  = .FALSE.
  ngribout_d     =  1_iintegers
  ncenter_d      = 78_iintegers
  num_gribtabs_d = 11
  lst_gribtabs_d(1:max_gribtabs) = (/2, 201, 202, 203, 204, 205, 241, 242, 243, 244, &
                                   250,   0,   0,   0,   0,   0,   0,   0,   0,   0/)

  nhour_restart_d(1)      = 12_iintegers
  nhour_restart_d(2)      =  0_iintegers  ! this means: no restart at all
  nhour_restart_d(3)      = 12_iintegers
  ydir_restart_d          = ''
  ytunit_restart_d        = 'f'
  yform_read_d            = 'grb1'
  lbdclim_d               = .FALSE.


  yncglob_institution_d   = '-'
  yncglob_title_d         = '-'
  yncglob_source_d        = '-'
  yncglob_project_id_d    = '-'
  yncglob_experiment_id_d = '-'
  yncglob_contact_d       = '-'
  yncglob_references_d    = '-'
  ncglob_realization_d    = 1

!-------------------------------------------------------------------------------
!- Section 2: Initialize variables with defaults
!-------------------------------------------------------------------------------

  nsma_stat     = nsma_stat_d
  nincwait      = nincwait_d
  nmaxwait      = nmaxwait_d
  nvers         = nvers_d
  ytrans_in     = ytrans_in_d
  ytrans_out    = ytrans_out_d
  ymode_read    = ymode_read_d
  ymode_write   = ymode_write_d
  lasync_io     = lasync_io_d
  ldwd_grib_use = ldwd_grib_use_d
  l_ke_in_gds   = l_ke_in_gds_d
  ngribout      = ngribout_d
  ncenter       = ncenter_d
  num_gribtabs  = num_gribtabs_d
  lst_gribtabs  = lst_gribtabs_d

  nhour_restart(:)        = nhour_restart_d(:)
  ydir_restart            = ydir_restart_d
  ytunit_restart          = ytunit_restart_d
  yform_read              = yform_read_d
  lbdclim                 = lbdclim_d

  yncglob_institution     = yncglob_institution_d
  yncglob_title           = yncglob_title_d
  yncglob_source          = yncglob_source_d
  yncglob_project_id      = yncglob_project_id_d
  yncglob_experiment_id   = yncglob_experiment_id_d
  yncglob_contact         = yncglob_contact_d
  yncglob_references      = yncglob_references_d
  ncglob_realization      = ncglob_realization_d

!-------------------------------------------------------------------------------
!- Section 3: Input of the namelist values
!-------------------------------------------------------------------------------

  READ (nuin, ioctl, IOSTAT=iz_err)
ENDIF

IF (nproc > 1) THEN
  ! distribute error status to all processors
  CALL distribute_values  (iz_err, 1, 0, imp_integers,  icomm_world, ierr)
ENDIF

IF (iz_err /= 0) THEN
  ierrstat = -1
  RETURN
ENDIF

IF (my_world_id == 0) THEN

!-------------------------------------------------------------------------------
!- Section 4: Check values for errors and consistency
!-------------------------------------------------------------------------------

  ! Check whether nprocio and lasync_io do fit
  IF ( (nprocio /= 0) .AND. (lasync_io .EQV. .FALSE.) ) THEN
    PRINT *,' WARNING  *** For non-asynchronous IO nprocio = 0 is needed ***'
    nprocio   = 0
    lasync_io = .FALSE.
  ENDIF
  IF ( (nprocio == 0) .AND. (lasync_io .EQV. .TRUE.) ) THEN
    PRINT *,' WARNING  *** For nprocio == 0 synchronous IO is needed ***'
    lasync_io = .FALSE.
  ENDIF

  ! Check nincwait and nmaxwait
  IF (nincwait < 0) THEN
    PRINT *, ' ERROR    *** nincwait < 0:  ', nincwait, '  *** '
    ierrstat = 1002
  ENDIF
  IF (nmaxwait < 0) THEN
    PRINT *, ' ERROR    *** nmaxwait < 0:  ', nmaxwait, '  *** '
    ierrstat = 1002
  ENDIF

  ! Check nsma_stat
  IF (nsma_stat < 0) THEN
    PRINT *, ' ERROR    *** nsma_stat < 0:  ', nsma_stat, '  *** '
    ierrstat = 1002
  ENDIF

  ! Check ngribout
  IF (ngribout < 0) THEN
    PRINT *, ' ERROR    *** ngribout < 0:  ', ngribout, '  *** '
    ierrstat = 1002
  ENDIF

  ! Restart files
  IF (nhour_restart(1) <= NINT (nstart*dt/3600.0_ireals) ) THEN
    DO WHILE (nhour_restart(1) <= NINT (nstart*dt/3600.0_ireals) )
      nhour_restart(1) = nhour_restart(1) + nhour_restart(3)
    ENDDO
    PRINT *, ' WARNING  *** restart possible only after ', nstart , ' steps ***'
    PRINT *, '          *** nhour_restart(1) set to ', nhour_restart(1)
  ENDIF

#ifdef NUDGING
  ! Check whether restart is started only after nudging ends
  IF (nhour_restart(2) > nhour_restart(1)) THEN
    IF (luseobs .AND. NINT (nudgend*dt/3600.0_ireals) >= nhour_restart(1)) THEN
      PRINT *, ' ERROR  *** restart files can only be written after nudging ***'
      ierrstat = 1002
    ENDIF
  ENDIF
#endif

  ! set ydir_restart
  nc = LEN_TRIM(ydir_restart)
  IF (nc > 0) THEN
    IF(ydir_restart /= '' .AND. ydir_restart(nc:nc)/='/') THEN
       ydir_restart(1:nc+1) = ydir_restart(1:nc)//'/'
    ENDIF
  ENDIF

  SELECT CASE (ytunit_restart)
  CASE('f','t','c','d')
    nc = 1
    ! that's o.k.
  CASE DEFAULT
    PRINT *, ytunit_restart,' is not a known parameter for YTUNIT_RESTART'
    ierrstat = 1002
  END SELECT

  ! Check format for reading/writing files
  IF ((yform_read /= 'grb1') .AND. (yform_read /= 'ncdf')) THEN
    PRINT *, ' ERROR  *** wrong format for input files *** ', yform_read
    ierrstat = 1002
  ENDIF

#ifndef GRIBDWD
  IF (yform_read == 'grb1') THEN
    PRINT *, ' ERROR  *** yform_read = grb1, but model is not compiled to use DWD GRIB library ***'
    ierrstat = 1002
  ENDIF
#endif

#ifndef NETCDF
  IF (yform_read == 'ncdf') THEN
    PRINT *, ' ERROR  *** yform_read = ncdf, but model is not compiled to use NetCDF ***'
    ierrstat = 1002
  ENDIF
#endif

ENDIF

!-------------------------------------------------------------------------------
!- Section 5: Distribute variables to all nodes
!-------------------------------------------------------------------------------

IF (nproc > 1) THEN

  IF (my_world_id == 0) THEN
    intbuf ( 1) = nincwait
    intbuf ( 2) = nmaxwait
    intbuf ( 3) = nsma_stat
    intbuf ( 4) = nvers  
    intbuf ( 5) = ncenter
    intbuf ( 6) = ngribout
    intbuf ( 7) = ncglob_realization
    intbuf ( 8) = nhour_restart(1)
    intbuf ( 9) = nhour_restart(2)
    intbuf (10) = nhour_restart(3)
    intbuf (11) = num_gribtabs
    intbuf (12:12+max_gribtabs-1) = lst_gribtabs(1:max_gribtabs)
    logbuf ( 2) = lasync_io
    logbuf ( 3) = lbdclim
    logbuf ( 4) = ldwd_grib_use
    logbuf ( 5) = l_ke_in_gds
    charbuf( 1) = ytrans_in
    charbuf( 2) = ytrans_out
    charbuf( 3) = ymode_read
    charbuf( 4) = ymode_write
    charbuf( 5) = yform_read
    charbuf( 7) = yncglob_institution
    charbuf( 8) = yncglob_title
    charbuf( 9) = yncglob_source
    charbuf(10) = yncglob_project_id
    charbuf(11) = yncglob_experiment_id
    charbuf(12) = yncglob_contact
    charbuf(13) = yncglob_references
    charbuf(14) = ydir_restart
    charbuf(15) = ytunit_restart
  ENDIF

  CALL distribute_values (intbuf,  11+max_gribtabs, 0, imp_integers,         &
                                                        icomm_world, ierr)
  CALL distribute_values (logbuf,  5, 0, imp_logical,   icomm_world, ierr)
  CALL distribute_values (charbuf,15, 0, imp_character, icomm_world, ierr)

  IF (my_world_id /= 0) THEN
    nincwait              = intbuf ( 1)
    nmaxwait              = intbuf ( 2)
    nsma_stat             = intbuf ( 3)
    nvers                 = intbuf ( 4)
    ncenter               = intbuf ( 5)
    ngribout              = intbuf ( 6)
    ncglob_realization    = intbuf ( 7)
    nhour_restart(1)      = intbuf ( 8)
    nhour_restart(2)      = intbuf ( 9)
    nhour_restart(3)      = intbuf (10)
    num_gribtabs          = intbuf (11)
    lst_gribtabs(1:max_gribtabs) = intbuf (12:12+max_gribtabs-1)
    lasync_io             = logbuf ( 2)
    lbdclim               = logbuf ( 3)
    ldwd_grib_use         = logbuf ( 4)
    l_ke_in_gds           = logbuf ( 5)
    ytrans_in             = charbuf( 1)
    ytrans_out            = charbuf( 2)
    ymode_read            = charbuf( 3)
    ymode_write           = charbuf( 4)
    yform_read            = charbuf( 5)
    yncglob_institution   = charbuf( 7)
    yncglob_title         = charbuf( 8)
    yncglob_source        = charbuf( 9)
    yncglob_project_id    = charbuf(10)
    yncglob_experiment_id = charbuf(11)
    yncglob_contact       = charbuf(12)
    yncglob_references    = charbuf(13)
    ydir_restart          = charbuf(14)
    ytunit_restart        = charbuf(15)
  ENDIF

ENDIF

! Set value num_io: for asnychronous IO num_io = nprocio,
!                   for  snychronous IO num_io = 1
IF (lasync_io) THEN
  num_io = nprocio
ELSE
  num_io = 1
ENDIF

!-------------------------------------------------------------------------------
!- Section 6: Output of the namelist variables and their default values
!-------------------------------------------------------------------------------

IF (my_world_id == 0) THEN

  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(A28)') '0     NAMELIST:  ioctl'
  WRITE (nuspecif, '(A28)') '      ----------------'
  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(T7,A,T21,A,T39,A,T58,A)') 'Variable', 'Actual Value',   &
                                               'Default Value', 'Format'

  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                 'nincwait   ',nincwait   ,nincwait_d, ' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                 'nmaxwait   ',nmaxwait   ,nmaxwait_d, ' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                'nsma_stat  ',nsma_stat  ,nsma_stat_d, ' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                'nvers      ',nvers      ,nvers_d    , ' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                         'ncenter', ncenter, ncenter_d,' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                                      'ngribout', ngribout, ngribout_d,' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                'start restart',    nhour_restart(1), nhour_restart(1),' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                'stop  restart',    nhour_restart(2), nhour_restart(2),' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                'incr. restart',    nhour_restart(3), nhour_restart(3),' I '
  WRITE (nuspecif, '(T8,A,T21,  A  ,T40,A    ,T59,A)')                       &
                   'ytunit_restart',ytunit_restart, ytunit_restart_d  ,'C*1'
  WRITE (nuspecif, '(T8,A,T30,A)')                                           &
                         'ydir_restart',          TRIM(ydir_restart)       
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                          'num_gribtabs', num_gribtabs, num_gribtabs_d,' I '

  DO nc = 1, max_gribtabs
    IF (lst_gribtabs(nc) /= 0) THEN
      WRITE (nuspecif, '(T8,A,I2,A,T21,I12  ,T40,I12  ,T59,A3)')             &
          'lst_gribtabs',nc,': ', lst_gribtabs(nc), lst_gribtabs_d(nc),' I '
    ENDIF
  ENDDO

  WRITE (nuspecif, '(T8,A,T21,  A  ,T40,  A  ,T59,A3)')                      &
                                'ytrans_in ' ,ytrans_in  ,ytrans_in_d, 'C*8'
  WRITE (nuspecif, '(T8,A,T21,  A  ,T40,  A  ,T59,A3)')                      &
                              'ytrans_out ' ,ytrans_out ,ytrans_out_d, 'C*8'
  WRITE (nuspecif, '(T8,A,T21,  A  ,T40,A    ,T59,A)')                       &
                              'ymode_read ',ymode_read, ymode_read_d  ,'C*3'
  WRITE (nuspecif, '(T8,A,T21,  A  ,T40,A    ,T59,A)')                       &
                              'ymode_write',ymode_write, ymode_write_d,'C*3'
  WRITE (nuspecif, '(T8,A,T21,  A  ,T40,A    ,T59,A)')                       &
                              'yform_read ',yform_read, yform_read_d  ,'C*4'

  WRITE (nuspecif, '(T8,A,T30,A)')                                           &
                         'yncglob_institution',   TRIM(yncglob_institution)
  WRITE (nuspecif, '(T8,A,T30,A)')                                           &
                         'yncglob_title',         TRIM(yncglob_title)
  WRITE (nuspecif, '(T8,A,T30,A)')                                           &
                         'yncglob_source',        TRIM(yncglob_source)
  WRITE (nuspecif, '(T8,A,T30,A)')                                           &
                         'yncglob_project_id',    TRIM(yncglob_project_id)
  WRITE (nuspecif, '(T8,A,T30,A)')                                           &
                         'yncglob_experiment_id', TRIM(yncglob_experiment_id)
  WRITE (nuspecif, '(T8,A,T30,A)')                                           &
                         'yncglob_contact',       TRIM(yncglob_contact)
  WRITE (nuspecif, '(T8,A,T30,A)')                                           &
                         'yncglob_references',    TRIM(yncglob_references)
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
            'ncglob_realization',ncglob_realization, ncglob_realization,' I '

  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                            'lasync_io',   lasync_io,   lasync_io_d,   ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                'ldwd_grib_use', ldwd_grib_use, ldwd_grib_use_d,       ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                'l_ke_in_gds', l_ke_in_gds, l_ke_in_gds_d,             ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A3)')                      &
                            'lbdclim',     lbdclim,   lbdclim_d,       ' L '
  WRITE (nuspecif, '(A2)')  '  '

ENDIF

!-------------------------------------------------------------------------------
!- End of the Subroutine
!-------------------------------------------------------------------------------

END SUBROUTINE input_ioctl

!==============================================================================
!+ Module procedure in "organize_data" for the input of NAMELIST dbase
!------------------------------------------------------------------------------

SUBROUTINE input_dbase  ( nuspecif, nuin, ierrstat )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the NAMELIST-group database.
!   The group gribin contains variables for controlling the input /output of
!   GRIB-fields from/to the database
!   Unlike the other routines that read a NAMELIST-group, this routine just
!   calls a routine of MPE_IO for reading the group dbase, because these
!   values must only be known by the IO-PEs.
!
! Method:
!   Directly call module procedure mpe_readnldb of => MPE_IO
!
!==============================================================================

! Parameters
  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    nuspecif,     & ! Unit number for protocolling the task
    nuin            ! Unit number for Namelist INPUT file

  INTEGER   (KIND=iintegers),   INTENT (OUT)     ::        &
    ierrstat        ! error status variable

  INTEGER   (KIND=iintegers)                     ::        &
    ierr

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE input_dbase
!------------------------------------------------------------------------------

  IF (my_world_id == 0) THEN
    CALL mpe_readnldb(nuin, nuspecif, ierrstat)
  ENDIF

  IF (nproc > 1) THEN
    ! distribute error status to all processors
    CALL distribute_values  (ierrstat, 1, 0, imp_integers,  icomm_world, ierr)
  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE input_dbase

!==============================================================================
!+ Internal procedure in "organize_data" for the input of NAMELIST gribin
!------------------------------------------------------------------------------

SUBROUTINE input_gribin (nuspecif, nuin, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the NAMELIST-group gribin.
!   The group gribin contains variables for controlling the input of
!   GRIB-fields for initial and boundary data.
!
! Method:
!   All variables are initialized with default values and then read in from
!   the file INPUT. The input values are checked for errors and for
!   consistency. If wrong input values are detected the program prints
!   an error message. The program is not stopped in this routine but an
!   error code is returned to the calling routine that aborts the program after
!   reading in all other namelists.
!   In parallel mode, the variables are distributed to all nodes with the
!   environment-routine distribute_values.   
!   Both, default and input values are written to the file YUSPECIF
!   (specification of the run).
!
!------------------------------------------------------------------------------

! Parameters
  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    nuspecif,     & ! Unit number for protocolling the task
    nuin            ! Unit number for Namelist INPUT file

  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable

! Local variables

! Variables for default values
  CHARACTER (LEN=250) :: ydirini_d     ! directory with initial data file
  CHARACTER (LEN=250) :: ydirbd_d      ! directory with boundary data file
  CHARACTER (LEN=  1) :: ytunitbd_d    ! Unit of time for the boundary data

  CHARACTER (LEN= 10)           ::             &
    yvarini_d(nzmxin), & !  names of the initial variables for input
    yvarbd_d (nzmxin)    !  names of the boundary variables for input

  INTEGER (KIND=iintegers)   ::       &
    nzmxini_d,    & ! number of input variables in the default list
    nzmxbd_d,     & ! number of input variables in the default list
    nlgw_ini_d,   & ! number of prognostic soil water levels in initial data
    nlgw_bd_d,    & ! number of prognostic soil water levels in boundary data
    newbcdt_d,    & ! time step increment of boundary update from analyses
    newbc_d,      & ! number of times that boundary update is analysis after 1 h
    nincboufac_d, & ! factor to 'nincbound' when new boundary update is analysis
    npstrframe_d, & ! width (number of points) of the strip around
                    ! the b.d. frame
    ilevbotnoframe_d! bottom model level with b.d. defined on the whole grid
                    ! (model levels below ilevbotnoframe are defined on a frame)

  REAL (KIND=ireals)         ::       &
    hincbound_d,  & ! hour increment of boundary update - default value
    hnewbcdt_d,   & ! hour increment of boundary update from analyses - default
    hnewbcdt,     & ! hour increment of boundary update from analyses
! UB>>
    htest
! UB<<

  LOGICAL                    ::       &
    lchkini_d,    & ! checking the initial data
    lchkbd_d,     & ! checking the boundary data
    lan_t_s_d,    & ! switch for selection of analysed (tri=0) sst, t_s
    lan_t_s,      &
    lan_t_so0_d,  & ! switch for selection of analysed (tri=0) sst, t_so(0)
    lan_t_so0,    &
    lan_t_snow_d, & ! switch for selection of analysed temp. of snow, t_surf
    lan_t_snow,   &
    lan_t_cl_d,   & ! switch for selection of anal. climate soil temp., t_cl
    lan_t_cl,     &
    lan_w_snow_d, & ! switch for selection of anal. water cont. of snow, w_snow
    lan_w_snow,   &
    lan_rho_snow_d,&! switch for selection of anal. snow density, rho_snow
    lan_rho_snow, &
    lan_w_i_d,    & ! switch for selection of anal. interception water, w_i
    lan_w_i,      &
    lan_w_cl_d,   & ! switch for selection of anal. clim. soil water cont., w_cl
    lan_w_cl,     &
    lan_vio3_d,   & ! switch for selection of anal. vert. integr. o3, vio3
    lan_vio3,     &
    lan_hmo3_d,   & ! switch for selection of anal. o3 maximum, hmo3
    lan_hmo3,     &
    lan_plcov_d,  & ! switch for selection of anal. plant cover, plcov
    lan_plcov,    &
    lan_lai_d,    & ! switch for selection of anal. index of plants, lai
    lan_lai,      &
    lan_rootdp_d, & ! switch for selection of anal. depth of roots, rootdp
    lan_rootdp

  LOGICAL                    ::       &
    lbdana_d,     & ! if boundary data are from analyses
    lana_qi_d,    & ! if .TRUE., take qi-values from analysis file
                    ! else, qi is set in the model
    llb_qi_d,     & ! if .TRUE., take qi_bd-values from lateral boundaries file
                    ! else, qi_bd is set in the model
    lana_qr_qs_d, & ! if .TRUE., take qr- and qs-values from analysis file
                    ! else, qr and qs are set in the model
    llb_qr_qs_d,  & ! if .TRUE., take qr_bd- and qs_bd-values from lateral
                    ! bound. file else, qr_bd and qs_bd are set in the model
    lana_qg_d,    & ! if .TRUE., take qg-values from analysis file
                    ! else, qg is set in the model
    llb_qg_d,     & ! if .TRUE., take qg_bd-values from lateral boundaries file
                    ! else, qg_bd is set in the model
    lana_rho_snow_d,&!if .TRUE., take rho-snow-values from analysis file
                    ! else, it is set in the model
    lbd_frame_d     ! if .TRUE., boundary data are on a frame

! Other Variables
  INTEGER (KIND=iintegers)   ::       &
    nzylen, ierr, invar, i, n2, n, nvartmp, iz_err

  LOGICAL                    ::       &
    lzcontained

! Define the namelist group
  NAMELIST /gribin/             ydirini   , lchkini   , nlgw_ini  ,        &
                                ydirbd    , lchkbd    , nlgw_bd   ,        &
                    ytunitbd  , lbdana    ,             hincbound ,        &
                    newbc     , newbcdt   , hnewbcdt  , nincboufac,        &
                    lan_t_s   , lan_t_so0 , lan_t_snow, lan_t_cl  ,        &
                    lan_w_snow, lan_w_i   , lan_w_cl  , lan_rootdp,        &
                    lan_vio3  , lan_hmo3  , lan_plcov , lan_lai   ,        &
                    lan_rho_snow,lbd_frame, lana_qi   , llb_qi    ,        &
                    lana_qr_qs, llb_qr_qs , lana_qg   , llb_qg    ,        &
                    lana_rho_snow, npstrframe, ilevbotnoframe

  ! NOTE: yvarini, yvarbd have been eliminated in Version 4.12

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE input_gribin
!------------------------------------------------------------------------------

iz_err = 0_iintegers

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 1: Initialize the default variables
!------------------------------------------------------------------------------

! Variables for initial data
  yvarini_d(:)   = '          '
  yvarini_d( 1)  = 'HSURF     '; yvarini_d( 2)  = 'FR_LAND   '
  yvarini_d( 3)  = 'Z0        '; yvarini_d( 4)  = 'SOILTYP   '
  yvarini_d( 5)  = 'PLCOV     '; yvarini_d( 6)  = 'LAI       '
  yvarini_d( 7)  = 'ROOTDP    '; yvarini_d( 8)  = 'VIO3      '
  yvarini_d( 9)  = 'HMO3      '; yvarini_d(10)  = 'U         '
  yvarini_d(11)  = 'V         '; yvarini_d(12)  = 'W         '
  yvarini_d(13)  = 'T         '; yvarini_d(14)  = 'QV        '
  yvarini_d(15)  = 'QC        '; yvarini_d(16)  = 'PP        '
  yvarini_d(17)  = 'T_SNOW    '; yvarini_d(18)  = 'W_I       '
  yvarini_d(19)  = 'QV_S      '; yvarini_d(20)  = 'W_SNOW    '
  yvarini_d(21)  = 'T_S       '

  IF (lmulti_layer) THEN
    yvarini_d(22)  = 'T_SO      '
    yvarini_d(23)  = 'W_SO      '
    yvarini_d(24)  = 'FRESHSNW  '
    nzmxini_d      = 24
  ELSE
    yvarini_d(22)  = 'T_M       '
    yvarini_d(23)  = 'T_CL      '
    yvarini_d(24)  = 'W_G1      '
    yvarini_d(25)  = 'W_G2      '
    yvarini_d(26)  = 'W_CL      '
    nzmxini_d      = 26
  ENDIF

  IF (lforest)THEN
    yvarini_d(nzmxini_d + 1)  = 'FOR_E     '
    yvarini_d(nzmxini_d + 2)  = 'FOR_D     '
    nzmxini_d      = nzmxini_d + 2
  ENDIF

  IF (lradtopo) THEN
    yvarini_d(nzmxini_d + 1)  = 'SKYVIEW   '
    yvarini_d(nzmxini_d + 2)  = 'SLO_ASP   '
    yvarini_d(nzmxini_d + 3)  = 'SLO_ANG   '
    yvarini_d(nzmxini_d + 4)  = 'HORIZON   '
    nzmxini_d      = nzmxini_d + 4
  ENDIF

  IF (lsso) THEN
    yvarini_d(nzmxini_d + 1)  = 'SSO_STDH  '
    yvarini_d(nzmxini_d + 2)  = 'SSO_GAMMA '
    yvarini_d(nzmxini_d + 3)  = 'SSO_THETA '
    yvarini_d(nzmxini_d + 4)  = 'SSO_SIGMA '
    nzmxini_d      = nzmxini_d + 4
  ENDIF

  IF (lemiss) THEN 
     yvarini_d(nzmxini_d + 1)  = 'EMIS_RAD  '
     nzmxini_d      = nzmxini_d + 1
  END IF
  
  IF (lstomata) THEN
     yvarini_d(nzmxini_d + 1)  = 'RSMIN     '
     nzmxini_d      = nzmxini_d + 1
  END IF
  
  IF ( itype_aerosol == 2 ) THEN
    yvarini_d(nzmxini_d + 1)  = 'AER_SO4   '
    yvarini_d(nzmxini_d + 2)  = 'AER_DUST  '
    yvarini_d(nzmxini_d + 3)  = 'AER_ORG   '
    yvarini_d(nzmxini_d + 4)  = 'AER_BC    '
    yvarini_d(nzmxini_d + 5)  = 'AER_SS    '
      nzmxini_d      = nzmxini_d + 5
  ENDIF

  ! Variables for the sea ice and/or the FLake model
  IF (lseaice .OR. llake) THEN
    yvarini_d(nzmxini_d + 1)  = 'T_ICE     '
    yvarini_d(nzmxini_d + 2)  = 'H_ICE     '
    nzmxini_d      = nzmxini_d + 2
  ENDIF

  IF (llake) THEN
    ! Warm start of FLake in LM.
    ! Values of FLake prognostic variables are read in from the INPUT file.
    yvarini_d(nzmxini_d + 1)  = 'FR_LAKE   '
    yvarini_d(nzmxini_d + 2)  = 'DEPTH_LK  '
    yvarini_d(nzmxini_d + 3)  = 'T_MNW_LK  '
    yvarini_d(nzmxini_d + 4)  = 'T_WML_LK  '
    yvarini_d(nzmxini_d + 5)  = 'T_BOT_LK  '
    yvarini_d(nzmxini_d + 6)  = 'C_T_LK    '
    yvarini_d(nzmxini_d + 7)  = 'H_ML_LK   '
    nzmxini_d      = nzmxini_d + 7
  ENDIF

  IF (lmulti_snow) THEN
    yvarini_d(nzmxini_d + 1)  = 'T_SNOW_MUL'
    yvarini_d(nzmxini_d + 2)  = 'DZH_SNOW  '
    yvarini_d(nzmxini_d + 3)  = 'WTOT_SNOW '
    yvarini_d(nzmxini_d + 4)  = 'WLIQ_SNOW '
    nzmxini_d      = nzmxini_d + 4
  ENDIF

  !_cdm The following FLake external parameters are assigned their default 
  !     values that are kept constant in space and time.
  !     These external-parameter fields are handled internally by LM
  !     and are not included into the LM IO lists
  !     (this should be done in the future).
  !       yvarini_d (nzmxini_d+NN)  = 'FETCH_LK  '
  !       yvarini_d (nzmxini_d+NN)  = 'DP_BS_LK  '
  !       yvarini_d (nzmxini_d+NN)  = 'T_BS_LK   '
  !       yvarini_d (nzmxini_d+NN)  = 'GAMSO_LK  '
  !     In the present configuration, 
  !     the bottom-sediment module of FLake is switched off.
  !     The respective fields are handled internally by LM
  !     and are not included into the LM IO lists.
  !       yvarini_d (nzmxini_d+NN)  = 'T_B1_LK   '
  !       yvarini_d (nzmxini_d+NN)  = 'H_B1_LK   '

#ifdef POLLEN
  IF (l_pollen) THEN
    IF (lpollenini) THEN
      yvarini_d(nzmxini_d +  1)  = 'CNC_BETU'
      yvarini_d(nzmxini_d +  2)  = 'QP0_BETU'
      yvarini_d(nzmxini_d +  3)  = 'TSS_BETU'
      yvarini_d(nzmxini_d +  4)  = 'TSE_BETU'
      nzmxini_d      = nzmxini_d + 4
    ENDIF
    yvarini_d(nzmxini_d +  1)  = 'FR_BETU'
    yvarini_d(nzmxini_d +  2)  = 'TTS_BETU'
    yvarini_d(nzmxini_d +  3)  = 'TTE_BETU'
    yvarini_d(nzmxini_d +  4)  = 'HCEM_BETU'
    yvarini_d(nzmxini_d +  5)  = 'SDES_BETU'
    nzmxini_d      = nzmxini_d + 5
  ENDIF
#endif

  ydirini_d      = ' '
  nlgw_ini_d     =   2
  hincbound_d    = 1.0_ireals
  lchkini_d      = .FALSE.
  lan_t_s_d      = .FALSE.
  lan_t_so0_d    = .FALSE.
  lan_t_snow_d   = .FALSE.
  lan_t_cl_d     = .FALSE.
  lan_w_snow_d   = .FALSE.
  lan_rho_snow_d = .FALSE.
  lan_w_i_d      = .FALSE.
  lan_w_cl_d     = .FALSE.
  lan_vio3_d     = .FALSE.
  lan_hmo3_d     = .FALSE.
  lan_plcov_d    = .FALSE.
  lan_lai_d      = .FALSE.
  lan_rootdp_d   = .FALSE.
  lbdana_d       = .FALSE.
  lana_qi_d      = .FALSE.
  lana_qr_qs_d   = .FALSE.
  lana_qg_d      = .FALSE.
  lana_rho_snow_d= .FALSE.

! Variables for boundary data
  yvarbd_d (:)   = '          '
  yvarbd_d ( 1)  = 'U         '; yvarbd_d ( 2)  = 'V         '
  yvarbd_d ( 3)  = 'T         '; yvarbd_d ( 4)  = 'QV        '
  yvarbd_d ( 5)  = 'QC        '; yvarbd_d ( 6)  = 'PP        '
  yvarbd_d ( 7)  = 'T_SNOW    '; yvarbd_d ( 8)  = 'W_SNOW    '
  yvarbd_d ( 9)  = 'QV_S      '
  nzmxbd_d       = 9

  IF (.NOT. lw_freeslip) THEN
    ! then we need W as a boundary variable
    yvarbd_d(nzmxbd_d +1)  = 'W         '
    nzmxbd_d               = nzmxbd_d + 1
  ENDIF

  IF (.NOT. lmulti_layer) THEN
    yvarbd_d(nzmxbd_d +1)  = 'T_S       '
    yvarbd_d(nzmxbd_d +2)  = 'T_M       '
    yvarbd_d(nzmxbd_d +3)  = 'W_G1      '
    yvarbd_d(nzmxbd_d +4)  = 'W_G2      '
    nzmxbd_d               = nzmxbd_d + 4
  ENDIF

  ! For climate runs, also most external parameters have to 
  ! be updated during the simulation
  IF (lbdclim) THEN
    yvarbd_d (nzmxbd_d +1)  = 'PLCOV   '
    yvarbd_d (nzmxbd_d +2)  = 'LAI     '
    yvarbd_d (nzmxbd_d +3)  = 'ROOTDP  '
    yvarbd_d (nzmxbd_d +4)  = 'VIO3    '
    yvarbd_d (nzmxbd_d +5)  = 'HMO3    '
    nzmxbd_d      = nzmxbd_d + 5
    IF (.NOT. lmulti_layer) THEN
      yvarbd_d (nzmxbd_d +1)  = 'T_CL    '
      yvarbd_d (nzmxbd_d +2)  = 'W_CL    '
      nzmxbd_d      = nzmxbd_d + 2
    ELSE
      yvarbd_d (nzmxbd_d +1)  = 'T_S     '
      nzmxbd_d      = nzmxbd_d + 1
    ENDIF
  ENDIF

#ifdef POLLEN
  IF (l_pollen) THEN
    IF (lpollenbd) THEN
      yvarbd_d(nzmxbd_d +  1)  = 'CNC_BETU '
      nzmxbd_d      = nzmxbd_d + 1
    ENDIF
  ENDIF
#endif

  ydirbd_d         = ' '
  ytunitbd_d       = 'f'
  nlgw_bd_d        =   2
  lchkbd_d         = .FALSE.
  llb_qi_d         = .FALSE.
  llb_qr_qs_d      = .FALSE.   ! in LMK preop binary this is .TRUE.
  llb_qg_d         = .FALSE.
  lbd_frame_d      = .FALSE.
  npstrframe_d     =   8
  ilevbotnoframe_d = 0
  newbcdt_d        = 540
  hnewbcdt_d       = 0.0_ireals
  newbc_d          = 0
  nincboufac_d     = 1

!------------------------------------------------------------------------------
!- Section 2: Initialize default variables and error status variable
!------------------------------------------------------------------------------

! Variables for initial data
  yvarini  (:)   = ''
  ydirini        = ydirini_d      
  nlgw_ini       = nlgw_ini_d   
  hincbound      = hincbound_d
  lchkini        = lchkini_d      
  lan_t_s        = lan_t_s_d
  lan_t_so0      = lan_t_so0_d
  lan_t_snow     = lan_t_snow_d
  lan_t_cl       = lan_t_cl_d
  lan_w_snow     = lan_w_snow_d
  lan_rho_snow   = lan_rho_snow_d
  lan_w_i        = lan_w_i_d
  lan_w_cl       = lan_w_cl_d
  lan_vio3       = lan_vio3_d
  lan_hmo3       = lan_hmo3_d
  lan_plcov      = lan_plcov_d
  lan_lai        = lan_lai_d
  lan_rootdp     = lan_rootdp_d
  lbdana         = lbdana_d
  lana_qi        = lana_qi_d
  lana_qr_qs     = lana_qr_qs_d
  lana_qg        = lana_qg_d
  lana_rho_snow  = lana_rho_snow_d

! Variables for boundary data
  yvarbd   (:)   = ''
  ydirbd         = ydirbd_d       
  ytunitbd       = ytunitbd_d
  nlgw_bd        = nlgw_bd_d    
  lchkbd         = lchkbd_d       
  llb_qi         = llb_qi_d
  llb_qr_qs      = llb_qr_qs_d
  llb_qg         = llb_qg_d
  lbd_frame      = lbd_frame_d
  npstrframe     = npstrframe_d
  ilevbotnoframe = ilevbotnoframe_d
  newbcdt        = newbcdt_d
  hnewbcdt       = hnewbcdt_d
  newbc          = newbc_d
  nincboufac     = nincboufac_d

!------------------------------------------------------------------------------
!- Section 3: Input of the namelist values
!------------------------------------------------------------------------------

  READ (nuin, gribin, IOSTAT=iz_err)
ENDIF

IF (nproc > 1) THEN
  ! distribute error status to all processors
  CALL distribute_values  (iz_err, 1, 0, imp_integers,  icomm_world, ierr)
ENDIF

IF (iz_err /= 0) THEN
  ierrstat = -1
  RETURN
ENDIF

IF (my_world_id == 0) THEN

!------------------------------------------------------------------------------
!- Section 4: Check values for errors and consistency
!------------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !- Section 4.1: Variables for initial data
  !----------------------------------------------------------------------------

  ! Adapt the default input lists (for initial and boundary fields) to the
  ! chosen parameters nlgw_ini/nlgw_bd.
  ! This is just for convenience that the user does not have to specify 
  ! corresponding input lists
  IF (nlgw_ini == 3) THEN
    yvarini_d(nzmxini_d + 1)  = 'W_G3      '
    nzmxini_d      = nzmxini_d + 1
  ENDIF
  IF (nlgw_bd  == 3) THEN
    yvarbd_d(nzmxbd_d + 1)    = 'W_G3      '
    nzmxbd_d       = nzmxbd_d  + 1
  ENDIF

  ! determination of the initial input variables
  invar = COUNT(yvarini(:) /= '')
  IF( invar == 0) THEN
    ! nothing has been read and the default list is used
    yvarini(1:nzmxini_d)  = yvarini_d(1:nzmxini_d)
    yvarini(nzmxini_d+1:) = ''
    nyvar_i             = nzmxini_d
  ELSEIF ( (invar >= 1) .AND. (invar <= nzmxin) ) THEN
    nyvar_i         = invar
  ELSEIF (invar > nzmxin) THEN
    PRINT *,' ERROR  *** number of variables for initial input is too big ***'
    ierrstat = 1002
  ELSE
    ! something went wrong: cannot determine number of variables for input
    PRINT *,' ERROR    *** cannot determine number of variables for input *** '
    ierrstat = 1002
  ENDIF

  ! Check the setting of lana_qi, llb_qi, if cloud ice is turned off
  IF (.NOT. lprog_qi) THEN
    ! no ice phase is computed, so no data should be read
    IF (lana_qi) THEN
      PRINT *,                                                          &
      ' WARNING: ** lana_qi is set .FALSE., because ice phase is turned off **'
      lana_qi = .FALSE.
    ENDIF
    IF (llb_qi) THEN
      PRINT *,                                                          &
      ' WARNING: ** llb_qi  is set .FALSE., because ice phase is turned off **'
      llb_qi = .FALSE.
    ENDIF
  ENDIF 
  
  IF (lana_qi) THEN
    ! check whether qi is contained in yvarini, else add it
    lzcontained = .FALSE.
    search_ini_qi: DO n=1, nyvar_i
      IF (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'QI') THEN
        lzcontained = .TRUE.
        EXIT search_ini_qi
      ENDIF
    ENDDO search_ini_qi
    IF (.NOT. lzcontained) THEN
      ! print a warning and add it
      IF (nyvar_i < nzmxin) THEN
        PRINT *,' WARNING: ** QI is added to initial state list **'
        PRINT *,'          ** because lana_qi is set .TRUE.     **'
        nyvar_i = nyvar_i + 1
        yvarini(nyvar_i) = 'QI        '
      ELSE
        PRINT *,' ERROR:   *** QI could not be added to initial state list ***'
        PRINT *,'          *** since this list would get too large         ***'
        ierrstat = 1002
      ENDIF
    ENDIF
  ELSE
    ! if qi is contained in varini, it must be deleted
    nvartmp = nyvar_i
    search_ini_qi_2: DO n=1, nvartmp
      IF (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'QI') THEN
        PRINT *, ' WARNING: ** QI deleted from initial state list **'
        PRINT *, '          ** because lana_qi is set .FALSE.     **'
        IF (n < nyvar_i) THEN
          DO n2 = n, nyvar_i-1
            yvarini(n2) = yvarini(n2+1)
          ENDDO
        ENDIF
        yvarini(nyvar_i) = ''
        nyvar_i = nyvar_i - 1
      ENDIF
      IF (n >= nyvar_i) EXIT search_ini_qi_2
    ENDDO search_ini_qi_2
  ENDIF

  ! Check the setting of lana_qr_qs and llb_qr_qs,
  ! if prognostic precip is turned off, ...
  IF (.NOT.lprogprec .OR. itype_gscp == 1) THEN
    ! no prognostic precipitation is computed, so no data should be read
    IF (lana_qr_qs) THEN
      PRINT *, ' WARNING: ** lana_qr_qs is set .FALSE., because '                &
             , '             prognostic precipitation is turned off **'          &
             , '             or only prog. Kessler scheme is used   **'
      lana_qr_qs = .FALSE.
    ENDIF
    IF (llb_qr_qs) THEN
      PRINT *, ' WARNING: ** llb_qr_qs is set .FALSE., because '                 &
             , '             prognostic precipitation is turned off **'          &
             , '             or only prog. Kessler scheme is used   **'
      llb_qr_qs = .FALSE.
    ENDIF
  ENDIF
 
  IF (lana_qr_qs) THEN
    ! check whether qr is contained in yvarini, else add it
    lzcontained = .FALSE.
    search_ini_qr: DO n=1, nyvar_i
      IF (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'QR') THEN
        lzcontained = .TRUE.
        EXIT search_ini_qr
      ENDIF
    ENDDO search_ini_qr
    IF (.NOT. lzcontained) THEN
      ! print a warning and add it
      IF (nyvar_i < nzmxin) THEN
        PRINT *,' WARNING: ** QR is added to initial state list **'
        PRINT *,'          ** because lana_qr_qs is set .TRUE.  **'
        nyvar_i = nyvar_i + 1
        yvarini(nyvar_i) = 'QR        '
      ELSE
        PRINT *,' ERROR:   *** QR could not be added to initial state list ***'
        PRINT *,'          *** since this list would get too large         ***'
        ierrstat = 1002
      ENDIF
    ENDIF
    ! check whether qs is contained in yvarini, else add it
    lzcontained = .FALSE.
    search_ini_qs: DO n=1, nyvar_i
      IF (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'QS') THEN
        lzcontained = .TRUE.
        EXIT search_ini_qs
      ENDIF
    ENDDO search_ini_qs
    IF (.NOT. lzcontained) THEN
      ! print a warning and add it
      IF (nyvar_i < nzmxin) THEN
        PRINT *,' WARNING: ** QS is added to initial state list **'
        PRINT *,'          ** because lana_qr_qs is set .TRUE.  **'
        nyvar_i = nyvar_i + 1
        yvarini(nyvar_i) = 'QS        '
      ELSE
        PRINT *,' ERROR:   *** QS could not be added to initial state list ***'
        PRINT *,'          *** since this list would get too large         ***'
        ierrstat = 1002
      ENDIF
    ENDIF
  ELSE
    ! if qr or qs is contained in varini, it must be deleted
    nvartmp = nyvar_i
    search_ini_qr_qs_2: DO n=1, nvartmp
      IF (     (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'QR')                    &
          .OR. (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'QS')) THEN
        IF (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'QR') THEN
          PRINT *, ' WARNING: ** QR deleted from initial state list **'
          PRINT *, '          ** because lana_qr_qs is set .FALSE.   **'
        ENDIF
        IF (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'QS') THEN
          PRINT *, ' WARNING: ** QS deleted from initial state list **'
          PRINT *, '          ** because lana_qr_qs is set .FALSE.   **'
        ENDIF
        IF (n < nyvar_i) THEN
          DO n2 = n, nyvar_i-1
            yvarini(n2) = yvarini(n2+1)
          ENDDO
        ENDIF
        yvarini(nyvar_i) = ''
        nyvar_i = nyvar_i - 1
      ENDIF
      IF (n >= nyvar_i) EXIT search_ini_qr_qs_2
    ENDDO search_ini_qr_qs_2
  ENDIF

  ! Check the setting of lana_qg, if prognostic precip is turned off, ...
  IF ( (.NOT.lprogprec .OR. itype_gscp /= 4 .OR. .NOT.lana_qr_qs)  &
       .AND. lana_qg ) THEN
    ! no prognostic precipitation is computed, so no data should be read
    ! or graupel scheme isn't used
    PRINT *, ' WARNING: ** lana_qg is set .FALSE., because '                &
           , '             prognostic precipitation is turned off,'         &
           , '             graupel scheme is not turned on'                 &
           , '             or lana_qr_qs is set .FALSE. **'
    lana_qg = .FALSE.
  ENDIF

  ! Check the setting of llb_qg, if prognostic precip is turned off, ...
  IF ( (.NOT.lprogprec .OR. itype_gscp /= 4 .OR. .NOT.llb_qr_qs)  &
       .AND. llb_qg ) THEN
    ! no prognostic precipitation is computed, so no data should be read
    ! or graupel scheme isn't used
    PRINT *, ' WARNING: ** llb_qg is set .FALSE., because '                 &
           , '             prognostic precipitation is turned off,'         &
           , '             graupel scheme is not turned on'                 &
           , '             or llb_qr_qs is set .FALSE. **'
    llb_qg = .FALSE.
  ENDIF

  IF (lana_qg) THEN
    ! check whether qg is contained in yvarini, else add it
    lzcontained = .FALSE.
    search_ini_qg: DO n=1, nyvar_i
      IF (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'QG') THEN
        lzcontained = .TRUE.
        EXIT search_ini_qg
      ENDIF
    ENDDO search_ini_qg
    IF (.NOT. lzcontained) THEN
      ! print a warning and add it
      IF (nyvar_i < nzmxin) THEN
        PRINT *,' WARNING: ** QG is added to initial state list **'
        PRINT *,'          ** because lana_qg is set .TRUE.     **'
        nyvar_i = nyvar_i + 1
        yvarini(nyvar_i) = 'QG        '
      ELSE
        PRINT *,' ERROR:   *** QG could not be added to initial state list ***'
        PRINT *,'          *** since this list would get too large         ***'
        ierrstat = 1002
      ENDIF
    ENDIF
  ELSE
    ! if qg is contained in varini, it must be deleted
    nvartmp = nyvar_i
    search_ini_qg_2: DO n=1, nvartmp
      IF (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'QG') THEN
        PRINT *, ' WARNING: ** QG deleted from initial state list **'
        PRINT *, '          ** because lana_qg is set .FALSE.     **'
        IF (n < nyvar_i) THEN
          DO n2 = n, nyvar_i-1
            yvarini(n2) = yvarini(n2+1)
          ENDDO
        ENDIF
        yvarini(nyvar_i) = ''
        nyvar_i = nyvar_i - 1
      ENDIF
      IF (n >= nyvar_i) EXIT search_ini_qg_2
    ENDDO search_ini_qg_2
  ENDIF

  IF (.NOT. lmulti_layer) THEN
    IF (lana_rho_snow) THEN
      PRINT *, ' WARNING: ** lana_rho_snow is set to .FALSE.,    ** '
      PRINT *, '          ** because lmulti_layer is turned off  ** '

      lana_rho_snow = .FALSE.
    ENDIF
  ENDIF

  IF (lana_rho_snow) THEN
    ! check whether rho_snow is contained in yvarini, else add it
    lzcontained = .FALSE.
    search_ini_rs: DO n=1, nyvar_i
      IF (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'RHO_SNOW') THEN
        lzcontained = .TRUE.
        EXIT search_ini_rs
      ENDIF
    ENDDO search_ini_rs
    IF (.NOT. lzcontained) THEN
      ! print a warning and add it
      IF (nyvar_i < nzmxin) THEN
        PRINT *,                                                          &
         ' WARNING:  ** RHO_SNOW is added to initial list because ',      &
                        'lana_rho_snow is set **'
        nyvar_i = nyvar_i + 1
        yvarini(nyvar_i) = 'RHO_SNOW  '
      ELSE
        PRINT *,                                                          &
         ' ERROR: *** RHO_SNOW could not be added to initial list because ***'
        PRINT *,                                                          &
         '        *** number of variables for initial input is too big    ***'
        ierrstat = 1002
      ENDIF
    ENDIF
  ELSE
    ! if rho_snow is contained in varini, it must be deleted
    lzcontained = .FALSE.
    search_ini_rs_2: DO n=1, nyvar_i
      IF (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'RHO_SNOW') THEN
        IF (n < nyvar_i) THEN
          IF (nyvar_i == nzmxin) THEN
            DO n2 = n, nyvar_i-1
              yvarini(n2) = yvarini(n2+1)
            ENDDO
            yvarini(nyvar_i) = ''
          ELSE
            DO n2 = n, nyvar_i
              yvarini(n2) = yvarini(n2+1)
            ENDDO
          ENDIF
        ELSE
          yvarini(n) = ''
        ENDIF
        PRINT *, ' WARNING:  ** RHO_SNOW is deleted from initial list ',     &
                                 'because lana_rho_snow is not set **'
        nyvar_i = nyvar_i - 1
        EXIT search_ini_rs_2
      ENDIF
    ENDDO search_ini_rs_2
  ENDIF

  ! Check for wrong names in case of lmulti_layer = TRUE
  IF (lmulti_layer) THEN
    DO n = 1, nyvar_i
      IF ( (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'T_M' )  .OR.     &
           (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'W_CL')  .OR.     &
           (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'T_CL')  .OR.     &
           (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'W_G1')  .OR.     &
           (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'W_G2')  .OR.     &
           (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'W_G3') ) THEN
        ! these fields are not allowed
        PRINT *,' ERROR    *** ', yvarini(n)(1:LEN_TRIM(yvarini(n))),    &
                ' is not allowed as initial field for lmulti_layer ***'
        ierrstat = 1002
      ENDIF
    ENDDO
  ENDIF

  ! Check for wrong names in case of lmulti_snow = FALSE
  IF (.NOT. lmulti_snow) THEN
    DO n = 1, nyvar_i
      IF ( (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'T_SNOW_MUL' )  .OR.     &
           (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'DZH_SNOW  ' )  .OR.     &
           (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'WTOT_SNOW ' )  .OR.     &
           (yvarini(n)(1:LEN_TRIM(yvarini(n))) == 'WLIQ_SNOW ' ) ) THEN
        ! these fields are not allowed
        PRINT *,' ERROR    *** ', yvarini(n)(1:LEN_TRIM(yvarini(n))),    &
                ' is not allowed as initial field for .NOT.lmulti_snow ***'
        ierrstat = 1002
      ENDIF
    ENDDO
  ENDIF

  ! determination of the boundary input variables
  invar = COUNT(yvarbd(:) /= '')
  IF( invar == 0) THEN
    ! nothing has been read and the default list is used
    yvarbd(1:nzmxbd_d)  = yvarbd_d(1:nzmxbd_d)
    yvarbd(nzmxbd_d+1:) = ''
    nyvar_b             = nzmxbd_d
  ELSEIF ( (invar >= 1) .AND. (invar <= nzmxin) ) THEN
    nyvar_b         = invar
  ELSEIF (invar > nzmxin) THEN
    PRINT *,' ERROR  *** number of variables for boundary input is too big ***'
    ierrstat = 1002
  ELSE
    ! something went wrong: cannot determine number of variables for input
    PRINT *,' ERROR  *** cannot determine number of variables for input ***'
    ierrstat = 1002
  ENDIF

  IF (llb_qi) THEN
    ! check whether qi is contained in yvarbd, else add it
    lzcontained = .FALSE.
    search_bd_qi: DO n=1, nyvar_b
      IF (yvarbd(n)(1:LEN_TRIM(yvarbd(n))) == 'QI') THEN
        lzcontained = .TRUE.
        EXIT search_bd_qi
      ENDIF
    ENDDO search_bd_qi
    IF (.NOT. lzcontained) THEN
      ! print a warning and add it
      IF (nyvar_b < nzmxin) THEN
        PRINT *,' WARNING: ** QI is added to boundary list    **'
        PRINT *,'          ** because llb_qi is set .TRUE.    **'
        nyvar_b = nyvar_b + 1
        yvarbd (nyvar_b) = 'QI        '
      ELSE
        PRINT *,' ERROR:   *** QI could not be added to boundary list ***'
        PRINT *,'          *** since this list would get too large    ***'
        ierrstat = 1002
      ENDIF
    ENDIF
  ELSE
    ! if qi is contained in varbd, it must be deleted
    nvartmp = nyvar_b
    search_bd_qi_2: DO n=1, nvartmp
      IF (yvarbd(n)(1:LEN_TRIM(yvarbd(n))) == 'QI') THEN
        PRINT *, ' WARNING: ** QI deleted from boundary list   **'
        PRINT *, '          ** because llb_qi is set .FALSE.   **'
        IF (n < nyvar_b) THEN
          DO n2 = n, nyvar_b-1
            yvarbd (n2) = yvarbd (n2+1)
          ENDDO
        ENDIF
        yvarbd (nyvar_b) = ''
        nyvar_b = nyvar_b - 1
      ENDIF
      IF (n >= nyvar_b) EXIT search_bd_qi_2
    ENDDO search_bd_qi_2
  ENDIF

  IF (llb_qr_qs) THEN
    ! check whether qr is contained in yvarbd, else add it
    lzcontained = .FALSE.
    search_bd_qr: DO n=1, nyvar_b
      IF (yvarbd(n)(1:LEN_TRIM(yvarbd(n))) == 'QR') THEN
        lzcontained = .TRUE.
        EXIT search_bd_qr
      ENDIF
    ENDDO search_bd_qr
    IF (.NOT. lzcontained) THEN
      ! print a warning and add it
      IF (nyvar_b < nzmxin) THEN
        PRINT *,' WARNING: ** QR is added to boundary list    **'
        PRINT *,'          ** because llb_qr_qs is set .TRUE. **'
        nyvar_b = nyvar_b + 1
        yvarbd (nyvar_b) = 'QR        '
      ELSE
        PRINT *,' ERROR:   *** QR could not be added to boundary list ***'
        PRINT *,'          *** since this list would get too large    ***'
        ierrstat = 1002
      ENDIF
    ENDIF
    ! check whether qs is contained in yvarbd, else add it
    lzcontained = .FALSE.
    search_bd_qs: DO n=1, nyvar_b
      IF (yvarbd(n)(1:LEN_TRIM(yvarbd(n))) == 'QS') THEN
        lzcontained = .TRUE.
        EXIT search_bd_qs
      ENDIF
    ENDDO search_bd_qs
    IF (.NOT. lzcontained) THEN
      ! print a warning and add it
      IF (nyvar_b < nzmxin) THEN
        PRINT *,' WARNING: ** QS is added to boundary list    **'
        PRINT *,'          ** because llb_qr_qs is set .TRUE. **'
        nyvar_b = nyvar_b + 1
        yvarbd(nyvar_b) = 'QS        '
      ELSE
        PRINT *,' ERROR:   *** QS could not be added to boundary list ***'
        PRINT *,'          *** since this list would get too large    ***'
        ierrstat = 1002
      ENDIF
    ENDIF
  ELSE
    ! if qr or qs is contained in varbd, it must be deleted
    nvartmp = nyvar_b
    search_bd_qr_qs_2: DO n=1, nvartmp
      IF (     (yvarbd(n)(1:LEN_TRIM(yvarbd(n))) == 'QR')                    &
          .OR. (yvarbd(n)(1:LEN_TRIM(yvarbd(n))) == 'QS')) THEN
        IF (yvarbd(n)(1:LEN_TRIM(yvarbd(n))) == 'QR') THEN
          PRINT *, ' WARNING: ** QR deleted from boundary list    **'
          PRINT *, '          ** because llb_qr_qs is set .FALSE. **'
        ENDIF
        IF (yvarbd(n)(1:LEN_TRIM(yvarbd(n))) == 'QS') THEN
          PRINT *, ' WARNING: ** QS deleted from boundary list    **'
          PRINT *, '          ** because llb_qr_qs is set .FALSE. **'
        ENDIF
        IF (n < nyvar_b) THEN
          DO n2 = n, nyvar_b-1
            yvarbd (n2) = yvarbd (n2+1)
          ENDDO
        ENDIF
        yvarbd (nyvar_b) = ''
        nyvar_b = nyvar_b - 1
      ENDIF
      IF (n >= nyvar_b) EXIT search_bd_qr_qs_2
    ENDDO search_bd_qr_qs_2
  ENDIF

  IF (llb_qg) THEN
    ! check whether qg is contained in yvarbd, else add it
    lzcontained = .FALSE.
    search_bd_qg: DO n=1, nyvar_b
      IF (yvarbd(n)(1:LEN_TRIM(yvarbd(n))) == 'QG') THEN
        lzcontained = .TRUE.
        EXIT search_bd_qg
      ENDIF
    ENDDO search_bd_qg
    IF (.NOT. lzcontained) THEN
      ! print a warning and add it
      IF (nyvar_b < nzmxin) THEN
        PRINT *,' WARNING: ** QG is added to boundary list    **'
        PRINT *,'          ** because llb_qg is set .TRUE.    **'
        nyvar_b = nyvar_b + 1
        yvarbd (nyvar_b) = 'QG        '
      ELSE
        PRINT *,' ERROR:   *** QG could not be added to boundary list ***'
        PRINT *,'          *** since this list would get too large    ***'
        ierrstat = 1002
      ENDIF
    ENDIF
  ELSE
    ! if qg is contained in varbd, it must be deleted
    nvartmp = nyvar_b
    search_bd_qg_2: DO n=1, nvartmp
      IF (yvarbd(n)(1:LEN_TRIM(yvarbd(n))) == 'QG') THEN
        PRINT *, ' WARNING: ** QG deleted from boundary list   **'
        PRINT *, '          ** because llb_qg is set .FALSE.   **'
        IF (n < nyvar_b) THEN
          DO n2 = n, nyvar_b-1
            yvarbd (n2) = yvarbd (n2+1)
          ENDDO
        ENDIF
        yvarbd (nyvar_b) = ''
        nyvar_b = nyvar_b - 1
      ENDIF
      IF (n >= nyvar_b) EXIT search_bd_qg_2
    ENDDO search_bd_qg_2
  ENDIF

  ! Check for wrong names in case of lmulti_layer = TRUE
  IF (lmulti_layer) THEN
    DO n = 1, nyvar_b
      IF ( (yvarbd (n)(1:LEN_TRIM(yvarbd (n))) == 'T_M' )  .OR.     &
           (yvarbd (n)(1:LEN_TRIM(yvarbd (n))) == 'W_G1')  .OR.     &
           (yvarbd (n)(1:LEN_TRIM(yvarbd (n))) == 'W_G2')  .OR.     &
           (yvarbd (n)(1:LEN_TRIM(yvarbd (n))) == 'W_G3') ) THEN
        ! these fields are not allowed
        PRINT *,' ERROR    *** ', yvarbd (n)(1:LEN_TRIM(yvarbd (n))),    &
                ' is not allowed as boundary field for lmulti_layer ***'
        ierrstat = 1002
      ENDIF
      IF ( (yvarbd (n)(1:LEN_TRIM(yvarbd (n))) == 'T_S' )  .AND.    &
           (.NOT. lbdclim) ) THEN
        ! this field is also not allowed
        PRINT *,' ERROR    *** ', yvarbd (n)(1:LEN_TRIM(yvarbd (n))),    &
                ' is not allowed as boundary field for lmulti_layer ***'
        ierrstat = 1002
      ENDIF
    ENDDO
  ENDIF

  ! Check length of the directory-names
  nzylen = LEN_TRIM (ydirini)
  IF( nzylen > 0 ) THEN
    IF( ydirini(nzylen:nzylen) /= '/') THEN
      IF( nzylen < LEN(ydirini) ) THEN
        ydirini = ydirini (1:nzylen)//'/'
      ELSE
        PRINT *,' ERROR    *** ydirini is too long *** '
        ierrstat = 1002
      ENDIF
    ENDIF
  ENDIF

  ! Maximal 3 soil water levels allowed
  IF (.NOT. lmulti_layer) THEN
    IF ( nlgw_ini > 3 ) THEN
      PRINT *,' ERROR    *** nlgw_ini = ',nlgw_ini,' > 3 *** '
      ierrstat = 1002
    ENDIF

    ! For restarts (nstart > 0): nlgw_ini = nlgw
    IF ( nstart > 0 ) THEN
      nlgw_ini = nlgw
    ENDIF
  ENDIF

  ! Check settings of lan_t_s and lan_t_so0
  IF (.NOT. lmulti_layer .AND. lan_t_so0) THEN
    PRINT *,                                                             &
     ' ERROR *** lmulti_layer=.FALSE. and lan_t_so0=.TRUE. is not possible ***'
    ierrstat = 1002
  ENDIF

  IF (lmulti_layer .AND. lan_t_s) THEN
    PRINT *,                                                             &
     ' ERROR *** lmulti_layer=.TRUE. and lan_t_s=.TRUE. is not possible ***'
    ierrstat = 1002
  ENDIF

  !----------------------------------------------------------------------------
  !- Section 4.2: Variables for boundary data
  !----------------------------------------------------------------------------

  ! Check length of the directory-names
  nzylen = LEN_TRIM (ydirbd)
  IF( nzylen > 0 ) THEN
    IF( ydirbd (nzylen:nzylen) /= '/') THEN
      IF( nzylen < LEN(ydirbd) ) THEN
        ydirbd  = ydirbd  (1:nzylen)//'/'
      ELSE
        PRINT *,' ERROR    *** ydirbd is too long *** '
        ierrstat = 1002
      ENDIF
    ENDIF
  ENDIF

  ! Maximal 3 soil water levels allowed
  IF (.NOT. lmulti_layer) THEN
    IF ( nlgw_bd > 3 ) THEN
      PRINT *,' ERROR    *** nlgw_bd = ',nlgw_bd,' > 3 *** '
      ierrstat = 1002
    ENDIF
  ENDIF

! UB>>
  ! Check whether the value for the increment of reading boundary data is
  ! given in full hours or is such a noninteger hour where the decimal part
  ! is a multiple of 5 minutes and is a divisor of 60 (ie., 5, 10, 15, 20, 30 minutes)
!!$ FOR NOW, SET TO 15 MINUTES INSTEAD OF 5, BUT 5 MINUTES MAY BE COMMENTED IN LATER ON!
  IF ( hincbound /= hincbound_d) THEN
    IF ( ABS(NINT(hincbound) - hincbound) > 1E-5) THEN
      ! then it is not a full hour, only allow multiples of 1/60 of the minutes part:
!!$     retrieve the decimal part of the time number and convert to 5 minutes intervals:
!!$       htest = MOD(hincbound, 1.0_ireals) * (60.0_ireals/5)  
      ! retrieve the decimal part of the time number and convert to 15 minutes intervals:
      htest = MOD(hincbound, 1.0_ireals) * (60.0_ireals/15)  
!      IF ( (hincbound /= 0.50_ireals) .AND. (hincbound /= 0.25_ireals) ) THEN
      IF ( ABS((htest) - NINT(htest)) > 1E-5 .OR. &
!!$           (htest > 0.5_ireals .AND. MOD(60_iintegers,NINT(htest)*5) /= 0) ) THEN
           (htest > 0.5_ireals .AND. MOD(60_iintegers,NINT(htest)*15) /= 0) ) THEN
        PRINT *, 'ERROR: *** This is not a valid hincbound: ', hincbound, ' ***'
!        PRINT *, '       *** only values = n.0 / 0.5 / 0.25 are allowed   ***'
!!$        PRINT *, '       *** only values = n.(5/60), i.e., whole 5 minutes '//&
!!$             'which are divisors of 60 (i.e., 5, 10, 15, 20, 30 minutes), are allowed for the decimal part  ***'
        PRINT *, '       *** only values = n.(15/60), i.e., whole 15 minutes '//&
             'which are divisors of 60 (i.e., 15 and 30 minutes), are allowed for the decimal part  ***'
        ierrstat = 1002
      ENDIF

      ! ytunitbd then cannot be 'd', because the filename can only resolve
      ! full hours:
      IF (ytunitbd == 'd') THEN
        PRINT *, 'ERROR: *** ytunitbd == d   can only be used,   ***'
        PRINT *, '       *** if hincbound is given in full hours ***'
        ierrstat = 1002
      ENDIF
    ENDIF
  ENDIF
! UB<<

  ! Compute an initial nincbound; but this is adjusted later on
  nincbound = NINT (3600.0_ireals * hincbound / dt)

  IF ( hnewbcdt /= hnewbcdt_d ) THEN
    newbcdt   = NINT(hnewbcdt  * 3600.0_ireals / dt)
  ENDIF
  newbcdt = ((newbcdt - 1) / nincbound + 1) * nincbound

  ! If boundary data are defined on a frame: npstrframe > 0
  IF ( lbd_frame .AND. npstrframe <= 0 ) THEN
    PRINT *,' ERROR    *** npstrframe = ',npstrframe,' <= 0 *** '
    ierrstat = 1002
  ENDIF

  IF ( lbd_frame .AND. ilevbotnoframe < 0 ) THEN
    PRINT *,' ERROR    *** ilevbotnoframe = ',ilevbotnoframe,' < 0 *** '
    ierrstat = 1002
  ENDIF
ENDIF

!------------------------------------------------------------------------------
!- Section 5: Distribute variables to all nodes
!------------------------------------------------------------------------------

IF (nproc > 1) THEN

  IF (my_world_id == 0) THEN
    intbuf  ( 1) = nlgw_ini
    intbuf  ( 2) = nlgw_bd
    intbuf  ( 3) = nincbound
    intbuf  ( 4) = newbcdt
    intbuf  ( 5) = newbc
    intbuf  ( 6) = nincboufac
    intbuf  ( 7) = nyvar_i
    intbuf  ( 8) = nyvar_b
    intbuf  ( 9) = npstrframe
    intbuf  (10) = ilevbotnoframe
    realbuf ( 1) = hincbound
    logbuf  ( 1) = lchkini
    logbuf  ( 2) = lchkbd
    logbuf  ( 3) = lan_t_s
    logbuf  ( 4) = lan_t_so0
    logbuf  ( 5) = lan_t_snow
    logbuf  ( 6) = lan_t_cl
    logbuf  ( 7) = lan_w_snow
    logbuf  ( 8) = lan_rho_snow
    logbuf  ( 9) = lan_w_i
    logbuf  (10) = lan_w_cl
    logbuf  (11) = lan_vio3
    logbuf  (12) = lan_hmo3
    logbuf  (13) = lan_plcov
    logbuf  (14) = lan_lai
    logbuf  (15) = lan_rootdp
    logbuf  (16) = lbdana
    logbuf  (17) = lana_qi
    logbuf  (18) = llb_qi
    logbuf  (19) = lana_qr_qs
    logbuf  (20) = llb_qr_qs
    logbuf  (21) = lana_qg
    logbuf  (22) = llb_qg
    logbuf  (23) = lana_rho_snow
    logbuf  (24) = lbd_frame

    ! ydirini, ydirbd have length 250, so they must be splitted
    ! because charbuf only has length 100
    charbuf ( 1)(  1:100) = ydirini(  1:100)
    charbuf ( 2)(  1:100) = ydirini(101:200)
    charbuf ( 3)(  1: 50) = ydirini(201:250)
    charbuf ( 4)(  1:100) = ydirbd (  1:100)
    charbuf ( 5)(  1:100) = ydirbd (101:200)
    charbuf ( 6)(  1: 50) = ydirbd (201:250)
    charbuf ( 7) = ytunitbd
    DO i=1,nzmxin
      charbuf(7+       i) = yvarini(i)
    ENDDO
    DO i=1,nzmxin
      charbuf(7+nzmxin+i) = yvarbd (i)
    ENDDO

  ENDIF

  CALL distribute_values (intbuf , 10, 0, imp_integers,  icomm_world, ierr)
  CALL distribute_values (realbuf,  1, 0, imp_reals,     icomm_world, ierr)
  CALL distribute_values (logbuf , 24, 0, imp_logical,   icomm_world, ierr)
  CALL distribute_values (charbuf,  7+2*nzmxin, 0, imp_character,    &
                                                         icomm_world, ierr)

  IF (my_world_id /= 0) THEN
    nlgw_ini       = intbuf  ( 1)
    nlgw_bd        = intbuf  ( 2)
    nincbound      = intbuf  ( 3)
    newbcdt        = intbuf  ( 4)
    newbc          = intbuf  ( 5)
    nincboufac     = intbuf  ( 6)
    nyvar_i        = intbuf  ( 7)
    nyvar_b        = intbuf  ( 8)
    npstrframe     = intbuf  ( 9)
    ilevbotnoframe = intbuf  (10)
    hincbound      = realbuf ( 1)
    lchkini        = logbuf  ( 1)
    lchkbd         = logbuf  ( 2)
    lan_t_s        = logbuf  ( 3)
    lan_t_so0      = logbuf  ( 4)
    lan_t_snow     = logbuf  ( 5)
    lan_t_cl       = logbuf  ( 6)
    lan_w_snow     = logbuf  ( 7)
    lan_rho_snow   = logbuf  ( 8)
    lan_w_i        = logbuf  ( 9)
    lan_w_cl       = logbuf  (10)
    lan_vio3       = logbuf  (11)
    lan_hmo3       = logbuf  (12)
    lan_plcov      = logbuf  (13)
    lan_lai        = logbuf  (14)
    lan_rootdp     = logbuf  (15)
    lbdana         = logbuf  (16)
    lana_qi        = logbuf  (17)
    llb_qi         = logbuf  (18)
    lana_qr_qs     = logbuf  (19)
    llb_qr_qs      = logbuf  (20)
    lana_qg        = logbuf  (21)
    llb_qg         = logbuf  (22)
    lana_rho_snow  = logbuf  (23)
    lbd_frame      = logbuf  (24)

    ydirini(  1:100) = charbuf ( 1)(  1:100)
    ydirini(101:200) = charbuf ( 2)(  1:100)
    ydirini(201:250) = charbuf ( 3)(  1: 50)
    ydirbd (  1:100) = charbuf ( 4)(  1:100)
    ydirbd (101:200) = charbuf ( 5)(  1:100)
    ydirbd (201:250) = charbuf ( 6)(  1: 50)
    ytunitbd       = charbuf ( 7)
    DO i=1,nzmxin
      yvarini(i)   = charbuf (7+i)
    ENDDO
    DO i=1,nzmxin
      yvarbd (i)   = charbuf (7+nzmxin+i)
    ENDDO
  ENDIF
ENDIF

! Set the array lanfld (has to be done in every processor)
IF (lmulti_layer) THEN
  lanfld( 1) = lan_t_so0
  ynaman( 1) = 'T_SO      '
ELSE
  lanfld( 1) = lan_t_s
ENDIF
lanfld( 2) = lan_t_snow
lanfld( 3) = lan_t_cl
lanfld( 4) = lan_w_snow
lanfld( 5) = lan_w_i
lanfld( 6) = lan_w_cl
lanfld( 7) = lan_vio3
lanfld( 8) = lan_hmo3
lanfld( 9) = lan_plcov
lanfld(10) = lan_lai
lanfld(11) = lan_rootdp
lanfld(12) = lan_rho_snow
IF (lmulti_layer) THEN
  lanfld(13) = lan_t_so0
  lanfld(14) = lan_t_so0
ELSE
  lanfld(13) = lan_t_s
  lanfld(14) = lan_t_s
ENDIF

!------------------------------------------------------------------------------
!- Section 6: Output of the namelist variables and their default values
!------------------------------------------------------------------------------

IF (my_world_id == 0) THEN

  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(A23)') '0     NAMELIST:  gribin'
  WRITE (nuspecif, '(A23)') '      -----------------'
  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(A26)')  'Variables for initial data'
  WRITE (nuspecif, '(T8,A,T21,A,T39,A,T58,A)') 'Variable', 'Actual Value',   &
                                               'Default Value', 'Format'

  WRITE(nuspecif, '(A)')  'Variables for initial input'
  DO i=1,nyvar_i
     WRITE (nuspecif, '(T8,A,I3,A,T21,A,T58,A)') 'yvarini(',i,')',     &
                                                            yvarini(i),'C10'
  ENDDO
  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(T8,A,T26,  A                  )')                       &
                             'ydirini (C*250)   ',ydirini(1:LEN_TRIM(ydirini))
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A)')                       &
                    'nlgw_ini      ',nlgw_ini      , nlgw_ini_d      , ' I '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lchkini       ',lchkini       , lchkini_d       , ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lbdana        ',lbdana        , lbdana_d        , ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lana_qi       ',lana_qi       , lana_qi_d       , ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'llb_qi        ',llb_qi        , llb_qi_d        , ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lana_qr_qs    ',lana_qr_qs    , lana_qr_qs_d    , ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'llb_qr_qs     ',llb_qr_qs     , llb_qr_qs_d     , ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lana_qg       ',lana_qg       , lana_qg_d       , ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'llb_qg        ',llb_qg        , llb_qg_d        , ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lana_rho_snow ',lana_rho_snow , lana_rho_snow_d , ' L '
IF(lmulti_layer) THEN
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lan_t_so0     ',lanfld( 1)    , lan_t_so0_d     , ' L '
ELSE
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lan_t_s       ',lanfld( 1)    , lan_t_s_d       , ' L '
ENDIF
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lan_t_snow    ',lanfld( 2)    , lan_t_snow_d    , ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lan_t_cl      ',lanfld( 3)    , lan_t_cl_d      , ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lan_w_snow    ',lanfld( 4)    , lan_w_snow_d    , ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lan_w_i       ',lanfld( 5)    , lan_w_i_d       , ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lan_w_cl      ',lanfld( 6)    , lan_w_cl_d      , ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lan_vio3      ',lanfld( 7)    , lan_vio3_d      , ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lan_hmo3      ',lanfld( 8)    , lan_hmo3_d      , ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lan_plcov     ',lanfld( 9)    , lan_plcov_d     , ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lan_lai       ',lanfld(10)    , lan_lai_d       , ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lan_rootdp    ',lanfld(11)    , lan_rootdp_d    , ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lan_rho_snow  ',lanfld(12)    , lan_rho_snow_d  , ' L '
  WRITE (nuspecif, '(A2)')  '  '

  WRITE(nuspecif, '(A)')  'Variables for boundary input'
  DO i=1,nyvar_b
     WRITE (nuspecif, '(T8,A,I3,A,T21,A,T58,A)') 'yvarbd(',i,')',     &
                                                             yvarbd(i),'C10'
  ENDDO
  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(T8,A,T26,  A                  )')                       &
                             'ydirbd  (C*250)   ',ydirbd(1:LEN_TRIM(ydirbd))
  WRITE (nuspecif, '(T8,A,T21,  A  ,T40,  A  ,T59,A)')                       &
                    'ytunitbd      ',ytunitbd      , ytunitbd_d      ,'C* 1'
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A)')                       &
                    'nlgw_bd       ',nlgw_bd       , nlgw_bd_d       , ' I '
  WRITE (nuspecif, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                      &
                              'hincbound   ',hincbound   ,hincbound_d ,' R '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                              'newbcdt     ',newbcdt     ,newbcdt_d   ,' I '
  WRITE (nuspecif, '(T8,A,T21,F12.4,T40,F12.4,T59,A3)')                      &
                              'hnewbcdt    ',hnewbcdt    ,hnewbcdt_d  ,' R '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                              'newbc       ',newbc       ,newbc_d     ,' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                              'nincboufac  ',nincboufac  ,nincboufac_d,' I '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lchkbd        ',lchkbd        , lchkbd_d        , ' L '
  WRITE (nuspecif, '(T8,A,T21,L12  ,T40,L12  ,T59,A)')                       &
                    'lbd_frame     ',lbd_frame     , lbd_frame_d     , ' L '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                              'npstrframe  ',npstrframe  ,npstrframe_d,' I '
  WRITE (nuspecif, '(T8,A,T21,I12  ,T40,I12  ,T59,A3)')                      &
                   'levbotnoframe  ',ilevbotnoframe  ,ilevbotnoframe_d,' I '
  WRITE (nuspecif, '(A2)')  '  '
ENDIF


! Check the setting of ldiniprec (!) if rain and snow data will be read
IF ((lana_qr_qs) .AND. (ldiniprec)) THEN
  ldiniprec = .FALSE.
  IF (my_world_id == 0) THEN
    PRINT *, ' WARNING: ** ldiniprec is set .FALSE., because '                 &
           , 'rain + snow data are to be read **'
    WRITE (nuspecif, '(A2)')  '  '
    WRITE (nuspecif, '(T8,A,T21,L12)')  'lana_qr_qs'  ,lana_qr_qs
    WRITE (nuspecif, '(T8,A        )')  '   ==>'
    WRITE (nuspecif, '(T8,A,T21,L12)')  'ldiniprec'   ,ldiniprec
    WRITE (nuspecif, '(A2)')  '  '
  ENDIF
ENDIF
  
!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE input_gribin

!==============================================================================
!+ Internal procedure in "organize_data" for the input of NAMELIST gribout
!------------------------------------------------------------------------------

SUBROUTINE input_gribout(outblock, nuspecif, nuin, ireadstat, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the NAMELIST-group gribout.
!   This group contains variables for controlling the output of GRIB-data.
!   The group gribout has to be at the end of the file INPUT and can appear
!   several times. The calling routine "organize_setup" creates a linked list
!   of all the different groups and "input_gribout" is called for every
!   group.
!
! Method:
!   All variables are initialized with default values and then read in from
!   the file INPUT. The input values are checked for errors and for
!   consistency. If wrong input values are detected the program prints
!   an error message. The program is not stopped in this routine but an
!   error code is returned to the calling routine that aborts the program after
!   reading in all other namelists.
!   In parallel mode, the variables are distributed to all nodes with the
!   environment-routine distribute_values.   
!   Both, default and input values are written to the file YUSPECIF
!   (specification of the run).
!
!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
  nuspecif,     & ! Unit number for protocolling the task
  nuin            ! Unit number for Namelist INPUT file

TYPE(pp_nl), POINTER                           ::        &
  outblock        ! pointer to the actual namelist group

! Scalar arguments with intent(inout):
INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
  ierrstat        ! error status variable

! Scalar arguments with intent(out):
INTEGER   (KIND=iintegers),   INTENT (OUT)     ::        &
  ireadstat       ! status of the read operation

!------------------------------------------------------------------------------

  INTEGER (KIND=iintegers)    ::  i, j, invar, i_huge, ndiff
  REAL (KIND=ireals)          ::  r_huge, ztti, ziv

! Inputvariable for the nest identification
  INTEGER (KIND=iintegers)     ::               &
     n_num,        n_num_d ! number of nest the output is organized for

!  Inputarrays for the output-triplets (start, end, increment)

   REAL (KIND=ireals)           ::               &
! UB>> 100 -> noutst_max
     hgrib(noutst_max),        & ! list of output steps given in hours
     hcomb(ntrip),      & ! triplets for building ngrib 
     hcomb_d(3),        & ! default values
     hgrib_d(1),        & ! default values
     htest                ! local help variable

   INTEGER (KIND=iintegers)     ::               &
! UB>> 100 -> noutst_max
     ngrib(noutst_max),        & ! list of output steps given in timesteps 
     ncomb(ntrip),      & ! triplets for building ngrib
     ngrib_d(1),        & ! default values
     ncomb_d(3)           ! default valuess

   INTEGER (KIND=iintegers)     ::               &
     nzmxml_d, nzmxpl_d, nzmxzl_d, nzmxsl_d, nzmxc_d, npls_d, nzls_d

!  Default for the name of the fields
   CHARACTER (LEN=10)           ::             &
     yvarml(nzmxml),        yvarml_d(nzmxml),  & !  name of the fields
     yvarpl(nzmxpl),        yvarpl_d(nzmxpl),  & !  name of fields for p-levels
     yvarzl(nzmxzl),        yvarzl_d(nzmxzl),  & !  name of fields for z-levels
     yvarsl(nzmxzl),        yvarsl_d(nzmxzl),  & !  name of fields for satellites
     yvarc (nzmxc ),        yvarc_d (nzmxc )     !  name of constant fields

!  Coordinates of the (sub-)domain
   REAL (KIND=ireals)           ::             &
     slon,                  slon_d,            & !  start longitude
     slat,                  slat_d,            & !  start latitude
     elon,                  elon_d,            & !  end longitude
     elat,                  elat_d,            & !  end latitude
     zlat1, zlat2, zlon1, zlon2,               & !  intermediate storage
     zendlon_tot, zendlat_tot

!  Inputvariables for the data handling
   CHARACTER (LEN=4)            ::             &
     yform_write,           yform_write_d        !  output format
   CHARACTER (LEN=4)            ::             &
     ysystem,               ysystem_d            !  datahandling system
   CHARACTER (LEN=5)            ::             &
     ydbtype,               ydbtype_d            !  database-typ
   CHARACTER (LEN=250)          ::             &
     ydir,                  ydir_d               !  directory
   CHARACTER (LEN=25)           ::             &
     ysuffix,               ysuffix_d            !  optional filename suffix
   CHARACTER                    ::             &
     ytunit,                ytunit_d             !  unit of timescale
   CHARACTER                    ::             &
     ydomain,               ydomain_d            !  sign of domain
   INTEGER (KIND=iintegers)     ::             &
     nrbit,                 nrbit_d,           & !  packrate
     nunit_of_time,         nunit_of_time_d,   & !  unit of time
     nprocess_ini_d,        nprocess_bd_d,     & !  generating process 
     nprocess_ini,          nprocess_bd          !  identification

! Inputvariable for the Checkflag
   LOGICAL                      ::             &
     lcheck,                lcheck_d,          & !  checkflag
     lwrite_const,          lwrite_const_d,    & !  write constant fields
     lanalysis,             lanalysis_d,       & !  whether analysis files 
                                                 !  have to be written
     luvmasspoint,          luvmasspoint_d       ! interpolate horizontal
                                                 ! winds to mass grid points


!  inputvariables for the pressure interpolation
   LOGICAL                      ::             &
     l_p_filter,            l_p_filter_d,      & !  filter for the p-levels
     l_z_filter,            l_z_filter_d,      & !  filter for the z-levels
     l_fi_ps_smooth,        l_fi_ps_smooth_d     !  smoothing of fi,pmsl

   REAL(KIND=ireals)            ::             &
     plev(nlevels+1),      plev_d(nlevels)       !  pressure for the p-levels

   REAL(KIND=ireals)            ::             &
     zlev(nlevels+1),      zlev_d(nlevels)       !  height of the z-levels

   REAL(KIND=ireals)            ::             &
     epsy = 1.0E-8_ireals

   INTEGER (KIND=iintegers)     :: ii, ierr, nc_hcomb, nc_hgrib, nz1hour

   CHARACTER (LEN=25)                   :: yroutine
   CHARACTER (LEN=80)                   :: yerrmsg

!  Input variables for specifying the type of vertical interpolation:

   INTEGER (KIND=iintegers)     :: itype_vertint,  itype_vertint_d

!------------------------------------------------------------------------------

!  Define the namelist group
   NAMELIST /gribout/ ngrib, hgrib, ncomb, hcomb, slon, slat, elon, elat,    &
                      yvarml, yvarpl, yvarzl, yvarsl, yform_write, ysystem,  &
                      ydir, ytunit, ydomain, ydbtype,                        &
                      nprocess_ini, nprocess_bd, n_num, nrbit, nunit_of_time,&
                      lcheck, lanalysis, lwrite_const, luvmasspoint,         &
                      l_p_filter, l_z_filter, plev, zlev, l_fi_ps_smooth,    &
                      ysuffix, itype_vertint

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
!- Section 1: Initialization with the defaults
!------------------------------------------------------------------------------

yroutine = 'input_gribout'
ireadstat = 0
npls_d    = 10
nzls_d    =  4
i_huge    = HUGE (1_iintegers)
r_huge    = HUGE (1.0_ireals)
!outblock%ngrib(:)  = i_huge
outblock%yvarml(:) = '          '
outblock%yvarpl(:) = '          '
outblock%yvarzl(:) = '          '
outblock%yvarsl(:) = '          '
outblock%yvarc (:) = '          '

IF (my_world_id == 0) THEN

   n_num_d       = 1
   hcomb_d(1)    = 0.0_ireals
   hcomb_d(2)    = REAL (INT (hstop), ireals) + 1.0_ireals
   hcomb_d(3)    = 1.0_ireals
   ncomb_d(:)    = i_huge         ! that means:  undefined
   ngrib_d(:)    = i_huge
   hgrib_d(:)    = r_huge

   yvarml_d(:)   = '          '
   yvarpl_d(:)   = '          '
   yvarzl_d(:)   = '          '
   yvarsl_d(:)   = '          '
   yvarc_d (:)   = '          '

   ! Defaults for model level output
   yvarml_d( 1)  = 'U         '; yvarml_d( 2)  = 'V         '
   yvarml_d( 3)  = 'W         '; yvarml_d( 4)  = 'T         '
   yvarml_d( 5)  = 'QV        '; yvarml_d( 6)  = 'QC        '
   yvarml_d( 7)  = 'P         '; yvarml_d( 8)  = 'PS        '
   yvarml_d( 9)  = 'T_SNOW    '; yvarml_d(10)  = 'T_S       '
   yvarml_d(11)  = 'T_G       '; yvarml_d(12)  = 'QV_S      '
   yvarml_d(13)  = 'W_SNOW    '; yvarml_d(14)  = 'W_I       '
   yvarml_d(15)  = 'RAIN_GSP  '; yvarml_d(16)  = 'SNOW_GSP  '
   yvarml_d(17)  = 'RAIN_CON  '; yvarml_d(18)  = 'SNOW_CON  '
   yvarml_d(19)  = 'U_10M     '; yvarml_d(20)  = 'V_10M     '
   yvarml_d(21)  = 'T_2M      '; yvarml_d(22)  = 'TD_2M     '
   yvarml_d(23)  = 'TMIN_2M   '; yvarml_d(24)  = 'TMAX_2M   '
   yvarml_d(25)  = 'VMAX_10M  '; yvarml_d(26)  = 'PMSL      '
   yvarml_d(27)  = 'TCM       '; yvarml_d(28)  = 'TCH       '
   yvarml_d(29)  = 'CLCT      '; yvarml_d(30)  = 'CLCH      '
   yvarml_d(31)  = 'CLCM      '; yvarml_d(32)  = 'CLCL      '
   yvarml_d(33)  = 'ALB_RAD   '; yvarml_d(34)  = 'ASOB_S    '
   yvarml_d(35)  = 'ATHB_S    '; yvarml_d(36)  = 'ASOB_T    '
   yvarml_d(37)  = 'ATHB_T    '; yvarml_d(38)  = 'APAB_S    '
   yvarml_d(39)  = 'TOT_PREC  '; yvarml_d(40)  = 'Z0        '
   yvarml_d(41)  = 'AUMFL_S   '; yvarml_d(42)  = 'AVMFL_S   '
   yvarml_d(43)  = 'ASHFL_S   '; yvarml_d(44)  = 'ALHFL_S   '
   yvarml_d(45)  = 'BAS_CON   '; yvarml_d(46)  = 'TOP_CON   '
   yvarml_d(47)  = 'HTOP_DC   '; yvarml_d(48)  = 'RUNOFF_S  '
   yvarml_d(49)  = 'RUNOFF_G  '; yvarml_d(50)  = 'HTOP_CON  '
   yvarml_d(51)  = 'HBAS_CON  '; yvarml_d(52)  = 'CLC       '

   IF (itype_conv == 3) THEN
     yvarml_d(50)  = 'HTOP_SC   '
     yvarml_d(51)  = 'HBAS_SC   '
   ENDIF

   IF (lprog_qi) THEN
     yvarml_d(53)  = 'QI        '
     nzmxml_d  = 53
   ELSE
     nzmxml_d  = 52
   ENDIF

   IF (.NOT. lmulti_layer) THEN
     yvarml_d(nzmxml_d+1)  = 'T_M       '
     yvarml_d(nzmxml_d+2)  = 'W_G1      '
     yvarml_d(nzmxml_d+3)  = 'W_G2      '
     nzmxml_d  = nzmxml_d+3
   ELSE
     yvarml_d(nzmxml_d+1)  = 'T_SO      '
     yvarml_d(nzmxml_d+2)  = 'W_SO      '
     yvarml_d(nzmxml_d+3)  = 'FRESHSNW  '
     nzmxml_d  = nzmxml_d+3
     IF (lmelt) THEN
       yvarml_d(nzmxml_d+1)  = 'W_SO_ICE  '
       nzmxml_d  = nzmxml_d+1
     ENDIF
   ENDIF

   IF (lseaice .OR. llake) THEN
     yvarml_d(nzmxml_d+1)  = 'T_ICE     '
     yvarml_d(nzmxml_d+2)  = 'H_ICE     '
     nzmxml_d              = nzmxml_d + 2
   ENDIF

   IF (llake) THEN
     yvarml_d (nzmxml_d+01)  = 'T_MNW_LK  '
     yvarml_d (nzmxml_d+02)  = 'T_WML_LK  '
     yvarml_d (nzmxml_d+03)  = 'T_BOT_LK  '
     yvarml_d (nzmxml_d+04)  = 'C_T_LK    '
     yvarml_d (nzmxml_d+05)  = 'H_ML_LK   '
     nzmxml_d                = nzmxml_d + 5
   END IF

   !_cdm In the present configuration, the bottom-sediment module of FLake 
   !     is switched off. The respective fields are handled internally by LM
   !     and are not included into the LM IO lists.
   !       yvarml_d (nzmxml_d+NN)  = 'T_B1_LK   '
   !       yvarml_d (nzmxml_d+NN)  = 'H_B1_LK   '

   IF (lmulti_snow) THEN
     yvarml_d (nzmxml_d+01)  = 'T_SNOW_MUL'
     yvarml_d (nzmxml_d+02)  = 'DZH_SNOW  '
     yvarml_d (nzmxml_d+03)  = 'WTOT_SNOW '
     yvarml_d (nzmxml_d+04)  = 'WLIQ_SNOW '
     nzmxml_d                = nzmxml_d + 4
   END IF

   IF(lbdclim) THEN
     yvarml_d (nzmxml_d+1)  = 'PLCOV     '
     yvarml_d (nzmxml_d+2)  = 'LAI       '
     yvarml_d (nzmxml_d+3)  = 'ROOTDP    '
     nzmxml_d = nzmxml_d + 3
     IF (.NOT. lmulti_layer) THEN
       yvarml_d (nzmxml_d+1)  = 'T_CL      '
       yvarml_d (nzmxml_d+2)  = 'W_CL      '
       nzmxml_d = nzmxml_d + 2
     ENDIF
   ENDIF

   ! Defaults for constant variables output
   yvarc_d ( 1)  = 'HHL       '; yvarc_d ( 2)  = 'HSURF     '
   yvarc_d ( 3)  = 'FIS       '; yvarc_d ( 4)  = 'FC        '
   yvarc_d ( 5)  = 'RLAT      '; yvarc_d ( 6)  = 'RLON      '
   yvarc_d ( 7)  = 'FR_LAND   '; yvarc_d ( 8)  = 'SOILTYP   '
   nzmxc_d = 8

   IF (.NOT. lbdclim) THEN
     yvarc_d ( 9)  = 'VIO3      '; yvarc_d (10)  = 'HMO3      '
     yvarc_d (11)  = 'PLCOV     '; yvarc_d (12)  = 'LAI       '
     yvarc_d (13)  = 'ROOTDP    '
     nzmxc_d = 13
     IF (.NOT. lmulti_layer) THEN
       yvarc_d (14)  = 'T_CL      '
       yvarc_d (15)  = 'W_CL      '
       nzmxc_d = 15
     ENDIF
   ENDIF

   IF (llake) THEN
     yvarc_d (nzmxc_d + 1)  = 'FR_LAKE   '
     yvarc_d (nzmxc_d + 2)  = 'DEPTH_LK  '
     nzmxc_d = nzmxc_d + 2
   END IF

   IF (lforest) THEN
     yvarc_d (nzmxc_d + 1)  = 'FOR_E     '
     yvarc_d (nzmxc_d + 2)  = 'FOR_D     '
     nzmxc_d = nzmxc_d + 2
   ENDIF

  IF (lradtopo) THEN
    yvarc_d(nzmxc_d + 1) = 'SKYVIEW   '
    yvarc_d(nzmxc_d + 2) = 'SLO_ASP   '
    yvarc_d(nzmxc_d + 3) = 'SLO_ANG   '
    yvarc_d(nzmxc_d + 4) = 'HORIZON   '
    nzmxc_d = nzmxc_d + 4
  ENDIF

   IF (lsso) THEN
     yvarc_d (nzmxc_d + 1)  = 'SSO_STDH  '; yvarc_d (nzmxc_d + 2)  = 'SSO_GAMMA '
     yvarc_d (nzmxc_d + 3)  = 'SSO_THETA '; yvarc_d (nzmxc_d + 4)  = 'SSO_SIGMA '
     nzmxc_d = nzmxc_d + 4
   ENDIF

   IF (lemiss) THEN
     yvarc_d(nzmxc_d+1)    = 'EMIS_RAD  '
     nzmxc_d  = nzmxc_d+1
   END IF

   IF (lstomata) THEN
     yvarc_d(nzmxc_d+1)    = 'RSMIN     '
     nzmxc_d  = nzmxc_d+1
   END IF
   
   IF ( itype_aerosol == 2 ) THEN
     yvarc_d(nzmxc_d+1)    = 'AER_SO4   '
     yvarc_d(nzmxc_d+2)    = 'AER_DUST  '
     yvarc_d(nzmxc_d+3)    = 'AER_ORG   '
     yvarc_d(nzmxc_d+4)    = 'AER_BC    '
     yvarc_d(nzmxc_d+5)    = 'AER_SS    '
     nzmxc_d  = nzmxc_d+5
   ENDIF

   !_cdm The following FLake external parameters are assigned their default 
   !     values that are kept constant in space and time.
   !    These external-parameter fields are handled internally by LM
   !    and are not included into the LM IO lists
   !    (this should be done in the future).
   !        yvarc_d (16)  = 'FETCH_LK  '
   !        yvarc_d (17)  = 'DP_BS_LK  '
   !        yvarc_d (18)  = 'T_BS_LK   '
   !        yvarc_d (19)  = 'GAMSO_LK  '

   slon_d = startlon_tot
   slat_d = startlat_tot

   elon_d = startlon_tot  + REAL(ie_tot-1, ireals) * dlon
   elat_d = startlat_tot  + REAL(je_tot-1, ireals) * dlat

   ! The end-longitude has to be limited to the range (-180.0,+180.0)
   IF (elon_d > 180.0_ireals) THEN
     elon_d = elon_d - 360.0_ireals
   ENDIF

   yform_write_d = 'grb1'
   ysystem_d     = 'FILE'
   ydbtype_d     = ''
   ydir_d(:)     = ' '
   ysuffix_d     = ''
   ytunit_d      = 'f'
   ydomain_d     = 'f'
   nrbit_d       = 16

   nunit_of_time_d = 1_iintegers
   nprocess_ini_d  = -999999_iintegers
   nprocess_bd_d   = -999999_iintegers
   lcheck_d        = .FALSE.
   lwrite_const_d  = .TRUE.
   lanalysis_d     = .FALSE.
   luvmasspoint_d  = .FALSE.

!  defaults of pressure interpolation
   l_p_filter_d     = .FALSE.
   l_z_filter_d     = .FALSE.
   l_fi_ps_smooth_d = .FALSE.

   plev_d(1)     = 200_ireals  ; plev_d(2)     = 250_ireals
   plev_d(3)     = 300_ireals  ; plev_d(4)     = 400_ireals
   plev_d(5)     = 500_ireals  ; plev_d(6)     = 600_ireals
   plev_d(7)     = 700_ireals  ; plev_d(8)     = 850_ireals
   plev_d(9)     = 950_ireals  ; plev_d(10)    =1000_ireals

   yvarpl_d( 1)= 'T         ' ; yvarpl_d( 2)= 'RELHUM    '
   yvarpl_d( 3)= 'U         ' ; yvarpl_d( 4)= 'V         '
   yvarpl_d( 5)= 'FI        ' ; yvarpl_d( 6)= 'OMEGA     '
   nzmxpl_d  =  6
   
   zlev_d(1)     = 1000_ireals  ; zlev_d(2)     = 2000_ireals
   zlev_d(3)     = 3000_ireals  ; zlev_d(4)     = 5000_ireals

   yvarzl_d( 1)= 'T         ' ; yvarzl_d( 2)= 'RELHUM    '
   yvarzl_d( 3)= 'U         ' ; yvarzl_d( 4)= 'V         '
   yvarzl_d( 5)= 'P         ' ; yvarzl_d( 6)= 'W         '
   nzmxzl_d  =  6

!  defaults for satellite output
   yvarsl_d( 1)= 'SYNME7    ' ; yvarsl_d( 2)= 'SYNMSG    '
   nzmxsl_d  =  2

   ! default for type of vertical interpolation to p- and z-levels
   ! 1 = cubic spline,   2 = linear interpolation
   itype_vertint_d = 1

!------------------------------------------------------------------------------
!- Section 2: Initialize variables with defaults
!------------------------------------------------------------------------------

   n_num             =  n_num_d
   yvarml(:)         =  ''
   yvarpl(:)         =  ''
   yvarzl(:)         =  ''
   yvarsl(:)         =  ''
   yvarc (:)         =  ''

   ! The variables for the hour- and/or time step lists are set after 
   ! reading the Namelist group; here everything is set to huge
   hgrib(:)          =  r_huge
   ngrib(:)          =  i_huge
   hcomb(:)          =  r_huge
   ncomb(:)          =  i_huge

   slon              =  slon_d
   slat              =  slat_d
   elon              =  elon_d
   elat              =  elat_d

   yform_write       =  yform_write_d
   ysystem           =  ysystem_d
   ydbtype           =  ydbtype_d
   ydir              =  ydir_d
   ysuffix           =  ysuffix_d
   ytunit            =  ytunit_d
   ydomain           =  ydomain_d
   nrbit             =  nrbit_d

   nunit_of_time     =  nunit_of_time_d
   nprocess_ini      =  nprocess_ini_d
   nprocess_bd       =  nprocess_bd_d
   lcheck            =  lcheck_d
   lwrite_const      =  lwrite_const_d
   lanalysis         =  lanalysis_d
   luvmasspoint      =  luvmasspoint_d

   l_p_filter        =  l_p_filter_d
   l_z_filter        =  l_z_filter_d
   l_fi_ps_smooth    =  l_fi_ps_smooth_d

   plev(:)           =  r_huge
   zlev(:)           =  r_huge

   itype_vertint     =  itype_vertint_d

!------------------------------------------------------------------------------
!- Section 3: Read namelist values
!------------------------------------------------------------------------------

  READ(nuin,gribout,IOSTAT=ireadstat)
ENDIF

IF (nproc > 1) THEN
  ! IOSTAT is distributed to all nodes
  CALL distribute_values (ireadstat, 1, 0, imp_integers, icomm_world, ierr)
ENDIF

IF (ireadstat /= 0) THEN
  ierrstat = -1
  RETURN
ENDIF

!------------------------------------------------------------------------------
!- Section 4: Store settings to the pointer in the linked list
!------------------------------------------------------------------------------

IF (my_world_id == 0 ) THEN

! UB>> In case of idealized runs, enable arbitrary output time specifications using hcomb 
!      by setting outblock%lhour = .false.:
  IF ( lartif_data ) THEN

    outblock%lhour = .FALSE.

  ELSE  
! UB<<

    ! If nothing is read for the hour- or time step lists, the defaults
    ! are set to hcomb_d
    IF     ( (hgrib(1) == r_huge) .AND. (hcomb(1) == r_huge) .AND.            &
         (ngrib(1) == i_huge) .AND. (ncomb(1) == i_huge) ) THEN
      hcomb(1:3)        =  hcomb_d(1:3)
      
      ! Indicate that all values have to be rounded to multiples of 0.25
      outblock%lhour = .TRUE.
    ELSEIF ( ((hgrib(1) /= r_huge) .OR. (hcomb(1) /= r_huge)) .AND.           &
         ((ngrib(1) /= i_huge) .OR. (ncomb(1) /= i_huge)) ) THEN
      ! values have been read for both: hour- and timestep-specification
      ! but only one kind is allowed
      PRINT *, '*** ERROR in output specification: ***'
      PRINT *, '*** You can either specify hour- or time-step values, but not both ***'
      ierrstat = 1001
      RETURN
    ELSEIF ( (hcomb(1) /= r_huge) .OR. (hgrib(1) /= r_huge) ) THEN
      ! hour-values have been specified: check that only full hour values,
      ! 0.5 or 0.25 are set
      
      nc_hcomb = COUNT(hcomb(:)  /= r_huge)
      DO i = 1, nc_hcomb
        ! this must be an integer value, to allow only multiples of 0.25
        htest = hcomb(i) / 0.25_ireals   
        IF ( ABS(NINT(htest) - htest) > 1E-5) THEN
          ! then it is not an integer value
          PRINT *, 'ERROR: *** This is not a valid value for hcomb: ', hcomb(i), ' ***'
          PRINT *, '       *** only values = n.0 / n.5 / n.25 are allowed   ***'
          ierrstat = 1002
        ENDIF
      ENDDO
      
      nc_hgrib = COUNT(hgrib(:)  /= r_huge)
      DO i = 1, nc_hgrib
        ! this must be an integer value, to allow only multiples of 0.25
        htest = hgrib(i) / 0.25_ireals   
        IF ( ABS(NINT(htest) - htest) > 1E-5) THEN
          ! then it is not an integer value
          PRINT *, 'ERROR: *** This is not a valid value for hgrib: ', hgrib(i), ' ***'
          PRINT *, '       *** only values = n.0 / n.5 / n.25 are allowed   ***'
          ierrstat = 1002
        ENDIF
      ENDDO
      
      ! Indicate that all values have to be rounded to multiples of 0.25
      outblock%lhour = .TRUE.
    ELSE
      ! Indicate that all values have to be rounded to multiples of 0.25
      outblock%lhour = .FALSE.
    END IF
    
! UB>>
  END IF
! UB<<

  IF ( (n_num < 1) .OR. (n_num > n_num_tot) ) THEN
    ! we have to go back, because with a wrong n_num, the program
    ! will abort!
    ierrstat =  2
    RETURN
  ENDIF

  ! calculation of the output-timesteps:  if output time steps are given
  ! in hcomb, hgrib or ncomb, the values for ngrib are calculated
  CALL calc_ngrib(hgrib, ngrib, hcomb, ncomb, outblock, grids_dt(n_num), &
                  ierrstat)

  ! adjustments for writing analysis files
  IF (lanalysis) THEN
    ! need pressure deviation for laf-files
    yvarml_d( 7)  = 'PP        '

    ! add additional variables to the default output list
    yvarml_d(nzmxml_d+ 1)  = 'VIO3      '
    yvarml_d(nzmxml_d+ 2)  = 'HMO3      '
    yvarml_d(nzmxml_d+ 3)  = 'PLCOV     '
    yvarml_d(nzmxml_d+ 4)  = 'LAI       '
    yvarml_d(nzmxml_d+ 5)  = 'ROOTDP    '
    yvarml_d(nzmxml_d+ 6)  = 'FR_LAND   '
    yvarml_d(nzmxml_d+ 7)  = 'SOILTYP   '
    yvarml_d(nzmxml_d+ 8)  = 'HHL       '
    yvarml_d(nzmxml_d+ 9)  = 'HSURF     '
    yvarml_d(nzmxml_d+10)  = 'FIS       '
    yvarml_d(nzmxml_d+11)  = 'RLAT      '
    yvarml_d(nzmxml_d+12)  = 'RLON      '
    yvarml_d(nzmxml_d+13)  = 'T_CL      '
    yvarml_d(nzmxml_d+14)  = 'W_CL      '
    nzmxml_d = nzmxml_d + 14

    ! Add external-parameter fields for sso scheme
    IF (lsso) THEN
      yvarml_d(nzmxml_d+ 1)  = 'SSO_STDH  '
      yvarml_d(nzmxml_d+ 2)  = 'SSO_GAMMA '
      yvarml_d(nzmxml_d+ 3)  = 'SSO_THETA '
      yvarml_d(nzmxml_d+ 4)  = 'SSO_SIGMA '
      nzmxml_d = nzmxml_d + 4
    ENDIF

    IF (lemiss) THEN
      yvarml_d(nzmxml_d+ 1)  = 'EMIS_RAD  '
      nzmxml_d = nzmxml_d + 1
    ENDIF

    IF (lstomata) THEN
      yvarml_d(nzmxml_d+ 1)  = 'RSMIN     '
      nzmxml_d = nzmxml_d + 1
    ENDIF

    IF ( itype_aerosol == 2 ) THEN
      yvarml_d(nzmxml_d + 1)  = 'AER_SO4   '
      yvarml_d(nzmxml_d + 2)  = 'AER_DUST  '
      yvarml_d(nzmxml_d + 3)  = 'AER_ORG   '
      yvarml_d(nzmxml_d + 4)  = 'AER_BC    '
      yvarml_d(nzmxml_d + 5)  = 'AER_SS    '
      nzmxml_d = nzmxml_d + 5
    ENDIF

    ! Add external-parameter fields for forest
    IF (lforest) THEN
      yvarml_d(nzmxml_d+ 1)  = 'FOR_E     '
      yvarml_d(nzmxml_d+ 2)  = 'FOR_D     '
      nzmxml_d = nzmxml_d + 2
    ENDIF

    ! Add external-parameter fields of the lake model FLake
    IF (llake) THEN
      yvarml_d(nzmxml_d+ 1)  = 'FR_LAKE   '
      yvarml_d(nzmxml_d+ 2)  = 'DEPTH_LK  '
      nzmxml_d = nzmxml_d + 2
    ENDIF

    ! only full hours > 0 can be output steps
    ! Additional check, if all hcomb- and hgrib-values are in full hours
    ! (most important: the increment hcomb(3)
    IF ( (hcomb(1) /= r_huge) .OR. (hgrib(1) /= r_huge) ) THEN
      ! hour-values have been specified: check that only full hour values,
      ! are set

      nc_hcomb = COUNT(hcomb(:)  /= r_huge)
      DO i = 1, nc_hcomb
        ! this must be an integer value, to allow only full hours
        htest = hcomb(i)
        IF ( ABS(NINT(htest) - htest) > 1E-5) THEN
          ! then it is not an integer value
          PRINT *, 'ERROR: *** This is not a valid value for hcomb: ', hcomb(i), ' ***'
          PRINT *, '       *** for lanalysis only full hours are allowed   ***'
          ierrstat = 1002
        ENDIF
      ENDDO

      nc_hgrib = COUNT(hgrib(:)  /= r_huge)
      DO i = 1, nc_hgrib
        ! this must be an integer value, to allow only full hours
        htest = hgrib(i)
        IF ( ABS(NINT(htest) - htest) > 1E-5) THEN
          ! then it is not an integer value
          PRINT *, 'ERROR: *** This is not a valid value for hgrib: ', hgrib(i), ' ***'
          PRINT *, '       *** for lanalysis only full hours are allowed   ***'
          ierrstat = 1002
        ENDIF
      ENDDO
    ENDIF

    ! and check output is done only from hour 1 onwards
    ! determine timestep for hour 1 = 3600.0 seconds:
    nz1hour = NINT (3600.0_ireals / grids_dt(n_num))

    i = 1
    DO WHILE (outblock%ngrib(i) < nz1hour)
      PRINT *, 'WARNING: *** Output for lanalysis can only start with hour 1!! ***'
      PRINT *, '         *** STEP ', outblock%ngrib(i),'  has been eliminated  ***'
      i = i+1
    ENDDO

    DO j = i, outblock%outsteps
      outblock%ngrib(j-i+1) = outblock%ngrib(j)
    ENDDO
    outblock%outsteps = outblock%outsteps - i + 1
  ENDIF

  ! determination of the output fields
  invar = COUNT(yvarml(:) /= '')
  IF((yvarml(1) == 'default').OR.(yvarml(1)=='DEFAULT'))THEN
    outblock%yvarml(1:nzmxml_d)  = yvarml_d(1:nzmxml_d)
    outblock%yvarml(nzmxml_d+1:) = ''
    outblock%nyvar_m             = nzmxml_d
  ELSE 
    IF( invar >= 1) THEN 
      DO ii = 1, invar
        outblock%yvarml(ii)(1:LEN_TRIM(yvarml(ii)))          &
                          = yvarml(ii)(1:LEN_TRIM(yvarml(ii)))
      ENDDO
      outblock%yvarml(invar+1:)= ''
      outblock%nyvar_m         = invar
    ELSE
      outblock%yvarml(:)       = ''
      outblock%nyvar_m         = 0
    ENDIF
  ENDIF

  invar = COUNT(yvarpl(:) /= '')
  IF((yvarpl(1) == 'default').OR.(yvarpl(1)=='DEFAULT'))THEN
    outblock%yvarpl(1:nzmxpl_d)  = yvarpl_d(1:nzmxpl_d)
    outblock%yvarpl(nzmxpl_d+1:) = ''
    outblock%nyvar_p               = nzmxpl_d
  ELSE
    IF( invar >= 1) THEN
      DO ii = 1, invar
        outblock%yvarpl(ii)(1:LEN_TRIM(yvarpl(ii)))          &
                          = yvarpl(ii)(1:LEN_TRIM(yvarpl(ii)))
      ENDDO
      outblock%yvarpl(invar+1:)= ''
      outblock%nyvar_p         = invar
    ELSE
      outblock%yvarpl(:)       = ''
      outblock%nyvar_p         = 0
    ENDIF
  ENDIF

  invar = COUNT(yvarzl(:) /= '')
  IF((yvarzl(1) == 'default').OR.(yvarzl(1)=='DEFAULT'))THEN
    outblock%yvarzl(1:nzmxzl_d)  = yvarzl_d(1:nzmxzl_d)
    outblock%yvarzl(nzmxzl_d+1:) = ''
    outblock%nyvar_z               = nzmxzl_d
  ELSE
    IF( invar >= 1) THEN
      DO ii = 1, invar
        outblock%yvarzl(ii)(1:LEN_TRIM(yvarzl(ii)))          &
                          = yvarzl(ii)(1:LEN_TRIM(yvarzl(ii)))
      ENDDO
      outblock%yvarzl(invar+1:)= ''
      outblock%nyvar_z         = invar
    ELSE
      outblock%yvarzl(:)       = ''
      outblock%nyvar_z         = 0
    ENDIF
  ENDIF

  invar = COUNT(yvarsl(:) /= '')
  IF((yvarsl(1) == 'default').OR.(yvarsl(1)=='DEFAULT'))THEN
    outblock%yvarsl(1:nzmxsl_d)  = yvarsl_d(1:nzmxsl_d)
    outblock%yvarsl(nzmxsl_d+1:) = ''
    outblock%nyvar_s             = nzmxsl_d
  ELSE
    IF( invar >= 1) THEN
      DO ii = 1, invar
        outblock%yvarsl(ii)(1:LEN_TRIM(yvarsl(ii)))          &
                          = yvarsl(ii)(1:LEN_TRIM(yvarsl(ii)))
      ENDDO
      outblock%yvarsl(invar+1:)= ''
      outblock%nyvar_s         = invar
    ELSE
      outblock%yvarsl(:)       = ''
      outblock%nyvar_s         = 0
    ENDIF
  ENDIF

  outblock%nyvar_c = nzmxc_d
  outblock%yvarc (1:outblock%nyvar_c) = yvarc_d (1:outblock%nyvar_c)

! assign the coordinates
  outblock%slon = slon
  outblock%slat = slat
  outblock%elon = elon
  outblock%elat = elat

! check the longitude coordinates
  zendlon_tot = startlon_tot + (ie_tot-1)*dlon
  IF (zendlon_tot > 180.0_ireals) THEN
    zendlon_tot = zendlon_tot - 360.0_ireals
  ENDIF

  IF (slon > elon) THEN
    IF (elon >= 0.0_ireals) THEN
      PRINT *, ' *** ERROR: end-lon is east of start-lon ', slon, elon
      ierrstat = 10
    ELSE
      ! domain is around the 180-Meridian
      IF (slon >= 0.0_ireals) THEN
        ! this is west of the 180-Meridian
        IF (slon < startlon_tot) THEN
          PRINT *, ' *** ERROR: start-lon is outside total domain: ', &
                                slon, startlon_tot
          ierrstat = 10
        ENDIF
      ELSE
        ! now it is: slon < 0, elon < 0 and slon > elon
        PRINT *, ' *** ERROR: start-lon is east of end-lon: ', slon, elon
        ierrstat = 10
      ENDIF
      IF (elon > zendlon_tot) THEN
        PRINT *, ' ERROR: end-lon is outside the domain: ', elon, zendlon_tot
        ierrstat = 10
      ENDIF
    ENDIF
  ELSE
    ! normal case: domain is not around the 180-Meridian
    IF (slon < startlon_tot) THEN
      PRINT *, ' ERROR: slon is outside total domain: ', slon, startlon_tot
      ierrstat = 10
    ENDIF
    IF (elon > zendlon_tot) THEN
      PRINT *, ' ERROR: elon is outside total domain: ', elon, zendlon_tot
      ierrstat = 10
    ENDIF
    IF (slon > elon) THEN
      PRINT *, ' ERROR: end-lon is west from start-lon: ', slon, elon
      ierrstat = 10
    ENDIF
  ENDIF

! check the latitude coordinates
  zendlat_tot = startlat_tot + (je_tot-1)*dlat
  IF (slat < startlat_tot) THEN
    PRINT *, ' ERROR: start-lat is outside total domain: ', slat, startlat_tot
    ierrstat = 10
  ENDIF
  IF (elat > zendlat_tot) THEN
    PRINT *, ' ERROR: end-lat is outside total domain: ', elat, zendlat_tot
    ierrstat = 10
  ENDIF
  IF (slat > elat) THEN
    PRINT *, ' ERROR: end-lat is south of start-lat: ', slat, elat
    ierrstat = 10
  ENDIF

! check and set of yform_write
#ifndef GRIBDWD
  IF (yform_write == 'grb1') THEN
    PRINT *, ' ERROR  *** yform_write = grb1, but model is not compiled to use DWD GRIB library ***'
    ierrstat = 1002
  ENDIF
#endif

#ifndef NETCDF
  IF (yform_write == 'ncdf') THEN
    PRINT *, ' ERROR  *** yform_write = ncdf, but model is not compiled to use NetCDF ***'
    ierrstat = 1002
  ENDIF
#endif

  SELECT CASE(yform_write)
  CASE('grb1','ncdf')
    outblock%yform_write(1:4) = yform_write(1:4)
  CASE('')
    PRINT *, 'Parameter yform_write is not set'
    ierrstat = 1
  CASE DEFAULT
    PRINT *, yform_write,' is not a known parameter for yform_write'
    ierrstat = 1
  END SELECT

! check and set of ysystem
  SELECT CASE(ysystem)
  CASE('file','FILE')
    outblock%ysystem = 'FILE'
  CASE('daba','DABA')
    PRINT *, 'Database not implemented'
    ierrstat = 1
  CASE('ecfs','ECFS')
    PRINT *, 'ECFS not implemented'
    ierrstat = 1
  CASE('')
    PRINT *, 'Parameter YSYSTEM not set'
    ierrstat = 1
  CASE DEFAULT
    PRINT *, ysystem,' is not a known parameter for YSYSTEM'
    ierrstat = 1
  END SELECT

! set of ydir
  IF(ydir(1:1) /= ' ' .AND. ydir(len_trim(ydir):len_trim(ydir))/='/') THEN
     ydir = ydir(1:len_trim(ydir))//'/'
  ENDIF
  outblock%ydir = ydir

! set of ysuffix
  outblock%ysuffix = ysuffix

! set of ytunit
  SELECT CASE (ytunit)
  CASE('f','t','c','d')
    outblock%ytunit = ytunit
  CASE DEFAULT
    PRINT *, ytunit,' is not a known parameter for YTUNIT'
    ierrstat = 3
  END SELECT

! set of ydomain
  IF ( (slat /= startlat_tot) .OR. (elat /= zendlat_tot)  .OR.   &
       (slon /= startlon_tot) .OR. (elon /= zendlon_tot) ) THEN
    outblock%ydomain = 's'
  ELSE
    outblock%ydomain = 'f'
  ENDIF

! set n_num, nrbit, nprocess*, ydbtype, lcheck,
! lanalysis, luvmasspoint
  outblock%n_num            = n_num
  outblock%nrbit            = nrbit
  outblock%nprocess_ini_out = nprocess_ini
  outblock%nprocess_bd_out  = nprocess_bd
  outblock%ydbtype          = ydbtype
  outblock%lcheck           = lcheck
  outblock%lwrite_const     = lwrite_const
  outblock%lanalysis        = lanalysis
  outblock%luvmasspoint     = luvmasspoint

! check and set the unit-of-time
  SELECT CASE (nunit_of_time)
  CASE (0, 1, 2, 10, 11, 12, 13, 14, 15)
    outblock%nunit_of_time    = nunit_of_time
  CASE DEFAULT
    PRINT *, 'unit-of-time = ', nunit_of_time,              &
                   ' is not a known value for unit-of-time in the LM'
    ierrstat = 4
  END SELECT

  invar = COUNT(plev(:) /= r_huge)
  IF (invar > nlevels) THEN
    PRINT *, ' *** ERROR:  number of P-levels greater than allowed ***'
    PRINT *, ' ***         set higher number of nlevels in data_io ***'
    ierrstat = 4
  ENDIF

  IF (invar >= 1) THEN
    DO i=1,invar-1
      IF(plev(i)>plev(i+1)) THEN
        yerrmsg = ' P-level are not in the right order '
        CALL model_abort (my_cart_id, 1000, yerrmsg, yroutine)
      ENDIF
    ENDDO
    outblock%plev(1:invar)    = plev(1:invar) * 100_ireals
    outblock%plev(invar+1:)   = r_huge
    outblock%kepin            = invar
  ELSE
    outblock%plev(1:npls_d)   = plev_d(1:npls_d) * 100.0_ireals
    outblock%plev(npls_d+1:)  = r_huge
    outblock%kepin            = npls_d
  ENDIF

  invar = COUNT(zlev(:) /= r_huge)
  IF (invar > nlevels) THEN
    PRINT *, ' *** ERROR:  number of Z-levels greater than allowed ***'
    PRINT *, ' ***         set higher number of nlevels in data_io ***'
    ierrstat = 4
  ENDIF

  IF (invar >= 1) THEN
    DO i=1,invar-1
      IF(zlev(i)>zlev(i+1)) THEN
        yerrmsg = ' Z-levels are not in the right order '
        CALL model_abort (my_cart_id, 1000, yerrmsg, yroutine)
      ENDIF
    ENDDO
    outblock%zlev(1:invar)    = zlev(1:invar) 
    outblock%zlev(invar+1:)   = r_huge
    outblock%kezin            = invar
  ELSE
    outblock%zlev(1:nzls_d)   = zlev_d(1:nzls_d)
    outblock%zlev(nzls_d+1:)  = r_huge
    outblock%kezin            = nzls_d
  ENDIF

  outblock%l_p_filter     = l_p_filter
  outblock%l_z_filter     = l_z_filter
  outblock%l_fi_ps_smooth = l_fi_ps_smooth

  outblock%itype_vertint  = itype_vertint

  IF (itype_vertint < 1 .OR. itype_vertint > 2) THEN
    PRINT *, ' *** ERROR:  wrong value for itype_vertint        ***'
    PRINT *, ' ***         only 1 or 2 is currently implemented ***'
    ierrstat = 4
  END IF

!------------------------------------------------------------------------------
!- Section 5: distribute record to the other processors
!------------------------------------------------------------------------------

  IF (nproc > 1) THEN
  ! broadcast namelist
    intbuf( 1)  = outblock%nyvar_m
    intbuf( 2)  = outblock%nyvar_p
    intbuf( 3)  = outblock%nyvar_z
    intbuf( 4)  = outblock%nyvar_s
    intbuf( 5)  = outblock%nyvar_c
    intbuf( 6)  = outblock%nrbit
    intbuf( 7)  = outblock%nprocess_ini_out
    intbuf( 8)  = outblock%nprocess_bd_out
    intbuf( 9)  = outblock%outsteps
    intbuf(10)  = outblock%kepin
    intbuf(11)  = outblock%kezin
    intbuf(12)  = outblock%n_num
    intbuf(13)  = outblock%nunit_of_time
    intbuf(14)  = outblock%itype_vertint

    realbuf(1) = outblock%slon
    realbuf(2) = outblock%slat
    realbuf(3) = outblock%elon
    realbuf(4) = outblock%elat
!   DO i=1,nlevels
!     realbuf(4+i) = outblock%plev(i)
!   ENDDO
!   DO i=1,nlevels
!     realbuf(4+nlevels+i) = outblock%zlev(i)
!   ENDDO

    charbuf(1) = outblock%ysystem
    charbuf(2) = outblock%ydbtype
    charbuf(3)(1:100) = outblock%ydir(  1:100)
    charbuf(4)(1:100) = outblock%ydir(101:200)
    charbuf(5)(1: 50) = outblock%ydir(201:250)
    charbuf(6) = outblock%ysuffix
    charbuf(7) = outblock%ytunit
    charbuf(8) = outblock%ydomain
    charbuf(9) = outblock%yform_write
    DO i=1,nzmxml
      charbuf(9+       i)(1:10) = outblock%yvarml(i)(1:10)
    ENDDO
    DO i=1,nzmxpl
      charbuf(9+nzmxml+i)(1:10) = outblock%yvarpl(i)(1:10)
    ENDDO
    DO i=1,nzmxzl
      charbuf(9+nzmxml+nzmxpl+i)(1:10) = outblock%yvarzl(i)(1:10)
    ENDDO
    DO i=1,nzmxzl
      charbuf(9+nzmxml+nzmxpl+nzmxzl+i)(1:10) = outblock%yvarsl(i)(1:10)
    ENDDO
    DO i=1,nzmxc
      charbuf(9+nzmxml+nzmxpl+2*nzmxzl+i)(1:10) = outblock%yvarc(i)(1:10)
    ENDDO

    logbuf(1)  = outblock%lcheck
    logbuf(2)  = outblock%lwrite_const
    logbuf(3)  = outblock%lanalysis
    logbuf(4)  = outblock%l_z_filter
    logbuf(5)  = outblock%l_p_filter
    logbuf(6)  = outblock%l_fi_ps_smooth
    logbuf(7)  = outblock%luvmasspoint
    logbuf(8)  = outblock%lhour
  ENDIF
ENDIF

IF (nproc > 1) THEN
  CALL distribute_values (intbuf , 14, 0, imp_integers, icomm_world, ierr)
  CALL distribute_values (realbuf,  4, 0, imp_reals,    icomm_world, ierr)
  CALL distribute_values (logbuf ,  8, 0, imp_logical,  icomm_world, ierr)
  CALL distribute_values (charbuf, 9+nzmxml+nzmxpl+2*nzmxzl+nzmxc, 0,     &
                                        imp_character, icomm_world, ierr)
  CALL distribute_values (outblock%plev, nlevels, 0, imp_reals, icomm_world, ierr)
  CALL distribute_values (outblock%zlev, nlevels, 0, imp_reals, icomm_world, ierr)

  IF (my_world_id /= 0) THEN
    outblock%nyvar_m          = intbuf( 1)
    outblock%nyvar_p          = intbuf( 2)
    outblock%nyvar_z          = intbuf( 3)
    outblock%nyvar_s          = intbuf( 4)
    outblock%nyvar_c          = intbuf( 5)
    outblock%nrbit            = intbuf( 6)
    outblock%nprocess_ini_out = intbuf( 7)
    outblock%nprocess_bd_out  = intbuf( 8)
    outblock%outsteps         = intbuf( 9)
    outblock%kepin            = intbuf(10)
    outblock%kezin            = intbuf(11)
    outblock%n_num            = intbuf(12)
    outblock%nunit_of_time    = intbuf(13)
    outblock%itype_vertint    = intbuf(14)
   
    outblock%slon     = realbuf(1)
    outblock%slat     = realbuf(2)
    outblock%elon     = realbuf(3)
    outblock%elat     = realbuf(4)
!   DO i=1,nlevels
!     outblock%plev(i)  = realbuf(4+i)
!   ENDDO
!   DO i=1,nlevels
!     outblock%zlev(i)  = realbuf(4+nlevels+i)
!   ENDDO


    outblock%ysystem(1:LEN_TRIM(outblock%ysystem))  = charbuf(1)
    outblock%ydbtype(1:LEN_TRIM(outblock%ydbtype))  = charbuf(2)
    outblock%ydir   (  1:100)                       = charbuf(3)(1:100)
    outblock%ydir   (101:200)                       = charbuf(4)(1:100)
    outblock%ydir   (201:250)                       = charbuf(5)(1: 50)
    outblock%ysuffix(1:LEN_TRIM(outblock%ysuffix))  = charbuf(6)
    outblock%ytunit (1:LEN_TRIM(outblock%ytunit))   = charbuf(7)
    outblock%ydomain(1:LEN_TRIM(outblock%ydomain))  = charbuf(8)
    outblock%yform_write(1:4)                       = charbuf(9)(1:4)
    DO i=1,nzmxml
      outblock%yvarml(i)(1:10) = charbuf(9+i)(1:10)
    ENDDO
    DO i=1,nzmxpl
      outblock%yvarpl(i)(1:10) = charbuf(9+nzmxml+i)(1:10)
    ENDDO
    DO i=1,nzmxzl
      outblock%yvarzl(i)(1:10) = charbuf(9+nzmxml+nzmxpl+i)(1:10)
    ENDDO
    DO i=1,nzmxzl
      outblock%yvarsl(i)(1:10) = charbuf(9+nzmxml+nzmxpl+nzmxzl+i)(1:10)
    ENDDO
    DO i=1,nzmxc
      outblock%yvarc (i)(1:10) = charbuf(9+nzmxml+nzmxpl+2*nzmxzl+i)(1:10)
    ENDDO

    outblock%lcheck         = logbuf(1)
    outblock%lwrite_const   = logbuf(2)
    outblock%lanalysis      = logbuf(3)
    outblock%l_z_filter     = logbuf(4)
    outblock%l_p_filter     = logbuf(5)
    outblock%l_fi_ps_smooth = logbuf(6)
    outblock%luvmasspoint   = logbuf(7)
    outblock%lhour          = logbuf(8)
  ENDIF
ENDIF

! look for special output variables:
DO n = 1, outblock%nyvar_m

  SELECT CASE (TRIM(outblock%yvarml(n)))
  CASE ('VORTIC1')
    l_dzeta_d_needed = .TRUE.
  CASE ('VORTIC2')
    l_dzeta_d_needed = .TRUE.
  CASE ('VORTIC3')
    l_dzeta_d_needed = .TRUE.
  CASE ('POT_VORTIC')
    l_dzeta_d_needed = .TRUE.
  END SELECT

ENDDO

! the list of output steps has to be treated separately
IF (nproc > 1) THEN
  IF (my_world_id /= 0) THEN
    ! Allocate the vector for the output steps in the namelist structure
    ALLOCATE(outblock%ngrib(outblock%outsteps), STAT=istat)
    IF(istat /= 0) THEN
      yerrmsg= 'Error allocating ngrib for namelist structure'
      CALL model_abort(my_world_id, istat, yerrmsg, yroutine)
    ENDIF
  ENDIF

  CALL distribute_values (outblock%ngrib, outblock%outsteps, 0, imp_integers,&
                                        icomm_world, ierr)
ENDIF

! determine noutst as maximum of all outblock%outsteps:
noutst = MAX (noutst, outblock%outsteps)

! determine the first output step
IF (nstart == 0) THEN
  outblock%nextstep = 1
ELSEIF (nstart >= outblock%ngrib(outblock%outsteps)) THEN
  ! then nothing has to be written for this namelist group in the restart run
  outblock%nextstep = nstop + 1
ELSE
  outblock%nextstep = 1
  ! compare with "< nstart" now; this was "<=" before, but the last output 
  ! step will now be written again after a restart
  DO WHILE (outblock%ngrib(outblock%nextstep) <  nstart)
    outblock%nextstep = outblock%nextstep + 1
  ENDDO
ENDIF

!------------------------------------------------------------------------------
!- Section 6: Calculation of specifications for output domain
!------------------------------------------------------------------------------

! From the slon/slat-elon/elat values the grid point specifications are 
! derived and stored in the namelist structure 
! This has to be done by every processor

! Initialize these variables to a reasonable value to avoid problems
! in case of small numerical discrepancies
outblock%i_out_start = 1
outblock%i_out_end = ie_tot
outblock%j_out_start = 1
outblock%j_out_end = je_tot

lon_loop: DO i=1,ie_tot-1
  zlon1 = startlon_tot + (i - 1) * dlon
  IF (zlon1 > 180.0_ireals) THEN
    zlon1 = zlon1 - 360.0_ireals
  ENDIF

  zlon2 = startlon_tot + (i    ) * dlon
  IF (zlon2 > 180.0_ireals) THEN
    zlon2 = zlon2 - 360.0_ireals
  ENDIF

  IF (zlon1 > zlon2) THEN
    ! this is around the 180-Meridian (could happen on southern hemisphere)
    IF (zlon1 <= outblock%slon) THEN
      ! Then we take zlon1 as starting point
      outblock%i_out_start = i
    ELSE
      IF ( (outblock%slon < zlon2) .AND. (outblock%slon < 0) ) THEN
        ! then we also take zlon1 as starting point
        outblock%i_out_start = i
      ENDIF
    ENDIF

    IF (zlon1 < outblock%elon) THEN
      ! Then we take zlon2 as ending point
      outblock%i_out_end   = i+1
      EXIT lon_loop
    ELSE
      IF (outblock%elon <= zlon2) THEN
        ! then we also take zlon2 as ending point
        outblock%i_out_end   = i+1
        EXIT lon_loop
      ENDIF
    ENDIF
  ELSE

    ! this is the normal case
    IF ( (zlon1 <= outblock%slon) .AND. (outblock%slon <  zlon2) ) THEN
      ! slon is in this interval; take zlon1 as starting point
      outblock%i_out_start = i
    ENDIF

    IF ( (zlon1 <  outblock%elon) .AND. (outblock%elon <= zlon2) ) THEN
      ! elon is in this interval; take zlon2 as ending point
      outblock%i_out_end   = i+1
      EXIT lon_loop
    ENDIF
  ENDIF
ENDDO lon_loop

lat_loop: DO j=1,je_tot-1
  zlat1 = startlat_tot + (j - 1) * dlat
  zlat2 = startlat_tot + (j    ) * dlat

  IF ( (zlat1 <= outblock%slat) .AND. (outblock%slat <  zlat2) ) THEN
    ! slat is in this interval; take zlat1 as starting point
    outblock%j_out_start = j
  ENDIF

  IF ( (zlat1 <  outblock%elat) .AND. (outblock%elat <= zlat2) ) THEN
    ! elat is in this interval; take zlat2 as ending point
    outblock%j_out_end   = j+1
    EXIT lat_loop
  ENDIF
ENDDO lat_loop


outblock%ie_out_tot = outblock%i_out_end - outblock%i_out_start + 1
outblock%je_out_tot = outblock%j_out_end - outblock%j_out_start + 1

!------------------------------------------------------------------------------
!- Section 7: Output of the namelist variables and their default values
!------------------------------------------------------------------------------

IF (my_world_id == 0) THEN

  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(A23)') '0     NAMELIST: gribout'
  WRITE (nuspecif, '(A23)') '      -----------------'
  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(T7,A,T21,A,T39,A,T58,A)') 'Variable', 'Actual Value',   &
                                               'Default Value', 'Format'

  WRITE(nuspecif,*)

  WRITE(nuspecif,'(T7,A)') 'hcomb triplets 1'
  WRITE(nuspecif,'(T7,A,T21,F6.2,T42,F6.2,T58,A)') 'begin',hcomb(1),hcomb(1),&
                                                   'R'
  WRITE(nuspecif,'(T7,A,T15,F12.2,T36,F12.2,T58,A)') 'end',hcomb(2),hcomb(2),&
                                                   'R'
  WRITE(nuspecif,'(T7,A,T21,F6.2,T42,F6.2,T58,A)') 'incr', hcomb(3),hcomb(3),&
                                                   'R'
  DO i=4,ntrip,3
     IF( hcomb(i) /= r_huge) THEN
        WRITE(nuspecif,*)
        WRITE(nuspecif,'(T7,A,I2)') 'hcomb triplets', i/3+1
        WRITE(nuspecif,'(T7,A,T21,F6.2,T58,A)') 'begin',hcomb(i),  'R'
        WRITE(nuspecif,'(T7,A,T21,F6.2,T58,A)') 'end',  hcomb(i+1),'R'
        WRITE(nuspecif,'(T7,A,T21,F6.2,T58,A)') 'incr', hcomb(i+2),'R'
     ENDIF
  ENDDO

  WRITE(nuspecif,*)

  DO i=1,ntrip,3
     IF( ncomb(i) /= i_huge) THEN
        WRITE(nuspecif, '(T7,A,I2)') 'ncomb triplets',i/3+1
        WRITE(nuspecif, '(T7,A,T21,I6,T58,A)') 'begin', ncomb(i),  'I'
        WRITE(nuspecif, '(T7,A,T21,I6,T58,A)') 'end'  , ncomb(i+1),'I'
        WRITE(nuspecif, '(T7,A,T21,I6,T58,A)') 'incr' , ncomb(i+2),'I'
     ENDIF
  ENDDO

  WRITE(nuspecif,*)
  WRITE(nuspecif,*)

! UB>> 100 -> noutst_max
  DO i=1, noutst_max
     IF( hgrib(i) /= r_huge) THEN
       WRITE(nuspecif, '(T7,A,I2,A,F6.2)') 'hgrib(',i,') = ',hgrib(i)
     ENDIF
  ENDDO

  WRITE(nuspecif,*)

! UB>> 100 -> noutst_max
  DO i=1, noutst_max
     IF( ngrib(i) /= i_huge) THEN
       WRITE(nuspecif, '(T7,A,I2,A,I6)') 'ngrib(',i,') = ',ngrib(i)
     ENDIF
  ENDDO

  WRITE(nuspecif,*)

  WRITE(nuspecif,'(T7,A,T21,F10.4,T39,F10.4,T58,A)') 'slon',outblock%slon,&
                                          startlon_tot,'R'
  WRITE(nuspecif,'(T7,A,T21,F10.4,T39,F10.4,T58,A)') 'slat',outblock%slat,&
                                          startlat_tot,'R'
  WRITE(nuspecif,'(T7,A,T21,F10.4,T39,F10.4,T58,A)') 'elon',outblock%elon,&
                                                                elon_d, 'R'
  WRITE(nuspecif,'(T7,A,T21,F10.4,T39,F10.4,T58,A)') 'elat',outblock%elat,&
                                                                elat_d, 'R'
  WRITE(nuspecif,*)

  WRITE(nuspecif,'(T7,A,T21,A,T39,A,T58,A)') 'yform_write',outblock%yform_write,  &
                                              yform_write_d,'C4'
  WRITE(nuspecif,'(T7,A,T21,A,T39,A,T58,A)') 'ysystem',outblock%ysystem,  &
                                              ysystem_d,'C4'
  WRITE(nuspecif,'(T7,A,T21,A,T39,A,T58,A)') 'ydbtype',outblock%ydbtype,  &
                                              ydbtype_d,'C5'
  WRITE(nuspecif,'(T7,A,T21,A,T39,A,T58,A)') 'ytunit' ,outblock%ytunit,   &
                                              ytunit_d,'C'
  WRITE(nuspecif,'(T7,A,T21,A,T39,A,T58,A)') 'ydomain',outblock%ydomain,  &
                                              ydomain_d,'C'
  WRITE(nuspecif,*)

  WRITE(nuspecif,'(T7,A,T25,I2,T46,I2,T58,A)') 'n_num',outblock%n_num,    &
                                                n_num_d,'I'
  WRITE(nuspecif,'(T7,A,T23,I6,T39,I6,T58,A)')                            &
                                       'nrbit',outblock%nrbit,nrbit_d,'I'
  WRITE(nuspecif,'(T7,A,T23,I7,T39,I7,T58,A)')                            &
          'nprocess_ini_out',outblock%nprocess_ini_out,nprocess_ini_d,'I'
  WRITE(nuspecif,'(T7,A,T23,I7,T39,I7,T58,A)')                            &
             'nprocess_bd_out',outblock%nprocess_bd_out,nprocess_bd_d,'I'
  WRITE(nuspecif,'(T7,A,T23,I7,T39,I7,T58,A)')                            &
             'nunit_of_time', outblock%nunit_of_time, nunit_of_time_d,'I'
  WRITE(nuspecif,'(T7,A,T23,L6,T39,L6,T58,A)')                            &
                                  'lcheck' ,outblock%lcheck, lcheck_d,'L'
  WRITE(nuspecif,'(T7,A,T21,L6,T39,L6,T58,A)')                            &
                'lwrite_const' ,outblock%lwrite_const, lwrite_const_d,'L'
  WRITE(nuspecif,'(T7,A,T21,L6,T39,L6,T58,A)')                            &
                         'lanalysis' ,outblock%lanalysis, lanalysis_d,'L'
  WRITE(nuspecif,'(T7,A,T21,L6,T39,L6,T58,A)')                            &
                         'l_z_filter' ,outblock%l_z_filter, l_z_filter_d,'L'
  WRITE(nuspecif,'(T7,A,T21,L6,T39,L6,T58,A)')                            &
                         'l_p_filter' ,outblock%l_p_filter, l_p_filter_d,'L'
  WRITE(nuspecif,'(T7,A,T21,L6,T39,L6,T58,A)')                            &
             'l_fi_ps_smooth' ,outblock%l_fi_ps_smooth, l_fi_ps_smooth_d,'L'
  WRITE(nuspecif,'(T7,A,T21,L6,T39,L6,T58,A)')                            &
                   'luvmasspoint' ,outblock%luvmasspoint, luvmasspoint_d,'L'
  WRITE(nuspecif,'(T7,A,T23,I6,T39,I6,T58,A)')                            &
                  'itype_vertint',outblock%itype_vertint,itype_vertint_d,'I'
  WRITE(nuspecif,*)

  WRITE(nuspecif,'(T8,A,T26,  A                  )')                      &
              'ydir    (C*250)   ',outblock%ydir(1:LEN_TRIM(outblock%ydir))
  WRITE(nuspecif,'(T8,A,T26,  A                  )')                      &
        'ysuffix (C*25)    ',outblock%ysuffix(1:LEN_TRIM(outblock%ysuffix))
  IF(outblock%nyvar_m >= 1) THEN
     WRITE(nuspecif,*)
     WRITE(nuspecif, '(A)')  'Variables for result data'
     DO i=1,outblock%nyvar_m
        WRITE (nuspecif, '(T7,A,I3,A,T21,A10,T58,A)') 'yvarml(',i,')',&
                                                  outblock%yvarml(i),'C10'
     ENDDO
  ENDIF

  IF(outblock%nyvar_c >= 1) THEN
     WRITE(nuspecif,*)
     WRITE(nuspecif, '(A)')  'Constant variables'
     DO i=1,outblock%nyvar_c
        WRITE (nuspecif, '(T7,A,I3,A,T21,A10,T58,A)') 'yvarc (',i,')',&
                                                  outblock%yvarc (i),'C10'
     ENDDO
  ENDIF

  IF(outblock%nyvar_p >= 1) THEN
     WRITE(nuspecif,*)
     WRITE(nuspecif, '(A)')  'Variables for p-interpolation'
     DO i=1,outblock%nyvar_p
        WRITE (nuspecif, '(T7,A,I3,A,T21,A10,T58,A)') 'yvarpl(',i,')',&
                                                  outblock%yvarpl(i),'C10'
     ENDDO
     WRITE(nuspecif,*)
     WRITE(nuspecif, '(A)')  'Pressure levels'
     DO i=1,outblock%kepin
        WRITE (nuspecif, '(T7,A,I3,A,T21,F8.1,T58,A)') 'plev(',i,')',&
                                                  outblock%plev(i)/100,'R'
     ENDDO
     WRITE(nuspecif,*)
  ENDIF

  IF(outblock%nyvar_z >= 1) THEN
     WRITE(nuspecif,*)
     WRITE(nuspecif, '(A)')  'Variables for z-interpolation'
     DO i=1,outblock%nyvar_z
        WRITE (nuspecif, '(T7,A,I3,A,T21,A10,T58,A)') 'yvarzl(',i,')',&
                                                  outblock%yvarzl(i),'C10'
     ENDDO
     WRITE(nuspecif,*)
     WRITE(nuspecif, '(A)')  'Geom. Height levels'
     DO i=1,outblock%kezin
        WRITE (nuspecif, '(T7,A,I3,A,T21,F8.1,T58,A)') 'zlev(',i,')',&
                                                  outblock%zlev(i),'R'
     ENDDO
     WRITE(nuspecif,*)
  ENDIF

  IF(outblock%nyvar_s >= 1) THEN
     WRITE(nuspecif,*)
     WRITE(nuspecif, '(A)')  'Variables for Satellites'
     DO i=1,outblock%nyvar_s
        WRITE (nuspecif, '(T7,A,I3,A,T21,A10,T58,A)') 'yvarsl(',i,')',&
                                                  outblock%yvarsl(i),'C10'
     ENDDO
     WRITE(nuspecif,*)
  ENDIF
ENDIF

END SUBROUTINE input_gribout

!==============================================================================
!+ Calculates the list of time steps with output of this special group
!------------------------------------------------------------------------------

SUBROUTINE calc_ngrib (hour, nstep, hour_comb, nstep_comb, outblock,     &
                       grid_dt, ierr)

!------------------------------------------------------------------------------
!
! Description:
!   This routine calculates the list of time steps when Grib output should
!   be performed for every namelist group "gribout". Input parameters are
!   the list of time steps (in hours or time steps) and the list of triplets
!   (in hours or time steps), whichever the user has set. calc_ngrib builds
!   the union of all four lists and determines the list of time steps ngrib
!   for the special pointer.
!
! Method:
!   The time steps when output of this group has to be performed are
!   determined  and put into a temporary vector. The steps are sorted
!   using the bubblesort algorithm, entries that appear more than once
!   are eliminated. The final list is put to ngrib of the pointer variable.
!
!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
! UB>> 100 -> noutst_max
REAL    (KIND=ireals)   , INTENT(IN)    :: hour(noutst_max),     hour_comb(ntrip)
INTEGER (KIND=iintegers), INTENT(IN)    :: nstep(noutst_max),    nstep_comb(ntrip)
REAL    (KIND=ireals)   , INTENT(IN)    :: grid_dt

! pointer to the actual namelist group
TYPE(pp_nl), POINTER                    :: outblock        

! Scalar arguments with intent(inout):
INTEGER (KIND=iintegers), INTENT(INOUT) :: ierr

!------------------------------------------------------------------------------
!
! Local scalars:
REAL    (KIND=ireals)                 :: r_huge, hgo
INTEGER (KIND=iintegers)              :: i, ii, tmp, i_huge, ncount, istat
INTEGER (KIND=iintegers)              :: nc_hour, nc_hour_comb, &
                                         nc_step, nc_step_comb
LOGICAL                               :: lchange, ldouble

CHARACTER(LEN=25),PARAMETER           :: yroutine='calc_ngrib'
CHARACTER(LEN=80)                     :: yerrmsg

!
! Local arrays:
INTEGER (KIND=iintegers), ALLOCATABLE :: ngrib(:)
! UB>> 100 -> noutst_max
INTEGER (KIND=iintegers)              :: nhtos(noutst_max), nhtos_comb(ntrip)
INTEGER (KIND=iintegers)              :: ninc_h, ninc_2, ninc_4, noutsteps_nl

!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
!- Section 1: Initialization with the defaults
!------------------------------------------------------------------------------

  i_huge = HUGE (1_iintegers)
  r_huge = HUGE (1.0_ireals)

  ! Determine, how many entries are in each list
  nc_hour      = COUNT(hour(:)       /= r_huge)
  nc_hour_comb = COUNT(hour_comb(:)  /= r_huge)
  nc_step      = COUNT(nstep(:)      /= i_huge)
  nc_step_comb = COUNT(nstep_comb(:) /= i_huge)

!------------------------------------------------------------------------------
!- Section 2: Transform values given in hours into time steps
!------------------------------------------------------------------------------

  WHERE(hour_comb(:) /= r_huge)
    nhtos_comb(:) =  NINT(hour_comb(:) * 3600.0_ireals / grid_dt)
  ENDWHERE

  WHERE(hour(:) /= r_huge)
    nhtos(:) =  NINT(hour(:) * 3600.0_ireals / grid_dt)
  ENDWHERE

!------------------------------------------------------------------------------
!- Section 3: Calculate the total number of output steps
!------------------------------------------------------------------------------

  ncount = 0

  ! For calculating the total number of output steps, we can still use
  ! the time step increment
  DO i = 1,nc_hour_comb,3
    IF (nhtos_comb(i+1) < nhtos_comb(i)) THEN
      PRINT *, 'end timestep < begin timestep in hour triplet ',i/3+1
      ierr = 1
      RETURN
    ENDIF
    IF (nhtos_comb(i) + nhtos_comb(i+2) > nhtos_comb(i+1)) THEN
      PRINT *, 'error in increment in hour triplet ',i/3+1
      ierr = 2
      RETURN
    ENDIF
    IF (nhtos_comb(i+2) == 0) THEN
      PRINT *, 'error in increment in hour triplet ',i/3+1,': inc = 0'
      ierr = 3
      RETURN
    ENDIF
! UB<<
    IF (hour_comb(i+2) <=  1.0_ireals) THEN
      ncount = ncount + CEILING((hour_comb(i+1)-hour_comb(i)) / hour_comb(i+2)) + 1
    ELSE
      ncount = ncount + CEILING((hour_comb(i+1)-hour_comb(i)))
    ENDIF
! UB<<    
!!$    IF     (hour_comb(i+2) ==  0.5_ireals) THEN
!!$      ncount = ncount + 2 * (nhtos_comb(i+1)-nhtos_comb(i)) / (3600.0_ireals / grid_dt) + 1
!!$    ELSEIF (hour_comb(i+2) == 0.25_ireals) THEN
!!$      ncount = ncount + 4 * (nhtos_comb(i+1)-nhtos_comb(i)) / (3600.0_ireals / grid_dt) + 1
!!$    ELSE
!!$      ncount = ncount + NINT((nhtos_comb(i+1)-nhtos_comb(i)) / (3600.0_ireals / grid_dt) + 1)
!!$    ENDIF
  ENDDO

  DO i = 1,nc_step_comb,3
    IF (nstep_comb(i+1) < nstep_comb(i)) THEN
      PRINT *, 'end timestep < begin timestep in step triplet ',i/3+1
      ierr = 1
      RETURN
    ENDIF
    IF (nstep_comb(i) + nstep_comb(i+2) > nstep_comb(i+1)) THEN
      PRINT *, 'error in increment in step triplet ',i/3+1
      ierr = 2
      RETURN
    ENDIF
    IF (nstep_comb(i+2) == 0) THEN
      PRINT *, 'error in increment in step triplet ',i/3+1,': inc = 0'
      ierr = 3
      RETURN
    ENDIF
    ncount = ncount + (nstep_comb(i+1)-nstep_comb(i)) / nstep_comb(i+2) + 1
  ENDDO
  ncount = ncount + nc_hour + nc_step

!------------------------------------------------------------------------------
!- Section 4: Allocate ngrib and put all output steps into this vector
!------------------------------------------------------------------------------

  ! Allocate the intermediate vector for the output steps
  ALLOCATE(ngrib(ncount), STAT=istat)
  IF(istat /= 0) THEN
    yerrmsg= 'Error allocating ngrib'
    CALL model_abort(my_world_id, istat, yerrmsg, yroutine)
  ENDIF

  ! Put all output steps into the vector ngrib
  noutsteps_nl = 0

  ! With a flexible time step, that does not fit into a full hour, 
  ! the time step increments can not be used any more.

  DO i = 1,nc_hour_comb,3
    hgo = hour_comb(i)
    DO WHILE (hgo <= hour_comb(i+1))
      noutsteps_nl = noutsteps_nl + 1
      ngrib(noutsteps_nl) = NINT (hgo * 3600.0_ireals / grid_dt)
      hgo = hgo + hour_comb(i+2)
    ENDDO
  ENDDO

  DO i = 1,nc_step_comb,3
    DO ii = nstep_comb(i), nstep_comb(i+1), nstep_comb(i+2)
      noutsteps_nl = noutsteps_nl + 1
      ngrib(noutsteps_nl) = ii
    ENDDO
  ENDDO
  DO i = 1,nc_hour
    noutsteps_nl = noutsteps_nl + 1
    ngrib(noutsteps_nl) = nhtos(i)
  ENDDO
  DO i = 1,nc_step
    noutsteps_nl = noutsteps_nl + 1
    ngrib(noutsteps_nl) = nstep(i)
  ENDDO

!------------------------------------------------------------------------------
!- Section 5: Sort the vector with bubblesort and eliminate identical entries
!------------------------------------------------------------------------------

  DO   ! bubblesort
    lchange = .FALSE.
    DO i=1,noutsteps_nl-1
      IF(ngrib(i) > ngrib(i+1)) THEN
        tmp        = ngrib(i+1)
        ngrib(i+1) = ngrib(i)
        ngrib(i)   = tmp
        lchange = .TRUE.
      ENDIF
    ENDDO
    IF (.NOT. lchange) EXIT
  ENDDO

  ! eliminate identical entries
  DO
    ldouble = .FALSE.
    DO i=1,noutsteps_nl-1
      IF(ngrib(i) == ngrib(i+1)) THEN
        ngrib(i:noutsteps_nl-1) = ngrib(i+1:noutsteps_nl)
        ldouble = .TRUE.
        EXIT
      ENDIF
    ENDDO
    IF (.NOT. ldouble) EXIT
    noutsteps_nl = noutsteps_nl-1
  ENDDO

!------------------------------------------------------------------------------
!- Section 6: Put ngrib to ngrib_nl
!------------------------------------------------------------------------------

  ! Allocate the vector for the output steps in the namelist structure
  outblock%outsteps = noutsteps_nl
  ALLOCATE(outblock%ngrib(noutsteps_nl), STAT=istat)
  IF(istat /= 0) THEN
    yerrmsg= 'Error allocating ngrib for namelist structure'
    CALL model_abort(my_world_id, istat, yerrmsg, yroutine)
  ENDIF

  outblock%ngrib(1:noutsteps_nl) = ngrib(1:noutsteps_nl)

  IF (izdebug >= 10) THEN
    PRINT *, '     Output will be done in steps:  ', ngrib(1:noutsteps_nl)
  ENDIF

  ! deallocate intermediate vector
  DEALLOCATE(ngrib, STAT=istat)
  IF(istat /= 0) THEN
    yerrmsg= 'Error de-allocating ngrib for namelist structure'
    CALL model_abort(my_world_id, istat, yerrmsg, yroutine)
  ENDIF

END SUBROUTINE calc_ngrib

!==============================================================================

!------------------------------------------------------------------------------
! End of external subroutine organize_data
!------------------------------------------------------------------------------

END SUBROUTINE organize_data
