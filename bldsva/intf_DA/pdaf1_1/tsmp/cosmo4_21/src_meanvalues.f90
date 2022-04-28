!+ Source module for computations of meanvalues for printout
!------------------------------------------------------------------------------

MODULE src_meanvalues 

!------------------------------------------------------------------------------
!
! Description:
!   This module provides routines for computing mean values of some
!   variables of the LM. The mass and energy related variables are printed
!   to a file YUPRMASS, the humidity related variables are printed to a file
!   YUPRHUMI.
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
! 1.5        1998/06/29 Guenther Doms
!  Additional printout of maximum absolute horizontal wind speed to YUPRMEAN
! 1.8        1998/08/03 Ulrich Schaettler
!  Eliminated dependency from data_io by introducing a parameterlist for
!  init_meanvalues.
! 1.10       1998/09/29 Ulrich Schaettler
!  Depending on llm, zwa300 and zwa500 are set to 0.0 in routine meanvalues.
! 1.20       1999/01/07 Guenther Doms
!  Renaming of some global variables
! 1.29       1999/05/11 Ulrich Schaettler
!  Adapted interfaces to utility-modules and prepared use of MPE_IO
! 1.34       1999/12/10 Ulrich Schaettler
!  Changed calls to timing routines
! 1.39       2000/05/03 Ulrich Schaettler
!  Included declaration of variables from former module data_diagnostics
!  or module data_runcontrol, resp.
! 2.8        2001/07/06 Ulrich Schaettler
!  Use new variables nvers for ASCII output
! 3.7        2004/02/18 Ulrich Schaettler
!  Changed treatment of unit-numbers for ASCII-file handling
! 3.13       2004/12/03 Thorsten Reinhardt
!  Adaptations for new Graupel scheme
! 3.15       2005/03/03 Ulrich Schaettler
!  Open and close file YUPRMEAN for every output step (to save the results)
!  Replaced FLOAT by REAL
! 3.18       2006/03/03 Ulrich Schaettler
!  Changes to allow the reproduction of meanvalues for restart runs
!  (Save some initial meanvalues to the first restart file)
! 3.21       2006/12/04 Klaus Stephan, Michael Baldauf, Ulrich Schaettler
!  Mean values for qi, qr, qs, qg
!  Splitted output for PRMEAN to PRMASS and PRHUMI
! V3_23        2007/03/30 Ulrich Schaettler
!  Introduced switch ldump_ascii to flush the ASCII files
! V4_3         2008/02/25 Ulrich Schaettler
!  Eliminated 'LM' from ASCII output
! V4_4         2008/07/16 Ulrich Schaettler
!  Changed NL parameter lyear_360 to itype_calendar, to have several options
! V4_9         2009/07/16 Michael Baldauf, Hans-Juergen Panitz
!  Editorial changes
!  Increase length of character variable ychh for long simulations
!  Write timing information independent of ltime
!  Combined some loops and inserted compiler directives
! V4_10        2009/09/11 Christian Bollmann
!  Added compiler directive to use option _on_adb for NEC
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_19        2011/08/01 Oliver Fuhrer
!  Bug fix for closing files with meanvalues
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

USE data_parameters, ONLY :   &
    ireals,    & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
    ie_tot,       & ! number of grid points in zonal direction
    je_tot,       & ! number of grid points in meridional direction
    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction
    ke1,          & ! KE+1

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the 
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from 
!    the other ones because of the use of the staggered Arakawa-B-grid.
!    
!   zonal direction
    istart,       & ! start index for the forecast of w, t, qv, qc and pp
    iend,         & ! end index for the forecast of w, t, qv, qc and pp
    istartu,      & ! start index for the forecast of u
    iendu,        & ! end index for the forecast of u

!   meridional direction
    jstart,       & ! start index for the forecast of w, t, qv, qc and pp
    jend,         & ! end index for the forecast of w, t, qv, qc and pp
    jstartv,      & ! start index for the forecast of v
    jendv,        & ! end index for the forecast of v

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------
    dt,           & ! long time-step

! 7. Layer index corresponding to a specified pressure level
! ----------------------------------------------------------
    klv850,       & ! k index of the LM-mainlevel, on 850 HPa
    klv500,       & ! k index of the LM-mainlevel, on 500 HPa
    klv300          ! k index of the LM-mainlevel, on 300 HPa

!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 2. physical constants and related variables
! -------------------------------------------
    cp_d,         & ! specific heat of dry air at constant pressure
    lh_v,         & ! latent heat of vapourization
    day_len,      & ! mean length of the day (s)
    g               ! acceleration due to gravity

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
    p0 ,            & ! reference pressure at main levels             ( pa  )
    dp0,            & ! base state pressure thik of levels        ( pa  )
    hhl,            & ! height of model half levels                   ( m )

! 3. prognostic variables                                             (unit)
! -----------------------
    u,              & ! zonal wind speed                              ( m/s )
    v,              & ! meridional wind speed                         ( m/s )
    w,              & ! vertical wind speed (defined on half levels)  ( m/s )
    t,              & ! temperature                                   (  k  )
    qv,             & ! specific water vapor content                  (kg/kg)
    qc,             & ! specific cloud water content                  (kg/kg)
    qi,             & ! specific cloud ice content                    (kg/kg)
    qr,             & ! specific cloud rain content                   (kg/kg)
    qs,             & ! specific cloud snow content                   (kg/kg)
    qg,             & ! specific cloud graupel content                (kg/kg)
    pp,             & ! deviation from the reference pressure         ( pa  )

! 5. fields for surface values and soil model variables               (unit )
! -----------------------------------------------------
    ps ,            & ! surface pressure                              ( pa  )

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------
    prr_con    ,    & ! precipitation rate of rain, convective        (kg/m2*s)
    prs_con    ,    & ! precipitation rate of snow, convective        (kg/m2*s)
    prr_gsp    ,    & ! precipitation rate of rain, grid-scale        (kg/m2*s)
    prs_gsp    ,    & ! precipitation rate of snow, grid-scale        (kg/m2*s)
    prg_gsp    ,    & ! precipitation rate of graupel, grid-scale     (kg/m2*s)
    rain_gsp   ,    & ! amount of rain from grid-scale precip. (sum)  (kg/m2)
    snow_gsp   ,    & ! amount of snow from grid-scale precip. (sum)  (kg/m2)
    grau_gsp   ,    & ! amount of graupel from grid-scale prec. (sum) (kg/m2)
    rain_con   ,    & ! amount of rain from convective precip. (sum)  (kg/m2)
    snow_con   ,    & ! amount of snow from convective precip. (sum)  (kg/m2)
    dpsdt             ! tendency of the surface pressure              ( Pa/s)

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    nstart,       & ! first time step of the forecast
    nstop,        & ! last time step of the forecast
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nnow,         & ! corresponds to ntstep
    nnew,         & ! corresponds to ntstep + 1

! 3. controlling the physics
! --------------------------   &
    itype_gscp,   & ! type of microphys. parametrization

! 7. additional control variables
! -------------------------------
    llm,          & ! if .TRUE., running with a lowered upper boundary
    nvers,        & ! version number of experiment for documentation
    ltime,        & ! the results are reproducible in parallel mode
    ldump_ascii,  & ! for flushing (close and re-open) the ASCII files

! 8. diagnostic calculations
! --------------------------
    psm0,   & ! initial value for mean surface pressure ps
    dsem0,  & ! initial value for mean dry static energy
    msem0,  & ! initial value for mean moist static energy
    kem0,   & ! initial value for mean kinetic energy
    qcm0,   & ! initial value for mean cloudwater content

! 9. Variables for Ascii file handling, time measuring, ...
! ---------------------------------------------------------
    itype_calendar,&! for specifying the calendar used

! 12. time measuring
! ------------------
    yakdat2         ! actual date (yadat+ntstep/dt) 
                    ! wd dd.mm.yy (weekday, day, month, year)

! end of data_runcontrol 
!------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
    num_compute,     & ! number of compute PEs
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    my_cart_id,      & ! rank of this subdomain in the global communicator
    icomm_cart,      & ! communicator that belongs to the cartesian grid
    imp_reals          ! determines the correct REAL type used in the model
                       ! for MPI

!------------------------------------------------------------------------------

USE utilities,          ONLY :  &
      elapsed_time,        & ! returns elapsed wall-clock time in seconds
      get_utc_date           ! actual date of the forecast in different forms

!------------------------------------------------------------------------------

USE environment,              ONLY :  &
      model_abort,            & ! aborts the program in case of errors
      comm_barrier              ! 

!------------------------------------------------------------------------------

USE parallel_utilities,       ONLY :  &
      global_values             ! collects values from all nodes (e.g. for
                                ! calculating mean values)

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Variables for module src_meanvalues
  REAL    (KIND=ireals)                  ::             &
    qim0,   & ! initial value for mean cloud ice  content
    qrm0,   & ! initial value for mean rain       content
    qsm0,   & ! initial value for mean snow       content
    qgm0      ! initial value for mean graupel    content

  REAL    (KIND=ireals)                  ::             &
    realtot_new,  & ! total real time since last print of the time
    realtot_old     ! total real time up to last print of the time

  INTEGER (KIND=iintegers)                 ::             &
    nlines    ! counts the printed lines (for adjusting output to a
              ! computer page)

  CHARACTER (LEN=22)    yinidate  ! date when the forecast begins

  INTEGER (KIND=iintegers)     ::    &
    nuprmass,   & ! mean-values for mass and energy fields
    nuprhumi      ! mean-values for humidity fields

  CHARACTER (LEN= 8)  ::   &
    yuprmass='YUPRMASS',   & ! mean-values for mass and energy fields
    yuprhumi='YUPRHUMI'      ! mean-values for humidity fields

  INTEGER   (KIND=iintegers)       ::           &
    n0meanval,    & ! time step of the first output
    nincmeanval,  & ! time step increment of outputs
    nextmeanval     ! next time step for mean value output

!==============================================================================

! Public and Private Subroutines

!==============================================================================
! Module procedures
!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure in "diagnostics" for initializing variables
!------------------------------------------------------------------------------

SUBROUTINE init_meanvalues (ydate_ini)

!------------------------------------------------------------------------------
!
! Description:
!   This routine performs the initializations of some variables to calculate  
!   the domain averages of prognostic and diagnostic variables.
!   Tasks are opening of files, allocation of space and computing initial 
!   values. Most of the actions are performed only in node 0
!
! Method:
!   
!------------------------------------------------------------------------------

! Parameterlist
CHARACTER (LEN=10), INTENT(IN)     ::   &
     ydate_ini  ! start of the forecast  yyyymmddhh (year, month, day, hour)
 
! Local variables

! support variables for computation
REAL (KIND=ireals)         ::       &
  zdsem, zmsem, zkem, zqcm, zke, zdse, zmse, zflar, zpsm, zrdzm, zanhh,   &
  zrealdiff, zqim, zqrm, zqsm, zqgm

REAL (KIND=ireals)         ::       &
  zrdz     (ie,je,ke),  & ! rho0*dz, factor for computing volume averages   
  realbuf  (10)           ! for communication

INTEGER (KIND=iintegers)   ::       &
  i,j,k,nzanjata,       & ! support variables
  nstat, nzerror

CHARACTER (LEN= 8)         :: ydatearg, ytime
CHARACTER (LEN=10)         :: ydate, ytimearg, yandat1
CHARACTER (LEN=10)         :: ychh
CHARACTER (LEN=20)         :: yroutine
CHARACTER (LEN=25)         :: yerrmsg 

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE init_meanvalues
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  Section 1: Initializations
!------------------------------------------------------------------------------

  yroutine = 'init_meanvalues'
  nzerror  = 0

!------------------------------------------------------------------------------
!  Section 2: Initializations for calculating mean values
!------------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !  Section 2.1: Definition and opening of the file
  !----------------------------------------------------------------------------

  IF (my_cart_id == 0) THEN
    ! file for mass and energy fields
    OPEN (nuprmass,FILE=yuprmass,FORM='FORMATTED',STATUS='REPLACE',IOSTAT=nstat)
    IF (nstat /= 0) THEN
        yerrmsg = 'opening of file YUPRMASS failed'
        CALL model_abort (my_cart_id, 7003, yerrmsg, yroutine)
    ENDIF
    REWIND nuprmass

    ! file for humidity fields
    OPEN (nuprhumi,FILE=yuprhumi,FORM='FORMATTED',STATUS='REPLACE',IOSTAT=nstat)
    IF (nstat /= 0) THEN
        yerrmsg = 'opening of file YUPRHUMI failed'
        CALL model_abort (my_cart_id, 7003, yerrmsg, yroutine)
    ENDIF
    REWIND nuprhumi
  ENDIF

  !----------------------------------------------------------------------------
  !  Section 2.2: Determining the initial values
  !----------------------------------------------------------------------------

  ! get actual date and time
  !-------------------------
  CALL DATE_AND_TIME (ydatearg, ytimearg)
  ydate = ydatearg(7:8)//'.'//ydatearg(5:6)//'.'//ydatearg(1:4)
  ytime = ytimearg(1:2)//'.'//ytimearg(3:4)//'.'//ytimearg(5:6)

  ! initialize line counter and total real time
  !----------------------------------------------
  nlines  = 0

  IF (nstart == 0) THEN

    ! Initializations of some variables
    !----------------------------------
    zflar  = 1.0 / (   ((ie_tot-nboundlines) - (1+nboundlines) + 1)         &
                     * ((je_tot-nboundlines) - (1+nboundlines) + 1) )
    zrdzm  = 0.0
    zpsm   = 0.0
    zdsem  = 0.0
    zmsem  = 0.0
    zkem   = 0.0
    zqcm   = 0.0
    zqim   = 0.0
    zqrm   = 0.0
    zqsm   = 0.0
    zqgm   = 0.0

    ! Summations for initial data
    !----------------------------
    DO k = 1,ke
      DO j = jstart,jend
        DO i = istart,iend
          zrdz(i,j,k) = dp0(i,j,k)/g
          zke   = 0.5 * (u(i,j,k,nnew)**2 + v(i,j,k,nnew)**2)
          zdse  = cp_d * t(i,j,k,nnew) + 0.5*g*( hhl(i,j,k) + hhl(i,j,k+1) )
          zmse  = zdse + lh_v * qv(i,j,k,nnew)
          zdsem = zdsem + zdse * zrdz(i,j,k)
          zmsem = zmsem + zmse * zrdz(i,j,k)
          zkem  = zkem  + zke  * zrdz(i,j,k)
          zqcm  = zqcm  + qc(i,j,k,nnow) * zrdz(i,j,k)
          zrdzm = zrdzm + zrdz(i,j,k)
        ENDDO
      ENDDO
    ENDDO

    IF (ALLOCATED(qi)) THEN
      DO k = 1,ke
        DO j = jstart,jend
          DO i = istart,iend
            zqim  = zqim  + qi(i,j,k,nnow) * zrdz(i,j,k)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF (ALLOCATED(qr)) THEN
      DO k = 1,ke
        DO j = jstart,jend
          DO i = istart,iend
            zqrm  = zqrm  + qr(i,j,k,nnow) * zrdz(i,j,k)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF (ALLOCATED(qs)) THEN
      DO k = 1,ke
        DO j = jstart,jend
          DO i = istart,iend
            zqsm  = zqsm  + qs(i,j,k,nnow) * zrdz(i,j,k)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF (ALLOCATED(qg)) THEN
      DO k = 1,ke
        DO j = jstart,jend
          DO i = istart,iend
            zqgm  = zqgm  + qg(i,j,k,nnow) * zrdz(i,j,k)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    DO j = jstart,jend
      DO i = istart,iend
        zpsm  = zpsm  + ps(i,j,nnew)
      ENDDO
    ENDDO

    ! Collecting the values from all nodes in the parallel mode
    IF (num_compute > 1) THEN
      realbuf(1) = zdsem
      realbuf(2) = zmsem
      realbuf(3) = zkem
      realbuf(4) = zqcm
      realbuf(5) = zpsm
      realbuf(6) = zrdzm
      realbuf(7) = zqim
      realbuf(8) = zqrm
      realbuf(9) = zqsm
      realbuf(10) = zqgm

      CALL global_values (realbuf, 10, 'SUM', imp_reals, icomm_cart, 0,   &
                           yerrmsg, nzerror)

      zdsem = realbuf(1)
      zmsem = realbuf(2)
      zkem  = realbuf(3)
      zqcm  = realbuf(4)
      zpsm  = realbuf(5)
      zrdzm = realbuf(6)
      zqim  = realbuf(7)
      zqrm  = realbuf(8)
      zqsm  = realbuf(9)
      zqgm  = realbuf(10)
    ENDIF

    IF (my_cart_id == 0) THEN
      ! Computing the mean values
      dsem0 = zdsem / zrdzm
      msem0 = zmsem / zrdzm
      kem0  = zkem  / zrdzm
      qcm0  = zqcm  / zrdzm
      qim0  = zqim  / zrdzm
      qrm0  = zqrm  / zrdzm
      qsm0  = zqsm  / zrdzm
      qgm0  = zqgm  / zrdzm
      psm0  = zpsm  * zflar
    ENDIF
  ENDIF

  ! get elapsed time since last call
  !---------------------------------
! IF (ltime .EQV. .TRUE.) THEN
    CALL elapsed_time (zrealdiff)
    realtot_new = 0.0
! ELSE
!   realtot_new = 0.0
! ENDIF

  !----------------------------------------------------------------------------
  !  Section 2.3: Output of the initial values
  !----------------------------------------------------------------------------

  IF (my_cart_id == 0) THEN
    ! get actual UTC-time
    !--------------------
    CALL get_utc_date (0_iintegers, ydate_ini, dt, itype_calendar,          &
                       yandat1, yinidate, nzanjata, zanhh)
    WRITE (ychh,'(I10)') NINT (ntstep * dt / 3600.0)

    ! output of initial values for mass and energy fields
    !----------------------------------------------------
    WRITE (nuprmass,'(A         ,I6,A25,A3,A10,A5,A22,A1,A15,A11)')          &
        '      Experiment:  COSMO-Model        Number: ',nvers, yinidate,   &
        ' + ',ychh,' H  (',yakdat2,')',ydate,ytime
    WRITE (nuprmass,'(A61,A16,F8.2)')                                       &
        '      Elapsed time for providing initial and boundary values:',    &
        '     REAL (S):  ',zrealdiff
    WRITE (nuprmass,'(A3)') '   '
    WRITE (nuprmass,'(A)')                                                  &
        '      Initial mean values for nstart = 0 for several variables:'

    WRITE (nuprmass,'(A54,A54)')                                            &
        '      surface pressure (hPa)  dry static energy (J/kg)',           &
        '  moist static energy (J/kg)   kinetic energy (J/kg)  '
    WRITE (nuprmass,'(A,F20.3,3F25.2)')                                     &
        '# ', psm0*0.01, dsem0, msem0, kem0
    WRITE (nuprmass,'(A3)') '   '

    WRITE (nuprmass,'(A         ,I6,A25,A3,A10,A5,A22,A1,A15,A11)')         &
        '      Experiment:  COSMO-Model        Number: ',nvers, yinidate,   &
        ' + ',ychh,' H  (',yakdat2,')',ydate,ytime

    WRITE (nuprmass,'(A6,11A10)')                                           &
        "ntstep",  "Real",    "dpsdt",   "ps",      "dse",     "mse",       &
        "ke",      "vhmax",   "wmax",    "wa850",   "wa500",   "wa300"
    WRITE (nuprmass,'(A6,11A10)')                                           &
        "  ",      "s",       "Pa/s",    "hPa",     "J/kg",    "J/kg",      &
        "J/kg",    "m/s",     "m/s",     "cm/s",    "cm/s",    "cm/s"
    WRITE (nuprmass,'(A6,11A10)')                                           &
        " ",       " ",       "E-2",     " ",       "E+3",     "E+3",       &
        " ",       " ",       " ",       " ",       " ",       " "
    WRITE (nuprmass,'(A3)') "   "


    ! output of initial values for humidity fields
    !---------------------------------------------
    WRITE (nuprhumi,'(A         ,I6,A25,A3,A10,A5,A22,A1,A15,A11)')          &
        '      Experiment:  COSMO-Model        Number: ',nvers, yinidate,   &
        ' + ',ychh,' H  (',yakdat2,')',ydate,ytime
    WRITE (nuprhumi,'(A61,A16,F8.2)')                                       &
        '      Elapsed time for providing initial and boundary values:',    &
        '     REAL (S):  ',zrealdiff
    WRITE (nuprhumi,'(A3)') '   '
    WRITE (nuprhumi,'(A39,I5,A25)')                                         &
        '      Initial mean values for nstart = 0 for several variables:'
    WRITE (nuprhumi,'(A,A)')                                                &
        '      cloud water  cloud ice     rain        snow       graupel',  &
        '      (all: g/kg)'

    WRITE (nuprhumi,'(A,F15.3,4F12.3)')                                     &
        '#  ', qcm0*1.0E3, qim0*1.0E3, qrm0*1.0E3, qsm0*1.0E3, qgm0*1.0E3

    WRITE (nuprhumi,'(A         ,I6,A25,A3,A10,A5,A22,A1,A15,A11)')         &
        '      Experiment:  COSMO-Model        Number: ',nvers, yinidate,   &
        ' + ',ychh,' H  (',yakdat2,')',ydate,ytime

    WRITE (nuprhumi,'(A6,11A10)')                                           &
        "ntstep",  "qc",      "qi",      "qr",      "qs",      "qg",        &
        "prrs",    "prss",    "prrk",    "prsk",    "rrn",     "rsn"
    WRITE (nuprhumi,'(A6,11A10)')                                           &
        " ",       "kg/kg",   "kg/kg",   "kg/kg",   "kg/kg",   "kg/kg",     &
        "mm/D",    "mm/D",    "mm/D",    "mm/D",    "mm",      "mm"
    WRITE (nuprhumi,'(A6,11A10)')                                           &
        " ",       "E-5",     "E-5",     "E-5",     "E-5",     "E-5",       &
        " ",       " ",       " ",       " ",       " ",       " "
    IF (itype_gscp == 4 ) THEN
      WRITE (nuprhumi,'(A)')                                                &
        '                                 (prss and rsn including graupel)'
    ENDIF
    WRITE (nuprhumi,'(A3)') '   '
  ENDIF

  IF ( (my_cart_id == 0) .AND. (ldump_ascii) ) THEN
    CLOSE (nuprmass, STATUS='KEEP')
    CLOSE (nuprhumi, STATUS='KEEP')
  ENDIF

!------------------------------------------------------------------------------
!  End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE init_meanvalues

!==============================================================================
!+ Module procedure computing mean values
!------------------------------------------------------------------------------

!option! -pvctl ifopt
!option! -pvctl _on_adb
SUBROUTINE meanvalues

!------------------------------------------------------------------------------
!
! Description:
!   This routine calculates mean values and control variables and prints them
!   on the unit 'yuprmass' (for mass and energy related variables) and 
!   'yuprhumi' (for humidity related variables), resp. 
!   The values are initialized in the routine init_diagnosis with the initial 
!   data. Then, meanvalues is called every 'nincmeanval' steps.
!   The variables n0meanval and nincmeanval can be set via the NAMELIST 
!   /RUNCTL/ and are 0 and 1 per default.
!
! Method:
!   Computing sums, mean values and maxima of certain fields
!
!------------------------------------------------------------------------------

! Local variables
! Mean values
REAL (KIND=ireals)         ::       &
  zpsm,    & ! mean surface pressure ps
  zrdz(ie,je,ke),   & ! rho0*dz
  zrdzm,   & ! 
  zdsem,   & ! mean dry static energy
  zmsem,   & ! mean moist static energy
  zkem,    & ! mean kinetic energy
  zdpsdtm, & ! mean absolute tendencies of the surface pressure
  zwa850,  & ! mean absolute vertical velocity in about 850 hpa
  zwa500,  & ! mean absolute vertical velocity in about 500 hpa
  zwa300,  & ! mean absolute vertical velocity in about 300 hpa
  zqcm,    & ! mean cloudwater content
  zqim,    & ! mean cloud ice content
  zqrm,    & ! mean rain content
  zqsm,    & ! mean snow content
  zqgm,    & ! mean graupel content
  zprrsm,  & ! mean precipitation rate for rain (scale)
  zprssm,  & ! mean precipitation rate for snow (scale)
  zprrcm,  & ! mean precipitation rate for rain (convective)
  zprscm,  & ! mean precipitation rate for snow (convective)
  zrrnm,   & ! mean sum of precipitaion for rain (scale+convective)
  zrsnm      ! mean sum of precipitaion for snow (scale+convective)

! Maximum values of the velocities
REAL (KIND=ireals)         ::       &
  zvamax,  & ! maximal horizontal velocity
  zwamax     ! maximal absolute vertical velocity

! Memory for sending variables
REAL (KIND=ireals)         ::  realbuf(20)

! support variables for computation
! (same meaning as above)
REAL (KIND=ireals), SAVE   ::       &
  zdse, zmse, zke, zvb, zwb, zflar, zrealdiff, ztot

INTEGER (KIND=iintegers)   ::       &
  i,j,k,nhg, izerror      ! support variables

CHARACTER (LEN= 8)         :: ydatearg, ytime
CHARACTER (LEN=10)         :: ydate, ytimearg
CHARACTER (LEN=10)         :: ychh
CHARACTER (LEN=80)         :: yzerrmsg
CHARACTER (LEN=20)         :: yzroutine

LOGICAL                    :: lqi, lqg, lqr, lqs

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE meanvalues
!------------------------------------------------------------------------------

  yzroutine = 'meanvalues'
  IF ( (my_cart_id == 0) .AND. (ldump_ascii) ) THEN
    ! file for mass and energy fields
    OPEN (nuprmass,FILE=yuprmass,FORM='FORMATTED',STATUS='OLD',   &
                   POSITION='APPEND', IOSTAT=izerror)
    IF (izerror /= 0) THEN
        yzerrmsg = 'opening of file YUPRMASS failed'
        CALL model_abort (my_cart_id, 7003, yzerrmsg, yzroutine)
    ENDIF

    ! file for humidity fields
    OPEN (nuprhumi,FILE=yuprhumi,FORM='FORMATTED',STATUS='OLD',   &
                   POSITION='APPEND', IOSTAT=izerror)
    IF (izerror /= 0) THEN
        yzerrmsg = 'opening of file YUPRHUMI failed'
        CALL model_abort (my_cart_id, 7003, yzerrmsg, yzroutine)
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
!  Section 1: date and time measurement
!------------------------------------------------------------------------------

! get actual date and time
  CALL DATE_AND_TIME (ydatearg, ytimearg)
  ydate = ydatearg(7:8)//'.'//ydatearg(5:6)//'.'//ydatearg(1:4)
  ytime = ytimearg(1:2)//'.'//ytimearg(3:4)//'.'//ytimearg(5:6)


!------------------------------------------------------------------------------
!  Section 2: Summations for 3-D arrays (energies and water contents)
!------------------------------------------------------------------------------

! Initializations of some variables
  zflar  = 1.0 / REAL (   ((ie_tot-nboundlines) - (1+nboundlines) + 1)      &
                        * ((je_tot-nboundlines) - (1+nboundlines) + 1), ireals)
  zrdzm = 0.0
  zdsem = 0.0
  zmsem = 0.0
  zkem  = 0.0
  zqcm  = 0.0
  zqim  = 0.0
  zqrm  = 0.0
  zqsm  = 0.0
  zqgm  = 0.0

! Summations
  lqg=ALLOCATED(qg)
  lqi=ALLOCATED(qi)
  lqr=ALLOCATED(qr)
  lqs=ALLOCATED(qs)

  DO k = 1,ke
!CDIR OUTERUNROLL=4
    DO j = jstart,jend
      DO i = istart,iend
        zrdz(i,j,k) = dp0(i,j,k)/g
        zrdzm = zrdzm + zrdz(i,j,k) 
        zke   = 0.5 * (u(i,j,k,nnew)**2 + v(i,j,k,nnew)**2)
        zdse  = cp_d * t(i,j,k,nnew) + g*0.5*(hhl(i,j,k)+hhl(i,j,k+1))
        zmse  = zdse + lh_v * qv(i,j,k,nnew)
        zdsem = zdsem + zdse * zrdz(i,j,k)
        zmsem = zmsem + zmse * zrdz(i,j,k)
        zkem  = zkem  + zke  * zrdz(i,j,k)
        zqcm  = zqcm  + qc(i,j,k,nnow) * zrdz(i,j,k)

        IF (lqi) zqim  = zqim  + qi(i,j,k,nnow) * zrdz(i,j,k)
        IF (lqr) zqrm  = zqrm  + qr(i,j,k,nnow) * zrdz(i,j,k)
        IF (lqs) zqsm  = zqsm  + qs(i,j,k,nnow) * zrdz(i,j,k)
        IF (lqg) zqgm  = zqgm  + qg(i,j,k,nnow) * zrdz(i,j,k)

      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
!  Section 3: Summations for 2-D arrays ( absolute tendency of the surface 
!             pressure, mean vertical velocities, precipitation rates) and
!             calculation the maximal velocities
!------------------------------------------------------------------------------

  ! initialize mean and maximal values
  zpsm    = 0.0
  zdpsdtm = 0.0
  zvamax  = 0.0
  zwamax  = 0.0
  zwa850  = 0.0
  zwa500  = 0.0
  zwa300  = 0.0
  zprrsm  = 0.0
  zprssm  = 0.0
  zprrcm  = 0.0
  zprscm  = 0.0
  zrrnm   = 0.0
  zrsnm   = 0.0

  ! Summations
  DO j = jstart,jend
    DO i = istart,iend
      zpsm    = zpsm    + ps(i,j,nnew)
      zdpsdtm = zdpsdtm + ABS ( dpsdt(i,j) )
      zwa850  = zwa850  + ABS ( w(i,j,klv850,nnow) )
      zwa500  = zwa500  + ABS ( w(i,j,klv500,nnow) )
      zwa300  = zwa300  + ABS ( w(i,j,klv300,nnow) )
      zprrsm  = zprrsm  + prr_gsp(i,j)
      IF (itype_gscp == 4 ) THEN
        zprssm  = zprssm  + prs_gsp(i,j) + prg_gsp(i,j)
      ELSE
        zprssm  = zprssm  + prs_gsp(i,j)
      ENDIF
      zprrcm  = zprrcm  + prr_con(i,j)
      zprscm  = zprscm  + prs_con(i,j)
      zrrnm   = zrrnm   + rain_gsp(i,j) + rain_con(i,j)
      IF (itype_gscp == 4 ) THEN
        zrsnm   = zrsnm   + snow_gsp(i,j) + snow_con(i,j) + grau_gsp(i,j)
      ELSE
        zrsnm   = zrsnm   + snow_gsp(i,j) + snow_con(i,j)
      ENDIF
    ENDDO
  ENDDO

  IF (llm) THEN
    ! zwa500 and zwa300 remain 0.0
    zwa500 = 0.0_ireals
    zwa300 = 0.0_ireals
  ENDIF

  ! maximal values of the velocities
  DO k = 1,ke
!CDIR OUTERUNROLL=4
    DO j = jstartv,jendv
      DO i = istartu,iendu
        zwb    = ABS ( w(i,j,k,nnew) )
!       zvb    = SQRT( u(i,j,k,nnew)**2 + v(i,j,k,nnew)**2 )
        zvb    =       u(i,j,k,nnew)**2 + v(i,j,k,nnew)**2
        zvamax = MAX ( zvamax, zvb)
        zwamax = MAX ( zwamax, zwb)
      ENDDO
    ENDDO
  ENDDO
  zvamax = SQRT(zvamax)

!------------------------------------------------------------------------------
!  Section 4: Collecting the data from all nodes in parallel mode and
!             computing the mean values and maximas
!------------------------------------------------------------------------------

  ! Collecting the values from all nodes in the parallel mode
  IF (num_compute > 1) THEN
    ! Collecting the values for summation
    realbuf( 1) = zdsem
    realbuf( 2) = zmsem
    realbuf( 3) = zkem
    realbuf( 4) = zqcm
    realbuf( 5) = zpsm
    realbuf( 6) = zrdzm
    realbuf( 7) = zdpsdtm
    realbuf( 8) = zprrsm
    realbuf( 9) = zprssm
    realbuf(10) = zprrcm
    realbuf(11) = zprscm
    realbuf(12) = zrrnm
    realbuf(13) = zrsnm
    realbuf(14) = zwa300
    realbuf(15) = zwa500
    realbuf(16) = zwa850
    realbuf(17) = zqim
    realbuf(18) = zqrm
    realbuf(19) = zqsm
    realbuf(20) = zqgm

    CALL global_values (realbuf, 20, 'SUM', imp_reals, icomm_cart, 0,  &
                         yzerrmsg, izerror)

    zdsem   = realbuf( 1)
    zmsem   = realbuf( 2)
    zkem    = realbuf( 3)
    zqcm    = realbuf( 4)
    zpsm    = realbuf( 5)
    zrdzm   = realbuf( 6)
    zdpsdtm = realbuf( 7)
    zprrsm  = realbuf( 8)
    zprssm  = realbuf( 9)
    zprrcm  = realbuf(10)
    zprscm  = realbuf(11)
    zrrnm   = realbuf(12)
    zrsnm   = realbuf(13)
    zwa300  = realbuf(14)
    zwa500  = realbuf(15)
    zwa850  = realbuf(16)
    zqim    = realbuf(17)
    zqrm    = realbuf(18)
    zqsm    = realbuf(19)
    zqgm    = realbuf(20)

    ! Collecting the values for maximation
    realbuf(1) = zvamax
    realbuf(2) = zwamax

    CALL global_values (realbuf, 2, 'MAX', imp_reals, icomm_cart, 0,  &
                         yzerrmsg, izerror)

    zvamax  = realbuf(1)
    zwamax  = realbuf(2)

  ENDIF

  IF (my_cart_id == 0) THEN
    ! mean values
    zdsem   = zdsem   / zrdzm
    zmsem   = zmsem   / zrdzm
    zkem    = zkem    / zrdzm
    zqcm    = zqcm    / zrdzm  
    zpsm    = zpsm    * zflar
    zqim    = zqim    / zrdzm
    zqrm    = zqrm    / zrdzm
    zqsm    = zqsm    / zrdzm
    zqgm    = zqgm    / zrdzm

    zdpsdtm = zdpsdtm * zflar
    zprrsm  = zprrsm  * zflar * day_len
    zprssm  = zprssm  * zflar * day_len
    zprrcm  = zprrcm  * zflar * day_len
    zprscm  = zprscm  * zflar * day_len
    zrrnm   = zrrnm   * zflar
    zrsnm   = zrsnm   * zflar
    zwa300  = zwa300  * zflar
    zwa500  = zwa500  * zflar
    zwa850  = zwa850  * zflar
  ENDIF

! get elapsed time for this time step
! IF (ltime .EQV. .TRUE.) THEN
    CALL elapsed_time (zrealdiff)
    realtot_new = realtot_new + zrealdiff
! ELSE
!   realtot_new = 0.0
! ENDIF

  ! increase line counter
  nlines   = nlines  + 1

!------------------------------------------------------------------------------
!  Section 5: Output of the results
!------------------------------------------------------------------------------

  IF (my_cart_id == 0) THEN
    WRITE (nuprmass,'(I8,11F10.3)')                                         &
        ntstep , zrealdiff, zdpsdtm*100.0, zpsm*0.01, (zdsem- dsem0)*0.001, &
        (zmsem- msem0)*0.001, (zkem- kem0), zvamax, zwamax,                 &
        zwa850*100.0, zwa500*100.0, zwa300*100.0

    WRITE (nuprhumi,'(I8,11F10.3)')                                         &
        ntstep , zqcm*1.0E5, zqim*1.0E5, zqrm*1.0E5, zqsm*1.0E5, zqgm*1.0E5,&
        zprrsm, zprssm, zprrcm, zprscm, zrrnm, zrsnm

    ! if necessary, start a new page
    IF ( (nlines >= 48) .OR. (ntstep == nstop) ) THEN
      IF ( (MOD(NINT(ntstep*dt),3600) == 0) .OR. (ntstep == nstop) ) THEN
        nhg = NINT ( nlines*nincmeanval * dt / 3600.0 )
        WRITE (nuprmass,'(A3)') '  '
        WRITE (nuprmass,'(A3)') '  '
        WRITE (nuprmass,'(A29,I2,A19,F10.2)')                               &
            '      Elapsed real time for ',nhg,' hour(s) forecast: ',       &
              realtot_new
        WRITE (nuprhumi,'(A3)') '  '
        WRITE (nuprhumi,'(A3)') '  '
        WRITE (nuprhumi,'(A29,I2,A19,F10.2)')                               &
            '      Elapsed real time for ',nhg,' hour(s) forecast: ',       &
              realtot_new
        nlines  = 0
        realtot_new = 0.0
        WRITE (ychh,'(I10)') NINT ( ntstep * dt / 3600.0 )

        IF (ntstep /= nstop) THEN
          WRITE (nuprmass,'(A         ,I6,A25,A3,A10,A5,A22,A1,A15,A11)')   &
              '#1    Experiment:  COSMO-Model        Number: ',nvers,       &
              yinidate,' + ',ychh,' H  (',yakdat2,')',ydate,ytime

          WRITE (nuprmass,'(A6,11A10)')                                       &
              "ntstep",  "Real",    "dpsdt",   "ps",      "dse",     "mse",   &
              "ke",      "vhmax",   "wmax",    "wa850",   "wa500",   "wa300"
          WRITE (nuprmass,'(A6,11A10)')                                       &
              "  ",      "s",       "Pa/s",    "hPa",     "J/kg",    "J/kg",  &
              "J/kg",    "m/s",     "m/s",     "cm/s",    "cm/s",    "cm/s"
          WRITE (nuprmass,'(A6,11A10)')                                       &
              " ",       " ",       "E-2",     " ",       "E+3",     "E+3",   &
              " ",       " ",       " ",       " ",       " ",       " "
          WRITE (nuprmass,'(A3)') '   '

          WRITE (nuprhumi,'(A         ,I6,A25,A3,A10,A5,A22,A1,A15,A11)')   &
              '      Experiment:  COSMO-Model        Number: ',nvers,       &
              yinidate, ' + ',ychh,' H  (',yakdat2,')',ydate,ytime

          WRITE (nuprhumi,'(A6,11A10)')                                       &
              "ntstep",  "qc",      "qi",      "qr",      "qs",      "qg",    &
              "prrs",    "prss",    "prrk",    "prsk",    "rrn",     "rsn"
          WRITE (nuprhumi,'(A6,11A10)')                                       &
              " ",       "kg/kg",   "kg/kg",   "kg/kg",   "kg/kg",   "kg/kg", &
              "mm/D",    "mm/D",    "mm/D",    "mm/D",    "mm",      "mm"
          WRITE (nuprhumi,'(A6,11A10)')                                       &
              " ",       "E-5",     "E-5",     "E-5",     "E-5",     "E-5",   &
              " ",       " ",       " ",       " ",       " ",       " "

          IF (itype_gscp == 4 ) THEN
            WRITE (nuprhumi,'(A)')                                          &
          '                                 (prss and rsn including graupel)'
          ENDIF
          WRITE (nuprhumi,'(A3)') '   '
        ENDIF
      ENDIF
    ENDIF
  ENDIF

IF (  (my_cart_id == 0) .AND. ((ldump_ascii) .OR. (ntstep == nstop)) ) THEN
  CLOSE (nuprmass, STATUS='KEEP')
  CLOSE (nuprhumi, STATUS='KEEP')
ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE meanvalues

!==============================================================================

END MODULE src_meanvalues 
