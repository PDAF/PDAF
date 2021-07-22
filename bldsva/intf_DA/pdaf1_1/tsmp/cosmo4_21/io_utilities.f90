!+ Source Module providing utility routines for grib I/O
!==============================================================================

MODULE io_utilities

!==============================================================================
!
! Description:
!  This module provides all routines that have to deal with input and output
!  of grib files. Some of these routines use the message-passing library MPI. 
!  If possible the routines are written plug-compatible (with the exception
!  of setup_io).
!
!  Routines (module procedures) currently contained:
!
!    - make_fn
!       builds the path and the name of I/O files.
!
!    - create_db_order
!       builds the database input order
!
!       
!    The following routines contain MPI calls:
!    
!    - open_file
!       opens a file with routines from special libraries (Grib, NetCDF, ...)
!
!    - close_file
!       closes a file with routines from special libraries (Grib, NetCDF, ...)
!
!    - read_grib
!       reads records from a grib file and distributes them to the PEs
!
!    - read_netcdf
!       reads records from a NetCDF file and distributes them to the PEs
!
!    - read_restart
!       reads records from a restart file and distributes them to the PEs
!
!    - write_grib
!       collects gribed records from the processors and writes them to disk.
!
!    - write_netcdf
!       collects records from the processors and writes them to netcdf files.
!
!    - write_restart
!       collects records from the processors and writes them to a binary
!       restart file
!
!    - check_ec_grid
!       checks the NAMELIST parameters against the EC values from the input.
!
!    - check_lm_grid 
!       checks the NAMELIST parameters against the LM values from the input.
!
!    - check_gme_grid 
!       checks the NAMELIST parameters against the GME values from the input.
!
!    - check_record
!       prints MIN, MAX and meanvalues of all records to the file YUCHKDAT
!
!    - difmin_360
!       computes the difference in minutes between days for a climatological
!       year with 360 days
!
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
! 1.4        1998/05/22 Guenther Doms
!  Adaptions for the two time-level integration scheme
! 1.7        1998/07/16 Guenther Doms
!  Removal of global array 'rrssk'.
! 1.8        1998/08/03 Ulrich Schaettler
!  Use grib parameters from module data_io.f90
! 1.10       1998/09/29 Ulrich Schaettler
!  Updated Use lists from the modules
! 1.17       1998/11/17 Ulrich Schaettler
!  Correct some ANSI violations.
! 1.20       1999/01/07 Guenhter Doms
!  Renaming of some global variables.
! 1.29       1999/05/11 Ulrich Schaettler
!  Adaptations to use this module also in GME2LM
! 1.32       1999/08/24 Guenther Doms
!  Correct a REAL-declaration
! 1.34       1999/12/10 Ulrich Schaettler
!  Changed routine check_record to use it also in GME2LM properly.
! 1.39       2000/05/03 Ulrich Schaettler
!  Adaptations to include asynchronous IO.
! 2.8        2001/07/06 Ulrich Schaettler
!  Checked and corrected the interfaces to mpe_io;
!  Introduced ymode_read and ymode_write;
!  In check_lm_grid and check_gme_grid now the actual date of the forecast is
!  checked instead of the initial date (or reference date).
! 2.9        2001/07/16 Ulrich Schaettler
!  Correction for reading Grib data.
! 2.12       2001/11/07 Ulrich Schaettler
!  Corrected the check for the pole of the rotated grid; adaptation to GME2LM
! 2.14       2002/02/15 Ulrich Schaettler
!  Modifications to SR check_lm_grid to allow the use of the SLEVE coordinate
! 2.15       2002/03/20 Ulrich Schaettler
!  Bug correction in check_lm_grid
! 2.17       2002/05/08 Ulrich Schaettler
!  Adaptations to perform communications for I/O in irealgrib-format;
!  in SR check_lm_grid, only bit 7 of igds(9) is tested for 0
! 3.6        2003/12/11 Ulrich Schaettler
!  Additional checks for "time range indicator" and "unit of time range"
!  in the Subroutines check_lm_grid and check_gme_grid
! 3.7        2004/02/18 Ulrich Schaettler
!  Adaptations in routine make_fn for checking the extension of file names
! 3.13       2004/12/03 Ulrich Schaettler
!  Adaptations to new version of INT2LM
! 3.15       2005/03/03 Ulrich Schaettler
!  Replaced FLOAT by REAL; Implemented unit_of_time in check-routines
! 3.16       2005/07/22 Ulrich Schaettler
!  Eliminated time check for T_CL and W_CL in case of LM2LM
! 3.18       2006/03/03 Ulrich Schaettler
!  Modifications to implement NetCDF and Restart I/O
!    (renamed and changed: open_grib, close_grib -> open_file, close_file
!     new: read_netcdf, read_restart, write_netcdf, write_restart)
!  New ytunit 'd' in SR make_fn (for CLM: file name with date as in laf-files)
!  New SR difmin_360: computes difference for 2 dates using a year with 360 days
!  Added parameter lyear_360 in SRs check_*_grid to choose proper SR difmin
!  Changes in SR check_lm_grid for new type of coding the vertical coordinates
!    (introduced by MeteoSwiss for SLEVE vertical coordinate)
! 3.19       2006/04/25 Ulrich Schaettler
!  In read_netcdf: skip variables in the initial list, that are not present
! 3.21       2006/12/04 Ulrich Schaettler
!  In case of ymode_write=append, restart files are opened as NEW
!  Adapt check_record to use of undefined values in netcdf
!  Bug corrections in read_restart, write_restart
! V3_23        2007/03/30 Ulrich Schaettler
!  Introduced idbg_level for verbosity of output
!  Changed SR make_fn to adapt the file name to quarterly hours, if necessary
! V3_26        2007/05/29 Ulrich Schaettler
!  Adaptations in SR make_fn for writing analysis files with flexible dt
! V4_1         2007/12/04 Ulrich Schaettler
!  Some editorial changes (print outs for debugging)
! V4_4         2008/07/16 Ulrich Schaettler, Burkhardt Rockel
!  Changed NL parameter lyear_360 to itype_calendar, to have several options
!  Allow 2D fields in netcdf input
! V4_5         2008/09/10 Guenther Zaengl
!  Adaptations for new reference atmosphere
! V4_8         2009/02/16 Guenther Zaengl
!  Change grib encoding for new reference atmosphere to facilitate future
!  extensions
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  In subroutine make_fn, introduced parameter ztgranul to make the
!  allowed time distances of output files more flexible in case of lhour=.true.
!  Before, the allowed time distance was set fixed to a quarter hour.
! V4_18        2011/05/26 Ulrich Schaettler
!  Adapted SR read_netcdf to identify the 3D external parameter field for
!  topographical corrections (Anne Roches)
! V4_19        2011/08/01 Ulrich Schaettler
!  Introduced conditional compilation for NetCDF and GRIBDWD
! V4_20        2011/08/31 Ulrich Schaettler
!  Safer communication between IO and compute PEs for little endian problem with MPI_BYTE
!  (send 4 additional bytes)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:
USE data_parameters , ONLY :   &
  ireals,    & ! KIND-type parameters for real variables
  iintegers, & ! KIND-type parameter for standard integer variables
  irealgrib, & ! KIND-type parameter for the real variables in the grib library
  intgribf,  & ! KIND-type parameter for the fortran files of the grib library
  intgribc,  & ! KIND-type parameter for the c files of the grib library
  iwlength     ! length of a integer word of the grib library in byte

!==============================================================================

! declare routines used from the parallel IO interface MPE_IO
USE mpe_io,           ONLY :   &
    mpe_io_open, mpe_io_read, mpe_io_write, mpe_io_close, mpe_io_complete,   &
    mpe_get_db_tab, mpe_db_tab

!==============================================================================

#ifdef NETCDF
! declare netcdf variables and functions
USE netcdf,           ONLY :   &
  NF90_clobber,            &
  NF90_close ,             &
  NF90_create,             &
  NF90_noerr,              &
  NF90_nofill,             &
  NF90_nowrite,            &
  NF90_open,               &
  NF90_inquire_variable,   &
  NF90_inquire_dimension,  &
  NF90_put_var,            &
  NF90_get_var,            &
  NF90_set_fill,           &
  NF90_strerror,           &
  NF90_write
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Include statements
  INCLUDE "mpif.h"

!==============================================================================

CONTAINS

!==============================================================================
!+ Creates the file name for grib-files or for ready-files
!------------------------------------------------------------------------------

SUBROUTINE make_fn (yhead, ydate, ytunit, yexten, nstep, dt, lhour,     &
                    ydir, ypathname, idbglv, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   *make_fn* creates the file names for the grib files or for ready-files
!   depending on yhead and ytunit.
!     Analysis:                 yhead//ydate_fn//yexten
!     Forecast and Boundary:    yhead//ytunit//yforrange//yexten
!
!   File name for forecasts ('f') or restart files ('r')
!   The form of the forecast range depends on "ytunit":
!     For ytunit = 't':  yforrange is the number of timesteps
!     For ytunit = 'f':  yforrange is given in the form "ddhhmmss" where
!                        dd: days, hh: hours, mm: minutes, ss: seconds
!     For ytunit = 'c':  yforrange is given in the form "yyydddhh" where
!                        yyy: years, ddd: days, hh: hours
!     For ytunit = 'd':  yforrange is given as ydate (for Climate LM Version)
!
! Method:
!
!------------------------------------------------------------------------------

!     Input
CHARACTER (LEN= 3),       INTENT(IN)  ::    &
  yhead        ! 3-Character header of the file, e.g.
               !   'gaf': analysis of GME for the full domain
               !   'laf': analysis of LM for the full domain
               !   'lbf': boundary for LM on the full domain
               !   'GME' or 'LMA': for ready files

CHARACTER (LEN=10),       INTENT(IN)  ::    &
  ydate        ! Initial date of the forecast in the form
               !   yyyymmddhh (yyyy: year, mm: month, dd: day, hh: time)

CHARACTER (LEN= 1),       INTENT(IN)  ::    &
  ytunit       ! 1-Character time unit for forecast range, e.g.
               !   't': timesteps, 'f': forecast mode, 
               !   'c': climate mode

CHARACTER (LEN= 1),       INTENT(IN)  ::    &
  yexten       ! 1-Character extension of file name, e.g.
               !   'c': constant files
               !   'z': data on z-levels
               !   'p': data on pressure levels
               !   's': data on "satellite" levels

CHARACTER (LEN= *),       INTENT(IN)  ::    &
  ydir         ! directory of the file

INTEGER (KIND=iintegers), INTENT(IN)  ::    &
  nstep,     & ! Current timestep number
  idbglv       ! debug level for verbosity of output

REAL (KIND=ireals),       INTENT(IN)  ::    &
  dt           ! Timestep of the model (in seconds)

LOGICAL,                  INTENT(IN)  ::    &
  lhour        ! if .TRUE., the filename has to be constructed as multiples
               ! of quarterly hours (only for ytunit = 'f' or 'c'

!------------------------------------------------------------------------------

! Output
CHARACTER (LEN= *),       INTENT(OUT) ::    &
  ypathname    ! Name of the file with full path name for access

INTEGER (KIND=iintegers), INTENT(OUT) ::    &
  ierror       ! Error flag; set to 0 if no error occured

!------------------------------------------------------------------------------

! Local variables
REAL (KIND=ireals)                    ::    &
  zforrange, & ! Forecast range in seconds
  z2, zhour, zdiff

INTEGER (KIND=iintegers)              ::    &
  mfor_s,    & ! Forecast range (seconds)
  mfor_m,    & ! Forecast range (minutes)
  mfor_h,    & ! Forecast range (hours)
  mfor_d,    & ! Forecast range (days)
  mfor_y       ! Forecast range (years)

INTEGER (KIND=iintegers)              ::    &
  izh

CHARACTER (LEN=10)                    ::    &
! yforrange, & ! Forecast range as character
  yzlocdate    ! if ydate has to be modified

CHARACTER (LEN=14)                    ::    &
  yforrange ! Forecast range as character

CHARACTER (LEN=4)                     ::    &
  yforranges

CHARACTER (LEN=20)                    ::    &
  yfname       ! name of the file

! ub>>
! for lhour=.true. : granularity of the output time (round to nearest multiple of ztgranul)   [s]
!                    NOTE: ztgranul has to fit exactly into 3600 s, i.e., ztgranul might be, 
!                          for example, something like 60.0, 300.0, 600.0, 900.0 or 1800.0
REAL(KIND=ireals), PARAMETER :: ztgranul = 900.0    ! in s, 15 min for now.
! ub<<

! End of header
!------------------------------------------------------------------------------

  ierror = 0
  yforrange = '              '
  yzlocdate = ydate

  zforrange = REAL(nstep, ireals)*dt        ! forecast time in seconds
  z2        = zforrange / 3600.0_ireals     ! forecast time in hours

! ub>>
  IF (lhour) THEN
    zhour = REAL (INT (z2, iintegers), ireals) * 3600.0_ireals
                               ! this is the last full hour in seconds
    zdiff = zforrange - zhour  ! difference to the last full hour in seconds

    IF (zdiff >= 3600.0 - ztgranul*0.5) THEN

      ! clip zforrange to one hour:

      zforrange = zhour + 3600.0_ireals ! add another hour

      ! then also the date has to be adapted!!!
      READ (ydate(9:10), '(I2.2)') izh
      WRITE(yzlocdate(9:10), '(I2.2)') izh+1
      IF (idbglv > 50) THEN
        PRINT *, '     in make_fn:  add one hour to the date !!! ', izh
      ENDIF

    ELSE

      ! make zforrange the next full interval of ztgranul:

      zforrange = zhour + NINT(zdiff/ztgranul, iintegers) * ztgranul
      
    ENDIF

  ENDIF

!!$  IF (lhour) THEN
!!$    zhour = REAL (INT (z2, iintegers), ireals) * 3600.0_ireals
!!$                               ! this is the last full hour in seconds
!!$    zdiff = zforrange - zhour  ! difference to the last full hour in seconds
!!$    IF       (zdiff < 450.0_ireals) THEN
!!$      zforrange = zhour                 ! nothing to add
!!$    ELSEIF ( ( 450.0_ireals <= zdiff) .AND. (zdiff < 1350.0_ireals) ) THEN
!!$      zforrange = zhour + 900.0_ireals  ! add a quarter
!!$    ELSEIF ( (1350.0_ireals <= zdiff) .AND. (zdiff < 2250.0_ireals) ) THEN
!!$      zforrange = zhour + 1800.0_ireals ! add half an hour
!!$    ELSEIF ( (2250.0_ireals <= zdiff) .AND. (zdiff < 3150.0_ireals) ) THEN
!!$      zforrange = zhour + 2700.0_ireals ! add 3 quarters
!!$    ELSE
!!$      zforrange = zhour + 3600.0_ireals ! add another hour
!!$
!!$      ! then also the date has to be adapted!!!
!!$      READ (ydate(9:10), '(I2.2)') izh
!!$      WRITE(yzlocdate(9:10), '(I2.2)') izh+1
!!$      IF (idbglv > 50) THEN
!!$        PRINT *, '     in make_fn:  add one hour to the date !!! ', izh
!!$      ENDIF
!!$    ENDIF
!!$  ENDIF
! ub<<

  IF (idbglv > 50) THEN
    PRINT *, '     in make_fn:  ', REAL(nstep, ireals)*dt, z2, zhour, zdiff, &
                                    zforrange, '  ', TRIM(yzlocdate)
  ENDIF

  ! Compute forecast range as character string depending on ytunit
  SELECT CASE (ytunit)
  CASE ('t')       ! forecast range in timesteps
    WRITE (yforrange,'(I8.8)') nstep

  CASE ('f')       ! forecast range in forecast-mode

    mfor_d    = INT ( zforrange/86400.0_ireals, iintegers)
    mfor_h    = INT ((zforrange - REAL(mfor_d, ireals)*86400.0_ireals) /    &
                                                   3600.0_ireals, iintegers)
    mfor_m    = INT ((zforrange - REAL(mfor_d, ireals)*86400.0_ireals       &
                                - REAL(mfor_h, ireals)* 3600.0_ireals) /    &
                                                     60.0_ireals, iintegers)
    mfor_s    = NINT( zforrange - REAL(mfor_d, ireals)*86400.0_ireals       &
                                - REAL(mfor_h, ireals)* 3600.0_ireals       &
                                - REAL(mfor_m, ireals)*   60.0_ireals,      &
                                                                  iintegers)
    WRITE (yforrange,'(4(I2.2), I4.4)') mfor_d, mfor_h, mfor_m, mfor_s, 0

  CASE ('c')

    mfor_y    = INT ( zforrange / (365.0_ireals*86400.0_ireals), iintegers)
    mfor_d    = INT ((zforrange -                                           &
                        REAL (mfor_y,ireals)*365.0_ireals*86400.0_ireals)   &
                                                / 86400.0_ireals, iintegers)
    mfor_h    = NINT((zforrange -                                           &
                        REAL (mfor_y, ireals)*365.0_ireals*86400.0_ireals   &
                      - REAL (mfor_d, ireals)*86400.0_ireals) /             &
                                                   3600.0_ireals, iintegers)
    WRITE (yforrange,'(I3.3, I3.3, I2.2, I4.4)') mfor_y, mfor_d, mfor_h, 0

  CASE ('d')


    mfor_d    = INT ( zforrange/86400.0_ireals, iintegers)
    mfor_h    = INT ((zforrange - REAL(mfor_d, ireals)*86400.0_ireals) /    &    
                                                   3600.0_ireals, iintegers)
    mfor_m    = INT ((zforrange - REAL(mfor_d, ireals)*86400.0_ireals       &    
                                - REAL(mfor_h, ireals)* 3600.0_ireals) /    &    
                                                     60.0_ireals, iintegers)
    mfor_s    = NINT( zforrange - REAL(mfor_d, ireals)*86400.0_ireals       &    
                                - REAL(mfor_h, ireals)* 3600.0_ireals       &    
                                - REAL(mfor_m, ireals)*   60.0_ireals,      &    
                                                                  iintegers)

    WRITE (yforranges,'(2(I2.2))') mfor_m, mfor_s
    yforrange = yzlocdate(1:LEN_TRIM(yzlocdate))//yforranges

  CASE DEFAULT
    ! Unknown file type or time unit, error exit
    IF (idbglv > 0) THEN
      PRINT *,'  Error in *make_fn*, unknown file type,',' yhead:  ', yhead
      PRINT *,'  Error in *make_fn*, unknown time unit,',' ytunit: ', ytunit
    ENDIF
    ierror = 1
  END SELECT

  ! put file name together
  IF (yhead(2:2) .EQ. 'a') THEN
    ! File name for analysis ('a')
    WRITE (yforranges,'(2(I2.2))') 0,0
    IF (yexten == ' ') THEN
      yfname = yhead//yzlocdate//yforranges
    ELSE
      yfname = yhead//yzlocdate//yforranges//yexten
    ENDIF
  ELSEIF (yhead(3:3) .EQ. 'A') THEN
    ! file name for LMA ready-file
    yfname = yhead//'_'//yzlocdate
  ELSEIF (yhead(2:2) .EQ. 'M') THEN
    ! file name for ready-files
    yfname = yhead//'_'//TRIM(yforrange)
  ELSE
    ! all other file  names
    IF (yexten == ' ') THEN
      yfname = yhead//ytunit//TRIM(yforrange)
    ELSE
      yfname = yhead//ytunit//TRIM(yforrange)//yexten
    ENDIF
  ENDIF

  ! Add the directory
  ypathname = ydir(1:LEN_TRIM(ydir))//yfname

END SUBROUTINE make_fn

!==============================================================================
!+ Module procedure in "io_utilities" for creating database order
!------------------------------------------------------------------------------

SUBROUTINE create_db_order (ydir, ydate, ymode, ydomain, ytunit, yextension,  &
                            nstep, dt, ydborder)

!-------------------------------------------------------------------------------!
! Description:
!  Creates the right database input order, analogous to 'create_filename'
!
! Method:
!  The name is constructed according to the file name conventions from the
!  input parameters. The directory also is an input parameter (see below).
!
!  Leave the code close to "create filename", additionally set database order
!
! Pallas GmbH, January 1999
!
!===============================================================================!
! Subroutine / Function arguments
! Scalar arguments with intent(in):
CHARACTER (LEN= *)      , INTENT(IN)  ::  &
  ydir,           & ! directory of the file
  ydate             ! date of the initial data

CHARACTER (LEN= 1)      , INTENT(IN)  ::  &
  ! see also the file name convention about these variables and their meaning
  ymode,          & ! characterizes the special kind of the data
  ydomain,        & ! specifies the region covered by the data
  ytunit,         & ! time unit of forecast range
  yextension        ! extension of the file name

INTEGER (KIND=iintegers), INTENT(IN)  ::  &
  nstep             ! actual time step

REAL    (KIND=ireals)   , INTENT(IN)  ::  &
  dt                ! length of the time step

! Scalar arguments with intent(out):
CHARACTER (LEN=*)     , INTENT(OUT) ::  &
  ydborder          ! name of the db order

!------------------------------------------------------------------------------
!
! Local scalars:
CHARACTER (LEN=10)              ::  &
  ytime                                ! string with output time

INTEGER (KIND=iintegers)        ::  &
  izyear, izday, izhour, izmin, izsec  ! integer variables to determine the
                                       ! output time

INTEGER (KIND=iintegers)        :: iz_ytablen
!
!- End of header
!==============================================================================
  izsec   = nstep * dt

  IF ( (ymode /= 'a') .AND. (ymode /= 'i') ) THEN
    SELECT CASE(ytunit)
    CASE('t')  ! forecast range in timesteps
      WRITE(ytime,'(I8.8)') nstep
    CASE('f')  ! forecast range in ddhhmmss
      izday  =  izsec / 86400
      izhour = (izsec - izday*86400) / 3600
      izmin  = (izsec - izday*86400 - izhour*3600) / 60
      izsec  =  izsec - izday*86400 - izhour*3600 - izmin*60
      WRITE(ytime,'(4(I2.2))') izday, izhour, izmin, izsec
    CASE('c')  ! forecast range in yyydddhh
      izyear =  izsec /(86400*365)
      izday  = (izsec - izyear*86400*365) / 86400
      izhour = (izsec - izyear*86400*365 - izday*86400) /3600
      WRITE(ytime,'(I3.3, I3.3, I2.2)') izyear, izday, izhour
    END SELECT

!! create input order with tabtyp 'bd' for 'boundary'

    ydborder=' '
    IF( izhour<10 ) THEN
    write(ydborder,'(a,i1)') 'vv=',izhour
    ELSE IF( izhour<100 ) THEN
    write(ydborder,'(a,i2)') 'vv=',izhour
    ELSE IF( izhour<1000 ) THEN
    write(ydborder,'(a,i3)') 'vv=',izhour
    ELSE IF( izhour<10000 ) THEN
    write(ydborder,'(a,i4)') 'vv=',izhour
    ENDIF

    iz_ytablen = LEN_TRIM(ydborder)
    call mpe_db_tab('bd',ydate,ydborder,iz_ytablen)

  ELSE

!! create input order with tabtyp 'ana' for 'analysis'
    iz_ytablen = 0
    call mpe_db_tab('ana',ydate,ydborder,iz_ytablen)
  ENDIF

END SUBROUTINE create_db_order

!==============================================================================
!+ Module procedure in "io_utilities" for opening a file
!------------------------------------------------------------------------------

SUBROUTINE open_file (nudat, datname, ymode, yformat, icomm, my_id, npes,   &
                      lasync, idbglv, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Open a file for reading or writing. The calling routine has to specify
!   the mode of the file (read, write, append, additional (Grib) information),
!   the format of the file (Grib1, NetCDF, binary) and the filename.
!   The routine returns the file specifier (nudat). Only in case of binary
!   restart files, the file specifier is already determined before.
!   This is a collective routine, i.e. it must be called by all compute
!   processors.
!
! Method:
!   Calls to special libraries.
!
!------------------------------------------------------------------------------
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):

! Parameter list 
CHARACTER (LEN= *)       , INTENT(IN)    ::  &
  datname     ! filename to open

CHARACTER (LEN= 3)       , INTENT(IN)    ::  &
  ymode       ! access mode: read, write, append: 'r  ', 'w  ', ...

CHARACTER (LEN= 4)       , INTENT(IN)    ::  &
  yformat     ! format of the file (grb1, ncdf, bina)

INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
  icomm,    & ! MPI communicator
  my_id,    & ! id in communicator icomm
  npes,     & ! number of PEs
  idbglv      ! for verbosity of debug output

! Scalar arguments with intent(out):
INTEGER  (KIND=iintegers), INTENT(INOUT) ::  &
  nudat    ! internal file descriptor

LOGICAL                  , INTENT(IN)    ::  &
  lasync   ! indicates whether asynchronous IO is used or not

CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg     ! string for error messages

INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ierror      ! error status

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER  (KIND=iintegers)           :: &
  joldMode,                            & ! for netcdf routines
  implcode                               ! Error variable for MPI

INTEGER  (KIND=intgribc)            :: &
  nudatc,                              & ! file descriptor for C-routine
  ierrc                                  ! error code for C-routine

CHARACTER (LEN=300)                 :: &
  ydborder                               ! for specifying database job

LOGICAL                             :: &
  stop_dummy                             ! for calling mpe_io_open

!- End of header
!------------------------------------------------------------------------------

ierror   = 0
yerrmsg  = '         '
ierrc    = 0_intgribc
nudatc   = 0_intgribc

IF (my_id == 0 .AND. idbglv > 1)                                         &
   PRINT *,'OPEN: ',yformat,'-file: ', datname(1:LEN_TRIM(datname))

SELECT CASE (yformat)

! At the moment, asynchronous IO is possible only for grib1-files
! I/O for all other formats is handled by processor 0

#ifdef GRIBDWD
CASE ('grb1')

  IF (lasync .OR. (npes > 1) ) THEN
    ! Open file using the parallel IO interface
    IF (ymode(1:1).EQ.'r') THEN
      ! For reading files, only PE 0 has to call mpe_io_open
      IF (my_id == 0) THEN
        CALL mpe_get_db_tab ('in', ydborder)
        CALL mpe_io_open (nudat, datname(1:LEN_TRIM(datname)), ydborder,   &
                          ymode, stop_dummy, ierror)
        IF (ierror /= 0) THEN
          yerrmsg = 'Error opening '//datname
          ! this error message has to be broadcasted to all other PEs
          ierror  = 3
        ENDIF
      ENDIF
    ELSEIF ( (ymode(1:1).EQ.'w') .OR. (ymode(1:1).EQ.'a') ) THEN
      ! For writing files, all PEs have to call mpe_io_open
      CALL mpe_io_open (nudat, datname(1:LEN_TRIM(datname)), ' ',          &
                        ymode, stop_dummy, ierror)
    ENDIF
  ELSE
    ! Open file using griblib calls
    IF (my_id == 0) THEN
      CALL copen(nudatc,datname(1:LEN_TRIM(datname)),ymode,ierrc)
      nudat = nudatc
      IF (ierrc /= 0) THEN
        yerrmsg = 'Error opening '//datname
        ! this error message has to be broadcasted to all other PEs
        ierror  = 3
      ENDIF
    ENDIF
  ENDIF
#endif

#ifdef NETCDF
CASE ('ncdf')

  IF (my_id == 0) THEN
    IF (ymode(1:1) == 'r') THEN
  
      ierror = nf90_open(TRIM(datname), nf90_nowrite, nudat)
      IF (ierror /= nf90_noerr) THEN
        PRINT *, TRIM(NF90_strerror(ierror))
      ENDIF
  
    ELSEIF (ymode(1:1) == 'w') THEN
  
      ierror = nf90_create(TRIM(datname), nf90_clobber, nudat)
      ierror = nf90_set_fill(nudat, NF90_NOFILL, joldMode)
      IF (ierror /= nf90_noerr) THEN
        PRINT *, TRIM(NF90_strerror(ierror))
      ENDIF
  
    ELSEIF (ymode(1:1) == 'a') THEN
  
      ierror = nf90_open(TRIM(datname), nf90_write, nudat)
      ierror = nf90_set_fill(nudat, NF90_NOFILL, joldMode)
      IF (ierror /= nf90_noerr) THEN
        PRINT *, TRIM(NF90_strerror(ierror))
      ENDIF
  
    ELSE
      PRINT *, ' *** ERROR:  Wrong mode for netcdf-format: ', ymode(1:1), &
               ' *** '
      ierror = -1
    ENDIF
  ENDIF
#endif

CASE ('bina')
 
  IF (my_id == 0) THEN
    IF (ymode(1:1) == 'r') THEN

      OPEN (nudat, FILE=TRIM(datname), FORM='UNFORMATTED', STATUS='OLD', &
                   ACTION='READ', IOSTAT=ierror)
    ELSEIF ((ymode(1:1) == 'w') .OR. (ymode(1:1) == 'a')) THEN

      ! there is no append for restart-files. They are NEW
      OPEN (nudat, FILE=TRIM(datname), FORM='UNFORMATTED', STATUS='NEW', &
                   ACTION='WRITE', IOSTAT=ierror)

    ELSE
      ierror = -1
    ENDIF
  ENDIF

CASE DEFAULT

  IF (my_id == 0) THEN
    ierror = -10
    PRINT *, ' *** ERROR: File format ', yformat, ' is not known to LM! ***'
  ENDIF

END SELECT


! The error message has to be broadcasted to all other PEs
IF (npes > 1) THEN
  CALL MPI_BCAST (ierror, 1, MPI_INTEGER, 0, icomm, implcode)
  IF (ierror /= 0) RETURN
ENDIF

END SUBROUTINE open_file

!==============================================================================
!+ Module procedure in "io_utilities" for closing a file
!------------------------------------------------------------------------------

SUBROUTINE close_file (nudat, yformat, icomm, my_id, npes, lasync,    &
                       idbglv, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Close a file. The calling routine has to specify the format of the file.
!   Depending on this format, the necessary closing routines are called.
!   At the moment, asynchronous IO is implemented only for the Grib1 format.
!   This is a collective routine, i.e. it must be called by all processors
!
! Method:
!   Calls to special libraries.
!
!------------------------------------------------------------------------------
!
! Scalar arguments with intent(in):
INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
  nudat, & ! internal file descriptor
  icomm, & ! MPI communicator
  my_id, & ! id in communicator icomm
  npes,  & ! number of PEs
  idbglv   ! for verbosity of output

CHARACTER (LEN= 4)       , INTENT(IN)    ::  &
  yformat     ! format of the file (grib1, netcdf, binary)

LOGICAL                  , INTENT(IN)    ::  &
  lasync   ! indicates whether asynchronous IO is used or not

! Scalar arguments with intent(out):
CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg     ! string for error messages

INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ierror      ! error status

!------------------------------------------------------------------------------

! Local scalars:
INTEGER  (KIND=iintegers)           :: &
  implcode                               ! Error variable for MPI

INTEGER  (KIND=intgribc)            ::  &
  nudatc,                               & ! file descriptor for C-routine
  ierrc                                   ! error code for C-routine

! End of Header
!------------------------------------------------------------------------------

ierror   = 0
implcode = 0
yerrmsg  = '         '

IF (my_id == 0 .AND. idbglv > 1) print *,'CLOSING ', yformat, ' FILE'

SELECT CASE (yformat)

! At the moment, asynchronous IO is possible only for grib1-files
! I/O for all other formats is handled by processor 0

#ifdef GRIBDWD
CASE ('grb1')

  IF (lasync .OR. (npes > 1) ) THEN
    ! Close file using the parallel IO interface
    CALL mpe_io_close( nudat , implcode)
  ELSE
    ! Close file using griblib calls
    IF (my_id == 0) THEN
      nudatc = nudat
      CALL cclose(nudatc,'exi',ierrc)
      implcode = INT (ierrc, iintegers)
    ENDIF
  ENDIF
#endif

#ifdef NETCDF
CASE ('ncdf')

  IF (my_id == 0) THEN
    implcode = nf90_close(nudat)

    IF (implcode /= NF90_NOERR) THEN
      PRINT *, TRIM(NF90_strerror(implcode))
    ENDIF
  ENDIF
#endif

CASE ('bina')

  IF (my_id == 0) THEN
    CLOSE (nudat, IOSTAT=implcode)
  ENDIF

CASE DEFAULT

END SELECT

  IF (implcode /= 0) THEN
    yerrmsg = 'Error closing '//yformat//' file'
    ierror  = 3
    RETURN
  ENDIF

END SUBROUTINE close_file

!==============================================================================
!==============================================================================
!+ Module procedure in "io_utilities" for reading a grib record from a file
!------------------------------------------------------------------------------

SUBROUTINE read_grib (nudat, maxlen, ilfd, icomm, data, ilen, npes,         &
                      lasync, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Processor 0 reads all records and distributes them to the other processors. 
!   As many records are read from the grib file as there are processors 
!   present.
!
! Some Remarks:  (from R. Johanni, SGI/Cray):
!   "maxlen" and "ilen" are lengthes in BYTES !!!!
!   "data" should also be considered as a BYTE-array since all the I/O Routines
!   are dealing with BYTE-numbers only.
!   Also the GRIB1 standard explicitly only speaks about bytes and
!   never about integer words - GRIB1 files therefore should always
!   be considered as byte-streams!
!   "data" is declared as (KIND=intgribf) here, since more obvious data types
!   (as CHARACTER or INTEGER*1 would be) lead to numerous problems with
!   different compilers (esp. CRAY).
!
!------------------------------------------------------------------------------
!
! Scalar arguments with intent(in):
INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
  nudat,     & ! internal file descriptor
  maxlen,    & ! max. length of data array
  ilfd,      & ! dimension for grib routines
  icomm,     & ! MPI communicator
  npes         ! number of PE

! Scalar arguments with intent(out):
INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ilen         ! length of data actually read

! Array arguments with intent(out):
INTEGER  (KIND=intgribf ), INTENT(OUT)   ::  &
  data(*)      ! array to be filled

LOGICAL                  , INTENT(IN)    ::  &
  lasync   ! indicates whether asynchronous IO is used or not

CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg      ! string for error messages

INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ierror       ! error status

!------------------------------------------------------------------------------

! Local scalars
INTEGER  (KIND=iintegers)                ::  &
  my_comm_id,    & ! id in communicator icomm
  ipr,           & ! loop counter
  ilr,           & ! record length read
  iw,            & ! for reading controlword
  implcode         ! Error variable for MPI

INTEGER  (KIND=iintegers)                ::  &
  istatus(MPI_STATUS_SIZE)                     ! for status in MPI calls

INTEGER  (KIND=intgribf )                ::  &
  ibuf(ilfd)

INTEGER  (KIND=intgribc)                 ::  &
  nudatc, maxlenc, ilenc, ierrc     ! corresponding variables for C-routines

!- End of header
!------------------------------------------------------------------------------

#ifdef GRIBDWD
ierrc    = 0_intgribc
ierror   = 0
yerrmsg  = '         '

! Safety first
IF ( maxlen > iwlength*ilfd ) THEN
  yerrmsg = 'lfd too small'
  ierror  = 1
  RETURN
ENDIF

! This is in all cases a collective Routine, therefore:
IF (npes > 1) THEN
  CALL MPI_BARRIER(icomm,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'MPI_BARRIER failed'
    ierror  = 2
    RETURN
  ENDIF

  ! Get id in communicator comm
  CALL MPI_COMM_RANK(icomm,my_comm_id,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'MPI_COMM_RANK failed'
    ierror  = 3
    RETURN
  ENDIF
ELSE
  my_comm_id = 0
ENDIF

! Read using DWD library calls
nudatc  = INT (nudat, intgribc)
maxlenc = INT (maxlen, intgribc)

IF(my_comm_id == 0) THEN

  ! Processor 0 reads the data and sends it to all others
  ! Every processor gets a whole record (as long as there are enough
  ! records
  DO ipr = 0, npes-1
    IF (ipr == 0) THEN
      ! Read record and keep it for processing
      IF (lasync .OR. (npes > 1) ) THEN
        CALL mpe_io_read(nudat, data, ilen, maxlen, ierror)
      ELSE
        CALL cuegin (nudatc, maxlenc, data, ilenc, ierrc)
        ilen = INT (ilenc, iintegers)
        IF ( ierrc /= 0 ) THEN
          yerrmsg = 'Error in cuegin'
          ierror  = 5
          RETURN
        ENDIF
      ENDIF
      ilr = ilen
    ELSE
      IF (ilr > 0) THEN
        ! Read record and send it to other subdomains for processing
        IF (lasync .OR. (npes > 1) ) THEN
          CALL mpe_io_read(nudat, ibuf, ilr, maxlen, ierror)
        ELSE
          CALL cuegin (nudatc, maxlenc, ibuf, ilenc, ierrc)
          ilr = INT (ilenc, iintegers)
          IF ( ierrc /= 0 ) THEN
            yerrmsg = 'Error in cuegin'
            ierror  = 5
            RETURN
          ENDIF
        ENDIF
      ENDIF
      ! Because of the little endian problem, send 4 additional bytes
      CALL MPI_SEND(ilr,1,MPI_INTEGER ,ipr,1001,icomm,implcode)
      IF ( implcode /= 0 ) THEN
        yerrmsg = 'Error in MPI_SEND length of record'
        ierror  = 6
        RETURN
      ENDIF
      IF (ilr > 0) THEN  ! ilr == 0 means end of file
        CALL MPI_SEND(ibuf,ilr+4,MPI_BYTE,ipr,1002,icomm,implcode)
        IF ( implcode /= 0 ) THEN
          yerrmsg = 'Error in MPI_SEND record'
          ierror  = 7
          RETURN
        ENDIF
      ENDIF
    ENDIF
  ENDDO

ELSE

  ! All other Processors just receive the data
  CALL MPI_RECV(ilen,1,MPI_INTEGER ,0,1001,icomm,istatus,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'Error in MPI_RECV length of record'
    ierror  = 8
    RETURN
  ENDIF
  ! Because of the little endian problem, read 4 more bytes
  IF (ilen > 0) THEN
    CALL MPI_RECV(data,ilen+4,MPI_BYTE,0,1002,icomm,istatus,implcode)
    IF ( implcode /= 0 ) THEN
      yerrmsg = 'Error in MPI_RECV record'
      ierror  = 9
      RETURN
    ENDIF
  ENDIF

ENDIF

IF (npes > 1) THEN
  CALL MPI_BARRIER(icomm,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'Error in MPI_BARRIER'
    ierror  = 10
    RETURN
  ENDIF
ENDIF
#endif

END SUBROUTINE read_grib

!==============================================================================
!==============================================================================
!+ Module procedure in "io_utilities" for reading a grib record from a file
!------------------------------------------------------------------------------

SUBROUTINE read_netcdf (ncid, idim, jdim, ivar_count, ilev_count, ivar_id,  &
                        nvarin, idims, ndims, icomm, my_id, npes,  record,  &
                        myvar, mylev, mylevtot, irsize, idattyp,            &
                        lasync, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Processor 0 reads all records and distributes them to the other processors. 
!   As many records are read from the file as there are processors present.
!
!   In case of NetCDF IO, the reading can be done in an ordered way, because
!   it is known which variables are in the file and at which location. This
!   has been investigated by PE 0 in the subroutine read_nc_vdefs.
!   So the variables are read in the order in which they appear in the list
!   of the input variables (listin). ivar_count gives the index of the 
!   variable that has to be processed.
!
!------------------------------------------------------------------------------
!
! Scalar arguments with intent(in):
INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
  ncid,        & ! internal file descriptor
  idim, jdim,  & ! dimensions of the 2D record
  nvarin,      & ! number of elements in the input variable list
  ndims,       & ! dimension of idims
  icomm,       & ! MPI communicator
  my_id,       & ! ID of this PE in communicator icomm
  npes,        & ! number of PEs
  idattyp,     & ! type of data for MPI calls
  idims(ndims),& ! for the dimension IDs
  ivar_id(nvarin)! array with the variable IDs

! Scalar arguments with intent(inout):
INTEGER  (KIND=iintegers), INTENT(INOUT) ::  &
  ivar_count,  & ! actual index in the variable list
  ilev_count     ! actual level of a 3D variable

! Array arguments with intent(out):
INTEGER  (KIND=intgribf ), INTENT(OUT)   ::  &
  myvar,       & ! index of the variable this PE has to process (in listin)
  mylev,       & ! level of the 3D variable this PE has to process
  mylevtot,    & ! total levels of the variable this PE has to process
  irsize         ! length of data actually read

REAL     (KIND=irealgrib), INTENT(OUT)   ::  &
  record(idim*jdim) ! data to be read

LOGICAL                  , INTENT(IN)    ::  &
  lasync         ! indicates whether asynchronous IO is used or not

CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg        ! string for error messages

INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ierror         ! error status

!------------------------------------------------------------------------------

! Local scalars
INTEGER  (KIND=iintegers)                ::  &
  ipr,           & ! loop counter
  iedims, jedims,& ! dimensions read from the NetCDF file
  ilen,          & ! length of a record
  ilenbuf,       & ! length of the buffer to receive data
  ndimids,       & ! number of dimension id's
  dimsids(4),    & ! for NetCDF dimension id's
  klev,          & ! number of vertical levels of a variable
  implcode         ! Error variable for MPI

INTEGER  (KIND=iintegers)                ::  &
  istatus(MPI_STATUS_SIZE)                     ! for status in MPI calls

REAL (KIND=irealgrib)                    ::  &
  rbuf(idim*jdim+3), reof

!- End of header
!------------------------------------------------------------------------------

#ifdef NETCDF
ierror   = 0
yerrmsg  = '         '
reof     = -999999.0_irealgrib
ilen     = idim * jdim
ilenbuf  = idim * jdim + 3

IF(my_id == 0) THEN

  ! Processor 0 reads the data and sends it to all others
  ! Every processor gets a 2D record (as long as there are enough records)

  processor_loop: DO ipr = 0, npes-1

  IF (ivar_count <= nvarin) THEN
    increase_loop: DO WHILE (ivar_id(ivar_count) == -1)
      ivar_count = ivar_count + 1
      IF (ivar_count > nvarin) THEN
        ! we are at the end
        EXIT increase_loop
      ENDIF
    ENDDO increase_loop
  ENDIF

  IF (ivar_count <= nvarin) THEN

    ! there is still something to read
    ! get number of dimension id's for variable ivar_count
    !  ndimids = 2: means 3 dimensions; a 2D space array
    !  ndimids = 3: means 3 dimensions; a 2D space array + time
    !  ndimids = 4: means 4 dimensions; a 3D space array + time
    ierror = NF90_inquire_variable (ncid, ivar_id(ivar_count), ndims=ndimids)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      PRINT *, '*** ERROR in read_netcdf:  ', ierror, TRIM(yerrmsg)
      RETURN
    ENDIF

    dimsids(:) = 0_iintegers

    ! get the values of the dimension ids: iedims, jedims, klev
    ierror = NF90_inquire_variable (ncid, ivar_id(ivar_count),     &
                                                 dimids = dimsids(1:ndimids))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      PRINT *, '*** ERROR in read_netcdf:  ', ierror, TRIM(yerrmsg)
      RETURN
    ENDIF

    ierror = NF90_inquire_dimension (ncid, dimsids(1), len=iedims)
    ierror = NF90_inquire_dimension (ncid, dimsids(2), len=jedims)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      PRINT *, '*** ERROR in read_netcdf:  ', ierror, TRIM(yerrmsg)
      RETURN
    ENDIF
    IF ( (iedims /= idim) .OR. (jedims /= jdim) ) THEN
      ierror  = 10
      yerrmsg = '  got wrong dimensions from NetCDF file ***'
      PRINT *, '*** ERROR in read_netcdf:  ', ierror, TRIM(yerrmsg)
      PRINT *, '      iedim = ', iedims, ';    jedim = ', jedims
      RETURN
    ENDIF

    ! Determine the number of vertical levels, dependent on the
    ! 3rd value of the dimension IDs
    ! (has been extended to also read the 3D external parameter
    ! field for the topographical corrections (by Anne Roches)
    IF (ndimids == 2) THEN
      ! this is a 2D field
      klev = 1
    ELSEIF (  (ndimids < 4) .OR. ( (dimsids(3) /= idims( 3)) .AND.   &
                                   (dimsids(3) /= idims( 4)) .AND.   &
                                   (dimsids(3) /= idims( 7)) .AND.   &
                                   (dimsids(3) /= idims( 8)) .AND.   &
                                   (dimsids(3) /= idims(11)) ) ) THEN
      ! this is a 2D field with time levels or for boundaries
      klev = 1
    ELSEIF ( (dimsids(3) == idims( 3))  .OR.                     &
             (dimsids(3) == idims( 4))  .OR.                     &
             (dimsids(3) == idims( 7))  .OR.                     &
             (dimsids(3) == idims( 8))  .OR.                     &
             (dimsids(3) == idims(11)) ) THEN
      ! this is a 3D field for the atmosphere (3/4), for the soil (7/8)
      ! or for the topographical corrections (11)
      ierror = NF90_inquire_dimension (ncid, dimsids(3), len=klev)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        PRINT *, '*** ERROR in read_netcdf:  ', ierror, TRIM(yerrmsg)
        RETURN
      ENDIF
    ENDIF

    ! read the data, depending on its dimension in space
    IF (ndimids == 2) THEN
      ! this is a 2D array
      ierror = NF90_get_var (ncid, ivar_id(ivar_count), rbuf,         &
                             start=(/1,1/), count=(/iedims, jedims/))
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        PRINT *, '*** ERROR in read_netcdf:  ', ierror, TRIM(yerrmsg)
        RETURN
      ENDIF
    ELSEIF (ndimids == 3) THEN
      ! this is a 3D array (2D space + 1D time)
      ierror = NF90_get_var (ncid, ivar_id(ivar_count), rbuf,         &
                             start=(/1,1,1/), count=(/iedims, jedims, 1/))
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        PRINT *, '*** ERROR in read_netcdf:  ', ierror, TRIM(yerrmsg)
        RETURN
      ENDIF
    ELSEIF (ndimids == 4) THEN
      ! this is a 4D array (3D space + 1D time)
      ilev_count = ilev_count + 1
      ierror = NF90_get_var (ncid, ivar_id(ivar_count), rbuf,         &
                    start=(/1,1,ilev_count,1/), count=(/iedims,jedims,1,1/))
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        PRINT *, '*** ERROR in read_netcdf:  ', ierror, TRIM(yerrmsg)
        RETURN
      ENDIF
    ENDIF
    irsize = ilen

    ! Distribute the record to other processors (or keep it for yourself)
    IF (ipr == 0) THEN
      record(1:ilen) = rbuf(1:ilen)
      myvar    = ivar_count
      IF (ilev_count == 0) THEN
        ! this is a 2D array
        mylev    = 1
      ELSE
        ! this is a 3D array
        mylev    = ilev_count
      ENDIF
      mylevtot = klev
    ELSE
      rbuf(ilen + 1) = REAL(ivar_count, irealgrib)
      IF (ilev_count == 0) THEN
        ! this is a 2D array
        rbuf(ilen + 2) = 1
      ELSE
        ! this is a 3D array
        rbuf(ilen + 2) = REAL(ilev_count, irealgrib)
      ENDIF
      rbuf(ilen + 3) = REAL(klev      , irealgrib)
      CALL MPI_SEND(rbuf, ilen+3, idattyp, ipr, 1002, icomm, implcode)
      IF ( implcode /= 0 ) THEN
        yerrmsg = 'Error in MPI_SEND record'
        ierror  = 7
        RETURN
      ENDIF
    ENDIF

    ! reset ilev_count and increase ivar_count, if necessary
    IF (ilev_count == klev) THEN
      ! reset ilev_count
      ilev_count = 0
    ENDIF
    IF (ilev_count == 0) THEN
      ! either a 2D field has been read or all levels of a 3D field
      ivar_count = ivar_count + 1
    ENDIF

  ELSE  !  (ivar_count > nvarin)

    ! the end of the file is reached: send reof to the rest of the processors
    IF (ipr == 0) THEN
      irsize = 0
    ELSE
      CALL MPI_SEND (reof, 1, idattyp, ipr, 1002, icomm, implcode)
      IF ( implcode /= 0 ) THEN
        yerrmsg = 'Error in MPI_SEND record'
        ierror  = 7
        RETURN
      ENDIF
    ENDIF

  ENDIF

  ENDDO processor_loop

ELSE   ! (my_id /= 0)

  ! All other Processors just receive the data
  CALL MPI_RECV(rbuf, ilenbuf, idattyp, 0, 1002, icomm, istatus, implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'Error in MPI_RECV record'
    ierror  = 9
    RETURN
  ENDIF

  IF (rbuf(1) /= reof) THEN
    record(1:ilen) = rbuf(1:ilen)
    myvar          = rbuf(ilen + 1)
    mylev          = rbuf(ilen + 2)
    mylevtot       = rbuf(ilen + 3)
    irsize   = ilen
  ELSE
    irsize   = 0
  ENDIF

ENDIF
#endif

END SUBROUTINE read_netcdf

!==============================================================================
!==============================================================================
!+ Module procedure in "io_utilities" for reading a grib record from a file
!------------------------------------------------------------------------------

SUBROUTINE read_restart (nudat, idim, jdim, record, ipds, npds, igds, ngds, &
                         icomm, my_id, npes, irsize, irealtyp, igribtyp,    &
                         lasync, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!
!------------------------------------------------------------------------------

INTEGER (KIND=iintegers),  INTENT(IN)     :: &
  nudat,       & ! unit number
  idim, jdim,  & ! dimensions of the 2D record
  icomm,       & ! MPI communicator
  my_id,       & ! ID of this PE in communicator icomm
  npes,        & ! number of PEs
  irealtyp,    & ! type of real data for MPI calls
  igribtyp       ! type of real data for MPI calls

INTEGER (KIND=intgribf),   INTENT(IN)     :: &
  npds, ngds     ! dimension of pds and gds

INTEGER (KIND=intgribf),   INTENT(OUT)    :: &
  ipds(npds),  & ! pds: product definition section
  igds(ngds)     ! pds: product definition section

INTEGER (KIND=iintegers),  INTENT(OUT)    :: &
  irsize         ! length of data actually read

REAL     (KIND=ireals),    INTENT(OUT)   ::  &
  record(idim*jdim) ! data to be read

LOGICAL                  , INTENT(IN)    ::  &
  lasync         ! indicates whether asynchronous IO is used or not

CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg        ! string for error messages

INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ierror         ! error status

!------------------------------------------------------------------------------

! Local scalars
INTEGER  (KIND=iintegers)                ::  &
  ipr,           & ! loop counter
  ilen, ilr,     & ! length of a record
  ilenbuf,       & ! length of the buffer to receive data
  izerr,         & ! local error variable
  implcode         ! Error variable for MPI

INTEGER  (KIND=iintegers)                ::  &
  istatus(MPI_STATUS_SIZE)                     ! for status in MPI calls

REAL (KIND=ireals)                       ::  &
  rbuf(idim*jdim), reof

INTEGER(KIND=intgribf)                   ::  &
  ipdsbuf(npds), igdsbuf(ngds)

LOGICAL                                  ::  &
  lzeof            ! for end-of-file

!- End of header
!------------------------------------------------------------------------------

ierror   = 0
yerrmsg  = '         '
reof     = -999999.0_irealgrib
izerr    = 0_iintegers
implcode = 0_iintegers
ilen     = idim * jdim
ilenbuf  = idim * jdim

IF(my_id == 0) THEN

  ! Processor 0 reads the data and sends it to all others
  ! Every processor gets a 2D record (as long as there are enough records)

  lzeof = .FALSE.

  processor_loop: DO ipr = 0, npes-1

    IF (.NOT. lzeof) THEN
      READ (nudat, IOSTAT=izerr) ipdsbuf, igdsbuf, rbuf

      IF (izerr < 0) THEN
        ! this is the end-of-file condition
        lzeof  = .TRUE.
        ilr    = 0
        IF (ipr > 0) THEN
          ! send reof to processor ipr
          CALL MPI_SEND (reof, 1, irealtyp, ipr, 1002, icomm, implcode)
          IF ( implcode /= 0 ) THEN
            yerrmsg = 'Error in MPI_SEND record'
            ierror  = 2
            RETURN
          ENDIF
        ELSE
          ! set irsize for processor 0 to eof:
          irsize = ilr
        ENDIF
      ELSEIF (izerr > 0) THEN
        ! an error occured
        ierror  = 1
        yerrmsg = 'error while reading restart file'
        RETURN
      ELSE  
        ! no error and no end-of-file occured

        ! Distribute the record to other processors (or keep it for yourself)
        IF (ipr == 0) THEN
          record(1:ilen) = rbuf(1:ilen)
          ipds(1:npds)   = ipdsbuf(1:npds)
          igds(1:ngds)   = igdsbuf(1:ngds)
          irsize         = ilen
        ELSE
          CALL MPI_SEND(rbuf,    ilen, irealtyp, ipr, 1002, icomm, implcode)
          CALL MPI_SEND(ipdsbuf, npds, igribtyp, ipr, 1003, icomm, implcode)
          CALL MPI_SEND(igdsbuf, ngds, igribtyp, ipr, 1004, icomm, implcode)
          IF ( implcode /= 0 ) THEN
            yerrmsg = 'Error in MPI_SEND record'
            ierror  = 7
            RETURN
          ENDIF
        ENDIF
      ENDIF

    ELSE

      ! send reof to all other processors
      IF (ipr /= 0) THEN
        CALL MPI_SEND (reof, 1, irealtyp, ipr, 1002, icomm, implcode)
        IF ( implcode /= 0 ) THEN
          yerrmsg = 'Error in MPI_SEND record'
          ierror  = 3
          RETURN
        ENDIF
      ENDIF

    ENDIF

  ENDDO processor_loop

ELSE   ! (my_id /= 0)

  ! All other Processors just receive the data
  CALL MPI_RECV(record, ilenbuf, irealtyp, 0, 1002, icomm, istatus, implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'Error in MPI_RECV record'
    ierror  = 9
    RETURN
  ENDIF

  IF (record(1) /= reof) THEN
    irsize = ilen
    ! also receive the pds and gds
    CALL MPI_RECV(ipds, npds, igribtyp, 0, 1003, icomm, istatus, implcode)
    CALL MPI_RECV(igds, ngds, igribtyp, 0, 1004, icomm, istatus, implcode)
  ELSE
    irsize = 0
  ENDIF

ENDIF

END SUBROUTINE read_restart

!==============================================================================
!==============================================================================
!+ Module procedure in "io_utilities" for writing a grib file
!------------------------------------------------------------------------------

SUBROUTINE write_grib (nudat, data, ilen, ilfd, icomm, npes, lflush,     &
                       ydbtype, lasync, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Processor 0 gathers all fields from the other processors and writes them
!   to disk.
!
! Some Remarks:  (from R. Johanni, SGI/Cray):
!   "ilen" is in BYTES !!!!
!   "data" should also be considered as a BYTE-array since all the I/O Routines
!   are dealing with BYTE-numbers only.
!   Also the GRIB1 standard explicitly only speaks about bytes and
!   never about integer words - GRIB1 files therefore should allways
!   be considered as byte-streams!
!   "data" is declared as (KIND=intgribf) here, since more obvious data types
!   (as CHARACTER or INTEGER*1 would be) lead to numerous problems with
!   different compilers (esp. CRAY).
!
!------------------------------------------------------------------------------
!
! Subroutine arguments
! Scalar arguments with intent(in):
INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
  nudat,     & ! internal file descriptor
  ilfd,      & ! dimension for grib routines
  icomm,     & ! MPI communicator
  npes         ! number of PE

INTEGER  (KIND=iintegers), INTENT(INOUT) ::  &
  ilen         ! length of data actually read

! Array arguments with intent(in):
INTEGER  (KIND=intgribf ), INTENT(INOUT) ::  &
  data(*)      ! array to be written

LOGICAL                  , INTENT(IN)    ::  &
  lasync,    & ! indicates whether asynchronous IO is used or not
  lflush       ! at the end of the writing

CHARACTER (LEN= 5)       , INTENT(IN)    ::  &
  ydbtype      ! type of database

! Scalar arguments with intent(out):
CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg      ! string for error messages

INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ierror       ! error status

!------------------------------------------------------------------------------

! Local scalars
INTEGER  (KIND=iintegers)                ::  &
  my_comm_id,    & ! id in communicator icomm
  ipr,           & ! loop counter
  ilr,           & ! record length to write
  knudat,        & ! just as nudat, but this is INOUT in mpe_io_write
  implcode         ! Error variable for MPI

INTEGER  (KIND=iintegers)                ::  &
  istatus(MPI_STATUS_SIZE)                     ! for status in MPI calls

INTEGER  (KIND=intgribf )                ::  &
  ibuf(ilfd)

INTEGER  (KIND=intgribc)                 ::  &
  nudatc, ilenc, ierrc       ! corresponding variables for C-routines

!- End of header
!------------------------------------------------------------------------------

#ifdef GRIBDWD
ierror   = 0
yerrmsg  = '         '

! Safety first
IF ( ilen > iwlength*ilfd ) THEN
  yerrmsg = 'lfd too small'
  ierror  = 1
  RETURN
ENDIF

! This is in all cases a collective Routine, therefore:
IF (npes > 1) THEN
  CALL MPI_BARRIER(icomm,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'MPI_BARRIER failed'
    ierror  = 2
    RETURN
  ENDIF

  ! Get id in communicator comm and number of processors
  CALL MPI_COMM_RANK(icomm,my_comm_id,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'MPI_COMM_RANK failed'
    ierror  = 3
    RETURN
  ENDIF
ELSE
  my_comm_id = 0
ENDIF

knudat = nudat  ! Necessary as knudat is INOUT in mpe_io_write
nudatc = nudat

IF (lasync .OR. (npes > 1) ) THEN
  CALL mpe_io_write(knudat, data, ilen, ilfd, ydbtype, ierror)
  IF( lflush ) THEN
    CALL mpe_io_complete()
  ENDIF
ELSE
  ! this is a pure sequential program without any MPI
  ierrc = 0 ! should be set in cuegex, but is not !!!!!!
  ilenc = ilen
  CALL cuegex (nudatc, data, ilenc, ierrc)
  IF ( ierrc /= 0 ) THEN
    yerrmsg = 'Error in cuegex'
    ierror  = 5
    RETURN
  ENDIF
ENDIF

IF (npes > 1) THEN
  CALL MPI_Barrier(icomm,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'Error in MPI_BARRIER'
    ierror  = 11
    RETURN
  ENDIF
ENDIF
#endif

END SUBROUTINE write_grib

!==============================================================================
!+ Module procedure in "io_utilities" for writing a netcdf file
!------------------------------------------------------------------------------

SUBROUTINE write_netcdf (nudat, data, i_tot, j_tot, irec_len, ncorg, icomm,  &
                         my_id, npes, idattyp, lasync, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Processor 0 gathers all fields from the other processors and writes them
!   to netcdf files.
!
!------------------------------------------------------------------------------
!
! Subroutine arguments
! Scalar arguments with intent(in):
INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
  nudat,     & ! internal file descriptor
  i_tot,     & ! i-dimension for the array data
  j_tot,     & ! j-dimension for the array data
  irec_len,  & ! length of the record (0 if no record available)
  icomm,     & ! MPI communicator
  my_id,     & ! Processor ID in the communicator icomm
  npes,      & ! number of total PEs
  idattyp      ! data type for MPI interface

INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
  ncorg(3,0:npes-1) 
               ! organizational variables for netcdf. Only Processor 0 has 
               ! correct information, which is necessary for writing the data
               !   ncorg(1,:): actual number of level
               !   ncorg(2,:): maximum number of levels for this variable
               !   ncorg(3,:): netcdf ID of this variable

! Array arguments with intent(in):
REAL (KIND=irealgrib),     INTENT(IN)    ::  &
  data(i_tot*j_tot)   
               ! array to be written

LOGICAL                  , INTENT(IN)    ::  &
  lasync       ! indicates whether asynchronous IO is used or not
               ! (not used at the moment)

! Scalar arguments with intent(out):
CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg      ! string for error messages

INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ierror       ! error status

!------------------------------------------------------------------------------

! Local scalars
INTEGER  (KIND=iintegers)                ::  &
  ipr,           & ! loop counter
  ilen,          & ! i_tot*j_tot
  implcode,      & ! Error variable for MPI
  izerr            ! another error variable

INTEGER  (KIND=iintegers)                ::  &
  istatus(MPI_STATUS_SIZE)                     ! for status in MPI calls

REAL(KIND=irealgrib)                     ::  &
  rbuf(i_tot*j_tot), & ! buffer for receiving the data
  reof                 ! for sending end-of-file message

!- End of header
!------------------------------------------------------------------------------

#ifdef NETCDF
ilen     = i_tot * j_tot
izerr    = 0_iintegers
ierror   = 0
reof     = -999999.0_irealgrib
yerrmsg  = '         '

IF(my_id == 0) THEN

  ! Processor 0 receives the data and writes it out
  DO ipr = 0, npes-1

    IF (ipr == 0) THEN

     IF (irec_len > 0) THEN
       ! Processor 0 has the data already and can write it to disk
       IF     (ncorg(2,ipr) == 1) THEN
         ! this is a 2D variable
         izerr = nf90_put_var (nudat, ncorg(3,ipr), data,              &
                    start=(/ 1, 1, 1 /), count=(/ i_tot, j_tot, 1 /))
         IF (izerr /= nf90_noerr) THEN
           yerrmsg = 'Error writing netcdf 2D variable'
           ierror  = 1
           RETURN
         ENDIF
       ELSEIF (ncorg(2,ipr) >  1) THEN
         ! this is a 2D-slice of a 3D variable
         izerr = nf90_put_var (nudat, ncorg(3,ipr), data,              &
            start=(/ 1, 1, ncorg(1,ipr),1 /), count=(/ i_tot, j_tot, 1, 1 /))
         IF (izerr /= nf90_noerr) THEN
           yerrmsg = 'Error writing netcdf 3D variable'
           ierror  = 2
           RETURN
         ENDIF
       ENDIF
     ENDIF

    ELSE

      ! Processor 0 gets the data from the other processors
      CALL MPI_RECV(rbuf, ilen, idattyp, ipr, 1002, icomm, istatus, implcode)
      IF ( implcode /= 0 ) THEN
        yerrmsg = 'Error in MPI_RECV record'
        ierror  = 3
        RETURN
      ENDIF

      IF (rbuf(1) /= reof) THEN
        ! processor ipr really sent a field that has to be written
        ! (otherwise rbuf(1) = reof: then just an end message is sent)
        IF     (ncorg(2,ipr) == 1) THEN
          ! this is a 2D variable
          izerr = nf90_put_var (nudat, ncorg(3,ipr), rbuf,             &
                     start=(/ 1, 1, 1 /), count=(/ i_tot, j_tot, 1 /))
          IF (izerr /= nf90_noerr) THEN
            yerrmsg = 'Error writing netcdf 2D variable'
            ierror  = 4
            RETURN
          ENDIF
        ELSEIF (ncorg(2,ipr) >  1) THEN
          ! this is a 2D-slice of a 3D variable
          izerr = nf90_put_var (nudat, ncorg(3,ipr), rbuf,             &
             start=(/ 1, 1, ncorg(1,ipr),1 /), count=(/ i_tot, j_tot, 1, 1 /))
          IF (izerr /= nf90_noerr) THEN
            yerrmsg = 'Error writing netcdf 3D variable'
            ierror  = 5
            RETURN
          ENDIF
        ENDIF
      ENDIF

    ENDIF

  ENDDO

ELSE

  ! All other Processors just send the data
  IF (irec_len > 0) THEN
    CALL MPI_SEND(data, ilen, idattyp, 0, 1002, icomm, implcode)
    IF ( implcode /= 0 ) THEN
      yerrmsg = 'Error in MPI_SEND record'
      ierror  = 6
      RETURN
    ENDIF
  ELSE
    CALL MPI_SEND(reof, 1, idattyp, 0, 1002, icomm, implcode)
    IF ( implcode /= 0 ) THEN
      yerrmsg = 'Error in MPI_SEND end message'
      ierror  = 7
      RETURN
    ENDIF
  ENDIF

ENDIF
#endif


END SUBROUTINE write_netcdf

!==============================================================================
!+ Module procedure in "io_utilities" for writing a binary restart file
!------------------------------------------------------------------------------

SUBROUTINE write_restart (nudat, data, i_tot, j_tot, irec_len, ipds_arr,     &
                          npds_dim, igds_arr, ngds_dim, icomm, my_id, npes,  &
                          idattyp, lasync, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   Processor 0 gathers all fields from the other processors and writes them
!   to a binary restart file
!
!------------------------------------------------------------------------------
!
! Subroutine arguments
! Scalar arguments with intent(in):
INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
  nudat,     & ! internal file descriptor
  i_tot,     & ! i-dimension for the array data
  j_tot,     & ! j-dimension for the array data
  irec_len,  & ! length of the record (0 if no record available)
  icomm,     & ! MPI communicator
  my_id,     & ! Processor ID in the communicator icomm
  npes,      & ! number of total PEs
  idattyp      ! data type for MPI interface

INTEGER  (KIND=intgribf),  INTENT(IN)    ::  &
  npds_dim,  & ! dimension for the product definition section
  ngds_dim     ! dimension for the grid definition section

INTEGER  (KIND=intgribf),  INTENT(IN)    ::  &
  ipds_arr(npds_dim), & ! product definition section
  igds_arr(ngds_dim)    ! grid definition section

! Array arguments with intent(in):
REAL (KIND=ireals),        INTENT(IN)    ::  &
  data(i_tot*j_tot)   
               ! array to be written

LOGICAL                  , INTENT(IN)    ::  &
  lasync       ! indicates whether asynchronous IO is used or not
               ! (not used at the moment)

! Scalar arguments with intent(out):
CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
  yerrmsg      ! string for error messages

INTEGER  (KIND=iintegers), INTENT(OUT)   ::  &
  ierror       ! error status

!------------------------------------------------------------------------------

! Local scalars
INTEGER  (KIND=iintegers)                ::  &
  ipr,           & ! loop counter
  ilen,          & ! i_tot*j_tot
  implcode,      & ! Error variable for MPI
  izerr            ! another error variable

INTEGER  (KIND=iintegers)                ::  &
  istatus(MPI_STATUS_SIZE)                     ! for status in MPI calls

REAL(KIND=ireals)                        ::  &
  rbuf(i_tot*j_tot),  & ! buffer for receiving the data
  reof                  ! to indicate end-of-file

INTEGER  (KIND=intgribf)                 ::  &
  ibuf_pds(npds_dim), & ! buffer for product definition section
  ibuf_gds(ngds_dim)    ! buffer for grid definition section

!- End of header
!------------------------------------------------------------------------------

ilen     = i_tot * j_tot
izerr    = 0_iintegers
ierror   = 0
reof     = HUGE (1.0_ireals)
yerrmsg  = '         '

IF(my_id == 0) THEN

  ! Processor 0 receives the data and writes it out
  DO ipr = 0, npes-1

    IF (ipr == 0) THEN

      IF (irec_len > 0) THEN
        ! Processor 0 has the data already and can write it to disk
        WRITE (nudat,IOSTAT=izerr) ipds_arr, igds_arr, data

        IF (izerr /= 0_iintegers) THEN
          yerrmsg = 'Error writing restart variable'
          ierror  = 1
          RETURN
        ENDIF
      ENDIF

    ELSE

      ! Processor 0 gets the data from the other processors
      CALL MPI_RECV (ibuf_pds, npds_dim, MPI_INTEGER, ipr, 1003, icomm,    &
                     istatus, implcode)
      IF ( implcode /= 0 ) THEN
        yerrmsg = 'Error in MPI_RECV pds'
        ierror  = 2
        RETURN
      ENDIF

      CALL MPI_RECV (ibuf_gds, ngds_dim, MPI_INTEGER, ipr, 1004, icomm,    &
                     istatus, implcode)
      IF ( implcode /= 0 ) THEN
        yerrmsg = 'Error in MPI_RECV gds'
        ierror  = 3
        RETURN
      ENDIF

      CALL MPI_RECV (rbuf,    ilen,      idattyp,     ipr, 1005, icomm,    &
                     istatus, implcode)
      IF ( implcode /= 0 ) THEN
        yerrmsg = 'Error in MPI_RECV record'
        ierror  = 4
        RETURN
      ENDIF

      IF (rbuf(1) /= reof) THEN
        ! processor ipr really sent a field that has to be written
        WRITE (nudat, IOSTAT=izerr) ibuf_pds, ibuf_gds, rbuf

        IF (izerr /= 0_iintegers) THEN
          yerrmsg = 'Error writing restart variable'
          ierror  = 5
          RETURN
        ENDIF

      ENDIF

    ENDIF

  ENDDO

ELSE

  ! All other Processors just send the data
  CALL MPI_SEND (ipds_arr, npds_dim, MPI_INTEGER, 0, 1003, icomm, implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'Error in MPI_SEND pds'
    ierror  = 6
    RETURN
  ENDIF

  CALL MPI_SEND (igds_arr, ngds_dim, MPI_INTEGER, 0, 1004, icomm, implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'Error in MPI_SEND gds'
    ierror  = 7
    RETURN
  ENDIF

  IF (irec_len > 0) THEN
    CALL MPI_SEND (data, ilen, idattyp, 0, 1005, icomm, implcode)
    IF ( implcode /= 0 ) THEN
      yerrmsg = 'Error in MPI_SEND record'
      ierror  = 8
      RETURN
    ENDIF
  ELSE
    CALL MPI_SEND (reof, 1, idattyp, 0, 1005, icomm, implcode)
    IF ( implcode /= 0 ) THEN
      yerrmsg = 'Error in MPI_SEND end message'
      ierror  = 9
      RETURN
    ENDIF
  ENDIF

ENDIF


END SUBROUTINE write_restart

!==============================================================================
!+ checks the input and the output data
!------------------------------------------------------------------------------

SUBROUTINE check_record                                                   &
           (field, idim1s, idim1e, idim2s, idim2e, idim3s, idim3e,        &
                   icom1s, icom1e, icom2s, icom2e, icom3s, icom3e, undef, &
            yname, iee, ilevel, lwork, nout, npe, icomm, myid, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!  Check_record calculates the minimum, the maximum and the meanvalue of a 
!  two-dimensional total field and prints them to a file with unit number
!  nout. Three dimensions have to be given for the two-dimensional field,
!  because this routine is also used for GME-records, where the third dimension
!  represents the diamond. In case of LM, the parameters of the third 
!  dimensions are 1. 
!  
!  This routine has to be called by all compute PEs. Every processor writes 
!  its results to a line of characters. These lines are gathered by compute
!  PE 0, which prints them.
!
! Method:
!  Fortran 90 intrinsic functions are used to determine the minimum, maximum
!  and the meanvalue of ufeld. The data from all processors are gathered with
!  MPI_GATHER to processor 0. For the T3E implementation, the characters have
!  to be transformed to integers (at the moment).
!
! Output files:
!  Ascii file with unit number nout
!
!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Scalar arguments with intent(in):

INTEGER (KIND=iintegers), INTENT(IN) ::  &
  idim1s, idim1e, idim2s, idim2e, idim3s, idim3e, & 
            ! start- and end-indices for the dimensions of field
  icom1s, icom1e, icom2s, icom2e, icom3s, icom3e
            ! start- and end-indices for the computations

INTEGER (KIND=intgribf),  INTENT(IN) ::  &
  iee,             & ! grib number                        "
  ilevel             ! level of this field

INTEGER (KIND=iintegers), INTENT(IN) ::  &
  nout               ! unit number of the output file

LOGICAL,                  INTENT(IN) ::  &
  lwork              ! is .TRUE. if there is work to do

INTEGER (KIND=iintegers), INTENT(IN) ::  &
  npe,             & ! number of pes in communicator icomm
  icomm,           & ! communicator of compute PEs
  myid               ! ID of this PE in communicator icomm

CHARACTER (LEN= *),       INTENT(IN) ::  &
  yname              ! name of the field in variable table

! Array arguments with intent(in):
REAL   (KIND=irealgrib),  INTENT(IN) ::  &
  field (idim1s:idim1e, idim2s:idim2e, idim3s:idim3e),  &
  undef              ! for points that are not defined by the bitmap

! Variables for output
CHARACTER (LEN= *),      INTENT(OUT) ::  &
  yerrmsg            ! for error message

INTEGER (KIND=iintegers),INTENT(OUT) ::  &
  ierror             ! error status

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER (KIND=iintegers)   ::   &
  j1min, j1max, j2min, j2max, jdmin, jdmax, idiv
                                  ! location of min and max in the 2D field

INTEGER (KIND=iintegers)   ::   &
  i, n, j1, j2, jd, implcode      ! Additional variables

REAL    (KIND=irealgrib)   ::   &
  zmin, zmax, zsum, zdiv,       & ! values of min, max and meanvalue
  undef_eps                       ! slightly changed undef-value to cope with
                                  ! numerical rounding
! Local arrays:
CHARACTER (LEN=200)        ::   &
  yline,                        & ! line to be written for one processor
  ylines(npe)                     ! array of lines to collect all lines from 
                                  ! the other processors

INTEGER                    ::   & ! really the standard integers!!
  iline (200),                  & ! the same for the integer values
  ilines(200,npe)   
!
!- End of header
!------------------------------------------------------------------------------

  ierror  = 0
  yerrmsg = '   '
  undef_eps = undef + 100.0_irealgrib

  j1min   = 1
  j2min   = 1
  j1max   = 1
  j2max   = 1

  IF (lwork) THEN
    ! if this processor has a required record then
    ! determine maximum and minimum with their locations
    zsum = 0.0_irealgrib
    zmin =   ABS (undef)
    zmax = - ABS (undef)
    idiv = 0_iintegers

    DO jd = icom3s, icom3e
      DO j2 = icom2s, icom2e
        DO j1 = icom1s, icom1e
          IF (field(j1,j2,jd) > undef_eps) THEN
            idiv = idiv + 1
            zsum = zsum + field(j1,j2,jd)
            IF (field(j1,j2,jd) > zmax) THEN
              zmax  = field(j1,j2,jd)
              j1max = j1
              j2max = j2
              jdmax = jd
            ENDIF
            IF (field(j1,j2,jd) < zmin) THEN
              zmin  = field(j1,j2,jd)
              j1min = j1
              j2min = j2
              jdmin = jd
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    IF (idiv /= 0_iintegers) THEN
      zdiv = REAL (idiv, ireals)
      zsum = zsum / zdiv
    ELSE
      ! this means that the whole field is undefined
      zsum = undef
    ENDIF

    ! Adjust the indices, if idim* is not equal to icom*, so that the
    ! indices are related to icom* for the first and the second dimension
    IF (icom1s > idim1s) THEN
      j1max = j1max - (icom1s-idim1s)
      j1min = j1min - (icom1s-idim1s)
    ENDIF
    IF (icom2s > idim2s) THEN
      j2max = j2max - (icom2s-idim2s)
      j2min = j2min - (icom2s-idim2s)
    ENDIF
   
    ! print the values to yline
    IF (icom3e == 1) THEN
      ! likely it is an LM field
      IF (zsum == undef) THEN
        WRITE(yline,'(A,A,I4,I4,A)')  '  ', yname, iee, ilevel, &
             '                       WHOLE FIELD UNDEFINED    '
      ELSE
        WRITE(yline,'(A,A,I4,I4,F15.6,2I5,A,F16.6,2I5,A,F16.7)')           &
                    '  ', yname, iee, ilevel, zmin, j1min, j2min, '     ', &
                    zmax, j1max, j2max, '     ', zsum
      ENDIF
    ELSE
      ! this really is an GME-field
      IF (zsum == undef) THEN
        WRITE(yline,'(A,A,I4,I4,A)')  '  ', yname, iee, ilevel, &
             '                       WHOLE FIELD UNDEFINED    '
      ELSE
        WRITE(yline,'(A,A,I4,I4,F15.6,3I5,F16.6,3I5,F16.7)')               &
                    '  ', yname, iee, ilevel, zmin, j1min, j2min, jdmin,   &
                    zmax, j1max, j2max, jdmax, zsum
      ENDIF
    ENDIF
  ELSE
    ! print just one blank to yline
    yline = ' '
  ENDIF
  
  ! collect all yline to ylines in processor 0
  IF (npe > 1) THEN
    ! for T3E's sake, convert the characters to standard integers
    DO i=1,200
       iline(i) = ICHAR(yline(i:i))
    ENDDO
    CALL MPI_Gather(iline, 200, MPI_INTEGER, ilines, 200, MPI_INTEGER,    &
                    0, icomm, implcode)
  
  ! CALL MPI_Gather(yline, 200, MPI_CHARACTER, ylines, 200, MPI_CHARACTER,&
  !                 0, icomm, implcode)
    IF ( implcode /= 0 ) THEN
      ierror  = implcode
      yerrmsg = 'Error in MPI_GATHER'
      RETURN
    ENDIF
  ELSE
    ylines(1) = yline
  ENDIF
  
  ! processor 0 prints all ylines
  IF( myid == 0) THEN
    IF(npe > 1) THEN
      DO n=1,npe
        DO i=1,200
          ylines(n)(i:i) = CHAR(ilines(i,n))
        ENDDO
      ENDDO
    ENDIF
    DO i=1,npe
      IF(ylines(i) /= ' ') THEN
        WRITE(nout,'(A)') ylines(i)(1:LEN_TRIM(ylines(i)))
      ENDIF
    ENDDO
  ENDIF
  
END SUBROUTINE check_record

!==============================================================================
!+ Check EC grid description section against the Namelist parameters
!------------------------------------------------------------------------------

SUBROUTINE check_ec_grid (igds, idimgds, ipds, idimpds,                      &
                 ydate, iedim, jedim, kedim, startlat, startlon, dlon, dlat, &
                 pollon, pollat, lwork, npe, icomm, myid, itype_calendar,    &
                 yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!  This routine checks the EC-grid against the Namelist parameters.
!  It is a collective routine and has to be called by all compute PEs.
!  For LM the following 10 values are checked:
!    ie_tot, je_tot, ke_tot, startlat, startlon, pollat, pollon, dlon, dlat
!    and the actual date of the forecast.
!
! Method:
!  If the values in the grid description section and the Namelist parameters
!  are not equal, an error status variable is set (/= 0). The error status
!  variables from all processors are gathered and scattered to all. If an
!  error occurs, the processor with lowest rank prints its (wrong) values.
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER (KIND=intgribf), INTENT(IN)   ::  &
  idimgds,          & ! dimension of the grid description section
  idimpds             ! dimension of the product definition section

INTEGER (KIND=intgribf), INTENT(IN)   ::  &
  igds(idimgds),    & ! grid description section
  ipds(idimpds)       ! product definition section

CHARACTER (LEN=*),        INTENT(IN)  ::  &
  ydate               ! actual date from Namelist parameter

INTEGER (KIND=iintegers), INTENT(IN)  ::  &
  iedim, jedim, kedim ! dimensions of the fields

REAL    (KIND=ireals),    INTENT(IN)  ::  &
  startlat,startlon,& ! coordinates of the lower left grid point
  dlat,  dlon,      & ! grid resolution in lambda and phi direction
  pollat, pollon      ! coordinates of the rotated north pole

LOGICAL,                  INTENT(IN) ::   &
  lwork               ! is .TRUE. if there is work to do

INTEGER (KIND=iintegers), INTENT(IN)  ::  &
  itype_calendar,   & ! for specifying the calendar used
  npe,              & ! number of PEs in communicator icomm
  icomm,            & ! communicator for compute PEs
  myid                ! ID of this PE in communicator icomm

! arguments with intent(out):
CHARACTER (LEN=*),        INTENT(OUT) ::  &
  yerrmsg             ! error message

INTEGER (KIND=iintegers), INTENT(OUT) ::  &
  ierror              ! error status variable

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER (KIND=iintegers)  ::  &
  kegrib, i, izmplcode, irlas

REAL    (KIND=ireals)     ::  &
  zstartlat_g, zstartlon_g,   & ! coordinates of the lower left grid point
  zpollat_grib, zpollon_grib, & ! coordinates of the rotated north pole
  zdlon_grib, zdlat_grib,     & ! grid resolution in latitude and longitude
  zsouthp_lat, zsouthp_lon      ! south pole coordinates from NAMELIST input

CHARACTER  (LEN=10)      ::  &
  ydate_ref     ! reference date from the grib

! Local arrays:
INTEGER                   ::  & ! really the standard integer
  isenderror (0:npe-1),       & !
  irecverror (0:npe-1)          !

INTEGER (KIND=intgribf)   ::                                               &
  iyear_ref, imonth_ref, iday_ref, ihour_ref, imin_ref, isec_ref,          &
  iyear, imonth, iday, ihour, imin, isec, iref_act, ierrf, ifactor,        &
  itimepassed

! UB>> hidden local switch to allow to deactivate time checking of grib records, 
!      see below ...
LOGICAL, PARAMETER :: lcheckdate = .true.

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ! initialize error status variable
  yerrmsg       = '   '
  ierror        = 0
  isenderror(:) = 0
  irecverror(:) = 0
  ierrf         = 0_intgribf

!------------------------------------------------------------------------------
! Section 2: Compare the values
!------------------------------------------------------------------------------

  IF (lwork) THEN
    ! check the field dimensions
    IF (iedim /= igds(5)) ierror = 1
    IF (jedim /= igds(6)) ierror = 2
    IF (igds(2) /= 0) THEN
      kegrib = (igds(2) - 1) / 2
      IF (kedim /= kegrib ) ierror = 3
    ENDIF

    IF ( (igds(4) == 10) .OR. (igds(4) == 14) ) THEN
      zpollat_grib = REAL(igds(20), ireals)*0.001_ireals
      zpollon_grib = REAL(igds(21), ireals)*0.001_ireals
      zsouthp_lat  = - pollat
      zsouthp_lon  = pollon + 180.0_ireals
      IF (zsouthp_lon > 180.0_ireals) THEN
        zsouthp_lon = zsouthp_lon - 360.0_ireals
      ENDIF
    ENDIF

    ! Check, whether the resolutions are given or not
    IF (IBITS(igds(9),7,1) /= 0) THEN
      zdlat_grib   = REAL(igds(13), ireals)*0.001_ireals
      zdlon_grib   = REAL(igds(12), ireals)*0.001_ireals
    ELSE
      zdlat_grib   = REAL(igds(10)-igds( 7), ireals)*0.001_ireals / (jedim-1)
      IF (igds(11)-igds( 8) < 0.0_ireals) THEN
        ! If the area is located around the 180-Meridian, igds(11)-igds(8)
        ! will be negative and 360 degrees have to be added
        zdlon_grib   = (REAL(igds(11)-igds( 8), ireals)*0.001_ireals +    &
                       360.0_ireals) / (iedim-1)
      ELSE
        zdlon_grib   = REAL(igds(11)-igds( 8), ireals)*0.001_ireals / (iedim-1)
      ENDIF
    ENDIF

    ! Southern latitude is first or last according to scanning mode
    IF (igds(14) == 0) THEN
      irlas=igds(10) ! ECMWF scanning mode
    ELSE
      irlas=igds( 7) ! almost all the others'
    ENDIF

    ! determine the values from the grid description section
    IF ( (ipds(2) == 2) .AND. (ipds(8) == 110) ) THEN
      IF (ipds(7) == 33) THEN
        ! u - gridpoint
        zstartlat_g = REAL(irlas,    ireals)*0.001_ireals
        zstartlon_g = REAL(igds( 8), ireals)*0.001_ireals                &
                                 - zdlon_grib*0.5_ireals
      ELSEIF (ipds(7) == 34) THEN
        ! v - gridpoint
        zstartlat_g = REAL(irlas,    ireals)*0.001_ireals                &
                                 - zdlat_grib*0.5_ireals
        zstartlon_g = REAL(igds( 8), ireals)*0.001_ireals
      ELSE
        ! all other gridpoints with tabtyp=2 and levtyp=110
        zstartlat_g = REAL(irlas,    ireals)*0.001_ireals
        zstartlon_g = REAL(igds( 8), ireals)*0.001_ireals
      ENDIF
    ELSE
       ! all other levtyps and tabtyps
       zstartlat_g = REAL(irlas,    ireals)*0.001_ireals
       zstartlon_g = REAL(igds( 8), ireals)*0.001_ireals
    ENDIF

    ! compare the grib values with the Namelist parameters
    IF (ABS(zstartlat_g  - startlat) > 1.0E-3_ireals) ierror = 4
    IF (ABS(zstartlon_g  - startlon) > 1.0E-3_ireals) ierror = 5

    ! compare pole of rotation only for rotated grids
    IF ( (igds(4) == 10) .OR. (igds(4) == 14) ) THEN
      IF (ABS(zpollat_grib - zsouthp_lat) > 1.0E-5_ireals) ierror = 6
      IF (ABS(zpollon_grib - zsouthp_lon) > 1.0E-5_ireals) ierror = 7
    ENDIF

    IF (ABS(zdlat_grib - dlat) > 1.0E-3_ireals) ierror = 8
    IF (ABS(zdlon_grib - dlon) > 1.0E-3_ireals) ierror = 9

    ! Determine the reference date from grib values
    IF (ipds(11) == 100) THEN
      WRITE(ydate_ref,'(5I2.2)')                          &
                       ipds(22)  ,       0 , ipds(12), ipds(13), ipds(14)
      iyear_ref = ipds(22) * 100
    ELSE
      WRITE(ydate_ref,'(5I2.2)')                          &
                       ipds(22)-1, ipds(11), ipds(12), ipds(13), ipds(14)
      iyear_ref = (ipds(22)-1) * 100 + ipds(11)
    ENDIF
    imonth_ref = ipds(12)
    iday_ref   = ipds(13)
    ihour_ref  = ipds(14)
    imin_ref   = ipds(15)

  ! UB>> Introduced local hidden switch to turn off date checking, if a user has
  !      to use grib files with an out-of-bounds time range indicator
  !      (e.g., if 5-min input is desired out to more than 254 minutes
  !      forecast range and the unit of time (ipds(16)) is minutes).
  IF (lcheckdate) THEN

    ! Determine the corresponding integer values for the actual date "ydate"
    ! no minutes are given here
    READ(ydate    ,'(I4,5I2)') iyear, imonth, iday, ihour, imin, isec

    ! Determine the difference between reference date and actual date
    ! in minutes (with routine from dwdlib)
    IF     (itype_calendar == 0) THEN
      CALL difmin     (iyear_ref, imonth_ref, iday_ref, ihour_ref, imin_ref, &
                       iyear,     imonth,     iday,     ihour,     imin,     &
                       iref_act,  ierrf)
    ELSEIF (itype_calendar == 1) THEN
      CALL difmin_360 (iyear_ref, imonth_ref, iday_ref, ihour_ref, imin_ref, &
                       iyear,     imonth,     iday,     ihour,     imin,     &
                       iref_act,  ierrf)
    ENDIF

    ! Determine the factor for calculating the time passed (since the
    ! reference date) in minutes. For that purpose the unit of time range
    ! (ipds(16)) is needed.
    SELECT CASE (ipds(16))
    CASE (0)
      ! unit is minute
      ifactor = 1
    CASE (1)
      ! unit is hour
      ifactor = 60
    CASE (2)
      ! unit is day
      ifactor = 1440
    CASE (10)
      ! unit is 3 hours
      ifactor = 180
    CASE (11)
      ! unit is 6 hours
      ifactor = 360
    CASE (12)
      ! unit is 12 hours
      ifactor = 720
    CASE (13)
      ! unit is 15 minutes (DWD)
      ifactor = 15
    CASE (14)
      ! unit is 30 minutes (DWD)
      ifactor = 30
    CASE DEFAULT
      ! this is not implemented or even not defined
      ierror  = 11
    END SELECT

    IF (ierror == 0) THEN
      ! Now determine the time passed (since the reference date) in minutes.
      ! For that the time range indicator (ipds(19)) is needed to select the
      ! correct grib octet (ipds(17) or ipds(18)).
      SELECT CASE (ipds(19))
      CASE (0)
        ! product valid for time in ipds(17)
        itimepassed = ifactor * ipds(17)
      CASE (1)
        ! initialized analyses for time = 0
        ! for IFS constant data, this should not be checked!
        ! Therefore set itimepassed to iref_act
        itimepassed = iref_act
      CASE (10)
        itimepassed = ifactor * ipds(18)
      CASE (2,3,4,5)
        ! this is a product which is valid for a time period and cannot be
        ! used as input field
        ierror = 12
      CASE DEFAULT
        ! these time range indicators are not used in the LM-Package
        ierror = 13
      END SELECT
    ENDIF

    IF (ierror == 0) THEN
      ! The timepassed and the difference of the dates are sometimes equal,
      ! but with hincbound=0.5/0.25, there can be differences of a timestep
      ! So allow for a difference of 2 minutes
      IF ( ABS(iref_act - itimepassed) > 2_iintegers) THEN
        ierror = 10
      ENDIF
    ENDIF

  ENDIF  ! if (lcheckdate)

  ENDIF  ! if (lwork)

!------------------------------------------------------------------------------
! Section 3: Exchange error status variable
!------------------------------------------------------------------------------

  isenderror (myid) = ierror

  IF (npe > 1) THEN
    CALL MPI_ALLGATHER (isenderror(myid), 1, MPI_INTEGER,              &
                        irecverror,       1, MPI_INTEGER, icomm, izmplcode)
  ELSE
    irecverror(:) = isenderror(:)
  ENDIF

!------------------------------------------------------------------------------
! Section 4: If necessary, print errors
!------------------------------------------------------------------------------

  DO i = 0, npe-1
    IF (irecverror(i)/= 0) THEN
      ierror = irecverror(i)

      ! only one processor prints the error message
      IF (myid == i) THEN
        IF (ierror < 10) THEN
          PRINT *, '         data file          namelist input'

          PRINT *, 'ie_tot                ',igds(5)     ,'       ',iedim
          PRINT *, 'je_tot                ',igds(6)     ,'       ',jedim
          PRINT *, 'ke_tot                ',kegrib      ,'       ',kedim

          PRINT *, 'startlat_tot          ',zstartlat_g ,'       ',startlat
          PRINT *, 'startlon_tot          ',zstartlon_g ,'       ',startlon

          PRINT *, 'rot. south pole (lat) ',zpollat_grib,'       ',zsouthp_lat
          PRINT *, 'rot. south pole (lon) ',zpollon_grib,'       ',zsouthp_lon

          PRINT *, 'dlat                  ',zdlat_grib  ,'       ',dlat
          PRINT *, 'dlon                  ',zdlon_grib  ,'       ',dlon
        ELSEIF (ierror == 10) THEN
          PRINT *, ' The actual date           ', ydate
          PRINT *, ' and the reference date    ', ydate_ref, ' + ', &
                     itimepassed, 'minutes    do not match!!'
        ELSEIF (ierror == 11) THEN
          PRINT *, ' Wrong (or not implemented) unit of time (uot)', ipds(16)
          PRINT *, '    Implemented values are uot =  0 (minute)     '
          PRINT *, '                           uot =  1 (hour)       '
          PRINT *, '                           uot =  2 (day)        '
          PRINT *, '                           uot = 10 ( 3 hours)   '
          PRINT *, '                           uot = 11 ( 6 hours)   '
          PRINT *, '                           uot = 12 (12 hours)   '
          PRINT *, '                           uot = 13 (15 minutes) '
          PRINT *, '                           uot = 14 (30 minutes) '
        ELSEIF (ierror == 12) THEN
          PRINT *, ' Wrong time range indicator (tri)', ipds(19)
          PRINT *, '    The values tri = 2,3,4,5 mean that the product is '
          PRINT *, '    valid for a  time period and cannot  be used as a '
          PRINT *, '    input field'
        ELSEIF (ierror == 13) THEN
          PRINT *, ' This time range indicator (tri)', ipds(19)
          PRINT *, ' is not used in the LM Package'
        ELSE
          PRINT *, ' ??? This is an unknown error ??? '
        ENDIF
      ENDIF
      EXIT
    ENDIF
  ENDDO

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE check_ec_grid

!==============================================================================
!+ Check LM grid description section against the Namelist parameters
!------------------------------------------------------------------------------

SUBROUTINE check_lm_grid (igds, idimgds, ipds, idimpds,                      &
                 ydate, iedim, jedim, kedim, startlat, startlon, dlon, dlat, &
                 pollon, pollat, lwork, npe, icomm, myid, itype_calendar,    &
                 yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!  This routine checks the LM-grid against the Namelist parameters. 
!  It is a collective routine and has to be called by all compute PEs. 
!  For LM the following 10 values are checked:
!    ie_tot, je_tot, ke_tot, startlat, startlon, pollat, pollon, dlon, dlat
!    and the actual date of the forecast.
!
! Method:
!  If the values in the grid description section and the Namelist parameters
!  are not equal, an error status variable is set (/= 0). The error status
!  variables from all processors are gathered and scattered to all. If an 
!  error occurs, the processor with lowest rank prints its (wrong) values.
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER (KIND=intgribf), INTENT(IN)   ::  &
  idimgds,          & ! dimension of the grid description section
  idimpds             ! dimension of the product definition section

INTEGER (KIND=intgribf), INTENT(IN)   ::  &
  igds(idimgds),    & ! grid description section
  ipds(idimpds)       ! product definition section

CHARACTER (LEN=*),        INTENT(IN)  ::  &
  ydate               ! actual date from Namelist parameter

INTEGER (KIND=iintegers), INTENT(IN)  ::  &
  iedim, jedim, kedim ! dimensions of the fields
 
REAL    (KIND=ireals),    INTENT(IN)  ::  &
  startlat, startlon, & ! coordinates of the lower left grid point
  dlon,  dlat,        & ! grid resolution in lambda and phi direction
  pollat, pollon        ! coordinates of the rotated north pole

LOGICAL,                  INTENT(IN) ::   &
  lwork               ! is .TRUE. if there is work to do

INTEGER (KIND=iintegers), INTENT(IN)  ::  &
  itype_calendar,   & ! for specifying the calendar used
  npe,              & ! number of PEs in communicator icomm
  icomm,            & ! communicator for compute PEs
  myid                ! ID of this PE in communicator icomm

LOGICAL                   ::  &
  lzclimate     ! to indicate whether T_CL or W_CL are processed

! arguments with intent(out):
CHARACTER (LEN=*),        INTENT(OUT) ::  &
  yerrmsg             ! error message

INTEGER (KIND=iintegers), INTENT(OUT) ::  &
  ierror              ! error status variable

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER (KIND=iintegers)  ::  &
  ivctype,                    & ! type of the vertical coordinate used
  kegrib, i, izmplcode

REAL    (KIND=ireals)     ::  &
  zstartlat_g, zstartlon_g,   & ! coordinates of the lower left grid point
  zpollat_grib, zpollon_grib, & ! coordinates of the rotated north pole
  zdlon_grib, zdlat_grib,     & ! grid resolution in latitude and longitude
  zsouthp_lat, zsouthp_lon      ! south pole coordinates from NAMELIST input

REAL    (KIND=irealgrib)  :: refstf

CHARACTER  (LEN=10)       ::  &
  ydate_ref     ! reference date from the grib

! Local arrays: 
INTEGER                   ::  & ! really the standard integer
  isenderror (0:npe-1),       & !
  irecverror (0:npe-1)          !

INTEGER (KIND=intgribf)   ::                                               &
  iyear_ref, imonth_ref, iday_ref, ihour_ref, imin_ref, isec_ref,          &
  iyear, imonth, iday, ihour, imin, isec, iref_act, ierrf, ifactor,        &
  itimepassed

! UB>> hidden local switch to allow to deactivate time checking of grib records, 
!      see below ...
LOGICAL, PARAMETER :: lcheckdate = .false.

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ! initialize error status variable
  yerrmsg       = '   '
  ierror        = 0
  isenderror(:) = 0
  irecverror(:) = 0
  ierrf         = 0_intgribf

!------------------------------------------------------------------------------
! Section 2: Compare the values
!------------------------------------------------------------------------------

  IF (lwork) THEN
    ! check the field dimensions
    IF (iedim /= igds(5)) ierror = 1
    IF (jedim /= igds(6)) ierror = 2

#ifdef GRIBDWD
    IF ((ipds(8) == 109) .OR. (ipds(8) == 110)) THEN
      ! Get number of vertical levels
      ! Check for the type of the vertical coordinate parameters
      ivctype = NINT (REFSTF(igds(26)), iintegers)
      IF ((ivctype >= 1) .AND. (ivctype <= 200)) THEN
         ! This is the new grib GDS coding style introduced for SLEVE
         ! kegrib is new on place 27 of the GDS
         kegrib = NINT (REFSTF(igds(27)),iintegers)
      ELSE
         ! This is the old grib GDS coding style, kept for backwards
         ! compatibility.
         ivctype = INT(igds(29 + kedim+1 + 1), iintegers)
         IF (ivctype == -999999) THEN
            ! This is the old type of coding (ivctype = 1 or 2)
            kegrib = igds(2) - 4 - 1      ! ivctype=1,2 has 4 additional variables
         ELSEIF ( (ivctype == 1) .OR. (ivctype == 2) ) THEN
            kegrib = igds(2) - 5 - 1      ! ivctype=1,2 has 5 additional variables
         ELSEIF (ivctype == 3) THEN
            kegrib = igds(2) - 8 - 1      ! ivctype=3 has 8 additional variables
         ELSEIF ( (ivctype >= 101) .AND. (ivctype <= 200) ) THEN
            kegrib = igds(2) - 10 - 1     ! new reference atmosphere has 10 additional variables
         ENDIF
      ENDIF
      IF (kedim /= kegrib ) ierror = 3
    ELSE
       kegrib = -1 ! means: not available
    ENDIF
#endif

    ! Instead of recomputing the given grib coordinates of the rotated south 
    ! pole to the north pole, the namelist coordinates of the north pole are
    ! computed to the south pole in the same way as it is done in GME2LM.
    zpollat_grib =   REAL(igds(20), ireals)*0.001_ireals
    zpollon_grib =   REAL(igds(21), ireals)*0.001_ireals
    zsouthp_lat  =   - pollat
    zsouthp_lon  =   pollon + 180.0_ireals
    IF (zsouthp_lon > 180.0_ireals) THEN
      zsouthp_lon = zsouthp_lon - 360.0_ireals
    ENDIF

    ! Check, whether the resolutions are given or not
    IF (IBITS(igds(9),7,1) /= 0) THEN
      zdlat_grib   = REAL(igds(13), ireals)*0.001_ireals
      zdlon_grib   = REAL(igds(12), ireals)*0.001_ireals
    ELSE
      zdlat_grib   = REAL(igds(10)-igds( 7), ireals)*0.001_ireals / (jedim-1)
      IF (igds(11)-igds( 8) < 0.0_ireals) THEN
        ! If the area is located around the 180-Meridian, igds(11)-igds(8)
        ! will be negative and 360 degrees have to be added
        zdlon_grib   = (REAL(igds(11)-igds( 8), ireals)*0.001_ireals +    &
                       360.0_ireals) / (iedim-1)
      ELSE
        zdlon_grib   = REAL(igds(11)-igds( 8), ireals)*0.001_ireals / (iedim-1)
      ENDIF
    ENDIF

    ! determine the values from the grid description section
    IF     ( ((ipds(2) == 2).AND.(ipds(8) == 110).AND.(ipds(7) ==  33)) .OR. &
             ((ipds(2) == 2).AND.(ipds(8) ==   1).AND.(ipds(7) == 124)) )THEN
        ! u - gridpoint
        zstartlat_g = REAL(igds( 7), ireals)*0.001_ireals
        zstartlon_g = REAL(igds( 8), ireals)*0.001_ireals                    &
                                 - zdlon_grib*0.5_ireals
    ELSEIF ( ((ipds(2) == 2).AND.(ipds(8) == 110).AND.(ipds(7) ==  34)) .OR. &
             ((ipds(2) == 2).AND.(ipds(8) ==   1).AND.(ipds(7) == 125)) )THEN
        ! v - gridpoint
        zstartlat_g = REAL(igds( 7), ireals)*0.001_ireals                    &
                                 - zdlat_grib*0.5_ireals
        zstartlon_g = REAL(igds( 8), ireals)*0.001_ireals
    ELSE
       ! all other variables
       zstartlat_g = REAL(igds( 7), ireals)*0.001_ireals
       zstartlon_g = REAL(igds( 8), ireals)*0.001_ireals
    ENDIF
  
    ! compare the grib values with the Namelist parameters
    IF (ABS(zstartlat_g  - startlat) > 1.0E-3_ireals) ierror = 4
    IF (ABS(zstartlon_g  - startlon) > 1.0E-3_ireals) ierror = 5
    
    IF (ABS(zpollat_grib - zsouthp_lat) > 1.0E-5_ireals) ierror = 6
    IF (ABS(zpollon_grib - zsouthp_lon) > 1.0E-5_ireals) ierror = 7
  
    IF (ABS(zdlat_grib - dlat) > 1.0E-3_ireals) ierror = 8
    IF (ABS(zdlon_grib - dlon) > 1.0E-3_ireals) ierror = 9
  
  ! UB>> Introduced local hidden switch to turn off date checking, if a user has
  !      to use grib files with an out-of-bounds time range indicator
  !      (e.g., if 5-min input is desired out to more than 254 minutes
  !      forecast range and the unit of time (ipds(16)) is minutes).
  IF (lcheckdate) THEN

    ! check the date
    ! Before, only the reference date was checked here.
    ! But what really is needed is to check the time for which the product
    ! is valid. This has to be calculated using ipds (16-19): unit of time
    ! range, time stamps and time range indicator

    ! Only for the fields T_CL and W_CL (climatological soil values) it is
    ! possible that older fields are taken (if LM initial data are computed
    ! from an older GME-file, e.g. from gfff00120000 from the run 12 hours
    ! before). T_CL and W_CL would then be taken from the corresponding
    ! giff00000000-file, and thus would be 12 hours older).

    ! Check whether T_CL or W_CL are processed
    IF ( ((ipds(7) == 85) .AND. (ipds(2) == 2) .AND. (ipds(8) == 111)         &
                                               .AND. (ipds(10) > 10 ))   .OR. &
         ((ipds(7) == 86) .AND. (ipds(2) == 2) .AND. (ipds(8) == 112)         &
                                               .AND. (ipds(10) > 100)) ) THEN
      lzclimate = .TRUE.
    ELSE
      lzclimate = .FALSE.
    ENDIF

    ! Determine the reference date from grib values
    IF (ipds(11) == 100) THEN
      WRITE(ydate_ref,'(5I2.2)')                          &
                       ipds(22)  ,       0 , ipds(12), ipds(13), ipds(14)
      iyear_ref = ipds(22) * 100
    ELSE
      WRITE(ydate_ref,'(5I2.2)')                          &
                       ipds(22)-1, ipds(11), ipds(12), ipds(13), ipds(14)
      iyear_ref = (ipds(22)-1) * 100 + ipds(11)
    ENDIF
    imonth_ref = ipds(12)
    iday_ref   = ipds(13)
    ihour_ref  = ipds(14)
    imin_ref   = ipds(15)

    ! Determine the corresponding integer values for the actual date "ydate"
    ! no minutes are given here
    READ(ydate    ,'(I4,5I2)') iyear, imonth, iday, ihour, imin, isec

    ! Determine the difference between reference date and actual date
    ! in minutes (with routine from dwdlib)
    IF     (itype_calendar == 0) THEN
      CALL difmin     (iyear_ref, imonth_ref, iday_ref, ihour_ref, imin_ref, &
                       iyear,     imonth,     iday,     ihour,     imin,     &
                       iref_act,  ierrf)
    ELSEIF (itype_calendar == 1) THEN
      CALL difmin_360 (iyear_ref, imonth_ref, iday_ref, ihour_ref, imin_ref, &
                       iyear,     imonth,     iday,     ihour,     imin,     &
                       iref_act,  ierrf)
    ENDIF

    ! Determine the factor for calculating the time passed (since the
    ! reference date) in minutes. For that purpose the unit of time range
    ! (ipds(16)) is needed.
    SELECT CASE (ipds(16))
    CASE (0)
      ! unit is minute
      ifactor = 1
    CASE (1)
      ! unit is hour
      ifactor = 60
    CASE (2)
      ! unit is day
      ifactor = 1440
    CASE (10)
      ! unit is 3 hours
      ifactor = 180
    CASE (11)
      ! unit is 6 hours
      ifactor = 360
    CASE (12)
      ! unit is 12 hours
      ifactor = 720
    CASE (13)
      ! unit is 15 minutes (DWD)
      ifactor = 15
    CASE (14)
      ! unit is 30 minutes (DWD)
      ifactor = 30
    CASE DEFAULT
      ! this is not implemented or even not defined
      ierror  = 11
    END SELECT

    IF (ierror == 0) THEN
      ! Now determine the time passed (since the reference date) in minutes.
      ! For that the time range indicator (ipds(19)) is needed to select the
      ! correct grib octet (ipds(17) or ipds(18)).
      SELECT CASE (ipds(19))
      CASE (0,1,13)
        itimepassed = ifactor * ipds(17)
      CASE (10)
        itimepassed = ifactor * ipds(18)
      CASE (2,3,4,5)
        ! this is a product which is valid for a time period and cannot be
        ! used as input field, but for restart fields
        !US:  ierror = 12
        itimepassed = ifactor * ipds(18)
      CASE DEFAULT
        ! these time range indicators are not used in the LM-Package
        ierror = 13
      END SELECT
    ENDIF

    IF (ierror == 0) THEN
      ! The timepassed and the difference of the dates are sometimes equal,
      ! but with hincbound=0.5/0.25, there can be differences of a timestep
      ! So allow for a difference of 2 minutes
      IF ( ABS(iref_act - itimepassed) > 2_iintegers) THEN
        IF (.NOT. lzclimate) THEN
          ierror = 10
        ELSE
          PRINT *, 'Using older climate values (T_CL, W_CL) from ', ydate_ref
        ENDIF
      ENDIF
    ENDIF

  ENDIF  ! if (lcheckdate)

  ENDIF  ! if (lwork)

!------------------------------------------------------------------------------
! Section 3: Exchange error status variable
!------------------------------------------------------------------------------

  isenderror (myid) = ierror  

  IF (npe > 1) THEN
    CALL MPI_ALLGATHER (isenderror(myid), 1, MPI_INTEGER,              &
                        irecverror,       1, MPI_INTEGER, icomm, izmplcode)
  ELSE 
    irecverror(:) = isenderror(:)
  ENDIF

!------------------------------------------------------------------------------
! Section 4: If necessary, print errors
!------------------------------------------------------------------------------

  DO i = 0, npe-1
    IF (irecverror(i)/= 0) THEN
      ierror = irecverror(i)

      ! only one processor prints the error message
      IF (myid == i) THEN
        IF (ierror < 10) THEN
          PRINT *, '         data file          namelist input'

          PRINT *, 'ie_tot                ',igds(5)     ,'       ',iedim
          PRINT *, 'je_tot                ',igds(6)     ,'       ',jedim
          IF (kegrib /= -1) THEN
          PRINT *, 'ke_tot                ',kegrib      ,'       ',kedim
          ENDIF

          PRINT *, 'startlat_tot          ',zstartlat_g ,'       ',startlat
          PRINT *, 'startlon_tot          ',zstartlon_g ,'       ',startlon

          PRINT *, 'rot. south pole (lat) ',zpollat_grib,'       ',zsouthp_lat
          PRINT *, 'rot. south pole (lon) ',zpollon_grib,'       ',zsouthp_lon

          PRINT *, 'dlat                  ',zdlat_grib  ,'       ',dlat
          PRINT *, 'dlon                  ',zdlon_grib  ,'       ',dlon
        ELSEIF (ierror == 10) THEN
          PRINT *, ' The actual date           ', ydate
          PRINT *, ' and the reference date    ', ydate_ref, ' + ', &
                     itimepassed, 'minutes    do not match!!'
        ELSEIF (ierror == 11) THEN
          PRINT *, ' Wrong (or not implemented) unit of time (uot)', ipds(16)
          PRINT *, '    Implemented values are uot =  0 (minute)     '
          PRINT *, '                           uot =  1 (hour)       '
          PRINT *, '                           uot =  2 (day)        '
          PRINT *, '                           uot = 10 ( 3 hours)   '
          PRINT *, '                           uot = 11 ( 6 hours)   '
          PRINT *, '                           uot = 12 (12 hours)   '
          PRINT *, '                           uot = 13 (15 minutes) '
          PRINT *, '                           uot = 14 (30 minutes) '
        ELSEIF (ierror == 12) THEN
          PRINT *, ' Wrong time range indicator (tri)', ipds(19)
          PRINT *, '    The values tri = 2,3,4,5 mean that the product is '
          PRINT *, '    valid for a  time period and cannot  be used as a '
          PRINT *, '    input field'
        ELSEIF (ierror == 13) THEN
          PRINT *, ' This time range indicator (tri)', ipds(19)
          PRINT *, ' is not used in the LM Package'
        ELSE
          PRINT *, ' ??? This is an unknown error ??? '
        ENDIF
      ENDIF
      EXIT
    ENDIF
  ENDDO

!------------------------------------------------------------------------------
!- End of the Subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE check_lm_grid

!==============================================================================
!+ Check GME grid description section against the Namelist parameters
!------------------------------------------------------------------------------

SUBROUTINE check_gme_grid (igds, idimgds, ipds, idimpds, ydate, kegme, ni,   &
                  ni2, ni3, ndiam, lwork, npe, icomm, myid, itype_calendar,  &
                  yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!  This routine checks the GME-grid against the Namelist parameters. 
!  It is a collective routine and has to be called by all compute PEs. 
!  For GME the following 5 values are checked:
!    i3e_gme, ni_gme, ni2, ni3, nd
!
! Method:
!  If the values in the grid description section and the Namelist parameters
!  are not equal, an error status variable is set (/= 0). The error status
!  variables from all processors are gathered and scattered to all. If an 
!  error occurs, the processor with lowest rank prints its (wrong) values.
!
!------------------------------------------------------------------------------

! Subroutine arguments
INTEGER (KIND=intgribf), INTENT(IN)   ::  &
  idimgds,          & ! dimension of the grid description section
  idimpds             ! dimension of the product definition section

INTEGER (KIND=intgribf), INTENT(IN)   ::  &
  igds(idimgds),    & ! grid description section
  ipds(idimpds)       ! product definition section

CHARACTER (LEN=*),        INTENT(IN)  ::  &
  ydate               ! actual date from Namelist parameter

INTEGER (KIND=iintegers), INTENT(IN)  ::  &
  kegme, ni, ni2, ni3, ndiam   ! values for GME 
 
LOGICAL,                  INTENT(IN) ::   &
  lwork               ! is .TRUE. if there is work to do

INTEGER (KIND=iintegers), INTENT(IN)  ::  &
  itype_calendar,   & ! for specifying the calendar used
  npe,              & ! number of PEs in communicator icomm
  icomm,            & ! communicator for compute PEs
  myid                ! ID of this PE in communicator icomm
 
! arguments with intent(out):
CHARACTER (LEN=*),        INTENT(OUT) ::  &
  yerrmsg             ! error message

INTEGER (KIND=iintegers), INTENT(OUT) ::  &
  ierror              ! error status variable

! UB>> hidden local switch to allow to deactivate time checking of grib records, 
!      see below ...
LOGICAL, PARAMETER :: lcheckdate = .true.

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER (KIND=iintegers)  ::  &
  kegrib, i, izmplcode

LOGICAL                   ::  &
  lzclimate     ! to indicate whether T_CL or W_CL are processed

CHARACTER  (LEN=10)       ::  &
  ydate_ref     ! reference date from the grib

! Local arrays: 
INTEGER                   ::  & ! really the standart integer
  isenderror (0:npe-1),       & !
  irecverror (0:npe-1)          !

INTEGER (KIND=intgribf)   ::                                               &
  iyear_ref, imonth_ref, iday_ref, ihour_ref, imin_ref, isec_ref,          &
  iyear, imonth, iday, ihour, imin, isec, iref_act, ierrf, ifactor,        &
  itimepassed

!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ! initialize error status variable
  yerrmsg       = '   '
  ierror        = 0
  isenderror(:) = 0
  irecverror(:) = 0
  ierrf         = 0_intgribf

!------------------------------------------------------------------------------
! Section 2: Compare the values
!------------------------------------------------------------------------------

  IF (lwork) THEN
    ! check the field dimensions
    kegrib = igds(2)/2 - 1
    IF (kegme /= kegrib ) ierror = 1
    IF (ni    /= igds(8)) ierror = 2
    IF (ni2   /= igds(5)) ierror = 3
    IF (ni3   /= igds(6)) ierror = 4
    IF (ndiam /= igds(7)) ierror = 5
  
  ! UB>> Introduced local hidden switch to turn off date checking, if a user has
  !      to use grib files with an out-of-bounds time range indicator
  !      (e.g., if 5-min input is desired out to more than 254 minutes
  !      forecast range and the unit of time (ipds(16)) is minutes).
  IF (lcheckdate) THEN

    ! check the date
    ! Before, only the reference date was checked here.
    ! But what really is needed is to check the time for which the product
    ! is valid. This has to be calculated using ipds (16-19): unit of time
    ! range, time stamps and time range indicator
    ! Only for the fields T_CL and W_CL (climatological soil values) it is
    ! possible that older fields are taken (if LM initial data are computed
    ! from an older GME-file, e.g. from gfff00120000 from the run 12 hours
    ! before). T_CL and W_CL would then be taken from the corresponding
    ! giff00000000-file, and thus would be 12 hours older).

    ! Check whether T_CL or W_CL are processed
    IF ( ((ipds(7) == 85) .AND. (ipds(2) == 2) .AND. (ipds(8) == 111)         &
                                               .AND. (ipds(10) > 10 ))   .OR. &
         ((ipds(7) == 86) .AND. (ipds(2) == 2) .AND. (ipds(8) == 112)         &
                                               .AND. (ipds(10) > 100)) ) THEN
      lzclimate = .TRUE.
    ELSE
      lzclimate = .FALSE.
    ENDIF


    ! Determine the reference date from grib values
    IF (ipds(11) == 100) THEN
      WRITE(ydate_ref,'(5I2.2)')                          &
                       ipds(22)  ,       0 , ipds(12), ipds(13), ipds(14)
      iyear_ref = ipds(22) * 100
    ELSE
      WRITE(ydate_ref,'(5I2.2)')                          &
                       ipds(22)-1, ipds(11), ipds(12), ipds(13), ipds(14)
      iyear_ref = (ipds(22)-1) * 100 + ipds(11)
    ENDIF
    imonth_ref = ipds(12)
    iday_ref   = ipds(13)
    ihour_ref  = ipds(14)
    imin_ref   = ipds(15)


    ! Determine the corresponding integer values for the actual date "ydate"
    ! no minutes are given here
    READ(ydate    ,'(I4,5I2)') iyear, imonth, iday, ihour, imin, isec

    ! Determine the difference between reference date and actual date
    ! in minutes (with routine from dwdlib)
    IF     (itype_calendar == 0) THEN
      CALL difmin     (iyear_ref, imonth_ref, iday_ref, ihour_ref, imin_ref, &
                       iyear,     imonth,     iday,     ihour,     imin,     &
                       iref_act,  ierrf)
    ELSEIF (itype_calendar == 1) THEN
      CALL difmin_360 (iyear_ref, imonth_ref, iday_ref, ihour_ref, imin_ref, &
                       iyear,     imonth,     iday,     ihour,     imin,     &
                       iref_act,  ierrf)
    ENDIF

    ! Determine the factor for calculating the time passed (since the
    ! reference date) in minutes. For that purpose the unit of time range
    ! (ipds(16)) is needed.
    SELECT CASE (ipds(16))
    CASE (0)
      ! unit is minute
      ifactor = 1
    CASE (1)
      ! unit is hour
      ifactor = 60
    CASE (2)
      ! unit is day
      ifactor = 1440
    CASE (10)
      ! unit is 3 hours
      ifactor = 180
    CASE (11)
      ! unit is 6 hours
      ifactor = 360
    CASE (12)
      ! unit is 12 hours
      ifactor = 720
    CASE (13)
      ! unit is 15 minutes (DWD)
      ifactor = 15
    CASE (14)
      ! unit is 30 minutes (DWD)
      ifactor = 30
    CASE DEFAULT
      ! this is not implemented or even not defined
      ierror  = 11
    END SELECT

    IF (ierror == 0) THEN
      ! Now determine the time passed (since the reference date) in minutes.
      ! For that the time range indicator (ipds(19)) is needed to select the
      ! correct grib octet (ipds(17) or ipds(18)).
      SELECT CASE (ipds(19))
      CASE (0,1,13)
        itimepassed = ifactor * ipds(17)
      CASE (10)
        itimepassed = ifactor * ipds(18)
      CASE (2,3,4,5)
        ! this is a product which is valid for a time period and cannot be
        ! used as input field
        ierror = 12
      CASE DEFAULT
        ! these time range indicators are not used in the LM-Package
        ierror = 13
      END SELECT
    ENDIF

    IF (ierror == 0) THEN
      ! The timepassed and the difference of the dates are sometimes equal,
      ! but with hincbound=0.5/0.25, there can be differences of a timestep
      ! So allow for a difference of 2 minutes
      IF ( ABS(iref_act - itimepassed) > 2_iintegers) THEN
        IF (.NOT. lzclimate) THEN
          ierror = 10
        ELSE
          PRINT *, 'Using older climate values (T_CL, W_CL) from ', ydate_ref
        ENDIF
      ENDIF
    ENDIF

  ENDIF  ! if (lcheckdate)

  ENDIF  ! if (lwork)

!------------------------------------------------------------------------------
! Section 3: Exchange error status variable
!------------------------------------------------------------------------------

  isenderror (myid) = ierror  

  IF (npe > 1) THEN
    CALL MPI_ALLGATHER (isenderror(myid), 1, MPI_INTEGER,              &
                        irecverror,       1, MPI_INTEGER, icomm, izmplcode)
  ELSE 
    irecverror(:) = isenderror(:)
  ENDIF

!------------------------------------------------------------------------------
! Section 4: If necessary, print errors
!------------------------------------------------------------------------------

  DO i = 0, npe-1
    IF (irecverror(i)/= 0) THEN
      ierror  = irecverror(i)
      yerrmsg = 'Wrong grid specification or date for GME data'

      ! only one processor prints the error message
      IF (myid == i) THEN
        IF (ierror < 6) THEN
          PRINT *, '         data file          namelist input'

          PRINT *, 'ni_gme     ',igds(8),'       ',ni
          PRINT *, 'ni2        ',igds(5),'       ',ni2
          PRINT *, 'ni3        ',igds(6),'       ',ni3
          PRINT *, 'nd         ',igds(7),'       ',ndiam
          PRINT *, 'i3e_gme    ',kegrib ,'       ',kegme
        ELSEIF (ierror == 10) THEN
          PRINT *, ' The actual date           ', ydate
          PRINT *, ' and the reference date    ', ydate_ref, ' + ', &
                     itimepassed, 'minutes    do not match!!'
        ELSEIF (ierror == 11) THEN
          PRINT *, ' Wrong (or not implemented) unit of time (uot)', ipds(16)
          PRINT *, '    Implemented values are uot =  0 (minute)     '
          PRINT *, '                           uot =  1 (hour)       '
          PRINT *, '                           uot =  2 (day)        '
          PRINT *, '                           uot = 10 ( 3 hours)   '
          PRINT *, '                           uot = 11 ( 6 hours)   '
          PRINT *, '                           uot = 12 (12 hours)   '
          PRINT *, '                           uot = 13 (15 minutes) '
          PRINT *, '                           uot = 14 (30 minutes) '
        ELSEIF (ierror == 12) THEN
          PRINT *, ' Wrong time range indicator (tri)', ipds(19)
          PRINT *, '    The values tri = 2,3,4,5 mean that the product is '
          PRINT *, '    valid for a  time period and cannot  be used as a '
          PRINT *, '    input field'
        ELSEIF (ierror == 13) THEN
          PRINT *, ' This time range indicator (tri)', ipds(19)
          PRINT *, ' is not used in the LM Package'
        ELSE
          PRINT *, ' ??? This is an unknown error ??? '
        ENDIF
      ENDIF
      EXIT
    ENDIF
  ENDDO

!------------------------------------------------------------------------------
!- End of the Subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE check_gme_grid

!==============================================================================
!==============================================================================
!+ Defines a subroutine for computing differences in dates
!------------------------------------------------------------------------------

SUBROUTINE difmin_360 (iyear, imonth, iday, ihour, iminute,                 &
                       nyear, nmonth, nday, nhour, nminute, imdif, irm)

!------------------------------------------------------------------------------
!
! Description:
!   difmin_360 computes from the input data
!         dat1: iyear, imonth, iday, ihour, iminute,
!         dat3: nyear, nmonth, nday, nhour, nminute,
!   the difference in minutes: imdif = dat2-dat1
!
!------------------------------------------------------------------------------

INTEGER(KIND=iintegers), INTENT(IN)   ::                 &
  iyear, imonth, iday, ihour, iminute,                   &
  nyear, nmonth, nday, nhour, nminute

INTEGER(KIND=iintegers), INTENT(OUT)  ::                 &
  imdif, irm

! Local scalars
INTEGER(KIND=iintegers)               ::                 &
  itdif, ind1, ind2, isign, j1, kschal, izlocday, nzlocday

!------------------------------------------------------------------------------

! presettings
  irm   = 0
  imdif = 0

! some checks
  IF (ihour < 0 .OR. nhour < 0 .OR.               &
      ihour >= 24 .OR. nhour >= 24 .OR.           &
      iminute < 0 .OR. nminute < 0 .OR.           &
      iminute >= 60 .OR. nminute >= 60) THEN
    irm   = 1
    RETURN
  ENDIF

! difference in days
  izlocday = (imonth - 1) * 30 + iday
  nzlocday = (nmonth - 1) * 30 + nday
  itdif = nzlocday - izlocday

! for different years
  IF (iyear /= nyear) THEN

    ! difference in days for changing of the year
    IF (nyear > iyear) THEN
      ind1  = iyear
      ind2  = nyear - 1
      isign =  1
    ELSE
      ind1  = nyear
      ind2  = iyear - 1
      isign = -1
      itdif = - itdif
    ENDIF

    ! "loop" over several years
    itdif = itdif + (ind2 - ind1 + 1) * 360
    itdif = itdif * isign

  ENDIF

! and now the difference in minutes
  imdif = (itdif*24 + nhour - ihour) *60 + nminute - iminute

END SUBROUTINE difmin_360

!==============================================================================

END MODULE io_utilities
