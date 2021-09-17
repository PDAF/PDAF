!+ Data module for organizational variables for the parallel program
!------------------------------------------------------------------------------

MODULE data_parallel

!------------------------------------------------------------------------------
!
! Description:
!  This module contains all variables that are necessary for the organization
!  of the parallel program (number of processors, virtual topology, neighbors,
!  etc.). In a sequential environment some of these variables are also
!  necessary and are set accordingly (e.g. nproc = 1 and my_*_id = 0).
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
! 1.10       1998/09/29 Ulrich Schaettler
!  Changed igroup_diag and icomm_diag to allocatable arrays and
!  eliminated use of parameters for grid point and diagnostic output.
! 1.29       1999/05/11 Ulrich Schaettler
!  Added variables for later use of MPE_IO
! 1.39       2000/05/03 Ulrich Schaettler
!  Added declarations of buffers for distributing NAMELIST Input.
! 2.8        2001/07/06 Ulrich Schaettler
!  Added new variable lreorder for specifying whether the PEs may be reordered
!  while building the cartesian MPI-communicator
! 2.17       2002/05/08 Ulrich Schaettler
!  Added variable imp_grib for specifying the MPI type of variables with 
!  KIND parameter irealgrib
! 3.5        2003/09/02 Ulrich Schaettler
!  Added variables ncomm_type, ldatatypes, ltime_barrier for the optimization
!  of the boundary exchange (as Namelist parameteres)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  Eliminated my_peri_neigh
! V4_22        2012/01/31 Christoph Schraff
!  Target attribute added for 'isubpos'.
! V4_23        2012/05/10 Ulrich Schaettler
!  Editorial Changes
! V4_25        2012/09/28 Florian Prill, Hans-Juergen Panitz, Carlos Osuna
!  Introduced new namelist variable lprefetch_io in group /IOCTL/
!  Introduced new variable nexch_tag used for MPI boundary exchange instead of ntstep
!    (which caused troubles with MPI maximal tag in climate runs) (HJP)
!  Added variables that configure asynchronous NetCDF I/O behaviour. Add communicators
!   for I/O PE's in Netcdf I/O
! V4_28        2013/07/12 Ulrich Schaettler
!  Introduced new MPI data type for special grib_api integers kindOfSize
!   (which could be a 4 or a 8-byte integer)
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
!  Added parameters imp_single and imp_double (SP and DP REAL kind for MPI) (OF)
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
USE data_parameters , ONLY :   &
           wp,        & ! KIND-type parameters for real variables
           iintegers    ! kind-type parameter for "normal" integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Global (i.e. public) Declarations:

! 1. Information about the processor grid
! ---------------------------------------

  LOGICAL                           ::           &
    ldatatypes,      & ! if .TRUE.: use MPI-Datatypes for some communications
    ltime_barrier,   & ! if .TRUE.: use additional barriers for determining the
                       ! load-imbalance
    lasync_io,       & ! if .TRUE.: the model runs with extra PEs for
                       ! asynchronous IO
    lprefetch_io       ! if .TRUE.: boundary data are read in advance

  INTEGER   (KIND=iintegers)       ::           &
    nprocx,          & ! number of processors in x-direction
    nprocy,          & ! number of processors in y-direction
    nprocio,         & ! number of extra processors for doing asynchronous IO
    nc_asyn_io,      & ! number of asynchronous I/O PEs (netcdf)
    num_asynio_comm, & ! number of asynchronous I/O communicators (netcdf)
    num_iope_percomm,& ! number of asynchronous I/O PE per communicator (netcdf)
    nproc,           & ! total number of processors: nprocx * nprocy
    num_compute,     & ! number of compute PEs
    num_io,          & ! number of IO PEs

    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    ncomm_type,      & ! type of communication for boundary exchange

    my_world_id,     & ! rank of this subdomain in the global communicator
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    my_cart_pos(2),  & ! position of this subdomain in the cartesian grid 
                       ! in x- and y-direction
    my_cart_neigh(4)   ! neighbors of this subdomain in the cartesian grid

  INTEGER   (KIND=iintegers), ALLOCATABLE , TARGET   ::     &
    isubpos(:,:)       ! positions of the subdomains in the total domain. Given
                       ! are the i- and the j-indices of the lower left and the
                       ! upper right grid point in the order 
                       !                  i_ll, j_ll, i_ur, j_ur.
                       ! Only the interior of the domains are considered, not 
                       ! the boundary lines.
                       ! ('target' attribute is required because a pointer will
                       !  point to this in the assimilation (obs processing))

! 2. Further information for MPI
! ------------------------------

  INTEGER   (KIND=iintegers)       ::           &
    igroup_world,        & ! group that belongs to MPI_COMM_WORLD, i.e. all
                           ! processors
    icomm_world,         & ! communicator that belongs to igroup_world, i.e.
                           ! = MPI_COMM_WORLD
    icomm_compute,       & ! communicator for the group of compute PEs
    icomm_asynio,        & ! communicator for the group of I/O PEs (netcdf).
    intercomm_asynio,    & ! communicator for the group of I/O PEs + PE0 (netcdf)
    igroup_cart,         & ! group of the compute PEs
    icomm_cart,          & ! communicator for the virtual cartesian topology
    icomm_row,           & ! communicator for a east-west row of processors
    iexch_req(4),        & ! stores the sends requests for the neighbor-exchange
                           ! that can be used by MPI_WAIT to identify the send
    imp_reals,           & ! determines the correct REAL type used in the model
                           ! for MPI
    imp_single,          & ! single precision REAL type for MPI
    imp_double,          & ! double precision REAL type for MPI
    imp_grib,            & ! determines the REAL type used for the GRIB library
    imp_integers,        & ! determines the correct INTEGER type used in the
                           ! model for MPI
    imp_integ_ga,        & ! determines the correct INTEGER type used for grib_api
    imp_byte,            & ! determines the correct BYTE type used in the model
                           ! for MPI
    imp_character,       & ! determines the correct CHARACTER type used in the
                           ! model for MPI
    imp_logical,         & ! determines the correct LOGICAL   type used in the
                           ! model for MPI
    nexch_tag              ! tag to be used for MPI boundary exchange
                           !  (in calls to exchg_boundaries)

  LOGICAL                                       &
    lcompute_pe,         & ! indicates whether this is a compute PE or not
    lreorder               ! during the creation of the virtual topology the
                           ! ranking of the processors may be reordered

  INTEGER (KIND=iintegers), ALLOCATABLE  ::           &
    igroup_diag(:),      & ! groups for the diagnostic domains
    icomm_diag (:)         ! communicators for the diagnostic domains

! 3. Static Send buffers
! ----------------------

  REAL (KIND=wp),     ALLOCATABLE  ::           &
    sendbuf(:,:)       ! sending buffer for boundary exchange:
                       ! 1-4 are used for sending, 5-8 are used for receiving

  INTEGER (KIND=iintegers)         ::           &
    isendbuflen        ! length of one column of sendbuf 

  ! Buffers for distributing the Namelists
  INTEGER (KIND=iintegers), ALLOCATABLE ::   intbuf  (:)
  REAL    (KIND=wp)       , ALLOCATABLE ::   realbuf (:)
  LOGICAL                 , ALLOCATABLE ::   logbuf  (:)
  CHARACTER (LEN=100)     , ALLOCATABLE ::   charbuf (:)

!==============================================================================
END MODULE data_parallel
