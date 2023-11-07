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
           ireals,    & ! KIND-type parameters for real variables
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
    lasync_io          ! if .TRUE.: the model runs with extra PEs for 
                       ! asynchronous IO

  INTEGER   (KIND=iintegers)       ::           &
    nprocx,          & ! number of processors in x-direction
    nprocy,          & ! number of processors in y-direction
    nprocio,         & ! number of extra processors for doing asynchronous IO
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
! UB >>
!!$    my_peri_neigh(4)   ! periodic neighbors of this subdomain in the periodic 
!!$                       ! grid if it lies at the boundary of the total domain

  INTEGER   (KIND=iintegers), ALLOCATABLE       ::           &
    isubpos(:,:)       ! positions of the subdomains in the total domain. Given
                       ! are the i- and the j-indices of the lower left and the
                       ! upper right grid point in the order 
                       !                  i_ll, j_ll, i_ur, j_ur.
                       ! Only the interior of the domains are considered, not 
                       ! the boundary lines.

! 2. Further information for MPI
! ------------------------------

  INTEGER   (KIND=iintegers)       ::           &
    igroup_world,        & ! group that belongs to MPI_COMM_WORLD, i.e. all
                           ! processors
    icomm_world,         & ! communicator that belongs to igroup_world, i.e.
                           ! = MPI_COMM_WORLD
    icomm_compute,       & ! communicator for the group of compute PEs
    igroup_cart,         & ! group of the compute PEs
    icomm_cart,          & ! communicator for the virtual cartesian topology
    icomm_row,           & ! communicator for a east-west row of processors
    iexch_req(4),        & ! stores the sends requests for the neighbor-exchange
                           ! that can be used by MPI_WAIT to identify the send
    imp_reals,           & ! determines the correct REAL type used in the model
                           ! for MPI
    imp_grib,            & ! determines the REAL type used for the GRIB library
    imp_integers,        & ! determines the correct INTEGER type used in the
                           ! model for MPI
    imp_byte,            & ! determines the correct BYTE type used in the model
                           ! for MPI
    imp_character,       & ! determines the correct CHARACTER type used in the
                           ! model for MPI
    imp_logical            ! determines the correct LOGICAL   type used in the
                           ! model for MPI

  LOGICAL                                       &
    lcompute_pe,         & ! indicates whether this is a compute PE or not
    lreorder               ! during the creation of the virtual topology the
                           ! ranking of the processors may be reordered

  INTEGER (KIND=iintegers), ALLOCATABLE  ::           &
    igroup_diag(:),      & ! groups for the diagnostic domains
    icomm_diag (:)         ! communicators for the diagnostic domains

! 3. Static Send buffers
! ----------------------

  REAL (KIND=ireals), ALLOCATABLE  ::           &
    sendbuf(:,:)       ! sending buffer for boundary exchange:
                       ! 1-4 are used for sending, 5-8 are used for receiving

  INTEGER (KIND=iintegers)         ::           &
    isendbuflen        ! length of one column of sendbuf 

  ! Buffers for distributing the Namelists
  INTEGER (KIND=iintegers), ALLOCATABLE ::   intbuf  (:)
  REAL    (KIND=ireals)   , ALLOCATABLE ::   realbuf (:)
  LOGICAL                 , ALLOCATABLE ::   logbuf  (:)
  CHARACTER (LEN=100)     , ALLOCATABLE ::   charbuf (:)

!==============================================================================
! input file suffix
  integer :: cosmo_input_suffix

END MODULE data_parallel
