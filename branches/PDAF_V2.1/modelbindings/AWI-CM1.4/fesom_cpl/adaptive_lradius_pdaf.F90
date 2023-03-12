!$Id: get_adaptive_lrange_pdaf.F90 2196 2020-03-26 13:26:59Z lnerger $
!>  Routine to compute adaptive localization radius
!!
!! The routine computes an adaptive localization
!! radius following Kirchgessner et al., MWR (2014).
!!
!! The routine is called by each filter process.
!!
!! __Revision history:__
!! 2012 - Lars Nerger - Initial code for FESOM
!! * Later revisions - see repository log
!!
SUBROUTINE get_adaptive_lradius_pdaf(domain_p, lradius, loc_radius)

  USE mod_assim_pdaf, &           ! Variables for assimilation
       ONLY: locweight, loctype, loc_ratio, dim_ens, eff_dim_obs, pi
  USE mod_parallel_pdaf, &        ! Parallelization variables
       ONLY: mype_filter
  USE g_parfe, &
       ONLY: mydim_nod2d
  USE o_mesh, &
       ONLY: coord_nod2D
  USE o_param, &
       ONLY: r_earth
  USE g_rotate_grid, &
       ONLY: r2g
  USE obs_sst_cmems_pdafomi, &
       ONLY: obs => thisobs
 
  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: domain_p     !< Current local analysis domain
  REAL, INTENT(out) :: lradius         !< Uniform localization radius
  REAL, INTENT(inout) :: loc_radius(mydim_nod2d)  !< Variable localization radius

! *** Local variables ***
  INTEGER :: node, i, j         ! Counter
  REAL :: wc_coord(2)           ! Coordinates of current water column
  REAL :: o_coord(2)            ! Coordinates of observation
  REAL :: dist2d(2)             ! Distance vector between water column and observation
  REAL :: dist                  ! distance for adaptive localization
  REAL :: distance2, lradius2   ! squared distance and localization radius
  INTEGER,SAVE :: allocflag=0   ! Allocation flag
  INTEGER,SAVE :: first=1       ! Flag for first call
  INTEGER,SAVE :: domain_save=1 ! Stored domain index
  INTEGER :: wtype              ! Type of weight function
  INTEGER :: rtype              ! Type of weight regulation
  REAL :: l_range               ! Localization radius
  REAL :: weight                ! Weight for observation localization
  LOGICAL, SAVE :: firstround = .true.   ! Flag for first analysis time
  REAL :: l_ranges(2)           ! Minimum and maximum allowed localization radii
  REAL :: l_step                ! Step size for search of localization radius


! *** Variable localization radius for fixed effective observation dimension ***

  ! Set parameters for localization radius search
  l_ranges(1) = 1.0e4
  l_ranges(2) = 1.0e7
  l_step     = 5.0e3

  IF (mype_filter==0 .AND. (domain_p<=domain_save .OR. first==1)) THEN
     WRITE (*,'(8x,a)') '--- Choose localization radius according to effective obs. dimension'
     WRITE (*,'(8x,a,es10.2)') '--- ratio to effective obs. dimension', loc_ratio
     WRITE (*,'(8x,a,2es10.2)') '--- minimum and maximum allowed radii', l_ranges
     first=0
  END IF
  domain_save = domain_p

  IF (locweight == 0) THEN
     ! Uniform (unit) weighting
     wtype = 0
     rtype = 0
  else if (locweight == 1 .OR. locweight == 3) THEN
     ! Exponential weighting
     wtype = 1
     rtype = 0
  ELSE IF (locweight == 2 .OR. locweight == 4 .OR. locweight == 5 &
       .OR. locweight == 6 .OR. locweight == 7) THEN
     ! 5th-order polynomial (Gaspari&Cohn, 1999)
     wtype = 2
     rtype = 0 ! Always non-regulated here
  END IF

  IF (firstround) THEN

     ! *** At the first analysis time, compute the loclaization radius ***

     ! Get location of current water column (basis point)
     ! wc_coord is the rotated coordinates in FESOM
     CALL r2g(wc_coord(1), wc_coord(2), coord_nod2d(1, domain_p), coord_nod2d(2, domain_p))

     lrange: DO l_range = l_ranges(1), l_ranges(2), l_step

        eff_dim_obs(domain_p) = 0

        ! Scan through full domain to initialize dimension and index array
        scanpointsB: DO node = 1, obs%dim_obs_f

           ! location of observation point
           o_coord(1) = obs%ocoord_f(1, node)
           o_coord(2) = obs%ocoord_f(2, node)
              
           ! approximate distances in longitude and latitude
           dist2d(1) = r_earth * MIN( ABS(wc_coord(1) - o_coord(1))* COS(wc_coord(2)), &
                ABS(ABS(wc_coord(1) - o_coord(1)) - 2.0*pi) * COS(wc_coord(2)))
           dist2d(2) = r_earth * ABS(wc_coord(2) - o_coord(2))

           ! full distance in meters
           dist = SQRT(dist2d(1)**2 + dist2d(2)**2)
           
           ! If distance below limit increment effective observation dimension
           IF (dist <= l_range) THEN
              CALL PDAF_local_weight(wtype, rtype, l_range, l_range, dist, &
                   1, 1, 1.0, 1.0, weight, 0)

              eff_dim_obs(domain_p) = eff_dim_obs(domain_p) + weight
           END IF

           IF (eff_dim_obs(domain_p) >= loc_ratio * dim_ens) THEN
              loc_radius(domain_p) = l_range
              EXIT lrange
           END IF

        END DO scanpointsB
           
     END DO lrange

     IF (domain_p==mydim_nod2d) firstround = .false.

  END IF

  lradius = loc_radius(domain_p)

END SUBROUTINE get_adaptive_lradius_pdaf

!-------------------------------------------------------
!>  Routine to compute statistics about adaptive localization radius
!!
!! The routine computes the minimum, maximum,
!! and mean values for the adaptive localization
!! radius following Kirchgessner et al., MWR (2014).
!!
!! The routine is called by each filter process.
!!
!! __Revision history:__
!! 2012 - Lars Nerger - Initial code for FESOM
!! * Later revisions - see repository log
!!
SUBROUTINE adaptive_lradius_stats_pdaf() 

  USE mod_parallel_pdaf, &
       ONLY: mype_filter=>mype_filter_fesom, npes=>npes_filter_fesom, &
       comm => COMM_filter_fesom, &
       MPI_REAL8, MPIerr, &
       MPI_INTEGER, MPI_SUM, MPI_MAX, MPI_MIN 
  USE mod_assim_pdaf, &       ! Variables for assimilation
       ONLY: eff_dim_obs, loctype
  USE g_parfe, &
       ONLY: dim_p => mydim_nod2d
  USE o_mesh, &
       ONLY: dim => nod2D

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: i                                   ! Counters
  REAL :: min_eff_dim_obs, max_eff_dim_obs       ! Stats on effective observation dimensions
  REAL :: min_eff_dim_obs_g, max_eff_dim_obs_g   ! Stats on effective observation dimensions
  REAL :: sum_eff_dim_obs, avg_eff_dim_obs_g     ! Stats on effective observation dimensions


! ***************************************************************
! *** Compute statistics for effective observation dimensions ***
! ***************************************************************

  IF (loctype==1) THEN

     max_eff_dim_obs = 0.0
     min_eff_dim_obs = 1.0e16
     sum_eff_dim_obs = 0.0

     DO i = 1, dim_p
        IF (eff_dim_obs(i) > max_eff_dim_obs) max_eff_dim_obs = eff_dim_obs(i)
        IF (eff_dim_obs(i) < min_eff_dim_obs) min_eff_dim_obs = eff_dim_obs(i)
        sum_eff_dim_obs = sum_eff_dim_obs + eff_dim_obs(i)
     END DO
     IF (npes>1) THEN
        CALL MPI_Reduce(sum_eff_dim_obs, avg_eff_dim_obs_g, 1, MPI_REAL8, MPI_SUM, &
             0, COMM, MPIerr)
        CALL MPI_Reduce(max_eff_dim_obs, max_eff_dim_obs_g, 1, MPI_REAL8, MPI_MAX, &
             0, COMM, MPIerr)
        CALL MPI_Reduce(min_eff_dim_obs, min_eff_dim_obs_g, 1, MPI_REAL8, MPI_MIN, &
             0, COMM, MPIerr)
     ELSE
        ! This is a work around for working with nullmpi.F90
        avg_eff_dim_obs_g = sum_eff_dim_obs
        min_eff_dim_obs_g = min_eff_dim_obs
        max_eff_dim_obs_g = max_eff_dim_obs
     END IF

     IF (mype_filter==0) THEN
        avg_eff_dim_obs_g = avg_eff_dim_obs_g / REAL(dim)

        WRITE (*, '(8x, a)') &
             '--- Effective observation dimensions for local analysis:'
        WRITE (*, '(12x, a, f12.2)') &
             'min. effective observation dimension:       ', min_eff_dim_obs_g
        WRITE (*, '(12x, a, f12.2)') &
             'max. effective observation dimension:       ', max_eff_dim_obs_g
        WRITE (*, '(12x, a, f12.2)') &
             'avg. effective observation dimension:       ', avg_eff_dim_obs_g
     END IF
  END IF

END SUBROUTINE adaptive_lradius_stats_pdaf
