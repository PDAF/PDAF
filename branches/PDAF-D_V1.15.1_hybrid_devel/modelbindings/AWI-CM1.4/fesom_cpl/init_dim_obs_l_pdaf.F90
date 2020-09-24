!$Id: init_dim_obs_l_pdaf.F90 2136 2019-11-22 18:56:35Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_l_pdaf --- Set dimension of local observation vector
!
! !INTERFACE:
SUBROUTINE init_dim_obs_l_pdaf(domain_p, step, dim_obs_f, dim_obs_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the loop over
! all local analysis domains. It has to set 
! the dimension of the local observation vector 
! for the current local analysis domain.
!
! The routine is called by each filter process.
!
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mod_assim_pdaf, &
       ONLY: local_range, local_obs_nod2d, ocoord_n2d, locweight, &
       dim_ens, loc_radius, &
       ivariance_obs_l, mean_ice_p, obs_depth, distance
  USE mod_parallel_pdaf, &
       ONLY: mype_filter
  USE g_parfe, &
       ONLY: mydim_nod2d
  USE o_mesh, &
       ONLY: nod2D, coord_nod2D
  USE o_param, &
       ONLY: r_earth
  USE i_array, &
       ONLY: a_ice
  USE g_rotate_grid
 
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: domain_p   ! Current local analysis domain
  INTEGER, INTENT(in)  :: step       ! Current time step
  INTEGER, INTENT(in)  :: dim_obs_f  ! Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  ! Local dimension of observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_obs_l)
! Called by: PDAF_lestkf_update  (as U_init_dim_obs_l)
! Called by: PDAF_letkf_update   (as U_init_dim_obs_l)
!EOP

! *** Local variables ***
  INTEGER :: node, i, j         ! Counter
  REAL :: wc_coord(2)           ! Coordinates of current water column
  REAL :: o_coord(2)            ! Coordinates of observation
  REAL :: dist2d(2)             ! Distance vector between water column and observation
  REAL :: dist                  ! distance for adaptive localization
  REAL :: distance2, local_range2 ! squared distance and localization radius
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


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! *** Fixed localization radius ***

  IF (mype_filter==0 .AND. (domain_p<=domain_save .OR. first==1)) THEN
     WRITE (*,'(a,5x,a,es11.3,a)') 'FESOM-PDAF','--- Range limit for observation domain', &
          local_range,' m'
     first=0
  END IF
  domain_save = domain_p

  ! Initialize localization radius for current local domain
  ! loc_radius is initialized in init_dim_obs_f.f90
  local_range = loc_radius(domain_p)


! *** Set local observation dimension and indices ***

  ! Get location of current water column (basis point)
  CALL r2g(wc_coord(1), wc_coord(2), coord_nod2d(1, domain_p), coord_nod2d(2, domain_p))

  ! Initialize squared localization radius
  local_range2 = local_range*local_range

  ! Count local observations
  dim_obs_l = 0
  scanpoints_count: DO node = 1, dim_obs_f

     ! location of observation point
     o_coord(1) = ocoord_n2d(1, node)
     o_coord(2) = ocoord_n2d(2, node)

     ! approximate distances in longitude and latitude
     dist2d(1) = r_earth * ABS(wc_coord(1) - o_coord(1))* COS(wc_coord(2))
     dist2d(2) = r_earth * ABS(wc_coord(2) - o_coord(2))

     ! full squared distance in meters
     distance2 = dist2d(1)*dist2d(1) + dist2d(2)*dist2d(2)

     ! If distance below limit, add observation to local domain
     IF (distance2 <= local_range2) THEN
        dim_obs_l = dim_obs_l + 1             ! dimension
     END IF
  END DO scanpoints_count

  
  ! allocate arrays
  IF (ALLOCATED(local_obs_nod2d)) DEALLOCATE(local_obs_nod2d)
  ALLOCATE(local_obs_nod2d(dim_obs_l))
  IF (ALLOCATED(distance)) DEALLOCATE(distance)
  ALLOCATE(distance(dim_obs_l))

  dim_obs_l = 0

  ! Scan through full domain to initialize dimension and index array
  scanpoints: DO node = 1, dim_obs_f

     ! location of observation point
     o_coord(1) = ocoord_n2d(1, node)
     o_coord(2) = ocoord_n2d(2, node)

     ! approximate distances in longitude and latitude
     dist2d(1) = r_earth * ABS(wc_coord(1) - o_coord(1))* COS(wc_coord(2))
     dist2d(2) = r_earth * ABS(wc_coord(2) - o_coord(2))

     ! full square distance in meters
     distance2 = dist2d(1)*dist2d(1) + dist2d(2)*dist2d(2)

     ! If distance below limit, add observation to local domain
     IF (distance2 <= local_range2) THEN
        dim_obs_l = dim_obs_l + 1             ! dimension
        local_obs_nod2d(dim_obs_l) = node     ! node index
        distance(dim_obs_l) = SQRT(distance2) ! distance
     END IF

  END DO scanpoints

END SUBROUTINE init_dim_obs_l_pdaf

