!>  Interface routine to call initialization of PDAF
!!
!! This routine collects the initialization of variables for PDAF.
!! In addition, the initialization routine PDAF_init is called
!! to perform the internal initialization of PDAF.
!!
!! This variant is for the offline mode of PDAF.
!!
!! This routine is generic. However, it assumes a constant observation
!! error (rms_obs). Further, with parallelization the local state
!! dimension dim_state_p is used.
!!
!! __Revision history:__
!! * 2008-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_pdaf()

  USE PDAF                        ! PDAF
  USE mod_parallel_pdaf, &        ! Parallelization variables
       ONLY: mype_world, mype_filter
  USE mod_assimilation, &         ! Variables for assimilation
       ONLY: nx_p, nx, ny, ndim, dim_state_p, local_dims, coords_p, &
       screen, filtertype, subtype, dim_ens, &
       type_forget, forget, rank_ana_enkf, locweight, cradius, sradius, &
       type_trans, type_sqrt, pf_res_type, pf_noise_type, pf_noise_amp, &
       observe_ens, type_obs_init, do_omi_obsstats, &
       type_coords, coords_origin, coords_scale, deg2rad, &
       omi_search_type, omi_sort_dir
  USE obs_A_pdafomi, &            ! Variables for observation type A
       ONLY: assim_A, rms_obs_A
  USE obs_B_pdafomi, &            ! Variables for observation type B
       ONLY: assim_B, rms_obs_B

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: filter_param_i(2) ! Integer parameter array for filter
  REAL    :: filter_param_r(1) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag
  REAL    :: lim_coords(2,2)   ! limiting coordinates of process sub-domain
  INTEGER :: i, off_nx         ! Counters
  INTEGER :: off_p             ! Process-local offset in global state vector

! *** External subroutines ***
  EXTERNAL :: init_ens_offline  ! Ensemble initialization
  

! ***************************
! ***   Initialize PDAF   ***
! ***************************

  IF (mype_world == 0) THEN
     WRITE (*,'(/1x,a)') 'INITIALIZE PDAF - OFFLINE MODE'
  END IF


! **********************************************************
! ***   CONTROL OF PDAF - used in call to PDAF_init      ***
! **********************************************************

! *** IO options ***
  screen = 2         ! Write screen output (1) for output, (2) add timings

! *** Ensemble size ***
  dim_ens = 9        ! Size of ensemble for all ensemble filters

! *** Options for filter method

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++
  ! +++ For available options see MOD_ASSIMILATION +++
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++

  filtertype = 6     ! Type of filter
  subtype = 0        ! Subtype of filter

  forget = 1.0       ! Forgetting factor value for inflation
  type_forget = 0    ! Type of forgetting factor

  type_trans = 0     ! Type of ensemble transformation (deterministic or random)
  type_sqrt = 0      ! SEIK/LSEIK/ESTKF/LESTKF: Type of transform matrix square-root

  !EnKF
  rank_ana_enkf = 0  ! EnKF: rank to be considered for inversion of HPH in analysis step

  pf_res_type = 1    ! Resampling type for PF
  pf_noise_type = 0  ! Resampling type for PF
  pf_noise_amp = 0.0 ! Noise amplitude


! *********************************************************************
! ***   Settings for analysis steps  - used in call-back routines   ***
! *********************************************************************

! *** Which observation type to assimilate
  assim_A = .true.
  assim_B = .false.

! *** specifications for observations ***
  rms_obs_A = 0.5    ! Observation error standard deviation for observation A
  rms_obs_B = 0.5    ! Observation error standard deviation for observation B

! *** Localization settings
  locweight = 0     ! Type of localizating weighting
  cradius = 0.0     ! Cut-off radius in grid points for observation domain in local filters
  sradius = cradius ! Support radius for 5th-order polynomial
                    ! or radius for 1/e for exponential weighting


! ***********************************
! *** Some optional functionality ***
! ***********************************

! *** Parse command line options   ***
! *** This is optional, but useful ***

  call init_pdaf_parse()

! *** Activate PDAF-OMI observation statistics ***

  IF (do_omi_obsstats) CALL PDAFomi_set_obs_diag(1)

! *** Set search type for local observations ***

  CALL PDAFomi_set_searchtype(omi_search_type, omi_sort_dir)

! *** Initial Screen output ***
! *** This is optional      ***

  IF (mype_world == 0) call init_pdaf_info()


! *****************************************************
! *** Call PDAF initialization routine on all PEs.  ***
! ***                                               ***
! *** Here, the full selection of filters is        ***
! *** implemented. In a real implementation, one    ***
! *** reduces this to selected filters.             ***
! ***                                               ***
! *** For all filters, PDAF_init is first called    ***
! *** specifying only the required parameters.      ***
! *** Further settings are done afterwards using    ***
! *** calls to PDAF_set_iparam & PDAF_set_rparam.   ***
! *****************************************************

  ! *** Here we specify only the required integer and real parameters
  ! *** Other parameters are set using calls to PDAF_set_iparam/PDAF_set_rparam
  filter_param_i(1) = dim_state_p ! State dimension
  filter_param_i(2) = dim_ens     ! Size of ensemble
  filter_param_r(1) = forget      ! Forgetting factor

  CALL PDAF3_init(filtertype, subtype, 0, &
       filter_param_i, 2,&
       filter_param_r, 1, &
       init_ens_offline, screen, status_pdaf)

  ! *** Additional parameter specifications ***

  ! Generic settings for all filters
  CALL PDAF_set_iparam(5, type_forget, status_pdaf)
  CALL PDAF_set_iparam(6, type_trans, status_pdaf)
  CALL PDAF_set_iparam(7, type_sqrt, status_pdaf)
  CALL PDAF_set_iparam(8, observe_ens, status_pdaf)
  CALL PDAF_set_iparam(9, type_obs_init, status_pdaf)

  ! Specific settings
  IF (filtertype==PDAF_DA_ENKF) CALL PDAF_set_iparam(4, rank_ana_enkf, status_pdaf)
  if (filtertype==PDAF_DA_PF) CALL PDAF_set_iparam(6, pf_res_type, status_pdaf)
  if (filtertype==PDAF_DA_PF) CALL PDAF_set_rparam(3, pf_noise_amp, status_pdaf)

! *** Check whether initialization of PDAF was successful ***
  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in initialization of PDAF - stopping! (PE ', mype_world,')'
     CALL PDAF_abort(1)
  END IF


! ***************************************************
! *** Set coordinates of elements in state vector ***
! *** (used for localization in EnKF/ENSRF)       ***
! ***************************************************

  ALLOCATE(coords_p(ndim, dim_state_p))

  ! Global coordinates of local analysis domain
  ! We use grid point indices as coordinates, but could e.g. use meters
  ! The particular way to initializate coordinates is because in this
  ! offline example we do not split the model domain, but the state vector.
  off_p = 0
  DO i = 1, mype_filter
     off_p = off_p + local_dims(i)
  END DO

  DO i = 1, dim_state_p
     coords_p(1, i) = REAL(CEILING(REAL(i+off_p)/REAL(ny)))
     coords_p(2, i) = REAL(i+off_p) - (coords_p(1,i)-1)*REAL(ny)
  END DO

  IF (type_coords>1) THEN
     ! Geographic coordinates - scale and shift to origin
     coords_p(1, :) = deg2rad * (coords_origin(1) + coords_scale * (coords_p(1, :)-1.0))
     coords_p(2, :) = deg2rad * (coords_origin(2) + coords_scale * (coords_p(2, :)-1.0))
  END IF


! ************************************************************************
! *** Set domain coordinate limits (for use with OMI's use_global_obs) ***
! ************************************************************************
  
    ! Get offset of local domain in global domain in x-direction
    off_nx = 0
    DO i = 1, mype_filter
       off_nx = off_nx + nx_p
    END DO

    lim_coords(1,1) = REAL(off_nx + 1)     ! West
    lim_coords(1,2) = REAL(off_nx + nx_p)  ! East
    lim_coords(2,1) = REAL(ny)             ! North
    lim_coords(2,2) = 1.0                  ! South

    IF (type_coords>1) THEN
       ! Geographic coordinates - scale and shift to origin
       lim_coords(1,1) = deg2rad * (coords_origin(1) + coords_scale * (lim_coords(1,1)-1.0))  ! West
       lim_coords(1,2) = deg2rad * (coords_origin(1) + coords_scale * (lim_coords(1,2)-1.0))  ! East
       lim_coords(2,1) = deg2rad * (coords_origin(2) + coords_scale * (lim_coords(2,1)-1.0))  ! North
       lim_coords(2,2) = deg2rad * (coords_origin(2) + coords_scale * (lim_coords(2,2)-1.0))  ! South
    END IF

    CALL PDAFomi_set_domain_limits(lim_coords)

END SUBROUTINE init_pdaf
