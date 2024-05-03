! Copyright (c) 2004-2024 Lars Nerger
!
! This file is part of PDAF.
!
! PDAF is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version
! 3 of the License, or (at your option) any later version.
!
! PDAF is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with PDAF.  If not, see <http://www.gnu.org/licenses/>.
!
!$Id: PDAFomi_obs_f.F90 1147 2023-03-12 16:14:34Z lnerger $

!> PDAF-OMI routines for full observations
!!
!! This module contains subroutines to handle full observations. 
!! Further, it contains routines to restrict the global full vector of observations
!! to those observations that are relevant for a process-local model subdomain.
!! The routines are
!!
!! * PDAFomi_gather_obs \n
!!        Gather full observation information
!! * PDAFomi_gather_obsstate \n
!!        Gather a full observed state vector (used in observation operators)
!! * PDAFomi_init_obs_f \n
!!        Initialize full vector of observations for adaptive forgetting factor
!! * PDAFomi_init_obsvar_f \n
!!        Compute mean observation error variance for adaptive forgetting factor
!! * PDAFomi_prodRinvA \n
!!        Multiply an intermediate matrix of the global filter analysis
!!        with the inverse of the observation error covariance matrix
!! * PDAFomi_likelihood \n
!!        Compute likelihood for an ensemble member
!! * PDAFomi_add_obs_err \n
!!        Add observation error to some matrix
!! * PDAFomi_init_obscovar \n
!!        Initialize global observation error covariance matrix
!! * PDAFomi_set_domain_limits \n
!!        Set min/max coordinate locations of decomposed grid
!! * PDAFomi_get_domain_limits_unstr \n
!!        Find min/max coordinate locations in unstructured grid
!! * PDAFomi_get_local_ids_obs_f \n
!!        Find observations inside or close to process domain
!! * PDAFomi_limit_obs_f \n
!!        Reduce full observation vector to part relevant for local process domain
!!
!! The routine PDAFomi_get_domain_limits_unstr assumed geographic coordinates in radians
!! and within the range -pi to +pi for longitude (- is westward) and -pi/2 to +pi/2 for
!! latitude.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
MODULE PDAFomi_obs_f

  USE mpi
  USE PDAF_mod_filtermpi, &
       ONLY: mype, npes, COMM_FILTER, MPIerr
  USE PDAF_mod_filter, &
       ONLY: screen, obs_member, filterstr, dim_p

  IMPLICIT NONE
  SAVE

! *** Module internal variables
  INTEGER :: debug=0                    !< Debugging flag
  INTEGER :: error=0                    !< Error flag

  REAL, ALLOCATABLE :: domain_limits(:) !< Limiting coordinates (NSWE) for process domain
  REAL, PARAMETER :: r_earth=6.3675e6   !< Earth radius in meters
  REAL, PARAMETER :: pi=3.141592653589793   !< Pi

! *** Data type to define the full observations by internally shared variables of the module
  TYPE obs_f
     ! ---- Mandatory variables to be set in INIT_DIM_OBS ----
     INTEGER :: doassim=0                 !< Whether to assimilate this observation type
     INTEGER :: disttype                  !< Type of distance computation to use for localization
                                          !<  (0) Cartesian, (1) Cartesian periodic
                                          !<  (2) simplified geographic, (3) geographic haversine function
                                          !<  (10,11,12,13) factorized 2+1D localization with distance
                                          !<    calculation from (0)-(3); obs. weighting is only done with
                                          !<    horizontal distance, which vertical uses only cut-off radius
     INTEGER :: ncoord                    !< Number of coordinates use for distance computation
     INTEGER, ALLOCATABLE :: id_obs_p(:,:) !< Indices of process-local observed field in state vector

     ! ---- Optional variables - they can be set in INIT_DIM_OBS ----
     REAL, ALLOCATABLE :: icoeff_p(:,:)   !< Interpolation coefficients for obs. operator (optional)
     REAL, ALLOCATABLE :: domainsize(:)   !< Size of domain for periodicity (<=0 for no periodicity) (optional)

     ! ---- Variables with predefined values - they can be changed in INIT_DIM_OBS  ----
     INTEGER :: obs_err_type=0            !< Type of observation error: (0) Gauss, (1) Laplace
     INTEGER :: use_global_obs=1          !< Whether to use (1) global full obs. 
                                          !< or (0) obs. restricted to those relevant for a process domain
     REAL :: inno_omit=0.0                !< Omit obs. if squared innovation larger this factor times
                                          !<     observation variance (only active for >0)
     REAL :: inno_omit_ivar=1.0e-12       !< Value of inverse variance to omit observation
                                          !<     (should be much larger than actual observation error variance)

     ! ----  The following variables are set in the routine PDAFomi_gather_obs ---
     INTEGER :: dim_obs_p                 !< number of PE-local observations
     INTEGER :: dim_obs_f                 !< number of full observations
     INTEGER :: dim_obs_g                 !< global number of observations
     INTEGER :: off_obs_f                 !< Offset of this observation in overall full obs. vector
     INTEGER :: off_obs_g                 !< Offset of this observation in overall global obs. vector
     INTEGER :: obsid                     !< Index of observation over all assimilated observations
     REAL, ALLOCATABLE :: obs_f(:)        !< Full observed field
     REAL, ALLOCATABLE :: ocoord_f(:,:)   !< Coordinates of full observation vector
     REAL, ALLOCATABLE :: ivar_obs_f(:)   !< Inverse variance of full observations
     INTEGER, ALLOCATABLE :: id_obs_f_lim(:) !< Indices of domain-relevant full obs. in global vector of obs.

     ! ----  Other internal variables ---
     INTEGER :: locweight_v               !< Type of localization function in vertical direction (for disttype>=10)
  END TYPE obs_f

  INTEGER :: n_obstypes = 0               ! Number of observation types
  INTEGER :: obscnt = 0                   ! current ID of observation type
  INTEGER :: offset_obs = 0               ! offset of current observation in overall observation vector
  INTEGER :: offset_obs_g = 0             ! offset of current observation in global observation vector
  LOGICAL :: omit_obs = .FALSE.           ! Flag whether observations are omitted for large innovation
  LOGICAL :: omi_was_used = .FALSE.       ! Flag whether OMI was used 
  INTEGER, ALLOCATABLE :: obsdims(:,:)    ! Observation dimensions over all types and process sub-domains
  INTEGER, ALLOCATABLE :: map_obs_id(:)   ! Index array to map obstype-first index to domain-first index

  INTEGER :: ostats_omit(7)             ! PE-local statistics
  ! ostats_omit(1): Number of local domains with excluded observations
  ! ostats_omit(2): Number of local domains without excluded observations
  ! ostats_omit(3): Sum of all excluded observations for all domains
  ! ostats_omit(4): Sum of all used observations for all domains
  ! ostats_omit(5): Sum of all used observations for all domains
  ! ostats_omit(6): Maximum count of excluded observations over all domains
  ! ostats_omit(7): Maximum count of used observations over all domains


  TYPE obs_arr_f
     TYPE(obs_f), POINTER :: ptr
  END TYPE obs_arr_f

  TYPE(obs_arr_f), ALLOCATABLE :: obs_f_all(:)


!$OMP THREADPRIVATE(debug)


!-------------------------------------------------------------------------------
  
CONTAINS
!> Gather full observational information
!!
!! This routine uses different gather routines from PDAF to obtain
!! combined full observational information.
!!
!! The observation-type specific variables that are initialized here are
!! * thisobs\%dim_obs_p   - PE-local number of module-type observations
!! * thisobs\%dim_obs_f   - full number of module-type observations
!! * thisobs\%off_obs_f   - Offset of full module-type observation in overall full obs. vector
!! * thisobs\%off_obs_g   - Offset of global module-type observation in overall full obs. vector
!! * thisobs\%obs_f       - full vector of module-type observations
!! * thisobs\%ocoord_f    - coordinates of observations in OBS_MOD_F
!! * thisobs\%ivar_obs_f  - full vector of inverse obs. error variances of module-type
!!
!! If the full observations are restricted to those relevant for a process domain, 
!! there are further initialized 
!! * thisobs\%dim_obs_g  - Number of global observations
!! * thisobs\%id_obs_f_lim - Ids of full observations in global observations
!!
!! __Revision history:__
!! * 2020-03 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
       ncoord, lradius, dim_obs_f)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs   !< Data type with full observation
    INTEGER, INTENT(in) :: dim_obs_p        !< Number of process-local observation
    REAL, INTENT(in) :: obs_p(:)            !< Vector of process-local observations
    REAL, INTENT(in) :: ivar_obs_p(:)       !< Vector of process-local inverse observation error variance
    REAL, INTENT(in) :: ocoord_p(:,:)       !< Array of process-local observation coordinates
    INTEGER, INTENT(in) :: ncoord           !< Number of rows of coordinate array
    REAL, INTENT(in) :: lradius             !< Localization radius (the maximum radius used in this process domain) 
    INTEGER, INTENT(out) :: dim_obs_f       !< Full number of observations

! *** Local variables ***
    INTEGER :: i                            ! Counter
    REAL, ALLOCATABLE :: obs_g(:)           ! Global full observation vector (used in case of limited obs.)
    REAL, ALLOCATABLE :: ivar_obs_g(:)      ! Global full inverse variances (used in case of limited obs.)
    REAL, ALLOCATABLE :: ocoord_g(:,:)      ! Global full observation coordinates (used in case of limited obs.)
    INTEGER :: status                       ! Status flag for PDAF gather operation
    INTEGER :: localfilter                  ! Whether the filter is domain-localized
    INTEGER :: globalobs                    ! Whether the filter needs global observations
    INTEGER :: maxid                        ! maximum index in thisobs%id_obs_p


! **********************
! *** Initialization ***
! **********************

    ! Increment counter of observation types
    n_obstypes = n_obstypes + 1

    ! Set observation ID
    thisobs%obsid = n_obstypes

    ! Set flag that OMI was used
    omi_was_used = .TRUE.

    IF (mype == 0 .AND. screen > 0) &
         WRITE (*, '(a, 5x, a, 1x, i3)') 'PDAFomi', '--- Initialize observation type ID', thisobs%obsid


! **************************************
! *** Gather full observation arrays ***
! **************************************

    ! Consistency check
    maxid = MAXVAL(thisobs%id_obs_p)
    IF (maxid > dim_p .AND. dim_obs_p>0) THEN
       ! Maximum value of id_obs_p point to outside of state vector
       WRITE (*,'(a)') 'PDAFomi - ERROR: thisobs%id_obs_p too large - index points to outside of state vector !!!'
       error = 1
    END IF

    ! Check  whether the filter is domain-localized
    CALL PDAF_get_localfilter(localfilter)

    ! Check  whether the filter needs global observations
    CALL PDAF_get_globalobs(globalobs)

    ! Print debug information
    IF (debug>0) THEN
       WRITE (*,*) '++ OMI-debug: ', debug, &
            'PDAFomi_gather_obs -- START Gather full observation vector'
       IF (localfilter==1) THEN
          WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'domain localized filter'
       ELSE
          WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'filter without domain-localization'
       END IF
       IF (globalobs==1) THEN
          WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'filter uses global observations'
       END IF
    END IF

    lfilter: IF (localfilter==1 .OR. globalobs==1) THEN

       ! For domain-localized filters: gather full observations

       fullobs: IF (thisobs%use_global_obs==1 .OR. globalobs==1) THEN

          ! *** Use global full observations ***

          IF (mype == 0 .AND. screen > 0) &
               WRITE (*, '(a, 5x, a)') 'PDAFomi', '--- Use global full observations'

          ! *** Initialize global dimension of observation vector ***
          CALL PDAF_gather_dim_obs_f(dim_obs_p, dim_obs_f)

          ! Store full and PE-local observation dimensions in module variables
          thisobs%dim_obs_p = dim_obs_p
          thisobs%dim_obs_f = dim_obs_f

          IF (mype == 0 .AND. screen > 0) &
               WRITE (*, '(a, 8x, a, i7)') 'PDAFomi', &
               '--- Number of full observations ', dim_obs_f

          IF (debug>0) THEN
             WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'thisobs%dim_obs_p', thisobs%dim_obs_p
             WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'thisobs%dim_obs_f', thisobs%dim_obs_f
             WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'obs_p', obs_p
             WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'ocoord_p', ocoord_p
             WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'thisobs%disttype', thisobs%disttype
          END IF

          ! *** Gather full observation vector and corresponding coordinates ***

          ! Allocate full observation arrays
          ! The arrays are deallocated in PDAFomi_deallocate_obs in PDAFomi_obs_l
          IF (dim_obs_f > 0) THEN
             ALLOCATE(thisobs%obs_f(dim_obs_f))
             ALLOCATE(thisobs%ivar_obs_f(dim_obs_f))
             ALLOCATE(thisobs%ocoord_f(ncoord, dim_obs_f))
          ELSE
             ALLOCATE(thisobs%obs_f(1))
             ALLOCATE(thisobs%ivar_obs_f(1))
             ALLOCATE(thisobs%ocoord_f(ncoord, 1))
          END IF

          CALL PDAFomi_gather_obs_f_flex(dim_obs_p, obs_p, thisobs%obs_f, status)
          CALL PDAFomi_gather_obs_f_flex(dim_obs_p, ivar_obs_p, thisobs%ivar_obs_f, status)
          CALL PDAFomi_gather_obs_f2_flex(dim_obs_p, ocoord_p, thisobs%ocoord_f, ncoord, status)

       ELSE fullobs

          ! *** Use full observations limited to those relevant for a process domain ***
          ! *** This can be more efficient as in the local analysis loop less        ***
          ! *** observations have a be checked for each analysis domain              ***

          IF (mype == 0 .AND. screen > 0) &
               WRITE (*, '(a, 5x, a)') 'PDAFomi', '--- Use limited full observations'

          ! *** Initialize global dimension of observation vector ***
          CALL PDAF_gather_dim_obs_f(dim_obs_p, thisobs%dim_obs_g)

          ! *** First gather global observation vector and corresponding coordinates ***

          ! Allocate global observation arrays
          IF (thisobs%dim_obs_g > 0) THEN
             ALLOCATE(obs_g(thisobs%dim_obs_g))
             ALLOCATE(ivar_obs_g(thisobs%dim_obs_g))
             ALLOCATE(ocoord_g(ncoord, thisobs%dim_obs_g))
          ELSE
             ALLOCATE(obs_g(1))
             ALLOCATE(ivar_obs_g(1))
             ALLOCATE(ocoord_g(ncoord, 1))
          END IF

          CALL PDAFomi_gather_obs_f_flex(dim_obs_p, obs_p, obs_g, status)
          CALL PDAFomi_gather_obs_f_flex(dim_obs_p, ivar_obs_p, ivar_obs_g, status)
          CALL PDAFomi_gather_obs_f2_flex(dim_obs_p, ocoord_p, ocoord_g, ncoord, status)


          ! *** Now restrict the global observation arrays to the process-relevant parts ***

          ! Get number of full observation relevant for the process domain
          ! and corresponding indices in global observation vector
     
          IF (thisobs%dim_obs_g > 0) THEN
             ALLOCATE(thisobs%id_obs_f_lim(thisobs%dim_obs_g))
          ELSE
             ALLOCATE(thisobs%id_obs_f_lim(1))
          END IF
          IF (ALLOCATED(thisobs%domainsize)) THEN
             CALL PDAFomi_get_local_ids_obs_f(thisobs%dim_obs_g, lradius, ocoord_g, dim_obs_f, &
                  thisobs%id_obs_f_lim, thisobs%disttype, thisobs%domainsize)
          ELSE
             CALL PDAFomi_get_local_ids_obs_f(thisobs%dim_obs_g, lradius, ocoord_g, dim_obs_f, &
                  thisobs%id_obs_f_lim, thisobs%disttype)
          END IF

          ! Store full and PE-local observation dimensions in module variables
          thisobs%dim_obs_p = dim_obs_p
          thisobs%dim_obs_f = dim_obs_f

          IF (debug>0) THEN
             WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'thisobs%dim_obs_p', thisobs%dim_obs_p
             WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'thisobs%dim_obs_f', thisobs%dim_obs_f
             WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'obs_p', obs_p
             WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'ocoord_p', ocoord_p
          END IF

          ! Allocate global observation arrays
          ! The arrays are deallocated in PDAFomi_deallocate_obs in PDAFomi_obs_l
          IF (dim_obs_f > 0) THEN
             ALLOCATE(thisobs%obs_f(dim_obs_f))
             ALLOCATE(thisobs%ivar_obs_f(dim_obs_f))
             ALLOCATE(thisobs%ocoord_f(ncoord, dim_obs_f))

             ! Get process-relevant full observation arrays
             CALL PDAFomi_limit_obs_f(thisobs, 0, obs_g, thisobs%obs_f)
             CALL PDAFomi_limit_obs_f(thisobs, 0, ivar_obs_g, thisobs%ivar_obs_f)
             DO i = 1, ncoord
                CALL PDAFomi_limit_obs_f(thisobs, 0, ocoord_g(i,:), thisobs%ocoord_f(i,:))
             END DO
          ELSE
             ALLOCATE(thisobs%obs_f(1))
             ALLOCATE(thisobs%ivar_obs_f(1))
             ALLOCATE(thisobs%ocoord_f(ncoord, 1))
          END IF

          IF (debug>0) THEN
             WRITE (*,*) '++ OMI-debug: ', debug, '   PDAFomi_gather_obs -- Limited full observations'
             WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'thisobs%dim_obs_g', thisobs%dim_obs_g
             WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'obs_g', obs_g
          END IF
          DEALLOCATE(obs_g, ivar_obs_g, ocoord_g)

       END IF fullobs

    ELSE lfilter

       ! *** For global filters use process-local observations without gathering ***

       IF (mype == 0 .AND. screen > 0) &
            WRITE (*, '(a, 5x, a)') 'PDAFomi', '--- Use process-local observations for global filters'

       ! *** Initialize global dimension of observation vector ***
       dim_obs_f = dim_obs_p


       ! *** Initialize global dimension of observation vector ***
       CALL PDAF_gather_dim_obs_f(dim_obs_p, thisobs%dim_obs_g)

       IF (TRIM(filterstr)=='ENKF' .OR. TRIM(filterstr)=='LENKF') THEN
          IF (mype == 0 .AND. screen > 0) &
               WRITE (*, '(a, 8x, a, i7)') 'PDAFomi', &
               '--- Number of global observations ', thisobs%dim_obs_g
       ELSE
          IF (mype == 0 .AND. screen > 0) &
               WRITE (*, '(a, 8x, a, i7)') 'PDAFomi', &
               '--- Number of full observations ', dim_obs_f
       END IF

       ! *** Gather full observation vector and corresponding coordinates ***

       ! Allocate full observation arrays
       ! The arrays are deallocated in PDAFomi_deallocate_obs in PDAFomi_obs_l
       IF (dim_obs_f > 0) THEN
          ALLOCATE(thisobs%obs_f(dim_obs_f))
          IF (TRIM(filterstr)=='ENKF' .OR. TRIM(filterstr)=='LENKF') THEN
             ! The LEnKF needs the global array ivar_obs_f
             ALLOCATE(thisobs%ivar_obs_f(thisobs%dim_obs_g))
             ALLOCATE(thisobs%ocoord_f(ncoord, thisobs%dim_obs_g))
          ELSE
             ALLOCATE(thisobs%ivar_obs_f(dim_obs_f))
             ALLOCATE(thisobs%ocoord_f(ncoord, dim_obs_f))
          END IF
       ELSE
          ALLOCATE(thisobs%obs_f(1))
          IF (thisobs%dim_obs_g>0 .AND. (TRIM(filterstr)=='ENKF' .OR. TRIM(filterstr)=='LENKF')) THEN
             ! The LEnKF needs the global array ivar_obs_f
             ! Here dim_obs_f=0, but dim_obs_g>0 is possible in case of domain-decomposition
             IF (thisobs%dim_obs_g>0) THEN
                ALLOCATE(thisobs%ivar_obs_f(thisobs%dim_obs_g))
                ALLOCATE(thisobs%ocoord_f(ncoord, thisobs%dim_obs_g))
             ELSE
                ALLOCATE(thisobs%ivar_obs_f(1))
                ALLOCATE(thisobs%ocoord_f(ncoord, 1))
             END IF
          ELSE
             ALLOCATE(thisobs%ivar_obs_f(1))
             ALLOCATE(thisobs%ocoord_f(ncoord, 1))
          END IF
       END IF

       thisobs%obs_f = obs_p

       IF (TRIM(filterstr)=='ENKF' .OR. TRIM(filterstr)=='LENKF') THEN
          ! The EnKF and LEnKF need the global array ivar_obs_f
          CALL PDAFomi_gather_obs_f_flex(dim_obs_p, ivar_obs_p, &
               thisobs%ivar_obs_f, status)
          CALL PDAFomi_gather_obs_f2_flex(dim_obs_p, ocoord_p, &
               thisobs%ocoord_f, ncoord, status)
       ELSE
          thisobs%ivar_obs_f = ivar_obs_p
       END IF

       ! Store full and PE-local observation dimensions in module variables
       thisobs%dim_obs_p = dim_obs_p
       thisobs%dim_obs_f = dim_obs_f

       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'thisobs%dim_obs_p', thisobs%dim_obs_p
          WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'thisobs%dim_obs_f', thisobs%dim_obs_f
          IF (TRIM(filterstr)=='ENKF' .OR. TRIM(filterstr)=='LENKF') &
               WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'thisobs%dim_obs_g', thisobs%dim_obs_g
          WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'obs_p', obs_p
       END IF

    END IF lfilter

    ! set observation offset
    thisobs%off_obs_f = offset_obs
    offset_obs = offset_obs + thisobs%dim_obs_f

    ! set global observation offset for EnKF/LEnKF
    IF (TRIM(filterstr)=='ENKF' .OR. TRIM(filterstr)=='LENKF') THEN
       thisobs%off_obs_g = offset_obs_g
       offset_obs_g = offset_obs_g + thisobs%dim_obs_g
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'thisobs%off_obs_g', thisobs%off_obs_g
       END IF
    END IF

    ! Set general flag if observations are omitted for large innovation
    IF (thisobs%inno_omit > 0.0) omit_obs = .TRUE.

    ! Initialize statistics for observations omitted for large innovation
    ostats_omit = 0

    ! Screen output
    IF (mype == 0 .AND. screen > 0.AND. thisobs%inno_omit > 0.0) THEN
       WRITE (*, '(a, 5x, a, 1x, i3, 1x , a, f8.2,a)') &
            'PDAFomi', '--- Exclude obs. type ID', n_obstypes, ' if innovation^2 > ', &
            thisobs%inno_omit,' times obs. error variance'
    END IF

    ! Print debug information
    IF (debug>0) THEN
       WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'thisobs%obs_f', thisobs%obs_f
       WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'thisobs%ivar_obs_f', thisobs%ivar_obs_f
       WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'thisobs%ocoord_f', thisobs%ocoord_f
       WRITE (*,*) '++ OMI-debug gather_obs:      ', debug, 'initialized obs. ID', thisobs%obsid
       WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_gather_obs -- END'
    END IF

  END SUBROUTINE PDAFomi_gather_obs



!-------------------------------------------------------------------------------
!> Gather full observational information
!!
!! This routine uses PDAFomi_gather_obs_f_flex to obtain
!! a full observed state vector. The routine is usually called
!! the observation operators.
!!
!! __Revision history:__
!! * 2020-05 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_gather_obsstate(thisobs, obsstate_p, obsstate_f)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), TARGET, INTENT(inout) :: thisobs  !< Data type with full observation
    REAL, INTENT(in) :: obsstate_p(:)      !< Vector of process-local observed state
    REAL, INTENT(inout) :: obsstate_f(:)   !< Full observed vector for all types

! *** Local variables ***
    INTEGER :: status                      ! Status flag for PDAF gather operation
    INTEGER :: localfilter                 ! Whether the filter is domain-localized
    INTEGER :: globalobs                   ! Whether the filter needs global observations
    REAL, ALLOCATABLE :: obsstate_tmp(:)   ! Temporary vector of globally full observations


! **************************************
! *** Gather full observation arrays ***
! **************************************

    ! Check  whether the filter is domain-localized
    CALL PDAF_get_localfilter(localfilter)

    ! Check  whether the filter needs global observations
    CALL PDAF_get_globalobs(globalobs)

    ! Print debug information
    IF (debug>0) THEN
       IF (obs_member==0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, &
               '  PDAFomi_gather_obsstate -- START Gather full observed ensemble mean'
       ELSE
          WRITE (*,*) '++ OMI-debug: ', debug, &
               '  PDAFomi_gather_obsstate -- START Gather full observed ensemble state', obs_member
       END IF
       WRITE (*,*) '++ OMI-debug gather_obsstate: ', debug, 'observation ID', thisobs%obsid
       WRITE (*,*) '++ OMI-debug gather_obsstate: ', debug, 'thisobs%dim_obs_p', thisobs%dim_obs_p
       WRITE (*,*) '++ OMI-debug gather_obsstate: ', debug, 'thisobs%dim_obs_f', thisobs%dim_obs_f
       WRITE (*,*) '++ OMI-debug gather_obsstate: ', debug, 'thisobs%off_obs_f', thisobs%off_obs_f
       IF (thisobs%use_global_obs==0) &
            WRITE (*,*) '++ OMI-debug gather_obsstate: ', debug, 'thisobs%dim_obs_g', thisobs%dim_obs_g
       WRITE (*,*) '++ OMI-debug gather_obsstate: ', debug, 'obsstate_p', obsstate_p
    END IF

    lfilter: IF (localfilter==1 .OR. globalobs==1) THEN

       ! For domain-localized filters: gather full observations

       fullobs: IF (thisobs%use_global_obs==1 .OR. globalobs==1) THEN

          ! *** Gather global full observation vector ***

          CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, obsstate_p, &
               obsstate_f(thisobs%off_obs_f+1 : thisobs%off_obs_f+thisobs%dim_obs_f), status)

       ELSE fullobs

          ! *** Use full observations limited to those relevant for a process domain ***
          ! *** This can be more efficient as in the local analysis loop less        ***
          ! *** observations have a be checked for each analysis domain              ***

          IF (thisobs%dim_obs_g>0) THEN
             ALLOCATE(obsstate_tmp(thisobs%dim_obs_g))
          ELSE
             ALLOCATE(obsstate_tmp(1))
          END IF

          ! *** Gather observation vector ***
          CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, obsstate_p, &
               obsstate_tmp, status)

          ! Now restrict observation vector to process-relevant part
          CALL PDAFomi_limit_obs_f(thisobs, thisobs%off_obs_f, obsstate_tmp, obsstate_f)

          DEALLOCATE(obsstate_tmp)

       END IF fullobs

    ELSE lfilter

       ! *** For global filters use process-local observations without gathering ***

       ! In case of a global filter store process-local observed state
       obsstate_f(thisobs%off_obs_f+1 : thisobs%off_obs_f+thisobs%dim_obs_p) &
            = obsstate_p(1:thisobs%dim_obs_p)

    END IF lfilter

    IF (debug>0) THEN
       WRITE (*,*) '++ OMI-debug gather_obsstate: ', debug, &
            'obsstate_f', obsstate_f(thisobs%off_obs_f+1 : thisobs%off_obs_f+thisobs%dim_obs_p)
    END IF

    ! Initialize pointer array
    IF (obscnt == 0) THEN
       IF (.NOT.ALLOCATED(obs_f_all)) ALLOCATE(obs_f_all(n_obstypes))
       obscnt = 1

       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug gather_obsstate: ', debug, &
               'initialize pointer array for ', n_obstypes, 'observation types'
       END IF
    END IF

    ! Set pointer to current observation
    obs_f_all(thisobs%obsid)%ptr => thisobs

    IF (debug>0) THEN
       WRITE (*,*) '++ OMI-debug: ', debug, '  PDAFomi_gather_obsstate -- END'
    END IF

  END SUBROUTINE PDAFomi_gather_obsstate



!-------------------------------------------------------------------------------
!> Initialize full vector of observations
!!
!! This routine initializes the part of the full vector of
!! observations for the current observation type.
!! It has to fill the observations to obsstate_f from
!! position OFFSET_OBS+1. For the return value OFFSET_OBS
!! has to be incremented by the number of added observations.
!! The routine will only be called if the adaptive forgetting
!! factor is used.
!!
!! __Revision history:__
!! * 2019-09 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_init_obs_f(thisobs, dim_obs_f, obsstate_f, offset)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    INTEGER, INTENT(in) :: dim_obs_f       !< Dimension of full observed state (all observed fields)
    REAL, INTENT(inout) :: obsstate_f(:)   !< Full observation vector (dim_obs_f)
    INTEGER, INTENT(inout) :: offset       !< input: offset of module-type observations in obsstate_f
                                           !< output: input + number of added observations


! ******************************************
! *** Initialize full observation vector ***
! ******************************************

    ! Consistency check
    IF (dim_obs_f < offset+thisobs%dim_obs_f) THEN
       WRITE (*,'(a)') 'PDAFomi - ERROR: PDAFomi_init_obs_f - dim_obs_f is too small !!!'
       error = 2
    END IF

    doassim: IF (thisobs%doassim == 1) THEN

       ! Print debug information
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, &
               'PDAFomi_init_obs_f -- START Initialize observation vector'
          WRITE (*,*) '++ OMI-debug init_obs_f:        ', debug, 'observation ID', thisobs%obsid
          WRITE (*,*) '++ OMI-debug init_obs_f:        ', debug, 'thisobs%dim_obs_f', thisobs%dim_obs_f
          WRITE (*,*) '++ OMI-debug init_obs_f:        ', debug, 'thisobs%obs_f', thisobs%obs_f
       END IF

       ! Fill part of full observation vector
       obsstate_f(offset+1 : offset+thisobs%dim_obs_f) = thisobs%obs_f(1 : thisobs%dim_obs_f)

       ! Increment offset
       offset = offset + thisobs%dim_obs_f

       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, &
               'PDAFomi_init_obs_f -- END'
       END IF

    END IF doassim

  END SUBROUTINE PDAFomi_init_obs_f



!-------------------------------------------------------------------------------
!> Compute mean observation error variance
!!
!! This routine will only be called, if the adaptive
!! forgetting factor feature is used. Please note that
!! this is an experimental feature.
!!
!! The routine is called in global filters (like ESTKF)
!! during the analysis or in local filters (e.g. LESTKF)
!! before the loop over local analysis domains 
!! by the routine PDAF_set_forget that estimates an 
!! adaptive forgetting factor.  The routine has to 
!! initialize the mean observation error variance.  
!! For global filters this should be the global mean,
!! while for local filters it should be the mean for the
!! PE-local  sub-domain. (init_obsvar_l_TYPE is the 
!! localized variant for local filters)
!!
!! The routine assumes a diagonal observation error
!! covariance matrix.
!!
!! If the observation counter is zero the computation
!! of the mean variance is initialized. The output is 
!! always the mean variance. If the observation counter
!! is >0 the computed was initialized for another 
!! observation. To proceed the computation, meanvar is
!! first multiplied with the observation counter to 
!! obtain the variance sum. Then the computation of the 
!! mean is continued.
!!
!! __Revision history:__
!! * 2019-09 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_init_obsvar_f(thisobs, meanvar, cnt_obs)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    REAL, INTENT(inout) :: meanvar         !< Mean variance
    INTEGER, INTENT(inout) :: cnt_obs      !< Observation counter

! Local variables ***
    INTEGER :: i        ! Counter


! ***********************************
! *** Compute local mean variance ***
! ***********************************

    doassim: IF (thisobs%doassim == 1) THEN

       IF (cnt_obs==0) THEN
          ! Reset mean variance
          meanvar = 0.0
       ELSE
          ! Compute sum of variances from mean variance
          meanvar = meanvar * REAL(cnt_obs)
       END IF

       ! Add observation error variances
       DO i = 1, thisobs%dim_obs_f
          meanvar = meanvar + 1.0 / thisobs%ivar_obs_f(i)
       END DO

       ! Increment observation count
       cnt_obs = cnt_obs + thisobs%dim_obs_f

       ! Compute updated mean variance
       meanvar = meanvar / REAL(cnt_obs)

    END IF doassim

  END SUBROUTINE PDAFomi_init_obsvar_f



!-------------------------------------------------------------------------------
!> Compute product of inverse of R with some matrix
!!
!! The routine is called during the analysis step
!! of the global square-root filters. It has to 
!! compute the product of the inverse of the
!! process-local observation error covariance matrix
!! with the matrix of process-local observed ensemble 
!! perturbations.
!!
!! This routine assumes a diagonal observation error
!! covariance matrix, but allows for varying observation
!! error variances.
!!
!! The routine can be applied with either all observations
!! of different types at once, or separately for each
!! observation type. The operation is done with all
!! process-local observations
!!
!! __Revision history:__
!! * 2019-12 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_prodRinvA(thisobs, ncols, A_p, C_p)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs !< Data type with full observation
    INTEGER, INTENT(in) :: ncols          !< Number of columns in A_p and C_p
    REAL, INTENT(in) :: A_p(:, :)         !< Input matrix (nobs_f, ncols)
    REAL, INTENT(out)   :: C_p(:, :)      !< Output matrix (nobs_f, ncols)


! *** local variables ***
    INTEGER :: i, j       ! index of observation component
    INTEGER :: off        ! row offset in A_l and C_l
    

! *************************************
! ***                -1             ***
! ***           C = R   A           ***
! ***                               ***
! *** The inverse observation error ***
! *** covariance matrix is not      ***
! *** computed explicitely.         ***
! *************************************

    doassim: IF (thisobs%doassim == 1) THEN

       ! Check process-local observation dimension

       IF (thisobs%dim_obs_p /= thisobs%dim_obs_f) THEN
          ! This error usually happens when localfilter=1
          WRITE (*,'(a)') 'PDAFomi - ERROR: PDAFomi_prodRinvA - INCONSISTENT value for DIM_OBS_P !!!'
          error = 3
       END IF

       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, &
               'PDAFomi_prodRinvA -- START Multiply with inverse observation variance'
          WRITE (*,*) '++ OMI-debug prodRinvA:         ', debug, 'observation ID', thisobs%obsid
          WRITE (*,*) '++ OMI-debug prodRinvA:         ', debug, 'thisobs%dim_obs_f', thisobs%dim_obs_f
          WRITE (*,*) '++ OMI-debug prodRinvA:         ', debug, 'thisobs%ivar_obs_f', thisobs%ivar_obs_f
          WRITE (*,*) '++ OMI-debug prodRinvA:         ', debug, 'Input matrix A_p', A_p
       END IF

       ! Initialize offset
       off = thisobs%off_obs_f

       DO j = 1, ncols
          DO i = 1, thisobs%dim_obs_f
             C_p(i+off, j) = thisobs%ivar_obs_f(i) * A_p(i+off, j)
          END DO
       END DO

       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_prodRinvA -- END'
       END IF

    END IF doassim

  END SUBROUTINE PDAFomi_prodRinvA



!-------------------------------------------------------------------------------
!> Compute likelihood for an ensemble member
!!
!! The routine is called during the analysis step
!! of the NETF or a particle filter.
!! It has to compute the likelihood of the
!! ensemble according to the difference from the
!! observation (residual) and the error distribution
!! of the observations.
!!
!! In general this routine is similar to the routine
!! prodRinvA used for ensemble square root Kalman
!! filters. As an addition to this routine, we here have
!! to evaluate the likelihood weight according the
!! assumed observation error statistics.
!!
!! __Revision history:__
!! * 2020-03 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_likelihood(thisobs, resid, lhood)

    USE PDAF_mod_filter, &
         ONLY: obs_member

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs   !< Data type with full observation
    REAL, INTENT(in)    :: resid(:)      ! Input vector of residuum
    REAL, INTENT(inout) :: lhood         ! Output vector - log likelihood

! *** local variables ***
    INTEGER :: i         ! index of observation component
    REAL, ALLOCATABLE :: Rinvresid(:) ! R^-1 times residual
    REAL :: lhood_one    ! Likelihood for this observation


    doassim: IF (thisobs%doassim == 1) THEN

! ****************************************
! *** First scale by observation error ***
! ***                   -1             ***
! ***      Rinvresid =  R  resid       ***
! ***                                  ***
! *** We assume a diagonal matrix R    ***
! ****************************************

       ! Screen output
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, &
               'PDAFomi_likelihood -- START Compute likelihood, member', obs_member
          IF (thisobs%obs_err_type==0) THEN
             WRITE (*,*) '++ OMI-debug likelihood:        ', debug, &
                  '  Assume Gaussian observation errors'
          ELSE
             WRITE (*,*) '++ OMI-debug likelihood:        ', debug, &
                  '  Assume double-exponential observation errors'
          END IF
       END IF

       ! Compute product of R^-1 with residuum
       ALLOCATE(Rinvresid(thisobs%dim_obs_f))

       DO i = 1, thisobs%dim_obs_f
          Rinvresid(i) = thisobs%ivar_obs_f(i) * resid(thisobs%off_obs_f+i)
       END DO


! ******************************
! *** Compute log likelihood ***
! ******************************

       IF (thisobs%obs_err_type==0) THEN

          ! Gaussian errors
          ! Calculate exp(-0.5*resid^T*R^-1*resid)

          ! Transform back to log likelihood to increment its values
          IF (lhood>0.0) lhood = - LOG(lhood)

          lhood_one = 0.0
          DO i = 1, thisobs%dim_obs_f
             lhood_one = lhood_one + 0.5*resid(thisobs%off_obs_f+i)*Rinvresid(i)
          END DO

          lhood = EXP(-(lhood + lhood_one))

       ELSE

          ! Double-exponential errors
          ! Calculate exp(-SUM(ABS(resid)))

          ! Transform pack to log likelihood to increment its values
          IF (lhood>0.0) lhood = - LOG(lhood)

          lhood_one = 0.0
          DO i = 1, thisobs%dim_obs_f
             lhood_one = lhood_one + ABS(Rinvresid(i))
          END DO

          lhood = EXP(-(lhood + lhood_one))

       END IF

       ! *** Clean up ***

       DEALLOCATE(Rinvresid)

       ! Screen output
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug likelihood:        ', debug, '  accumulated likelihood', lhood
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_likelihood -- END'
       END IF

    END IF doassim
    
  END SUBROUTINE PDAFomi_likelihood



!-------------------------------------------------------------------------------
!> Add observation error to some matrix
!!
!! The routine is called during the analysis step
!! of the stochastic EnKF. It it provided with a
!! matrix in observation space and has to add the 
!! observation error covariance matrix.
!!
!! This routine assumes a diagonal observation error
!! covariance matrix, but allows for varying observation
!! error variances.
!!
!! The routine can be applied with either all observations
!! of different types at once, or separately for each
!! observation type. The operation is done with all
!! process-local observations
!!
!! __Revision history:__
!! * 2020-03 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_add_obs_error(thisobs, nobs_all, matC)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(in) :: thisobs      !< Data type with full observation
    INTEGER, INTENT(in) :: nobs_all         !< Number of observations
    REAL, INTENT(inout) :: matC(:, :)       !< Input/Output matrix (nobs_f, rank)


! *** local variables ***
    INTEGER :: i, pe, cnt               ! Counters
    INTEGER :: idummy                   ! Dummy to access nobs_all
    INTEGER, ALLOCATABLE :: id_start(:) ! Start index of obs. type in global averall obs. vector
    INTEGER, ALLOCATABLE :: id_end(:)   ! End index of obs. type in global averall obs. vector


    doassim: IF (thisobs%doassim == 1) THEN

      ! Initialize dummy to prevent compiler warning
       idummy = nobs_all


! *************************************************
! *** Check process-local observation dimension ***
! *************************************************

       IF (thisobs%dim_obs_p /= thisobs%dim_obs_f) THEN
          ! This error usually happens when localfilter=1
          WRITE (*,'(a)') 'PDAFomi - ERROR: PDAFomi_add_obs_error - INCONSISTENT  VALUE for DIM_OBS_P !!!'
          error = 4
       END IF


! ********************************************************
! *** Initialize indices of observation type in global ***
! *** state vector over all observation types.         ***
! ********************************************************

       ALLOCATE(id_start(npes), id_end(npes))

       pe = 1
       id_start(1) = 1
       IF (thisobs%obsid>1) id_start(1) = id_start(1) + sum(obsdims(1, 1:thisobs%obsid-1))
       id_end(1)   = id_start(1) + obsdims(1,thisobs%obsid) - 1
       DO pe = 2, npes
          id_start(pe) = id_start(pe-1) + SUM(obsdims(pe-1,thisobs%obsid:))
          IF (thisobs%obsid>1) id_start(pe) = id_start(pe) + sum(obsdims(pe,1:thisobs%obsid-1))
          id_end(pe) = id_start(pe) + obsdims(pe,thisobs%obsid) - 1
       END DO


! *************************************
! ***   Add observation error       ***
! ***                               ***
! *** Measurements are uncorrelated ***
! *** here, thus R is diagonal      ***
! *************************************

       cnt = 1 
       DO pe = 1, npes
          DO i = id_start(pe), id_end(pe)
             matC(i, i) = matC(i, i) + 1.0/thisobs%ivar_obs_f(cnt)
             cnt = cnt + 1
          ENDDO
       ENDDO

       DEALLOCATE(id_start, id_end)

    END IF doassim

  END SUBROUTINE PDAFomi_add_obs_error



!-------------------------------------------------------------------------------
!> Initialize global observation error covariance matrix
!!
!! The routine is called during the analysis
!! step when an ensemble of observations is
!! generated by PDAF_enkf_obs_ensemble. 
!! It has to initialize the global observation 
!! error covariance matrix.
!!
!! This routine assumes a diagonal observation error
!! covariance matrix, but allows for varying observation
!! error variances.
!!
!! The routine can be applied with either all observations
!! of different types at once, or separately for each
!! observation type. The operation is done with all
!! process-local observations
!!
!! __Revision history:__
!! * 2020-03 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_init_obscovar(thisobs, nobs_all, covar, isdiag)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    INTEGER, INTENT(in) :: nobs_all        !< Number of observations
    REAL, INTENT(inout) :: covar(:, :)     !< Input/Output matrix (nobs_all, nobs_all)
                                           !< (needs to be set =0 before calling the routine)
    LOGICAL, INTENT(out) :: isdiag         !< Whether matrix R is diagonal

! *** local variables ***
    INTEGER :: i, pe, cnt               ! Counters
    INTEGER :: idummy                   ! Dummy to access nobs_all
    INTEGER, ALLOCATABLE :: id_start(:) ! Start index of obs. type in global averall obs. vector
    INTEGER, ALLOCATABLE :: id_end(:)   ! End index of obs. type in global averall obs. vector


    doassim: IF (thisobs%doassim == 1) THEN

       ! Initialize dummy to prevent compiler warning
       idummy = nobs_all


! *************************************************
! *** Initialize indices of observation type in ***
! *** global state vector over all obsservation *** 
! *** types and corresponding mapping vector.   ***
! *************************************************

       ALLOCATE(id_start(npes), id_end(npes))

       ! Initialize indices
       pe = 1
       id_start(1) = 1
       IF (thisobs%obsid>1) id_start(1) = id_start(1) + sum(obsdims(1, 1:thisobs%obsid-1))
       id_end(1)   = id_start(1) + obsdims(1,thisobs%obsid) - 1
       DO pe = 2, npes
          id_start(pe) = id_start(pe-1) + SUM(obsdims(pe-1,thisobs%obsid:))
          IF (thisobs%obsid>1) id_start(pe) = id_start(pe) + sum(obsdims(pe,1:thisobs%obsid-1))
          id_end(pe) = id_start(pe) + obsdims(pe,thisobs%obsid) - 1
       END DO

       ! Initialize mapping vector (to be used in PDAF_enkf_obs_ensemble)
       cnt = 1
       IF (thisobs%obsid-1 > 0) cnt = cnt+ SUM(obsdims(:,1:thisobs%obsid-1))
       DO pe = 1, npes
          DO i = id_start(pe), id_end(pe)
             map_obs_id(i) = cnt
             cnt = cnt + 1
          END DO
       END DO


! *************************************
! ***   Initialize covariances      ***
! ***                               ***
! *** Measurements are uncorrelated ***
! *** here, thus R is diagonal      ***
! *************************************

       cnt = 1 
       DO pe = 1, npes
          DO i = id_start(pe), id_end(pe)
             covar(i, i) = covar(i, i) + 1.0/thisobs%ivar_obs_f(cnt)
             cnt = cnt + 1
          ENDDO
       ENDDO

       ! The matrix is diagonal
       ! This setting avoids the computation of the SVD of COVAR
       ! in PDAF_enkf_obs_ensemble
       isdiag = .TRUE.

       DEALLOCATE(id_start, id_end)

    END IF doassim
    
  END SUBROUTINE PDAFomi_init_obscovar



!-------------------------------------------------------------------------------
!> Initialize vector of observation errors for generating synthetic obs.
!!
!! This routine initializes a vector of observation errors
!! used to perturb an observe dmodel state when generating 
!! synthetic observations.
!!
!! __Revision history:__
!! * 2020-05 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_init_obserr_f(thisobs, obserr_f)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(in) :: thisobs  !< Data type with full observation
    REAL, INTENT(inout) :: obserr_f(:)     !< Full vector of observation errors

! *** Local variables ***
    INTEGER :: i                           ! Counter


! *****************************************************************************
! *** Initialize vector of observation errors for generating synthetic obs. ***
! *****************************************************************************

    doassim: IF (thisobs%doassim == 1) THEN

       DO i = 1, thisobs%dim_obs_f
          IF (thisobs%ivar_obs_f(i)>0.0) THEN
             obserr_f(i + thisobs%off_obs_f) = 1.0/SQRT(thisobs%ivar_obs_f(i))
          ELSE    
             obserr_f(i + thisobs%off_obs_f) = 1.0e12
          END IF
       END DO

    END IF doassim

  END SUBROUTINE PDAFomi_init_obserr_f



!-------------------------------------------------------------------------------
!> Set min/max coordinate locations of a decomposed grid
!!
!! This routine sets the limiting coordinates of a 
!! process domain, i.e. the northern-, southern-,
!! eastern-, and western-most coordinate. The
!! information can be used to restrict the full
!! observations for PDAF to those that might be
!! used for the local analysis. 
!!
!! __Revision history:__
!! * 2020-03 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_set_domain_limits(lim_coords)

    IMPLICIT NONE

! *** Arguments ***
    REAL, INTENT(in) :: lim_coords(2,2)     !< geographic coordinate array (1: longitude, 2: latitude)
                                            !< ranges: longitude (-pi, pi), latitude (-pi/2, pi/2)

    ! Store domain limiting coordinates in module array
    IF (.NOT.ALLOCATED(domain_limits)) ALLOCATE(domain_limits(4))

    domain_limits(1) = lim_coords(2,1)  ! Northern edge
    domain_limits(2) = lim_coords(2,2)  ! Southern edge
    domain_limits(3) = lim_coords(1,1)  ! Western edge
    domain_limits(4) = lim_coords(1,2)  ! Eastern edge

    IF (debug>0) THEN
       WRITE (*,*) '++ OMI-debug: ', debug, &
            'PDAFomi_set_domain_limits -- START'
       WRITE (*,*) '++ OMI-debug set_domain_limits:      ', debug, 'domain limits', domain_limits
       WRITE (*,*) '++ OMI-debug: ', debug, &
            'PDAFomi_set_domain_limits -- END'
    END IF

  END SUBROUTINE PDAFomi_set_domain_limits



!-------------------------------------------------------------------------------
!> Find min/max coordinate locations in unstructured grid
!!
!! This routine finds the limiting coordinates of a 
!! process domain, i.e. the northern-, southern-,
!! eastern-, and western-most coordinate. The
!! information can be used to restrict the full
!! observations for PDAF to those that might be
!! used for the local analysis. Usually this is used
!! for unstructured grids, like in FESOM, because the
!! grid point indices do not contain information on 
!! coordinates in this case.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_get_domain_limits_unstr(npoints_p, coords_p)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: npoints_p        !< number of process-local grid points
    REAL, INTENT(in) :: coords_p(:,:)       !< geographic coordinate array (row 1: longitude, 2: latitude)
                                            !< ranges: longitude (-pi, pi), latitude (-pi/2, pi/2)

! *** Local variables ***
    INTEGER :: i                            ! Counter
    REAL :: nlimit, slimit, elimit, wlimit  ! Limiting coordinates
    REAL :: abslonmin                       ! absolute minimum longitude


! *** Determine limiting coordinates ***

    IF (debug>0) THEN
       WRITE (*,*) '++ OMI-debug: ', debug, &
            'PDAFomi_get_domain_limits_unstr -- START'
    END IF

    ! Initialize limiting values
    nlimit = -100.0
    slimit = 100.0
    wlimit = 100.0
    elimit = -100.0
    abslonmin = 100.0

    DO i=1, npoints_p
       ! Get North/South Limits
       IF (coords_p(2,i) < slimit) slimit = coords_p(2,i)
       IF (coords_p(2,i) > nlimit) nlimit = coords_p(2,i)

       ! Get East/West Limits
       IF (coords_p(1,i) < wlimit) wlimit = coords_p(1,i)
       IF (coords_p(1,i) > elimit) elimit = coords_p(1,i)
       IF (ABS(coords_p(1,i)) < abslonmin) THEN
          abslonmin = ABS(coords_p(1,i))
       END IF
    ENDDO

    IF (elimit*wlimit<0.0) THEN
       ! Domain crosses prime meridian or date line

       IF (wlimit<-3.1 .AND. elimit>3.1 .AND. abslonmin>0.5) THEN

          ! If the domain crosses the date line, we have to search the longitudinal limits differently
          elimit = -100.0
          wlimit = 100.0
          DO i=1, npoints_p
             IF (coords_p(1,i)<0.0 .AND. coords_p(1,i)>elimit) elimit = coords_p(1,i)
             IF (coords_p(1,i)>0.0 .AND. coords_p(1,i)<wlimit) wlimit = coords_p(1,i)
          END DO
          IF (debug>0) THEN
             WRITE (*,*) '++ OMI-debug get_domain_limits_unstr: ', debug, &
                  'limits (crossing date line) ', nlimit, slimit, wlimit, elimit
          END IF
       ELSE
          ! In this case the domain crosses the prime meridian
          IF (debug>0) THEN
             WRITE (*,*) '++ OMI-debug get_domain_limits_unstr: ', debug, &
                  'limits (cossing prime meridian) ', nlimit, slimit, wlimit, elimit
          END IF
       END IF
    ELSE
       ! Standard case
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug get_domain_limits_unstr: ', debug, &
               'limits (standard case) ', nlimit, slimit, wlimit, elimit
          END IF
    END IF

    ! Store domain limiting coordinates in module array

    IF (.NOT.ALLOCATED(domain_limits)) ALLOCATE(domain_limits(4))

    domain_limits(1) = nlimit
    domain_limits(2) = slimit
    domain_limits(3) = wlimit
    domain_limits(4) = elimit

    IF (debug>0) THEN
       WRITE (*,*) '++ OMI-debug: ', debug, &
            'PDAFomi_get_domain_limits_unstr -- END'
    END IF

  END SUBROUTINE PDAFomi_get_domain_limits_unstr


  
!-------------------------------------------------------------------------------
!> Find observations inside or close to process domain
!!
!! This routine finds observations that lie inside the 
!! local process sub-domain or within the distance
!! LRADIUS around it. The observations are counted and
!! an index array is initialized storing the indices
!! of the process-local relevant full observations in the
!! global full observation vector.
!!
!! The routine has to be called by all filter processes
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_get_local_ids_obs_f(dim_obs_g, lradius, oc_f, cnt_lim, id_lim, disttype, domainsize)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_obs_g       !< Global full number of observations
    REAL, INTENT(in) :: lradius            !< Localization radius (used is a constant one here)
    REAL, INTENT(in) :: oc_f(:,:)          !< observation coordinates (radians), row 1: lon, 2: lat
                                           !< ranges: longitude (-pi, pi), latitude (-pi/2, pi/2)
    INTEGER, INTENT(out) :: cnt_lim        !< Number of full observation for local process domain
    INTEGER, INTENT(out) :: id_lim(:)      !< Indices of process-local full obs. in global full vector
                                           !< It has to be allocated sufficiently large
    INTEGER, INTENT(in) :: disttype        !< type of distance computation
    REAL, INTENT(in), OPTIONAL :: domainsize(:)   !< Global size of model domain

! *** Local variables ***
    INTEGER :: i         ! Counter
    INTEGER :: flag      ! Counting flag
    REAL :: limdist      ! Limit distance normalized by r_earth
    INTEGER :: cnt_lim_max, cnt_lim_min  ! min/max number over all domains
    REAL :: maxlat       ! Highest latitude of a domain
    REAL :: period(2)    ! Size of domain in case of periodicity


! **********************
! *** Initialization ***
! **********************

    IF (debug>0) THEN
       WRITE (*,*) '++ OMI-debug: ', debug, &
            'PDAFomi_get_local_ids_obs_f -- START Limit observations to process sub-domains'
       WRITE (*,*) '++ OMI-debug get_local_ids_obs_f: ', debug, 'domain limiting coordinates', domain_limits
    END IF

    IF (.NOT.ALLOCATED(domain_limits)) THEN
       WRITE (*,'(a)') 'PDAFomi - ERROR: PDAFomi_get_local_ids_obs_f - DOMAIN_LIMITS is not initialized !!!'
       error = 5
    END IF

    IF ((disttype==1 .OR. disttype==11) .AND. .NOT.PRESENT(domainsize)) THEN
       WRITE (*,'(a)') 'PDAFomi - ERROR: PDAFomi_get_local_ids_obs_f - THISOBS%DOMAINSIZE is not initialized !!!'
       error = 6
    END IF

    ! initialize index array
    id_lim = 0


! ***************************************
! *** Find relevant full observations ***
! ***************************************

    cnt_lim = 0

    dtype: IF (disttype==2 .OR. disttype==3 .OR. disttype==12 .OR. disttype==13) THEN

       ! Limit distance around the domain
       limdist = lradius / r_earth

       IF (debug==1) THEN
          WRITE (*,*) '++ OMI-debug get_local_ids_obs_f: ', debug, 'limit for geographic coordinates'
          WRITE (*,*) '++ OMI-debug get_local_ids_obs_f: ', debug, 'limiting distance (m)', limdist
       END IF

       fullobsloop: DO i = 1, dim_obs_g

          ! Init flag for latitudinal check
          flag = 0

          ! First check in latitudinal direction
          checklat: IF (oc_f(2,i)<=domain_limits(1) .AND. oc_f(2,i)>=domain_limits(2)) THEN
             ! inside domain north-south extent
             flag=1
          ELSEIF (oc_f(2,i)>domain_limits(1)) THEN
             ! north of the domain
             IF (ABS(oc_f(2,i)-domain_limits(1)) <= limdist) flag=1
          ELSEIF (oc_f(2,i)<domain_limits(2)) THEN
             ! south of the domain
             IF (ABS(oc_f(2,i)-domain_limits(2)) <= limdist) flag=1
          END IF checklat

          ! Store highest latitude
          maxlat = MAX(ABS(domain_limits(1)), ABS(domain_limits(2)))

          ! if observation fits in the latitudinal direction check longitudinal direction
          lat_ok: IF (flag==1) THEN
             lontypes: IF (domain_limits(4)>=0.0 .OR. (domain_limits(4)<0.0 .AND. domain_limits(3)<0.0)) THEN

                IF (oc_f(1,i)>=domain_limits(3) .AND. oc_f(1,i)<=domain_limits(4)) THEN

                   ! fully inside domain extent
                   cnt_lim = cnt_lim+1
                   id_lim(cnt_lim) = i
                ELSEIF (oc_f(1,i)<domain_limits(3)) THEN

                   ! west of the domain
                   IF (ABS(COS(maxlat)*(oc_f(1,i)-domain_limits(3))) <= limdist .OR. &
                        (ABS(COS(maxlat)*(oc_f(1,i)-domain_limits(3)))-2.0*pi) <= limdist ) THEN
                      cnt_lim = cnt_lim+1
                      id_lim(cnt_lim) = i
                   END IF
                ELSEIF (oc_f(1,i)>domain_limits(4)) THEN

                   ! east of the domain
                   IF (ABS(COS(maxlat)*(oc_f(1,i)-domain_limits(4))) <= limdist .OR. &
                        (ABS(COS(maxlat)*(oc_f(1,i)-domain_limits(4)))-2.0*pi) <= limdist ) THEN
                      cnt_lim = cnt_lim+1
                      id_lim(cnt_lim) = i
                   END IF
                ENDIF
             ELSE lontypes
                IF ((oc_f(1,i)>=domain_limits(3) .AND. oc_f(1,i)<=pi) .OR. &
                     (oc_f(1,i)<=domain_limits(4)) .AND. oc_f(1,i)>=-pi) THEN

                   ! fully inside domain extent
                   cnt_lim = cnt_lim+1
                   id_lim(cnt_lim) = i

                ELSEIF (oc_f(1,i)<domain_limits(3) .AND. oc_f(1,i)>=0.0) THEN

                   ! east of the domain
                   IF (ABS(COS(maxlat)*(oc_f(1,i)-domain_limits(3))) <= limdist) THEN
                      cnt_lim = cnt_lim+1
                      id_lim(cnt_lim) = i
                   ENDIF
                ELSEIF (oc_f(1,i)>domain_limits(4)) THEN

                   ! west of the domain
                   IF (ABS(COS(maxlat)*(oc_f(1,i)-domain_limits(4))) <= limdist) THEN
                      cnt_lim = cnt_lim+1
                      id_lim(cnt_lim) = i
                   ENDIF
                ENDIF
             ENDIF lontypes
          ENDIF lat_ok
       END DO fullobsloop

    ELSE IF (disttype==0 .OR. disttype==10) THEN

       ! *** Check Cartesian coordinates without periodicity ***

       limdist = lradius

       IF (debug==1) THEN
          WRITE (*,*) '++ OMI-debug get_local_ids_obs_f: ', debug, 'limit for Cartesian coordinates'
          WRITE (*,*) '++ OMI-debug get_local_ids_obs_f: ', debug, 'limiting distance', limdist
       END IF

       fullobsloopB: DO i = 1, dim_obs_g

          ! Init flag for latitudinal check
          flag = 0

          ! First check in latitudinal direction
          checklatC: IF (oc_f(2,i)<=domain_limits(1) .AND. oc_f(2,i)>=domain_limits(2)) THEN
             ! inside domain north-south extent
             flag=1
          ELSEIF (oc_f(2,i)>domain_limits(1)) THEN
             ! north of the domain
             IF (ABS(oc_f(2,i)-domain_limits(1)) <= limdist) flag=1
          ELSEIF (oc_f(2,i)<domain_limits(2)) THEN
             ! south of the domain
             IF (ABS(oc_f(2,i)-domain_limits(2)) <= limdist) flag=1
          END IF checklatC

          ! if observation fits in the latitudinal direction check longitudinal direction
          lat_okB: IF (flag==1) THEN

             IF (oc_f(1,i)>=domain_limits(3) .AND. oc_f(1,i)<=domain_limits(4)) THEN

                ! fully inside domain extent
                cnt_lim = cnt_lim+1
                id_lim(cnt_lim) = i
             ELSEIF (oc_f(1,i)<domain_limits(3)) THEN

                ! West of the domain
                IF (ABS(oc_f(1,i)-domain_limits(3)) <= limdist) THEN
                   cnt_lim = cnt_lim+1
                   id_lim(cnt_lim) = i
                END IF
             ELSEIF (oc_f(1,i)>domain_limits(4)) THEN

                ! East of the domain
                IF (ABS(oc_f(1,i)-domain_limits(4)) <= limdist) THEN
                   cnt_lim = cnt_lim+1
                   id_lim(cnt_lim) = i
                END IF
             ENDIF
          END IF lat_okB
       END DO fullobsloopB
    ELSE IF (disttype==1 .OR. disttype==11) THEN

       ! *** Check Cartesian coordinates with periodicity ***

       limdist = lradius

       IF (debug==1) THEN
          WRITE (*,*) '++ OMI-debug get_local_ids_obs_f: ', debug, 'limit for periodic Cartesian coordinates'
          WRITE (*,*) '++ OMI-debug get_local_ids_obs_f: ', debug, 'limiting distance', limdist
          WRITE (*,*) '++ OMI-debug get_local_ids_obs_f: ', debug, 'thisobs%domainsize', domainsize
       END IF

       fullobsloopC: DO i = 1, dim_obs_g

          ! Init flag for latitudinal check
          flag = 0
          
          IF (domainsize(1)>0.0) THEN 
             period(1) = domainsize(1)
          ELSE
             period(1) = 0.0
          END IF
          IF (domainsize(2)>0.0) THEN 
             period(2) = domainsize(2)
          ELSE
             period(2) = 0.0
          END IF

          ! First check in latitudinal direction
          checklatB: IF (oc_f(2,i)<=domain_limits(1) .AND. oc_f(2,i)>=domain_limits(2)) THEN
             ! inside domain north-south extent
             flag=1
          ELSEIF (oc_f(2,i)>domain_limits(1)) THEN
             ! north of the domain
             IF ((ABS(oc_f(2,i)-domain_limits(1)) <= limdist) .OR. &
                  (ABS(oc_f(2,i)-domain_limits(1) - period(2)) <= limdist)) flag=1
          ELSEIF (oc_f(2,i)<domain_limits(2)) THEN
          ! south of the domain
             IF ((ABS(oc_f(2,i)-domain_limits(2)) <= limdist) .OR. &
                 (ABS(oc_f(2,i)-domain_limits(2) + period(2)) <= limdist)) flag=1
          END IF checklatB

          ! if observation fits in the latitudinal direction check longitudinal direction
          lat_okC: IF (flag==1) THEN

             IF (oc_f(1,i)>=domain_limits(3) .AND. oc_f(1,i)<=domain_limits(4)) THEN

                ! fully inside domain extent
                cnt_lim = cnt_lim+1
                id_lim(cnt_lim) = i
             ELSEIF (oc_f(1,i)<domain_limits(3)) THEN

                ! West of the domain
                IF (ABS(oc_f(1,i)-domain_limits(3)) <= limdist .OR. &
                     (ABS(oc_f(1,i)-domain_limits(3)) - period(1)) <= limdist) THEN
                   cnt_lim = cnt_lim+1
                   id_lim(cnt_lim) = i
                END IF
             ELSEIF (oc_f(1,i)>domain_limits(4)) THEN

                ! East of the domain
                IF (ABS(oc_f(1,i)-domain_limits(4)) <= limdist .OR. &
                     (ABS(oc_f(1,i)-domain_limits(4)) - period(1)) <= limdist) THEN
                   cnt_lim = cnt_lim+1
                   id_lim(cnt_lim) = i
                END IF
             ENDIF
          END IF lat_okC
       END DO fullobsloopC

    END IF dtype

    IF (debug>0) THEN
       WRITE (*,*) '++ OMI-debug get_local_ids_obs_f: ', debug, 'obs. ids for process domains', id_lim(1:cnt_lim)
       WRITE (*,*) '++ OMI-debug: ', debug, &
            'PDAFomi_get_local_ids_obs_f -- END'
    END IF

    ! Get number of min/max process-local full observation dimensions
    CALL MPI_Allreduce (cnt_lim, cnt_lim_max, 1, MPI_INTEGER, MPI_MAX, &
         COMM_filter, MPIerr)
    CALL MPI_Allreduce (cnt_lim, cnt_lim_min, 1, MPI_INTEGER, MPI_MIN, &
         COMM_filter, MPIerr)
  
    IF (mype == 0 .AND. screen > 0) THEN
       WRITE (*,'(a,8x,a,i8)') 'PDAFomi','--- global observation dimension', dim_obs_g
       WRITE (*,'(a,8x,a,i7,1x,i7)') 'PDAFomi','--- process-local min/max full obs. dimensions', &
            cnt_lim_min, cnt_lim_max
    END IF

  END SUBROUTINE PDAFomi_get_local_ids_obs_f


  
!-------------------------------------------------------------------------------
!> Reduce full observation vector to part relevant for local process domain
!!
!! This routine initializes a full vector of observations that only
!! contains those full observations that are relevant for a process
!! subdomain. The indices of these observations were determined
!! using get_local_ids_obs_f.
!!
!! __Revision history:__
!! * 2019-07 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_limit_obs_f(thisobs, offset, obs_f_one, obs_f_lim)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    REAL, INTENT(in) :: obs_f_one(:)       !< Global full observation vector (nobs_f)
    REAL, INTENT(out) :: obs_f_lim(:)      !< full observation vector for process domains (nobs_lim)
    INTEGER, INTENT(in) :: offset          !< offset of this observation in obs_f_lim

! *** Local variables ***
    INTEGER :: i         ! Counter


! ********************************************
! *** Initialize process-local full vector ***
! ********************************************

    IF (.NOT.ALLOCATED(thisobs%id_obs_f_lim)) THEN
       WRITE (*,'(a)') 'PDAFomi - ERROR: PDAFomi_limit_obs_f - thisobs%id_obs_f_lim is not allocated !!!'
       error = 7
    END IF

    DO i = 1, thisobs%dim_obs_f
       obs_f_lim(i+offset) = obs_f_one(thisobs%id_obs_f_lim(i))
    END DO

  END SUBROUTINE PDAFomi_limit_obs_f

SUBROUTINE PDAFomi_gather_obs_f_flex(dim_obs_p, obs_p, obs_f, status)

! !DESCRIPTION:
! If the local filter is used with a domain-decomposed model,
! the observational information from different sub-domains
! has to be combined into the full observation vector. 
! In this routine the process-local parts of the observation
! vector are gathered into a full observation vector. 
! The routine requires that PDAF_gather_dim_obs_f was executed
! before, because this routine initializes dimensions that are 
! used here. 
! The routine can also be used to gather full arrays of coordinates.
! It is however, only usable if the coordinates are stored row-
! wise, i.e. each row represents the set of coordinates for one
! observation point. It has to be called separately for each column. 
! A  better alternative is the row-wise storage of coordinates. In this
! case the routine PDAF_gather_dim_obs_f allows the gather the full
! coordinate array in one step.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2019-03 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_mod_filtermpi, &
       ONLY: COMM_filter, MPIerr, mype_filter, npes_filter

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_obs_p    ! PE-local observation dimension
  REAL, INTENT(in)  :: obs_p(:)  ! PE-local vector
  REAL, INTENT(out) :: obs_f(:)  ! Full gathered vector
  INTEGER, INTENT(out) :: status   ! Status flag: (0) no error

! !CALLING SEQUENCE:
! Called by: user code
! Calls: MPI_Allreduce
! Calls: MPI_Allgather
! Calls: MPI_AllgatherV
!EOP

! Local variables
  INTEGER :: i                              ! Counter
  INTEGER :: dimobs_f                       ! full dimension of observation vector obtained from allreduce
  INTEGER, ALLOCATABLE :: all_dim_obs_p(:)  ! PE-Local observation dimensions
  INTEGER, ALLOCATABLE :: all_dis_obs_p(:)  ! PE-Local observation displacements


! **********************************************************
! *** Compute global sum of local observation dimensions ***
! **********************************************************

  IF (npes_filter>1) THEN
     CALL MPI_Allreduce(dim_obs_p, dimobs_f, 1, MPI_INTEGER, MPI_SUM, &
          COMM_filter, MPIerr)
  ELSE
     dimobs_f = dim_obs_p
  END IF


! ****************************************************************************
! *** Gather and store array of process-local dimensions and displacements ***
! ****************************************************************************

  ALLOCATE(all_dim_obs_p(npes_filter))
  ALLOCATE(all_dis_obs_p(npes_filter))

  IF (npes_filter>1) THEN
     CALL MPI_Allgather(dim_obs_p, 1, MPI_INTEGER, all_dim_obs_p, 1, &
          MPI_INTEGER, COMM_filter, MPIerr)

     ! Init array of displacements for observation vector
     all_dis_obs_p(1) = 0
     DO i = 2, npes_filter
        all_dis_obs_p(i) = all_dis_obs_p(i-1) + all_dim_obs_p(i-1)
     END DO
  ELSE
     all_dim_obs_p = dim_obs_p
     all_dis_obs_p = 0
  END IF


! **********************************************************
! *** Gather full observation vector                     ***
! **********************************************************

  IF (npes_filter>1) THEN
     CALL MPI_AllGatherV(obs_p, all_dim_obs_p(mype_filter+1), MPI_REALTYPE, &
          obs_f, all_dim_obs_p, all_dis_obs_p, MPI_REALTYPE, &
          COMM_filter, MPIerr)
  
     status = MPIerr
  ELSE
     obs_f = obs_p

     status = 0
  END IF


! ****************
! *** Clean up ***
! ****************

  DEALLOCATE(all_dim_obs_p, all_dis_obs_p)

END SUBROUTINE PDAFomi_gather_obs_f_flex

SUBROUTINE PDAFomi_gather_obs_f2_flex(dim_obs_p, coords_p, coords_f, &
     nrows, status)

! !DESCRIPTION:
! If the local filter is used with a domain-decomposed model,
! the observational information from different sub-domains
! has to be combined into the full observation vector. 
! In this routine the process-local parts of a coordinate array
! accompanying the observation vector are gathered into a full
! array of coordinates. 
! The routine is for the case that the observation coordinates
! are stored column-wise, i.e. each column is the set of coordinates
! for one observation. This should be the usual case, as in this
! case the set of coordinates of one observations are stored
! next to each other in memory. If the coordinates are stored row-
! wise, the routine PDAF_gather_obs_f can be used, but has to be
! called separately for each column. 
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2019-03 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_mod_filtermpi, &
       ONLY: COMM_filter, MPIerr, mype_filter, npes_filter

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_obs_p    ! PE-local observation dimension
  INTEGER, INTENT(in) :: nrows        ! Number of rows in array
  REAL, INTENT(in)  :: coords_p(:,:)  ! PE-local array
  REAL, INTENT(out) :: coords_f(:,:)  ! Full gathered array
  INTEGER, INTENT(out) :: status   ! Status flag: (0) no error

! !CALLING SEQUENCE:
! Called by: user code
! Calls: MPI_Allreduce
! Calls: MPI_Allgather
! Calls: MPI_AllgatherV
!EOP

! local variables
  INTEGER :: i                              ! Counter
  INTEGER :: dimobs_f                       ! full dimension of observation vector obtained from allreduce
  INTEGER, ALLOCATABLE :: all_dim_obs_p(:)  ! PE-Local observation dimensions
  INTEGER, ALLOCATABLE :: all_dim_obs_p2(:) ! local-dims for multi-row array
  INTEGER, ALLOCATABLE :: all_dis_obs_p2(:) ! displacements to gather multi-row array


! **********************************************************
! *** Compute global sum of local observation dimensions ***
! **********************************************************

  IF (npes_filter>1) THEN
     CALL MPI_Allreduce(dim_obs_p, dimobs_f, 1, MPI_INTEGER, MPI_SUM, &
          COMM_filter, MPIerr)
  ELSE
     dimobs_f = dim_obs_p
  END IF


! ****************************************************************************
! *** Gather and store array of process-local dimensions and displacements ***
! ****************************************************************************

  ALLOCATE(all_dim_obs_p(npes_filter))

  IF (npes_filter>1) THEN
     CALL MPI_Allgather(dim_obs_p, 1, MPI_INTEGER, all_dim_obs_p, 1, &
          MPI_INTEGER, COMM_filter, MPIerr)
  ELSE
     all_dim_obs_p = dim_obs_p
  END IF


! **********************************************************
! *** Gather full observation coordinates array          ***
! **********************************************************

  IF (npes_filter>1) THEN
     ALLOCATE(all_dis_obs_p2(npes_filter))
     ALLOCATE(all_dim_obs_p2(npes_filter))

     ! Init array of local dimensions
     do i = 1, npes_filter
        all_dim_obs_p2(i) = nrows * all_dim_obs_p(i)
     end do

     ! Init array of displacements for observation vector
     all_dis_obs_p2(1) = 0
     DO i = 2, npes_filter
        all_dis_obs_p2(i) = all_dis_obs_p2(i-1) + all_dim_obs_p2(i-1)
     END DO

     CALL MPI_AllGatherV(coords_p, all_dim_obs_p2(mype_filter+1), MPI_REALTYPE, &
          coords_f, all_dim_obs_p2, all_dis_obs_p2, MPI_REALTYPE, &
          COMM_filter, MPIerr)

     DEALLOCATE(all_dim_obs_p2, all_dis_obs_p2)

     status = MPIerr
  ELSE
     coords_f = coords_p
     
     status = 0
  END IF

END SUBROUTINE PDAFomi_gather_obs_f2_flex



!-------------------------------------------------------------------------------
!> Exclude observations for too high innovation
!!
!! The routine is called during the analysis step
!! of a global filter. It checks the size of the 
!! innovation and sets the observation error to a
!! high value if the squared innovation exceeds a
!! limit relative to the observation error variance.
!!
!! __Revision history:__
!! * 2024-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_omit_by_inno(thisobs, inno_f, obs_f_all, obsid, cnt_all)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    REAL, INTENT(in)    :: inno_f(:)         !< Input vector of observation innovation
    REAL, INTENT(in)    :: obs_f_all(:)      !< Input vector of local observations
    INTEGER, INTENT(in) :: obsid             !< ID of observation type
    INTEGER, INTENT(inout) :: cnt_all        !< Count of omitted observation over all types


! *** local variables ***
    INTEGER :: i                    ! Index of observation component
    INTEGER :: cnt                  ! Counter
    REAL :: inno2                   ! Squared innovation
    REAL :: limit2                  ! Squared limit


! **********************
! *** INITIALIZATION ***
! **********************

    doassim: IF (thisobs%doassim == 1) THEN

       IF (thisobs%inno_omit > 0.0) THEN

          IF (debug>0) THEN
             WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_omit_by_inno -- START   obs-ID', obsid
             WRITE (*,*) '++ OMI-debug omit_by_inno:', debug, 'limit for innovation', &
                  thisobs%inno_omit
             WRITE (*,*) '++ OMI-debug omit_by_inno:', debug, 'inno_omit_ivar', &
                  thisobs%inno_omit_ivar
             WRITE (*,*) '++ OMI-debug omit_by_inno:', debug, 'innovation_f(1:10)', inno_f(1:10)
          ENDIF

          ! Squared limit factor
          limit2 = thisobs%inno_omit * thisobs%inno_omit

          ! Check for observations to be excluded
          cnt = 0
          DO i = 1, thisobs%dim_obs_f

             ! Squared innovation
             inno2 = inno_f(i + thisobs%off_obs_f)* inno_f(i + thisobs%off_obs_f)

             IF (inno2 > limit2 * 1.0/thisobs%ivar_obs_f(i)) THEN
                IF (debug>0) THEN
                   WRITE (*,*) '++ OMI-debug omit_by_inno:', debug, 'omit: innovation:', &
                        inno_f(i + thisobs%off_obs_f), 'observation:', obs_f_all(i + thisobs%off_obs_f)
                END IF

                ! Exclude observation by increased its observation error
                thisobs%ivar_obs_f(i) = thisobs%inno_omit_ivar

                ! Count excluded obs
                cnt = cnt + 1
             END IF
          ENDDO

          IF (debug>0 .and. cnt>0) THEN
             WRITE (*,*) '++ OMI-debug omit_by_inno:', debug, 'count of excluded obs.: ', cnt
             WRITE (*,*) '++ OMI-debug omit_by_inno:', debug, 'updated thisobs_f%ivar_obs_f ', &
                  thisobs%ivar_obs_f
          ENDIF

          IF (debug>0) &
               WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_omit_by_inno_f -- END   obs-ID', obsid

          cnt_all = cnt_all + cnt

       END IF

    ENDIF doassim

  END SUBROUTINE PDAFomi_omit_by_inno



!-------------------------------------------------------------------------------
!> Get statistics on local observations
!!
!! The routine is called in the update routine of
!! global filters and writes statistics on 
!! used and excluded observations.
!!
!! __Revision history:__
!! * 2023-03 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_obsstats(screen)

    USE MPI
    USE PDAF_mod_filtermpi, &
         ONLY: COMM_filter, MPIerr, npes_filter

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: screen            !< Verbosity flag

! *** Local variables ***
    INTEGER :: ostats_omit_g(7)

    IF (npes_filter>1) THEN
       CALL MPI_Reduce(ostats_omit, ostats_omit_g, 7, MPI_INTEGER, MPI_SUM, &
            0, COMM_filter, MPIerr)
    ELSE
       ! This is a work around for working with nullmpi.F90
       ostats_omit_g = ostats_omit
    END IF

    IF (mype == 0 .AND. screen > 0 .AND. ostats_omit_g(1)>0) THEN
       WRITE (*, '(a, 9x, a)') 'PDAFomi', 'Global statistics for omitted observations:'
       WRITE (*, '(a, 11x, a, i10)') &
            'PDAFomi', 'Global number of omitted observations: ', ostats_omit_g(6)
       WRITE (*, '(a, 11x, a, i10)') &
            'PDAFomi', 'Global number of used observations:    ', ostats_omit_g(7)
    ELSEIF (mype == 0 .AND. screen > 0) THEN
       WRITE (*, '(a, 9x, a)') 'PDAFomi', 'Global statistics for omitted observations:'
       WRITE (*, '(a, 11x, a)') &
            'PDAFomi', 'Zero observations omitted'
    END IF

  END SUBROUTINE PDAFomi_obsstats


!-------------------------------------------------------------------------------
!> Gather global observation dimension information
!!
!! This routine gathers the information about
!! the full dimension of each observation type
!! in each process-local subdomain.
!!
  SUBROUTINE PDAFomi_gather_obsdims()

    IMPLICIT NONE

! *** local variables ***
    INTEGER :: i                ! Loop counter
    INTEGER :: dim_obs_all      ! Full number of global observations


! *****************************************
! *** Gather all observation dimensions ***
! *****************************************

    ALLOCATE(obsdims(npes,n_obstypes))

    DO i = 1, n_obstypes
       CALL MPI_Allgather(obs_f_all(i)%ptr%dim_obs_f, 1, MPI_INTEGER, obsdims(:,i), 1, &
          MPI_INTEGER, COMM_filter, MPIerr)
    END DO

    ! Determine overall number of observations
    dim_obs_all = SUM(obsdims)

    ! Allocate mapping vector
    ALLOCATE(map_obs_id(dim_obs_all))

  END SUBROUTINE PDAFomi_gather_obsdims



!-------------------------------------------------------------------------------
!> Check error flag
!!
!! This routine returns the value of the PDAF-OMI internal error flag.
!!
!! __Revision history:__
!! * 2024-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_check_error(flag)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: flag            !< Error flag

    ! Set error flag
    flag = error

  END SUBROUTINE PDAFomi_check_error

END MODULE PDAFomi_obs_f
