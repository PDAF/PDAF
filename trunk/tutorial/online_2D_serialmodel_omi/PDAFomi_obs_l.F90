! Copyright (c) 2004-2019 Lars Nerger
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
!$Id$

!> PDAF-OMI routines for local observations
!!
!! This module contains generic routines for several observation-related
!! operations for local filters. The routines are
!!
!! * init_dim_obs_l \n
!!        Initialize dimension of local obs. vetor and arrays for
!!        local observations
!! * init_obs_l \n
!!        Initialize local observation vector and inverse 
!!        error variances
!! * cnt_dim_obs_l \n
!!        Set dimension of local obs. vector
!! * init_obsarrays_l \n
!!        Initialize arrays for the index of a local observation in 
!!        the full observation vector and its corresponding distance.
!! * init_obs_l \n
!!        Initialize the local vector of observations
!! * g2l_obs \n
!!        Initialize local observation vector from full observation vector
!! * prodRinvA_l \n
!!        Multiply an intermediate matrix of the local filter analysis
!!        with the inverse of the observation error covariance matrix
!!        and apply observation localization
!! * init_obsvar_l \n
!!        Compute mean observation error variance
!! *set_debug_flag \n
!!        Set or unset the debugging flag for PDAFomi routines
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE PDAFomi_obs_l

  IMPLICIT NONE
  SAVE

! *** Module internal variables
  REAL, PARAMETER :: r_earth=6.3675e6     !< Earth radius in meters
  INTEGER :: debug=0                      !< Debugging flag

  ! Data type to define the local observations by internally shared variables of the module
  type obs_l
     INTEGER :: dim_obs_l                 !< number of local observations
     INTEGER :: off_obs_l                 !< Offset of this observation in overall local obs. vector
     INTEGER, ALLOCATABLE :: id_obs_l(:)  !< Indices of local observations in full obs. vector 
     REAL, ALLOCATABLE :: distance_l(:)   !< Distances of local observations
     REAL, ALLOCATABLE :: ivar_obs_l(:)   !< Inverse variance of local observations
  end type obs_l

!$OMP THREADPRIVATE(debug)


!-------------------------------------------------------------------------------

CONTAINS


!> Set dimension of local obs. vector and local obs. arrays
!!
!! This routine sets the number of local observations for the
!! current observation type for the local analysis domain
!! with coordinates COORD_l and localization radius LRADIUS.
!! Further the routine initializes arrays for the index of a
!! local observation in the full observation vector and its 
!! corresponding distance.
!! The operation are performed by calling the routines 
!! cnt_dim_obs_l and init_obsarrays_l.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE init_dim_obs_l(thisobs, thisobs_l, coord_l, lradius, nobs_l_one, &
       off_obs_l_all, off_obs_f_all)

    USE PDAFomi_obs_f, &
         ONLY: obs_f

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    REAL, INTENT(in) :: coord_l(:)           !< Coordinates of current analysis domain
    REAL, INTENT(in) :: lradius              !< Localization radius in meters
    INTEGER, INTENT(out) :: nobs_l_one       !< Local dimension of current observation vector
    INTEGER, INTENT(inout) :: off_obs_l_all  !< input: offset of current obs. in local obs. vector
                                             !< output: input + nobs_l_one
    INTEGER, INTENT(inout) :: off_obs_f_all  !< input: offset of current obs. in full obs. vector
                                             !< output: input + nobs_f_one


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    CALL cnt_dim_obs_l(thisobs%disttype, thisobs%ncoord, coord_l, lradius, &
         thisobs%dim_obs_f, thisobs%ocoord_f, nobs_l_one)
    
    ! Store number of local module-type observations in module
    thisobs_l%dim_obs_l = nobs_l_one


! ************************************************************
! *** Initialize internal local arrays for local distances ***
! *** and indices of local obs. in full obs. vector        ***
! ************************************************************

    ! Allocate module-internal index array for indices in module-type observation vector
    IF (ALLOCATED(thisobs_l%id_obs_l)) DEALLOCATE(thisobs_l%id_obs_l)
    IF (ALLOCATED(thisobs_l%distance_l)) DEALLOCATE(thisobs_l%distance_l)
    IF (thisobs_l%dim_obs_l>0) THEN
       ALLOCATE(thisobs_l%id_obs_l(thisobs_l%dim_obs_l))
       ALLOCATE(thisobs_l%distance_l(thisobs_l%dim_obs_l))
    ELSE
       ALLOCATE(thisobs_l%id_obs_l(1))
       ALLOCATE(thisobs_l%distance_l(1))
    END IF

    ! Store offsets
    thisobs_l%off_obs_l = off_obs_l_all

    ! Initialize ID_OBS_L and DISTANCE_L and increment offsets
    CALL init_obsarrays_l(thisobs%disttype, thisobs%ncoord, coord_l, lradius, &
         thisobs_l%dim_obs_l, thisobs%dim_obs_f, thisobs%ocoord_f, &
         thisobs_l%distance_l, thisobs_l%id_obs_l, off_obs_l_all, off_obs_f_all)

  END SUBROUTINE init_dim_obs_l



!-------------------------------------------------------------------------------
!> Set dimension of local observation vector
!!
!! This routine sets the number of local observations for the
!! current observation type for the local analysis domain
!! with coordinates COORD_L and localization radius LRADIUS.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE cnt_dim_obs_l(disttype, ncoord, coord_l, lradius, nobs_f_one, ocoord_f_one, nobs_l_one)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: disttype          !< Type of distance computation
                                             !< 0: Cartesian
                                             !< 1: geographic/spherical 
    INTEGER, INTENT(in) :: ncoord            !< Number of coordinates to consider
    REAL, INTENT(in) :: coord_l(ncoord)      !< Coordinates of current analysis domain
    REAL, INTENT(in) :: lradius              !< Localization radius in meters
    INTEGER, INTENT(in) :: nobs_f_one        !< Full dimension of current observation vector
    REAL, INTENT(in) :: ocoord_f_one(:, :)   !< coordinate array for current observations
    INTEGER, INTENT(out) :: nobs_l_one       !< Local dimension of current observation vector

! *** Local variables ***
    INTEGER :: i, k         ! Counters
    REAL :: ocoord(ncoord)  ! Coordinates of observation
    REAL :: dists(ncoord)   ! Distance vector between analysis point and observation
    REAL :: lradius2        ! squared localization radius
    REAL :: distance2       ! squared distance


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    ! Initialize squared localization radius
    lradius2 = lradius*lradius

    ! Count local observations
    nobs_l_one = 0

    norm: IF (disttype==0) THEN

       ! *** Count with Cartesian distance ***

       scancount_Cart: DO i = 1, nobs_f_one

          ! location of observation point
          ocoord(1:ncoord) = ocoord_f_one(1:ncoord, i)

          ! approximate distances in longitude and latitude
          DO k = 1, ncoord
             dists(k) = ABS(coord_l(k) - ocoord(k))
          END DO

          ! full squared distance
          distance2 = 0.0
          DO k = 1, ncoord
             distance2 = distance2 + dists(k)*dists(k)
          END DO

          ! If distance below limit, add observation to local domain
          IF (distance2 <= lradius2) THEN
             nobs_l_one = nobs_l_one + 1
          END IF
       END DO scancount_Cart

    ELSEIF (disttype==1) THEN norm

       ! *** Count with distance from geographic coordinates ***

       scancount_Geo: DO i = 1, nobs_f_one

          ! location of observation point
          ocoord(1:ncoord) = ocoord_f_one(1:ncoord, i)

          ! approximate distances in longitude and latitude
          dists(1) = r_earth * ABS(coord_l(1) - ocoord(1))* COS(coord_l(2))
          dists(2) = r_earth * ABS(coord_l(2) - ocoord(2))
          IF (ncoord>2) dists(3) = ABS(coord_l(3) - ocoord(3))

          ! full squared distance in meters
          distance2 = 0.0
          DO k = 1, ncoord
             distance2 = distance2 + dists(k)*dists(k)
          END DO

          ! If distance below limit, add observation to local domain
          IF (distance2 <= lradius2) THEN
             nobs_l_one = nobs_l_one + 1
          END IF
       END DO scancount_Geo

    END IF norm

  END SUBROUTINE cnt_dim_obs_l



!-------------------------------------------------------------------------------
!> Initialize local arrays for an observation
!!
!! This routine has to initialize for the current 
!! observation type the indices of the local observations
!! in the full observation vector and the corresponding 
!! distances from the local analysis domain. The offset
!! of the observation type in the local onbservation 
!! vector is given by OFF_OBS_L_ALL. 
!!
!! The routine has also to return OFF_OBS_L_ALL incremented
!! by the number of initialized local observations. 
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE init_obsarrays_l(disttype, ncoord, coord_l, lradius, nobs_l_one, nobs_f_one, &
       ocoord_f_one, dist_l_one, id_obs_l_one, off_obs_l_all, off_obs_f_all)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: disttype          !< Type of distance computation
                                             !< 0: Cartesian
                                             !< 1: geographic/spherical 
    INTEGER, INTENT(in) :: ncoord            !< Number of coordinates to consider
    REAL, INTENT(in) :: coord_l(ncoord)      !< Coordinates of current water column
    REAL, INTENT(in) :: lradius              !< Localization radius in meters
    INTEGER, INTENT(in) :: nobs_l_one        !< Local dimension of current observation vector
    INTEGER, INTENT(in) :: nobs_f_one        !< Full dimension of current observation vector
    REAL, INTENT(in) :: ocoord_f_one(:, :)   !< coordinate array for current observations
    REAL, INTENT(inout) :: dist_l_one(nobs_l_one)        !< Distance of obs. from water column
    INTEGER, INTENT(inout) :: id_obs_l_one(nobs_l_one)   !< index of this local obs. in this full obs. vector
    INTEGER, INTENT(inout) :: off_obs_l_all  !< input: offset of current obs. in local obs. vector
                                             !< output: input + nobs_l_one
    INTEGER, INTENT(inout) :: off_obs_f_all  !< input: offset of current obs. in full obs. vector
                                             !< output: input + nobs_f_one

! *** Local variables ***
    INTEGER :: i, k, off_obs ! Counter
    REAL :: ocoord(ncoord)  ! Coordinates of observation
    REAL :: dists(ncoord)   ! Distance vector between analysis point and observation
    REAL :: lradius2        ! squared localization radius
    REAL :: distance2       ! squared distance


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    ! Initialize squared localization radius
    lradius2 = lradius*lradius

    off_obs = 0

    ! Count local observations
    IF (nobs_l_one>0) THEN

       norm: IF (disttype==0) THEN

          ! *** Compute Cartesian distance ***

          scancount_Cart: DO i = 1, nobs_f_one

             ! location of observation point
             ocoord(1:ncoord) = ocoord_f_one(1:ncoord, i)

             ! Distance for each dimension
             DO k = 1, ncoord
                dists(k) = ABS(coord_l(k) - ocoord(k))
             END DO

             ! full squared distance in unit of coordinates
             distance2 = 0.0
             DO k = 1, ncoord
                distance2 = distance2 + dists(k)*dists(k)
             END DO

             ! If distance below limit, add observation to local domain
             IF (distance2 <= lradius2) THEN
                ! Count overall local observations
                off_obs_l_all = off_obs_l_all + 1     ! dimension
          
                ! For internal storage (use in init_obs_prof_l)
                off_obs = off_obs + 1                 ! dimension
                id_obs_l_one(off_obs) = i             ! node index
                dist_l_one(off_obs) = SQRT(distance2) ! distance
             END IF
          END DO scancount_Cart

       ELSEIF (disttype==1) THEN norm

          ! *** Count with distance from geographic coordinates ***
          
          ! For this computation of distances the coordinate array
          ! OCOORD is assumed to be structured as
          ! OCOORD(1): longitude in radians
          ! OCOORD(2): latitude in radians
          ! OCOORD(3): depth in meters (optional)

          scancount_Geo: DO i = 1, nobs_f_one

             ! location of observation point
             ocoord(1:ncoord) = ocoord_f_one(1:ncoord, i)

             ! approximate distances in longitude and latitude
             dists(1) = r_earth * ABS(coord_l(1) - ocoord(1))* COS(coord_l(2))
             dists(2) = r_earth * ABS(coord_l(2) - ocoord(2))
             IF (ncoord>2) dists(3) = ABS(coord_l(3) - ocoord(3))

             ! full squared distance in meters
             distance2 = 0.0
             DO k = 1, ncoord
                distance2 = distance2 + dists(k)*dists(k)
             END DO

             ! If distance below limit, add observation to local domain
             IF (distance2 <= lradius2) THEN
                ! Count overall local observations
                off_obs_l_all = off_obs_l_all + 1     ! dimension
          
                ! For internal storage (use in init_obs_prof_l)
                off_obs = off_obs + 1                 ! dimension
                id_obs_l_one(off_obs) = i             ! node index
                dist_l_one(off_obs) = SQRT(distance2) ! distance
             END IF
          END DO scancount_Geo
       END IF norm
    END IF

    ! Increment offset for next observation type
    off_obs_f_all = off_obs_f_all + nobs_f_one

  END SUBROUTINE init_obsarrays_l



!-------------------------------------------------------------------------------
!> Initialize local observation vector and inverse error variance
!!
!! This routine has to initialize the part of the 
!! overall local observation vector corresponding
!! to the current observation type. The offset of
!! the current observation type in the local obs.
!! vector is given by OFFSET_OBS_l_ALL.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE init_obs_l(nobs_l, thisobs_l, thisobs, obs_l_all)

    USE PDAFomi_obs_f, &
         ONLY: obs_f

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: nobs_l             !< Local dimension of obs. vector for all variables
    TYPE(obs_l), INTENT(inout) :: thisobs_l   !< Data type with local observation
    TYPE(obs_f), INTENT(inout) :: thisobs     !< Data type with full observation
    REAL, INTENT(inout) :: obs_l_all(:)       !< Local observation vector for all variables


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

    ! Initialize local observations
    CALL g2l_obs(nobs_l, thisobs_l%dim_obs_l, thisobs%dim_obs_f, thisobs_l%id_obs_l, &
         thisobs%obs_f, thisobs_l%off_obs_l, obs_l_all)

    ! Initialize local inverse variances for current observation
    ! they will be used in prodRinva_l
    IF (ALLOCATED(thisobs_l%ivar_obs_l)) DEALLOCATE(thisobs_l%ivar_obs_l)
    IF (thisobs_l%dim_obs_l>0) THEN
       ALLOCATE(thisobs_l%ivar_obs_l(thisobs_l%dim_obs_l))
    ELSE
       ALLOCATE(thisobs_l%ivar_obs_l(1))
    END IF

    CALL g2l_obs(thisobs_l%dim_obs_l, thisobs_l%dim_obs_l, thisobs%dim_obs_f, thisobs_l%id_obs_l, &
         thisobs%ivar_obs_f, 0, thisobs_l%ivar_obs_l)

  END SUBROUTINE init_obs_l



!-------------------------------------------------------------------------------
!> Initialize local observation vector
!!
!! This routine has to initialize the part of the 
!! overall local observation vector corresponding
!! to the current observation type. The offset of
!! the current observation type in the local obs.
!! vector is given by OFFSET_OBS_l_ALL.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE g2l_obs(nobs_l_all, nobs_l_one, nobs_f_one, id_obs_l_one, &
       obs_f_one, offset_obs_l_all, obs_l_all)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: nobs_l_all          !< Local dimension of obs. vector for all variables
    INTEGER, INTENT(in) :: nobs_l_one          !< Local number of obs. of current obs. type
    INTEGER, INTENT(in) :: nobs_f_one          !< Full dimension of obs. vector for current obs. type
    REAL, INTENT(in) :: obs_f_one(nobs_f_one)  !< Full obs. vector of current obs. type
    INTEGER, INTENT(in) :: id_obs_l_one(nobs_l_one) !< local index of obs. for current obs. type
    REAL, INTENT(inout) :: obs_l_all(nobs_l_all)    !< Local observation vector for all variables
    INTEGER, INTENT(in) :: offset_obs_l_all    !< offset of current observation in obs_l_all and ivar_l_all

! *** Local variables ***
    INTEGER :: i  ! Counter


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

    ! Local observations
    DO i = 1, nobs_l_one
       obs_l_all(i+offset_obs_l_all) = obs_f_one(id_obs_l_one(i))
    ENDDO

  END SUBROUTINE g2l_obs




!-------------------------------------------------------------------------------
!> Compute product of inverse of R with some matrix
!!
!! The routine is called during the analysis step
!! on each local analysis domain. It has to 
!! compute the product of the inverse of the local
!! observation error covariance matrix with
!! the matrix of locally observed ensemble 
!! perturbations.
!!
!! Next to computing the product, a localizing 
!! weighting ("observation localization") can be
!! applied to matrix A.
!!
!! This implementation assumes a diagonal observation
!! error covariance matrix, and supports varying
!! observation error variances.
!!
!! The routine can be applied with either all observations
!! of different types at once, or separately for each
!! observation type.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE prodRinvA_l(verbose, nobs_l, rank, locweight, lradius, sradius, &
       ivar_obs_l, dist_l, A_l, C_l)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: verbose        !< Verbosity flag
    INTEGER, INTENT(in) :: nobs_l         !< Dimension of local obs. vector (one or all obs. types)
    INTEGER, INTENT(in) :: rank           !< Rank of initial covariance matrix
    INTEGER, INTENT(in) :: locweight      !< Localization weight type
    REAL, INTENT(in)    :: lradius        !< localization radius
    REAL, INTENT(in)    :: sradius        !< support radius for weight functions
    REAL, INTENT(in)    :: ivar_obs_l(:)  !< Local vector of inverse obs. variances (nobs_l)
    REAL, INTENT(in)    :: dist_l(:)      !< Local vector of obs. distances (nobs_l)
    REAL, INTENT(inout) :: A_l(:, :)      !< Input matrix (nobs_l, rank)
    REAL, INTENT(out)   :: C_l(:, :)      !< Output matrix (nobs_l, rank)


! *** local variables ***
    INTEGER :: i, j          ! Index of observation component
    INTEGER :: verbose_w     ! Verbosity flag for weight computation
    INTEGER :: wtype         ! Type of weight function
    INTEGER :: rtype         ! Type of weight regulation
    REAL    :: var_obs_l     ! Variance of observation error
    REAL, ALLOCATABLE :: weight(:)     ! Localization weights
    REAL, ALLOCATABLE :: A_obs(:,:)    ! Array for a single row of A_l
    

! **********************
! *** INITIALIZATION ***
! **********************

  ! Screen output
  IF (verbose == 1) THEN
     WRITE (*, '(a, 5x, a, 1x)') &
          'PDAFomi', '--- Domain localization'
     WRITE (*, '(a, 8x, a, 1x, es11.3)') &
          'PDAFomi', '--- Local influence radius', lradius

     IF (locweight == 5 .OR. locweight == 6 .OR. locweight == 7) THEN
        WRITE (*, '(a, 8x, a)') &
             'PDAFomi', '--- Use distance-dependent weight for observed ensemble'
     ELSE IF (locweight == 1 .OR. locweight == 2 .OR. locweight == 3 &
          .OR. locweight == 4) THEN
        WRITE (*, '(a, 8x, a)') &
             'PDAFomi', '--- Use distance-dependent weight for observation errors'

        IF (locweight == 3) THEN
           write (*, '(a, 8x, a)') &
                'PDAFomi', '--- Use regulated weight with mean error variance'
        ELSE IF (locweight == 4) THEN
           write (*, '(a, 8x, a)') &
                'PDAFomi', '--- Use regulated weight with single-point error variance'
        END IF
     END IF
  ENDIF


! ***********************************************
! *** Apply a weight matrix with correlations ***
! *** of compact support to matrix A or the   ***
! *** observation error covariance matrix.    ***
! ***********************************************

! *** Initialize weight array

    ! Allocate weight array
    ALLOCATE(weight(nobs_l))

    IF (locweight == 0) THEN
       ! Uniform (unit) weighting
       wtype = 0
       rtype = 0
    ELSE IF (locweight == 1 .OR. locweight == 5) THEN
       ! Exponential weighting
       wtype = 1
       rtype = 0
    ELSE IF (locweight == 2 .OR. locweight == 3 .OR. locweight == 4 &
         .OR. locweight == 6 .OR. locweight == 7) THEN
       ! 5th-order polynomial (Gaspari&Cohn, 1999)
       wtype = 2

       IF (locweight == 3 .OR. locweight == 4) THEN
          ! Use regulated weight
          rtype = 1
       ELSE   
          ! No regulated weight
          rtype = 0
       END IF

    END IF

    IF (locweight == 4) THEN
       ! Allocate array for single observation point
       ALLOCATE(A_obs(1, rank))
    END IF

    DO i = 1, nobs_l

       ! Control verbosity of PDAF_local_weight
       IF (verbose==1 .AND. i==1) THEN
          verbose_w = 1
       ELSE
          verbose_w = 0
       END IF

       ! set observation variance value
       var_obs_l = 1.0 / ivar_obs_l(i)

       IF (locweight /= 4) THEN
          ! All localizations except regulated weight based on variance at 
          ! single observation point
          CALL PDAF_local_weight(wtype, rtype, lradius, sradius, dist_l(i), &
               nobs_l, rank, A_l, var_obs_l, weight(i), verbose_w)
       ELSE
          ! Regulated weight using variance at single observation point
          A_obs(1,:) = A_l(i,:)
          CALL PDAF_local_weight(wtype, rtype, lradius, sradius, dist_l(i), &
               1, rank, A_obs, var_obs_l, weight(i), verbose_w)
       END IF
    END DO
  
    IF (locweight == 4) DEALLOCATE(A_obs)

! *** Handling of special weighting types ***

    lw2: IF (locweight ==6 ) THEN
       ! Use square-root of 5th-order polynomial on A

       DO i = 1, nobs_l
          ! Check if weight >0 (Could be <0 due to numerical precision)
          IF (weight(i) > 0.0) THEN
             weight(i) = SQRT(weight(i))
          ELSE
             weight(i) = 0.0
          END IF
       END DO
    END IF lw2


! *** Apply weight

    doweighting: IF (locweight >= 5) THEN

       ! *** Apply weight to matrix A
       DO j = 1, rank
          DO i = 1, nobs_l
             A_l(i, j) = weight(i) * A_l(i, j)
          END DO
       END DO

       ! ***       -1
       ! ***  C = R   A 
       DO j = 1, rank
          DO i = 1, nobs_l
             C_l(i, j) = ivar_obs_l(i) * A_l(i, j)
          END DO
       END DO
  
    ELSE doweighting

       ! *** Apply weight to matrix R only
       DO j = 1, rank
          DO i = 1, nobs_l
             C_l(i, j) = ivar_obs_l(i) * weight(i) * A_l(i, j)
          END DO
       END DO
     
    END IF doweighting

! *** Clean up ***

    DEALLOCATE(weight)

  END SUBROUTINE prodRinvA_l



!-------------------------------------------------------------------------------
!> Compute mean observation error variance
!!
!! This routine will only be called, if the adaptive
!! forgetting factor feature is used. Please note that
!! this is an experimental feature.
!!
!! The routine is called in the loop over all
!! local analysis domains during each analysis
!! by the routine PDAF_set_forget_local that 
!! estimates a local adaptive forgetting factor.
!! The routine has to initialize the mean observation 
!! error variance for the current local analysis 
!! domain.  (See init_obsvar_f for a global variant)
!!
!! The routine assumed a diagonal observation error
!! covariance matrix.
!!
!! If the observation counter is zero the computation
!! of the mean variance is initialized. The output is 
!! always the mean variance. If the observation counter
!! is >0 first the variance sum is computed by 
!! multiplying with the observation counter.
!!
!! __Revision history:__
!! * 2019-09 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE init_obsvar_l(thisobs_l, meanvar_l, cnt_obs_l)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    REAL, INTENT(inout) :: meanvar_l         !< Mean variance
    INTEGER, INTENT(inout) :: cnt_obs_l      !< Observation counter

! Local variables
    INTEGER :: i        ! Counter


! ***********************************
! *** Compute local mean variance ***
! ***********************************

    IF (cnt_obs_l==0) THEN
       ! Reset mean variance
       meanvar_l = 0.0
    ELSE
       ! Compute sum of variances from mean variance
       meanvar_l = meanvar_l * REAL(cnt_obs_l)
    END IF

    ! Add observation error variances
    DO i = 1, thisobs_l%dim_obs_l
       meanvar_l = meanvar_l + 1.0 / thisobs_l%ivar_obs_l(i)
    END DO

    ! Increment observation count
    cnt_obs_l = cnt_obs_l + thisobs_l%dim_obs_l

    ! Compute updated mean variance
    meanvar_l = meanvar_l / REAL(cnt_obs_l)

  END SUBROUTINE init_obsvar_l


  
!-------------------------------------------------------------------------------
!> Set debugging flag
!!
!! This routine sets the debug flag for PDAF-OMI.
!! One can set the flag dependent on the local analysis
!! domain, the MPI rank, or the OpenMP thread ID, or
!! and combination of them.
!!
!! For debugval>0 additional information is written by
!! the OMI routine to stdout. One should activate the 
!! debugging before calling some selected routine(s) and
!! deactivate it with debugval=0 afterwards. This allows 
!! for a targeted checking of the functionality.
!!
!! __Revision history:__
!! * 2019-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE set_debug_flag(debugval)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: debugval          !< Value for debugging flag

    debug = debugval

  END SUBROUTINE set_debug_flag

END MODULE PDAFomi_obs_l
