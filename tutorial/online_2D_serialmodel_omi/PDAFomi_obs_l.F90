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
!BOP
!
! !MODULE:
MODULE PDAFomi_obs_l
!
! !DESCRIPTION:
! This module contains generic routines for several observation-related
! operations for local filters. The routines are
!
! cnt_dim_obs_l
!        Set dimension of local obs. vector
! init_obsarrays_l
!        Initialize arrays for the index of a local observation in 
!        the full observation vector and its corresponding distance.
! init_obs_l
!        Initialize the local vector of observations
! prodRinvA_l
!        Multiply an intermediate matrix of the fitler analysis with
!        the inverse of the observation error covariance matrix
!
! !REVISION HISTORY:
! 2019-06 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  SAVE

  REAL, PARAMETER :: r_earth=6.3675e6  ! Earth radius in meters
  INTEGER :: debug     ! Debugging flag

  ! Data type to define the local observations by internally shared variables of the module
  type obs_l
     INTEGER :: dim_obs_l                 ! number of local observations
     INTEGER :: off_obs_l                 ! Offset of this observation in overall local obs. vector
     INTEGER, ALLOCATABLE :: id_obs_l(:)  ! Indices of local observations in full obs. vector 
     REAL, ALLOCATABLE :: distance_l(:)   ! Distances of local observations
     REAL, ALLOCATABLE :: ivar_obs_l(:)   ! Inverse variance of local observations
  end type obs_l

!$OMP THREADPRIVATE(debug)

!EOP
!-------------------------------------------------------------------------------

CONTAINS


!BOP
!
! !ROUTINE: cnt_dim_obs_l --- Set dimension of local obs. vector
!
! !INTERFACE:
  SUBROUTINE cnt_dim_obs_l(disttype, ncoord, wc_coord, lradius, nobs_f_one, ocoord_f_one, nobs_l_one)

! !DESCRIPTION:
! This routine sets the number of local observations for the
! current observation type for the local analysis domain
! with coordinates WC_COORD and localization radius LRADIUS.
!
! The routine is called by each filter process.
!
! !REVISION HISTORY:
! 2019-06 - Lars Nerger - Initial code from restructuring observation routines
! Later revisions - see svn log
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: disttype          ! Type of distance computation
                                             ! 0: Cartesian
                                             ! 1: geographic/spherical 
    INTEGER, INTENT(in) :: ncoord            ! Number of coordinates to consider
    REAL, INTENT(in) :: wc_coord(ncoord)     ! Coordinates of current analysis domain
    REAL, INTENT(in) :: lradius              ! Localization radius in meters
    INTEGER, INTENT(in) :: nobs_f_one        ! Full dimension of current observation vector
    REAL, INTENT(in) :: ocoord_f_one(:, :)   ! coordinate array for current observations
    INTEGER, INTENT(out) :: nobs_l_one       ! Local dimension of current observation vector
!EOP

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
             dists(k) = ABS(wc_coord(k) - ocoord(k))
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
          dists(1) = r_earth * ABS(wc_coord(1) - ocoord(1))* COS(wc_coord(2))
          dists(2) = r_earth * ABS(wc_coord(2) - ocoord(2))
          IF (ncoord>2) dists(3) = ABS(wc_coord(3) - ocoord(3))

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
!BOP
!
! !ROUTINE: init_obsarrays_l --- Init. local arrays for an observation
!
! !INTERFACE:
  SUBROUTINE init_obsarrays_l(disttype, ncoord, wc_coord, lradius, nobs_l_one, nobs_f_one, &
       ocoord_f_one, dist_l_one, id_obs_l_one, off_obs_l_all, off_obs_f_all)

! !DESCRIPTION:
! This routine has to initialize for the current 
! observation type the indices of the local observations
! in the full observation vector and the corresponding 
! distances from the local analysis domain. The offset
! of the observation type in the local onbservation 
! vector is given by OFF_OBS_L_ALL. 
!
! The routine has also to return OFF_OBS_L_ALL incremented
! by the number of initialized local observations. 
!
! The routine is called by each filter process.
!
! !REVISION HISTORY:
! 2019-06 - Lars Nerger - Initial code from restructuring observation routines
! Later revisions - see svn log
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: disttype          ! Type of distance computation
                                             ! 0: Cartesian
                                             ! 1: geographic/spherical 
    INTEGER, INTENT(in) :: ncoord            ! Number of coordinates to consider
    REAL, INTENT(in) :: wc_coord(ncoord)     ! Coordinates of current water column
    REAL, INTENT(in) :: lradius              ! Localization radius in meters
    INTEGER, INTENT(in) :: nobs_l_one        ! Local dimension of current observation vector
    INTEGER, INTENT(in) :: nobs_f_one        ! Full dimension of current observation vector
    REAL, INTENT(in) :: ocoord_f_one(:, :)   ! coordinate array for current observations
    REAL, INTENT(inout) :: dist_l_one(nobs_l_one)        ! Distance of obs. from water column
    INTEGER, INTENT(inout) :: id_obs_l_one(nobs_l_one)   ! index of current local obs. in current full obs. vector
    INTEGER, INTENT(inout) :: off_obs_l_all  ! input: offset of current obs. in local obs. vector
                                             ! output: input + nobs_l_one
    INTEGER, INTENT(inout) :: off_obs_f_all  ! input: offset of current obs. in full obs. vector
                                             ! output: input + nobs_f_one
!EOP

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

             ! approximate distances in longitude and latitude
             DO k = 1, ncoord
                dists(k) = ABS(wc_coord(k) - ocoord(k))
             END DO

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
          END DO scancount_Cart

       ELSEIF (disttype==1) THEN norm

          ! *** Count with distance from geographic coordinates ***

          scancount_Geo: DO i = 1, nobs_f_one

             ! location of observation point
             ocoord(1:ncoord) = ocoord_f_one(1:ncoord, i)

             ! approximate distances in longitude and latitude
             dists(1) = r_earth * ABS(wc_coord(1) - ocoord(1))* COS(wc_coord(2))
             dists(2) = r_earth * ABS(wc_coord(2) - ocoord(2))
             IF (ncoord>2) dists(3) = ABS(wc_coord(3) - ocoord(3))

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
!BOP
!
! !ROUTINE: init_obs_l --- Initialize local observation vector
!
! !INTERFACE:
  SUBROUTINE init_obs_l(nobs_l_all, nobs_l_one, nobs_f_one, id_obs_l_one, &
       obs_f_one, offset_obs_l_all, obs_l_all)

! !DESCRIPTION:
! This routine has to initialize the part of the 
! overall local observation vector corresponding
! to the current observation type. The offset of
! the current observation type in the local obs.
! vector is given by OFFSET_OBS_l_ALL.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2019-06 - Lars Nerger - Initial code from restructuring observation routines
! Later revisions - see svn log
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: nobs_l_all          ! Local dimension of obs. vector for all variables
    INTEGER, INTENT(in) :: nobs_l_one          ! Local number of obs. of current obs. type
    INTEGER, INTENT(in) :: nobs_f_one          ! Full dimension of obs. vector for current obs. type
    REAL, INTENT(in) :: obs_f_one(nobs_f_one)  ! Full obs. vector of current obs. type
    INTEGER, INTENT(in) :: id_obs_l_one(nobs_l_one) ! local index of obs. for current obs. type
    REAL, INTENT(inout) :: obs_l_all(nobs_l_all)    ! Local observation vector for all variables
    INTEGER, INTENT(in) :: offset_obs_l_all    ! offset of current observation in obs_l_all and ivar_l_all
!EOP

! *** Local variables ***
    INTEGER :: i  ! Counter


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

    ! Local observations
    DO i = 1, nobs_l_one
       obs_l_all(i+offset_obs_l_all) = obs_f_one(id_obs_l_one(i))
    ENDDO

  END SUBROUTINE init_obs_l




!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: prodRinvA_l --- Compute product of inverse of R with some matrix
!
! !INTERFACE:
  SUBROUTINE prodRinvA_l(verbose, nobs_l_all, rank, locweight, lradius, sradius, &
       ivar_obs_l_all, dist_l_all, A_l, C_l)

! !DESCRIPTION:
! The routine is called during the analysis step
! on each local analysis domain. It has to 
! compute the product of the inverse of the local
! observation error covariance matrix with
! the matrix of locally observed ensemble 
! perturbations.
! Next to computing the product, a localizing 
! weighting ("observation localization") can be
! applied to matrix A.
! This implementation assumes a diagonal observation
! error covariance matrix, and supports varying
! observation error variances.
!
! This routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2019-06 - Lars Nerger - Initial code from restructuring observation routines
! Later revisions - see svn log
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: verbose           ! Verbosity flag
    INTEGER, INTENT(in) :: nobs_l_all        ! Dimension of local observation vector
    INTEGER, INTENT(in) :: rank              ! Rank of initial covariance matrix
    INTEGER, INTENT(in) :: locweight         ! Localization weight type
    REAL, INTENT(in)    :: lradius           ! localization radius
    REAL, INTENT(in)    :: sradius           ! support radius for weight functions
    REAL, INTENT(in)    :: ivar_obs_l_all(nobs_l_all) ! Local vector of inverse obs. variances
    REAL, INTENT(in)    :: dist_l_all(nobs_l_all)     ! Local vector of obs. distances
    REAL, INTENT(inout) :: A_l(nobs_l_all, rank)      ! Input matrix
    REAL, INTENT(out)   :: C_l(nobs_l_all, rank)      ! Output matrix
!EOP


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
          'PDAF_MOD_OBS_L', '--- Domain localization'
     WRITE (*, '(a, 8x, a, 1x, es11.3)') &
          'PDAF_MOD_OBS_L', '--- Local influence radius', lradius

     IF (locweight == 5 .OR. locweight == 6 .OR. locweight == 7) THEN
        WRITE (*, '(a, 8x, a)') &
             'PDAF_MOD_OBS_L', '--- Use distance-dependent weight for observed ensemble'
     ELSE IF (locweight == 1 .OR. locweight == 2 .OR. locweight == 3 &
          .OR. locweight == 4) THEN
        WRITE (*, '(a, 8x, a)') &
             'PDAF_MOD_OBS_L', '--- Use distance-dependent weight for observation errors'

        IF (locweight == 3) THEN
           write (*, '(a, 8x, a)') &
                'PDAF_MOD_OBS_L', '--- Use regulated weight with mean error variance'
        ELSE IF (locweight == 4) THEN
           write (*, '(a, 8x, a)') &
                'PDAF_MOD_OBS_L', '--- Use regulated weight with single-point error variance'
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
    ALLOCATE(weight(nobs_l_all))

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

    DO i = 1, nobs_l_all

       ! Control verbosity of PDAF_local_weight
       IF (verbose==1 .AND. i==1) THEN
          verbose_w = 1
       ELSE
          verbose_w = 0
       END IF

       ! set observation variance value
       var_obs_l = 1.0 / ivar_obs_l_all(i)

       IF (locweight /= 4) THEN
          ! All localizations except regulated weight based on variance at 
          ! single observation point
          CALL PDAF_local_weight(wtype, rtype, lradius, sradius, dist_l_all(i), &
               nobs_l_all, rank, A_l, var_obs_l, weight(i), verbose_w)
       ELSE
          ! Regulated weight using variance at single observation point
          A_obs(1,:) = A_l(i,:)
          CALL PDAF_local_weight(wtype, rtype, lradius, sradius, dist_l_all(i), &
               1, rank, A_obs, var_obs_l, weight(i), verbose_w)
       END IF
    END DO
  
    IF (locweight == 4) DEALLOCATE(A_obs)

! *** Handling of special weighting types ***

    lw2: IF (locweight ==6 ) THEN
       ! Use square-root of 5th-order polynomial on A

       DO i = 1, nobs_l_all
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
          DO i = 1, nobs_l_all
             A_l(i, j) = weight(i) * A_l(i, j)
          END DO
       END DO

       ! ***       -1
       ! ***  C = R   A 
       DO j = 1, rank
          DO i = 1, nobs_l_all
             C_l(i, j) = ivar_obs_l_all(i) * A_l(i, j)
          END DO
       END DO
  
    ELSE doweighting

       ! *** Apply weight to matrix R only
       DO j = 1, rank
          DO i = 1, nobs_l_all
             C_l(i, j) = ivar_obs_l_all(i) * weight(i) * A_l(i, j)
          END DO
       END DO
     
    END IF doweighting

! *** Clean up ***

    DEALLOCATE(weight)

  END SUBROUTINE prodRinvA_l


  
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: set_debug_flag - Set debugging flag
!
! !INTERFACE:
  SUBROUTINE set_debug_flag(debugval)


! !DESCRIPTION:
! This routine sets the debug flag for PDAF-OMI.
! One can set the flag dependent on the local analysis
! domain, the MPI rank, or the OpenMP thread ID, or
! and combination of them.
! For debugval>0 additional information is written by
! the OMI routine to stdout. One should activate the 
! debugging before calling some selected routine(s) and
! deactivate it with debugval=0 afterwards. This allows 
! for a targeted checking of the functionality.
!
! !REVISION HISTORY:
! 2019-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: debugval          ! Value for debugging flag

    debug = debugval

  END SUBROUTINE set_debug_flag

END MODULE PDAFomi_obs_l
