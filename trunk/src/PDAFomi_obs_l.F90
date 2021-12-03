! Copyright (c) 2004-2021 Lars Nerger
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
!! * PDAFomi_set_debug_flag \n
!!        Set or unset the debugging flag for PDAFomi routines
!! * PDAFomi_init_dim_obs_l \n
!!        Initialize dimension of local obs. vetor and arrays for
!!        local observations
!! * PDAFomi_cnt_dim_obs_l \n
!!        Set dimension of local obs. vector
!! * PDAFomi_init_obsarrays_l \n
!!        Initialize arrays for the index of a local observation in 
!!        the full observation vector and its corresponding distance.
!! * PDAFomi_g2l_obs \n
!!        Initialize local observation vector from full observation vector
!! * PDAFomi_init_obs_l \n
!!        Initialize the local vector of observations
!! * PDAFomi_prodRinvA_l \n
!!        Multiply an intermediate matrix of the local filter analysis
!!        with the inverse of the observation error covariance matrix
!!        and apply observation localization
!! * PDAFomi_init_obsvar_l \n
!!        Compute mean observation error variance
!! * PDAFomi_likelihood_l \n
!!        Compute local likelihood for an ensemble member
!! * PDAFomi_localize_covar \n
!!        Apply covariance localization in LEnKF
!! * PDAFomi_g2l_obs_internal \n
!!        Internal routine to initialize local observation vector from full
!!        observation vector (used by PDAFomi_init_obs_l and PDAFomi_g2l_obs)
!! * PDAFomi_comp_dist2 \n
!!        Compute squared distance
!! * PDAFomi_weights_l \n
!!        Compute a vector of localization weights
!! * PDAFomi_deallocate_obs \n
!!        Deallocate arrays in observation type
!! * PDAFomi_dealloc \n
!!        Deallocate arrays in all observation types
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE PDAFomi_obs_l

  USE PDAFomi_obs_f, ONLY: obs_f, r_earth, pi, debug, n_obstypes
  USE PDAF_mod_filter, ONLY: screen, obs_member
  USE PDAF_mod_filtermpi, ONLY: mype

  IMPLICIT NONE
  SAVE

! *** Module internal variables

  ! Data type to define the local observations by internally shared variables of the module
  TYPE obs_l
     INTEGER :: dim_obs_l                 !< number of local observations
     INTEGER :: off_obs_l                 !< Offset of this observation in overall local obs. vector
     INTEGER, ALLOCATABLE :: id_obs_l(:)  !< Indices of local observations in full obs. vector 
     REAL, ALLOCATABLE :: distance_l(:)   !< Distances of local observations
     REAL, ALLOCATABLE :: ivar_obs_l(:)   !< Inverse variance of local observations
     INTEGER :: locweight                 !< Specify localization function
     REAL :: lradius                      !< localization radius
     REAL :: sradius                      !< support radius for localization function
  END TYPE obs_l

  TYPE obs_arr_l                          ! Type for pointer array over all observation types
     TYPE(obs_l), POINTER :: ptr
  END TYPE obs_arr_l

  TYPE(obs_arr_l), ALLOCATABLE :: obs_l_all(:) ! Declare pointer array

  INTEGER :: firstobs = 0                 ! Flag for very first call to init_dim_obs_l
  INTEGER :: offset_obs_l = 0             ! offset of current observation in overall local obs. vector

!$OMP THREADPRIVATE(obs_l_all, firstobs, offset_obs_l)


!-------------------------------------------------------------------------------

CONTAINS


!!> Set debugging flag
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
  SUBROUTINE PDAFomi_set_debug_flag(debugval)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: debugval          !< Value for debugging flag

    debug = debugval

    ! Print debug information
    IF (debug>0) THEN
       WRITE (*,*) '++ OMI-debug set_debug_flag: mype_filter', mype, 'activate', debug
    END IF

  END SUBROUTINE PDAFomi_set_debug_flag




!-------------------------------------------------------------------------------
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
  SUBROUTINE PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, locweight, lradius, &
       sradius, cnt_obs_l)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    TYPE(obs_l), TARGET, INTENT(inout) :: thisobs_l  !< Data type with local observation
    REAL, INTENT(in) :: coords_l(:)          !< Coordinates of current analysis domain
    INTEGER, INTENT(in) :: locweight         !< Type of localization function
    REAL, INTENT(in) :: lradius              !< Localization radius
    REAL, INTENT(in) :: sradius              !< Support radius of localization function
    INTEGER, INTENT(inout) :: cnt_obs_l      !< Local dimension of current observation vector


    doassim: IF (thisobs%doassim == 1) THEN


! ***********************************************
! *** Check offset in full observation vector ***
! ***********************************************

       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_init_dim_obs_l -- START'

       ! Store ID of first observation type that call the routine
       ! This is reset in PDAFomi_deallocate_obs
       IF (firstobs == 0) THEN
          firstobs = thisobs%obsid
       END IF

       ! Reset offset of currrent observation in overall local obs. vector
       IF (thisobs%obsid == firstobs) THEN
          offset_obs_l = 0
          cnt_obs_l = 0
       END IF


! **************************************
! *** Store localization information ***
! **************************************

       thisobs_l%locweight = locweight
       thisobs_l%lradius = lradius
       thisobs_l%sradius = sradius


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, &
               '   PDAFomi_init_dim_obs_l -- count local observations'
          IF (thisobs%obsid == firstobs) THEN
             WRITE (*,*) '++ OMI-debug init_dim_obs_l:', debug, '  Re-init dim_obs_l=0'
          END IF
          WRITE (*,*) '++ OMI-debug init_dim_obs_l:', debug, '  coords_l', coords_l
       END IF

       CALL PDAFomi_cnt_dim_obs_l(thisobs_l, thisobs, coords_l, lradius)

       ! Store number of local module-type observations for output
       cnt_obs_l = cnt_obs_l + thisobs_l%dim_obs_l


! **************************************************
! *** Initialize local observation pointer array ***
! **************************************************

       ! Initialize pointer array
       IF (thisobs%obsid == firstobs) THEN
          IF (ALLOCATED(obs_l_all)) DEALLOCATE(obs_l_all)
          ALLOCATE(obs_l_all(n_obstypes))
       END IF

       ! Set pointer to current observation
       obs_l_all(thisobs%obsid)%ptr => thisobs_l


! ************************************************************
! *** Initialize internal local arrays for local distances ***
! *** and indices of local obs. in full obs. vector        ***
! ************************************************************

       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, &
            '   PDAFomi_init_dim_obs_l -- initialize local observation arrays'

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

       ! Store offset
       thisobs_l%off_obs_l = offset_obs_l

       ! Initialize ID_OBS_L and DISTANCE_L and increment offsets
       CALL PDAFomi_init_obsarrays_l(thisobs_l, thisobs, coords_l, lradius, &
            offset_obs_l)

       ! Print debug information
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug init_dim_obs_l:', debug, '  thisobs_l%dim_obs_l', thisobs_l%dim_obs_l
          IF (thisobs_l%dim_obs_l>0) THEN
             WRITE (*,*) '++ OMI-debug init_dim_obs_l:', debug, '  thisobs_l%id_obs_l', thisobs_l%id_obs_l
             WRITE (*,*) '++ OMI-debug init_dim_obs_l:', debug, '  thisobs_l%distance_l', thisobs_l%distance_l
          END IF
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_init_dim_obs_l -- END'
       END IF

    ELSE doassim

       cnt_obs_l = cnt_obs_l + 0

    END IF doassim


  END SUBROUTINE PDAFomi_init_dim_obs_l




!-------------------------------------------------------------------------------
!> Set dimension of local observation vector
!!
!! This routine sets the number of local observations for the
!! current observation type for the local analysis domain
!! with coordinates COORDS_L and localization radius LRADIUS.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_cnt_dim_obs_l(thisobs_l, thisobs, coords_l, lradius)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    REAL, INTENT(in) :: coords_l(:)          !< Coordinates of current analysis domain (thisobs%ncoord)
    REAL, INTENT(in) :: lradius              !< Localization radius in meters

! *** Local variables ***
    INTEGER :: i            ! Counters
    REAL :: ocoord(thisobs%ncoord)  ! Coordinates of observation
    REAL :: lradius2        ! squared localization radius
    REAL :: distance2       ! squared distance


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    ! Initialize squared localization radius
    lradius2 = lradius*lradius

    ! Count local observations
    thisobs_l%dim_obs_l = 0

    IF (debug>0) THEN
       WRITE (*,*) '++ OMI-debug cnt_dim_obs_l: ', debug, '  thisobs%ncoord', thisobs%ncoord
       WRITE (*,*) '++ OMI-debug cnt_dim_obs_l: ', debug, '  thisobs_l%lradius', thisobs_l%lradius
       WRITE (*,*) '++ OMI-debug cnt_dim_obs_l: ', debug, '  Check for observations within radius'
    END IF

    scancount: DO i = 1, thisobs%dim_obs_f

       ! location of observation point
       ocoord(1:thisobs%ncoord) = thisobs%ocoord_f(1:thisobs%ncoord, i)

       CALL PDAFomi_comp_dist2(thisobs, coords_l, ocoord, distance2, i-1)

       ! If distance below limit, add observation to local domain
       IF (distance2 <= lradius2) THEN
          IF (debug>0) THEN
             WRITE (*,*) '++ OMI-debug cnt_dim_obs_l: ', debug, &
                  '  valid observation with coordinates', ocoord(1:thisobs%ncoord)
          END IF

          thisobs_l%dim_obs_l = thisobs_l%dim_obs_l + 1
       END IF
    END DO scancount

  END SUBROUTINE PDAFomi_cnt_dim_obs_l



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
  SUBROUTINE PDAFomi_init_obsarrays_l(thisobs_l, thisobs, coords_l, lradius, &
       off_obs_l_all)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    REAL, INTENT(in) :: coords_l(:)          !< Coordinates of current water column (thisobs%ncoord)
    REAL, INTENT(in) :: lradius              !< Localization radius in meters
    INTEGER, INTENT(inout) :: off_obs_l_all  !< input: offset of current obs. in local obs. vector
                                             !< output: input + thisobs_l%dim_obs_l

! *** Local variables ***
    INTEGER :: i, off_obs   ! Counters
    REAL :: ocoord(thisobs%ncoord)  ! Coordinates of observation
    REAL :: lradius2        ! squared localization radius
    REAL :: distance2       ! squared distance


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    ! Initialize squared localization radius
    lradius2 = lradius*lradius

    off_obs = 0

    ! Count local observations
    IF (thisobs_l%dim_obs_l>0) THEN

       scancount: DO i = 1, thisobs%dim_obs_f

          ! location of observation point
          ocoord(1:thisobs%ncoord) = thisobs%ocoord_f(1:thisobs%ncoord, i)

          CALL PDAFomi_comp_dist2(thisobs, coords_l, ocoord, distance2, i-1)

          ! If distance below limit, add observation to local domain
          IF (distance2 <= lradius2) THEN
             ! Count overall local observations
             off_obs_l_all = off_obs_l_all + 1     ! dimension
          
             ! For internal storage (use in init_obs_prof_l)
             off_obs = off_obs + 1                 ! dimension
             thisobs_l%id_obs_l(off_obs) = i             ! node index
             thisobs_l%distance_l(off_obs) = SQRT(distance2) ! distance
          END IF
       END DO scancount
    END IF

  END SUBROUTINE PDAFomi_init_obsarrays_l



!-------------------------------------------------------------------------------
!> Initialize local observation vector
!!
!! This routine has to initialize the part of the 
!! overall local observation vector corresponding
!! to the current observation type. The offset of
!! the current observation type in the local obs.
!! vector is given by OFFSET_OBS_l_ALL. 
!! This routine uses a shortened interface and just
!! passed the operation to the actually routine.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_g2l_obs(thisobs_l, thisobs, obs_f_all, obs_l_all)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    REAL, INTENT(in) :: obs_f_all(:)         !< Full obs. vector of current obs. for all variables
    REAL, INTENT(inout) :: obs_l_all(:)      !< Local observation vector for all variables


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

    doassim: IF (thisobs%doassim == 1) THEN

       IF (debug>0) THEN
          If (obs_member==0) THEN
             WRITE (*,*) '++ OMI-debug: ', debug, &
                  'PDAFomi_g2l_obs -- START Get local observed ensemble mean'
          ELSE
             WRITE (*,*) '++ OMI-debug: ', debug, &
                  'PDAFomi_g2l_obs -- START Get local observed ensemble member', obs_member
          END If
       END IF

       CALL PDAFomi_g2l_obs_internal(thisobs_l, &
            obs_f_all(thisobs%off_obs_f+1:thisobs%off_obs_f+thisobs%dim_obs_f), &
            thisobs_l%off_obs_l, obs_l_all)

       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_g2l_obs -- END'

    END IF doassim

  END SUBROUTINE PDAFomi_g2l_obs



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
  SUBROUTINE PDAFomi_init_obs_l(thisobs_l, thisobs, obs_l_all)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l   !< Data type with local observation
    TYPE(obs_f), INTENT(inout) :: thisobs     !< Data type with full observation
    REAL, INTENT(inout) :: obs_l_all(:)       !< Local observation vector for all variables


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

    doassim: IF (thisobs%doassim == 1) THEN

       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_init_obs_l -- START'
          WRITE (*,*) '++ OMI-debug init_obs_l:    ', debug, '  thisobs_l%dim_obs_l', thisobs_l%dim_obs_l
          WRITE (*,*) '++ OMI-debug init_obs_l:    ', debug, '  thisobs_l%off_obs_l', thisobs_l%off_obs_l
          WRITE (*,*) '++ OMI-debug: ', debug, &
               '   PDAFomi_init_obs_l -- Get local vector of observations'
       END IF

       ! Initialize local observations
       CALL PDAFomi_g2l_obs_internal(thisobs_l, thisobs%obs_f, thisobs_l%off_obs_l, obs_l_all)

       ! Initialize local inverse variances for current observation
       ! they will be used in prodRinva_l
       IF (ALLOCATED(thisobs_l%ivar_obs_l)) DEALLOCATE(thisobs_l%ivar_obs_l)
       IF (thisobs_l%dim_obs_l>0) THEN
          ALLOCATE(thisobs_l%ivar_obs_l(thisobs_l%dim_obs_l))
       ELSE
          ALLOCATE(thisobs_l%ivar_obs_l(1))
       END IF

       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, &
               '   PDAFomi_init_obs_l -- Get local vector of inverse obs. variances'
       END IF

       CALL PDAFomi_g2l_obs_internal(thisobs_l, thisobs%ivar_obs_f, 0, thisobs_l%ivar_obs_l)

       ! Print debug information
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_init_obs_l -- END'
       END IF

    END IF doassim

  END SUBROUTINE PDAFomi_init_obs_l



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
  SUBROUTINE PDAFomi_init_obsvar_l(thisobs_l, thisobs, meanvar_l, cnt_obs_l)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    REAL, INTENT(inout) :: meanvar_l         !< Mean variance
    INTEGER, INTENT(inout) :: cnt_obs_l      !< Observation counter

! Local variables
    INTEGER :: i        ! Counter


! ***********************************
! *** Compute local mean variance ***
! ***********************************

    doassim: IF (thisobs%doassim == 1) THEN

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

       ! Print debug information
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug init_obsvar_l: ', debug, 'thisobs_l%dim_obs_l', thisobs_l%dim_obs_l
          WRITE (*,*) '++ OMI-debug init_obsvar_l: ', debug, 'cnt_obs_l', cnt_obs_l
          WRITE (*,*) '++ OMI-debug init_obsvar_l: ', debug, 'thisobs_l%ivar_obs_l', thisobs_l%ivar_obs_l(:)
          WRITE (*,*) '++ OMI-debug init_obsvar_l: ', debug, 'meanvar_l', meanvar_l
       END IF

    END IF doassim

  END SUBROUTINE PDAFomi_init_obsvar_l



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
  SUBROUTINE PDAFomi_prodRinvA_l(thisobs_l, thisobs, nobs_all, ncols, &
       A_l, C_l, verbose)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    INTEGER, INTENT(in) :: nobs_all          !< Dimension of local obs. vector (all obs. types)
    INTEGER, INTENT(in) :: ncols             !< Rank of initial covariance matrix
    REAL, INTENT(inout) :: A_l(:, :)         !< Input matrix (thisobs_l%dim_obs_l, ncols)
    REAL, INTENT(out)   :: C_l(:, :)         !< Output matrix (thisobs_l%dim_obs_l, ncols)
    INTEGER, INTENT(in) :: verbose           !< Verbosity flag


! *** local variables ***
    INTEGER :: i, j                    ! Index of observation component
    REAL, ALLOCATABLE :: weight(:)     ! Localization weights
    INTEGER :: idummy                  ! Dummy to access nobs_all
    INTEGER :: off                     ! row offset in A_l and C_l


! **********************
! *** INITIALIZATION ***
! **********************

    doassim: IF (thisobs%doassim == 1) THEN

       ! Initialize dummy to prevent compiler warning
       idummy = nobs_all

       ! Initialize offset
       off = thisobs_l%off_obs_l

       ! Screen output
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, &
               'PDAFomi_prodrinva_l -- START Multiply with inverse R and and apply localization'
          WRITE (*,*) '++ OMI-debug prodrinva_l:    ', debug, '  thisobs_l%locweight', thisobs_l%locweight
          WRITE (*,*) '++ OMI-debug prodRinvA_l:    ', debug, 'thisobs%dim_obs_f', thisobs_l%dim_obs_l
          WRITE (*,*) '++ OMI-debug prodRinvA_l:    ', debug, 'thisobs%ivar_obs_f', thisobs_l%ivar_obs_l
          WRITE (*,*) '++ OMI-debug prodRinvA_l:    ', debug, 'Input matrix A_l', A_l
       END IF

       IF (verbose == 1) THEN
          WRITE (*, '(a, 5x, a, 1x)') &
               'PDAFomi', '--- Domain localization'
          WRITE (*, '(a, 8x, a, 1x, es11.3)') &
               'PDAFomi', '--- Support radius', thisobs_l%lradius
       ENDIF


! ***********************************************
! *** Apply a weight matrix with correlations ***
! *** of compact support to matrix A or the   ***
! *** observation error covariance matrix.    ***
! ***********************************************

       ! *** Initialize weight array

       ALLOCATE(weight(thisobs_l%dim_obs_l))

       CALL PDAFomi_weights_l(verbose, thisobs_l%dim_obs_l, ncols, thisobs_l%locweight, &
            thisobs_l%lradius, thisobs_l%sradius, &
            A_l, thisobs_l%ivar_obs_l, thisobs_l%distance_l, weight)


       ! *** Handling of special weighting types ***

       lw2: IF (thisobs_l%locweight ==6 ) THEN
          ! Use square-root of 5th-order polynomial on A

          DO i = 1, thisobs_l%dim_obs_l
             ! Check if weight >0 (Could be <0 due to numerical precision)
             IF (weight(i) > 0.0) THEN
                weight(i) = SQRT(weight(i))
             ELSE
                weight(i) = 0.0
             END IF
          END DO
       END IF lw2


       ! *** Apply weight

       doweighting: IF (thisobs_l%locweight >= 5) THEN

          ! *** Apply weight to matrix A
          DO j = 1, ncols
             DO i = 1, thisobs_l%dim_obs_l
                A_l(i+off, j) = weight(i) * A_l(i+off, j)
             END DO
          END DO

          ! ***       -1
          ! ***  C = R   A 
          DO j = 1, ncols
             DO i = 1, thisobs_l%dim_obs_l
                C_l(i+off, j) = thisobs_l%ivar_obs_l(i) * A_l(i+off, j)
             END DO
          END DO
  
       ELSE doweighting

          ! *** Apply weight to matrix R only
          DO j = 1, ncols
             DO i = 1, thisobs_l%dim_obs_l
                C_l(i+off, j) = thisobs_l%ivar_obs_l(i) * weight(i) * A_l(i+off, j)
             END DO
          END DO
     
       END IF doweighting

       ! *** Clean up ***

       DEALLOCATE(weight)

       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_prodrinva_l -- END'

    ENDIF doassim

  END SUBROUTINE PDAFomi_prodRinvA_l



!-------------------------------------------------------------------------------
!> Compute local likelihood for an ensemble member
!!
!! The routine is called during the analysis step
!! of the localized NETF.
!! It has to compute the likelihood of the
!! ensemble according to the difference from the
!! observation (residual) and the error distribution
!! of the observations.
!!
!! In addition, a localizing weighting of the 
!! inverse of R by expotential decrease or a 5-th order 
!! polynomial of compact support can be applied. This is 
!! defined by the variables 'locweight', 'local_range, 
!! 'local_range2' and 'srange' in the main program.
!!
!! In general this routine is similar to the routine
!! prodRinvA_l used for ensemble square root Kalman
!! filters. As an addition to this routine, we here have
!! to evaluate the likelihood weight according the
!! assumed observation error statistics.
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
!! * 2020-03 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_likelihood_l(thisobs_l, thisobs, resid_l, lhood_l, verbose)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    REAL, INTENT(inout) :: resid_l(:)        !< Input vector of residuum
    REAL, INTENT(inout) :: lhood_l           !< Output vector - log likelihood
    INTEGER, INTENT(in) :: verbose           !< Verbosity flag


! *** local variables ***
    INTEGER :: i                          ! Index of observation component
    REAL, ALLOCATABLE :: weight(:)        ! Localization weights
    REAL, ALLOCATABLE :: resid_obs(:,:)   ! Array for a single row of resid_l
    REAL, ALLOCATABLE :: Rinvresid_l(:)   ! R^-1 times residual
    REAL :: lhood_one                     ! Likelihood for this observation


    doassim: IF (thisobs%doassim == 1) THEN

! **********************
! *** INITIALIZATION ***
! **********************

       ! Screen output
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, &
               'PDAFomi_likelihood_l -- START localization and likelihood, member', obs_member
          WRITE (*,*) '++ OMI-debug likelihood_l:  ', debug, '  thisobs_l%locweight', thisobs_l%locweight
       END IF

       ! Screen output
       IF (verbose == 1) THEN
          IF (thisobs%obs_err_type==0) THEN
             WRITE (*, '(a, 8x, a)') &
                  'PDAFomi', '--- Assume Gaussian observation errors'
          ELSE
             WRITE (*, '(a, 8x, a)') &
                  'PDAFomi', '--- Assume double-exponential observation errors'
          END IF
          WRITE (*, '(a, 8x, a, 1x)') &
               'PDAFomi', '--- Domain localization'
          WRITE (*, '(a, 12x, a, 1x, f12.2)') &
               'PDAFomi', '--- Local influence radius', thisobs_l%lradius
       ENDIF


! ***********************************************
! *** Before computing the likelihood, apply  ***
! *** the localization weight and scale by    ***
! *** observation error variance              ***
! ***                           -1            ***
! ***       Rinvresid = Weight R  resid       ***
! ***                                         ***
! *** Apply a weight matrix with correlations ***
! *** of compact support to residual or the   ***
! *** observation error variance.             ***
! ***********************************************

       ! *** Initialize weight array

       ALLOCATE(weight(thisobs_l%dim_obs_l))
       ALLOCATE(resid_obs(thisobs_l%dim_obs_l,1))

       resid_obs(:,1) = resid_l(:)

       CALL PDAFomi_weights_l(verbose, thisobs_l%dim_obs_l, 1, thisobs_l%locweight, &
            thisobs_l%lradius, thisobs_l%sradius, &
            resid_obs, thisobs_l%ivar_obs_l, thisobs_l%distance_l, weight)

       DEALLOCATE(resid_obs)


       ! *** Handling of special weighting types ***

       lw2: IF (thisobs_l%locweight ==6 ) THEN
          ! Use square-root of 5th-order polynomial on A

          DO i = 1, thisobs_l%dim_obs_l
             ! Check if weight >0 (Could be <0 due to numerical precision)
             IF (weight(i) > 0.0) THEN
                weight(i) = SQRT(weight(i))
             ELSE
                weight(i) = 0.0
             END IF
          END DO
       END IF lw2


       ! *** Apply weight

       ALLOCATE(Rinvresid_l(thisobs_l%dim_obs_l))

       DO i = 1, thisobs_l%dim_obs_l
          Rinvresid_l(i) = thisobs_l%ivar_obs_l(i) * weight(i) * resid_l(i)
       END DO


! ********************************
! *** Compute local likelihood ***
! ********************************

       IF (thisobs%obs_err_type == 0) THEN

          ! Gaussian errors
          ! Calculate exp(-0.5*resid^T*R^-1*resid)

          ! Transform pack to log likelihood to increment its values
          IF (lhood_l>0.0) lhood_l = - LOG(lhood_l)

          CALL dgemv('t', thisobs_l%dim_obs_l, 1, 0.5, resid_l, &
               thisobs_l%dim_obs_l, Rinvresid_l, 1, 0.0, lhood_one, 1)

          lhood_l = EXP(-(lhood_l + lhood_one))

       ELSE

          ! Double-exponential errors
          ! Calculate exp(-SUM(ABS(resid)))

          ! Transform pack to log likelihood to increment its values
          IF (lhood_l>0.0) lhood_l = - LOG(lhood_l)

          lhood_one = 0.0
          DO i = 1, thisobs_l%dim_obs_l
             lhood_one = lhood_one + ABS(Rinvresid_l(i))
          END DO

          lhood_l = EXP(-(lhood_l + lhood_one))

       END IF

       ! *** Clean up ***

       DEALLOCATE(weight, Rinvresid_l)

       ! Screen output
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug likelihood_l:  ', debug, '  accumulated likelihood', lhood_l
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_likelihood_l -- END'
       END IF

    END IF doassim

  END SUBROUTINE PDAFomi_likelihood_l




!-------------------------------------------------------------------------------
!> Apply covariance localization
!!
!! This routine applies a localization matrix B
!! to the matrices HP and HPH^T of the localized EnKF.
!!
!! __Revision history:__
!! * 2020-03 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_localize_covar(thisobs, dim,  locweight, lradius, sradius, &
       coords, HP, HPH)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(in) :: thisobs    !< Data type with full observation
    INTEGER, INTENT(in) :: dim            !< State dimension
    INTEGER, INTENT(in) :: locweight      !< Localization weight type
    REAL, INTENT(in)    :: lradius        !< localization radius
    REAL, INTENT(in)    :: sradius        !< support radius for weight functions
    REAL, INTENT(in)    :: coords(:,:)    !< Coordinates of state vector elements
    REAL, INTENT(inout) :: HP(:, :)       !< Matrix HP, dimension (nobs, dim)
    REAL, INTENT(inout) :: HPH(:, :)      !< Matrix HPH, dimension (nobs, nobs)

! *** local variables ***
    INTEGER :: i, j          ! Index of observation component
    INTEGER :: ncoord        ! Number of coordinates
    REAL    :: distance      ! Distance between points in the domain 
    REAL    :: weight        ! Localization weight
    REAL, ALLOCATABLE :: weights(:) ! Localization weights array
    REAL    :: tmp(1,1)= 1.0 ! Temporary, but unused array
    INTEGER :: wtype         ! Type of weight function
    INTEGER :: rtype         ! Type of weight regulation
    REAL, ALLOCATABLE :: co(:), oc(:)


    doassim: IF (thisobs%doassim == 1) THEN

! **********************
! *** INITIALIZATION ***
! **********************

       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_localize_covar -- START'
          WRITE (*,*) '++ OMI-debug localize_covar:', debug, 'thisobs%off_obs_f', thisobs%off_obs_f
       END IF

       ! Screen output
       IF (screen > 0 .and. mype==0) THEN
          WRITE (*,'(a, 8x, a)') &
               'PDAFomi', '--- Apply covariance localization'
          WRITE (*, '(a, 12x, a, 1x, f12.2)') &
               'PDAFomi', '--- Local influence radius', lradius

          IF (locweight == 0) THEN
             WRITE (*, '(a, 12x, a)') &
                  'PDAFomi', '--- Use uniform weight'
          ELSE IF (locweight == 1) THEN
             WRITE (*, '(a, 12x, a)') &
                  'PDAFomi', '--- Use exponential distance-dependent weight'
          ELSE IF (locweight == 2) THEN
             WRITE (*, '(a, 12x, a)') &
                  'PDAFomi', '--- Use distance-dependent weight by 5th-order polynomial'
          END IF
       ENDIF

       ! Set ncoord locally for compact code
       ncoord = thisobs%ncoord


! **************************
! *** Apply localization ***
! **************************


       ! Set parameters for weight calculation
       IF (locweight == 0) THEN
          ! Uniform (unit) weighting
          wtype = 0
          rtype = 0
       ELSE IF (locweight == 1) THEN
          ! Exponential weighting
          wtype = 1
          rtype = 0
       ELSE IF (locweight == 2) THEN
          ! 5th-order polynomial (Gaspari&Cohn, 1999)
          wtype = 2
          rtype = 0
       END IF

       ALLOCATE(oc(ncoord))
       ALLOCATE(co(ncoord))


       ! *** Localize HP ***

       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug localize_covar:', debug, '  localize matrix HP'
       END IF

       ALLOCATE(weights(thisobs%dim_obs_f))

       DO i = 1, dim

          ! Initialize coordinate
          co(1:ncoord) = coords(1:thisobs%ncoord, i)

          DO j = 1, thisobs%dim_obs_f

             ! Initialize coordinate
             oc(1:ncoord) = thisobs%ocoord_f(1:thisobs%ncoord, j)

             ! Compute distance
             CALL PDAFomi_comp_dist2(thisobs, co, oc, distance, (i*j)-1)
             distance = SQRT(distance)

             ! Compute weight
             CALL PDAF_local_weight(wtype, rtype, lradius, sradius, distance, &
                  1, 1, tmp, 1.0, weights(j), 0)
          END DO

          IF (debug==i) THEN
             WRITE (*,*) '++ OMI-debug localize_covar:  ', debug, 'weights for row in HP', weights
          END IF

          DO j = 1, thisobs%dim_obs_f

             ! Apply localization
             HP(j + thisobs%off_obs_f, i) = weights(j) * HP(j + thisobs%off_obs_f, i)

          END DO
       END DO

       DEALLOCATE(weights)


       ! *** Localize HPH^T ***

       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug localize_covar:', debug, '  localize matrix HPH^T'
       END IF

       DO i = 1, thisobs%dim_obs_f

          ! Initialize coordinate
          co(1:ncoord) = thisobs%ocoord_f(1:thisobs%ncoord, i)

          DO j = 1, thisobs%dim_obs_f

             ! Initialize coordinate
             oc(1:ncoord) = thisobs%ocoord_f(1:thisobs%ncoord, j)

             ! Compute distance
             CALL PDAFomi_comp_dist2(thisobs, co, oc, distance, (i*j)-1)
             distance = SQRT(distance)

             ! Compute weight
             CALL PDAF_local_weight(wtype, rtype, lradius, sradius, distance, &
                  1, 1, tmp, 1.0, weight, 0)

             ! Apply localization
             HPH(j + thisobs%off_obs_f, i + thisobs%off_obs_f) = weight * HPH(j + thisobs%off_obs_f, i + thisobs%off_obs_f)

          END DO
       END DO

       ! clean up
       DEALLOCATE(co, oc)

       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_localize_covar -- END'

    END IF doassim

  END SUBROUTINE PDAFomi_localize_covar



!-------------------------------------------------------------------------------
!> Initialize local observation vector
!!
!! This routine has to initialize the part of the 
!! overall local observation vector corresponding
!! to the current observation type. The offset of
!! the current observation type in the local obs.
!! vector is given by OFFSET_OBS_l_ALL.
!! This routine is both used directly to initialize
!! the local part of a vector in observation space
!! (PDAFomi_g2l_obs) and it's called from 
!! PDAFomi_init_obs_l to initialize the local part
!! of the vector of observations and the corresponding
!! vector of variances.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_g2l_obs_internal(thisobs_l, obs_f_one, offset_obs_l_all, obs_l_all)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    REAL, INTENT(in) :: obs_f_one(:)         !< Full obs. vector of current obs. type (nobs_f_one)
    REAL, INTENT(inout) :: obs_l_all(:)      !< Local observation vector for all variables (nobs_l_all)
    INTEGER, INTENT(in) :: offset_obs_l_all  !< Offset of current observation in obs_l_all and ivar_l_all

! *** Local variables ***
    INTEGER :: i  ! Counter


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

    ! Local observations
    DO i = 1, thisobs_l%dim_obs_l
       obs_l_all(i+offset_obs_l_all) = obs_f_one(thisobs_l%id_obs_l(i))
    ENDDO

    ! Print debug information
    IF (debug>0) THEN
       WRITE (*,*) '++ OMI-debug g2l_obs:       ', debug, '  thisobs%id_obs_l', thisobs_l%id_obs_l
       WRITE (*,*) '++ OMI-debug g2l_obs:       ', debug, '  obs_l', &
            obs_l_all(1+offset_obs_l_all:offset_obs_l_all+thisobs_l%dim_obs_l)
    END IF

  END SUBROUTINE PDAFomi_g2l_obs_internal




!-------------------------------------------------------------------------------
!> Compute square distance between two locations
!!
!! This routine computes the distance between two locations.
!! The computation can be for Cartesian grids with and without
!! periodicity and for geographic coordinates. For Cartesian
!! grids, the coordinats can be in any unit, while geographic
!! coordinates must be provided in radians and the resulting
!! distance will be in meters.
!!
!! Choices for distance computation - disttype:
!! 0: Cartesian distance in ncoord dimensions
!! 1: Cartesian distance in ncoord dimensions with periodicity
!!    (Needs specification of domsize(ncoord))
!! 2: Aproximate geographic distance with horizontal coordinates in radians (-pi/+pi)
!! 3: Geographic distance computation using haversine formula
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_comp_dist2(thisobs, coordsA, coordsB, distance2, verbose)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(in) :: thisobs    !< Data type with full observation
    REAL, INTENT(in) :: coordsA(:)           !< Coordinates of current analysis domain (ncoord)
    REAL, INTENT(in) :: coordsB(:)           !< Coordinates of observation (ncoord)
    REAL, INTENT(out) :: distance2           !< Squared distance
    INTEGER, INTENT(in) :: verbose           ! Control screen output

! *** Local variables ***
    INTEGER :: k                    ! Counters
    REAL :: dists(thisobs%ncoord)   ! Distance vector between analysis point and observation
    REAL :: slon, slat              ! sine of distance in longitude or latitude
    INTEGER :: domsize              ! Flag whether domainsize is set


! ************************
! *** Compute distance ***
! ************************

    IF (.NOT.ALLOCATED(thisobs%domainsize)) THEN
       domsize = 0
    ELSE
       domsize = 1
    END IF

    norm: IF ((thisobs%disttype==0) .OR.(thisobs%disttype==1 .AND. domsize==0)) THEN

       ! *** Compute Cartesian distance ***

       IF (debug>0 .AND. verbose==0) THEN
          WRITE (*,*) '++ OMI-debug comp_dist2:    ', debug, '  compute Cartesian distance'
       END IF

       ! Distance per direction
       DO k = 1, thisobs%ncoord
          dists(k) = ABS(coordsA(k) - coordsB(k))
       END DO

       ! full squared distance
       distance2 = 0.0
       DO k = 1, thisobs%ncoord
          distance2 = distance2 + dists(k)*dists(k)
       END DO

    ELSEIF (thisobs%disttype==1 .AND. domsize==1) THEN norm

       ! *** Compute periodic Cartesian distance ***

       IF (debug>0 .AND. verbose==0) THEN
          WRITE (*,*) '++ OMI-debug comp_dist2:    ', debug, '  compute periodic Cartesian distance'
       END IF

       ! Distance per direction
       DO k = 1, thisobs%ncoord
          IF (thisobs%domainsize(k)<=0.0) THEN 
             dists(k) = ABS(coordsA(k) - coordsB(k))
          ELSE
             dists(k) = MIN(ABS(coordsA(k) - coordsB(k)), &
                  ABS(ABS(coordsA(k) - coordsB(k))-thisobs%domainsize(k)))
          END IF
       END DO

       ! full squared distance
       distance2 = 0.0
       DO k = 1, thisobs%ncoord
          distance2 = distance2 + dists(k)*dists(k)
       END DO

    ELSEIF (thisobs%disttype==2) THEN norm

       ! *** Compute distance from geographic coordinates ***

       IF (debug>0 .AND. verbose==0) THEN
          WRITE (*,*) '++ OMI-debug comp_dist2:    ', debug, '  compute geographic distance'
       END IF

       ! approximate distances in longitude and latitude
       dists(1) = r_earth * MIN( ABS(coordsA(1) - coordsB(1))* COS(coordsA(2)), &
            ABS(ABS(coordsA(1) - coordsB(1)) - 2.0*pi) * COS(coordsA(2)))
       dists(2) = r_earth * ABS(coordsA(2) - coordsB(2))
       IF (thisobs%ncoord>2) dists(3) = ABS(coordsA(3) - coordsB(3))

       ! full squared distance in meters
       distance2 = 0.0
       DO k = 1, thisobs%ncoord
          distance2 = distance2 + dists(k)*dists(k)
       END DO

    ELSEIF (thisobs%disttype==3) THEN norm

       ! *** Compute distance from geographic coordinates with haversine formula ***

       IF (debug>0 .AND. verbose==0) THEN
          WRITE (*,*) '++ OMI-debug comp_dist2:    ', debug, &
               '  compute geographic distance using haversine function'
       END IF

       slon = SIN((coordsA(1) - coordsB(1))/2)
       slat = SIN((coordsA(2) - coordsB(2))/2)

       dists(2) = SQRT(slat*slat + COS(coordsA(2))*COS(coordsB(2))*slon*slon)
       IF (dists(2)<=1.0) THEN
          dists(2) = 2.0 * r_earth* ASIN(dists(2))
       ELSE
          dists(2) = r_earth* pi
       END IF

       IF (thisobs%ncoord>2) dists(3) = ABS(coordsA(3) - coordsB(3))

       ! full squared distance in meters
       distance2 = 0.0
       DO k = 2, thisobs%ncoord
          distance2 = distance2 + dists(k)*dists(k)
       END DO

    END IF norm

  END SUBROUTINE PDAFomi_comp_dist2




!-------------------------------------------------------------------------------
!> Compute weight vector for localization
!!
!! The routine computes a weight vector according to the
!! distances of observations from the local analysis
!! domain.
!!
!! __Revision history:__
!! * 2020-03 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_weights_l(verbose, nobs_l, ncols, locweight, lradius, sradius, &
       matA, ivar_obs_l, dist_l, weight_l)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: verbose        !< Verbosity flag
    INTEGER, INTENT(in) :: nobs_l         !< Number of local observations
    INTEGER, INTENT(in) :: ncols          !< 
    INTEGER, INTENT(in) :: locweight      !< Localization weight type
    REAL, INTENT(in)    :: lradius        !< localization radius
    REAL, INTENT(in)    :: sradius        !< support radius for weight functions
    REAL, INTENT(in)    :: matA(:,:)      !< 
    REAL, INTENT(in)    :: ivar_obs_l(:)  !< Local vector of inverse obs. variances (nobs_l)
    REAL, INTENT(in)    :: dist_l(:)      !< Local vector of obs. distances (nobs_l)
    REAL, INTENT(out) :: weight_l(:)      !< Output: vector of weights 

! *** local variables ***
    INTEGER :: i             ! Index of observation component
    INTEGER :: verbose_w     ! Verbosity flag for weight computation
    INTEGER :: wtype         ! Type of weight function
    INTEGER :: rtype         ! Type of weight regulation
    REAL    :: var_obs_l     ! Variance of observation error
    REAL, ALLOCATABLE :: A_obs(:,:)    ! Array for a single row of matA
    

! **********************
! *** INITIALIZATION ***
! **********************

    ! Screen output
    IF (verbose == 1) THEN
       IF (locweight == 1 .OR. locweight == 2 .OR. locweight == 3 &
            .OR. locweight == 4) THEN
          WRITE (*, '(a, 8x, a)') &
               'PDAFomi', '--- Use distance-dependent weight for observation errors'

          IF (locweight == 3) THEN
             WRITE (*, '(a, 8x, a)') &
                  'PDAFomi', '--- Use regulated weight with mean error variance'
          ELSE IF (locweight == 4) THEN
             WRITE (*, '(a, 8x, a)') &
                  'PDAFomi', '--- Use regulated weight with single-point error variance'
          END IF
       ELSE IF (locweight == 5 .OR. locweight == 6 .OR. locweight == 7) THEN
          WRITE (*, '(a, 8x, a)') &
               'PDAFomi', '--- Use distance-dependent weight for observed ensemble'
       END IF
    ENDIF


! *** Initialize weight vector

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
       ALLOCATE(A_obs(1, ncols))
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
               nobs_l, ncols, matA, var_obs_l, weight_l(i), verbose_w)
       ELSE
          ! Regulated weight using variance at single observation point
          A_obs(1,:) = matA(i,:)
          CALL PDAF_local_weight(wtype, rtype, lradius, sradius, dist_l(i), &
               1, ncols, A_obs, var_obs_l, weight_l(i), verbose_w)
       END IF
    END DO

    ! Print debug information
    IF (debug>0) THEN
       WRITE (*,*) '++ OMI-debug weights_l:     ', debug, 'thisobs_l%dim_obs_l', nobs_l
       WRITE (*,*) '++ OMI-debug weights_l:     ', debug, 'thisobs%distance_l', dist_l(1:nobs_l)
       WRITE (*,*) '++ OMI-debug weights_l:     ', debug, 'weight_l', weight_l(1:nobs_l)
       WRITE (*,*) '++ OMI-debug weights_l:     ', debug, 'thisobs_l%ivar_obs_l', ivar_obs_l(1:nobs_l)
    END IF
  
    IF (locweight == 4) DEALLOCATE(A_obs)

  END SUBROUTINE PDAFomi_weights_l



!-------------------------------------------------------------------------------
!> Deallocate arrays in observation type
!!
!! This routine deallocates arrays in the data type THISOBS.
!! The routine mainly operates on the full observation type. 
!! It is included here to avoid cross-dependences between
!! PDAFomi_obs_f and PDAFomi_obs_l.
!!
!! The routine is called by all filter processes.
!!
!! __Revision history:__
!! * 2019-10 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_deallocate_obs(thisobs)

    USE PDAFomi_obs_f, &
         ONLY: obs_f, n_obstypes, obscnt, offset_obs, obs_f_all, &
         offset_obs_g

    IMPLICIT NONE

! *** Arguments
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation

   ! *** Perform deallocation ***

    IF (ALLOCATED(thisobs%obs_f)) DEALLOCATE(thisobs%obs_f)
    IF (ALLOCATED(thisobs%ocoord_f)) DEALLOCATE(thisobs%ocoord_f)
    IF (ALLOCATED(thisobs%id_obs_p)) DEALLOCATE(thisobs%id_obs_p)
    IF (ALLOCATED(thisobs%ivar_obs_f)) DEALLOCATE(thisobs%ivar_obs_f)
    IF (ALLOCATED(thisobs%icoeff_p)) DEALLOCATE(thisobs%icoeff_p)
    IF (ALLOCATED(thisobs%domainsize)) DEALLOCATE(thisobs%domainsize)
    IF (ALLOCATED(thisobs%id_obs_f_lim)) DEALLOCATE(thisobs%id_obs_f_lim)
    IF (ALLOCATED(obs_f_all)) DEALLOCATE(obs_f_all)

    ! Reset counters over all observation types
    n_obstypes = 0
    obscnt = 0
    offset_obs = 0
    offset_obs_g = 0

    ! Reset flag for first call to local observations
    firstobs = 0

  END SUBROUTINE PDAFomi_deallocate_obs



!-------------------------------------------------------------------------------
!> Deallocate arrays in all observation types
!!
!! This routine deallocates arrays in all observation types.
!! The routine is only called internally in PDAF. 
!!
!! The routine is called by all filter processes.
!!
!! __Revision history:__
!! * 2021-04 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_dealloc()

    USE PDAFomi_obs_f, &
         ONLY: obs_f, n_obstypes, obscnt, offset_obs, obs_f_all, &
         offset_obs_g

! *** Local variables
    INTEGER :: i

    ! *** Perform deallocation ***

    IF (n_obstypes>0) THEN
       DO i=1, n_obstypes
          IF (ALLOCATED(obs_f_all(i)%ptr%obs_f)) DEALLOCATE(obs_f_all(i)%ptr%obs_f)
          IF (ALLOCATED(obs_f_all(i)%ptr%ocoord_f)) DEALLOCATE(obs_f_all(i)%ptr%ocoord_f)
          IF (ALLOCATED(obs_f_all(i)%ptr%id_obs_p)) DEALLOCATE(obs_f_all(i)%ptr%id_obs_p)
          IF (ALLOCATED(obs_f_all(i)%ptr%ivar_obs_f)) DEALLOCATE(obs_f_all(i)%ptr%ivar_obs_f)
          IF (ALLOCATED(obs_f_all(i)%ptr%icoeff_p)) DEALLOCATE(obs_f_all(i)%ptr%icoeff_p)
          IF (ALLOCATED(obs_f_all(i)%ptr%ivar_obs_f)) DEALLOCATE(obs_f_all(i)%ptr%ivar_obs_f)
          IF (ALLOCATED(obs_f_all(i)%ptr%icoeff_p)) DEALLOCATE(obs_f_all(i)%ptr%icoeff_p)
          IF (ALLOCATED(obs_f_all(i)%ptr%domainsize)) DEALLOCATE(obs_f_all(i)%ptr%domainsize)
          IF (ALLOCATED(obs_f_all(i)%ptr%id_obs_f_lim)) DEALLOCATE(obs_f_all(i)%ptr%id_obs_f_lim)
       END DO
       IF (ALLOCATED(obs_f_all)) DEALLOCATE(obs_f_all)

       ! Reset counters over all observation types
       n_obstypes = 0
       obscnt = 0
       offset_obs = 0
       offset_obs_g = 0

       ! Reset flag for first call to local observations
       firstobs = 0
    END IF 
   
  END SUBROUTINE PDAFomi_dealloc

END MODULE PDAFomi_obs_l
