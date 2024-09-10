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
!$Id: PDAFomi_obs_l.F90 1147 2023-03-12 16:14:34Z lnerger $

!> PDAF-OMI routines for local observations
!!
!! This module contains generic routines for several observation-related
!! operations for local filters. The routines are
!!
!! * PDAFomi_set_debug_flag \n
!!        Set or unset the debugging flag for PDAFomi routines
!! * PDAFomi_init_dim_obs_l_old \n
!!        Initialize dimension of local obs. vetor and arrays for
!!        local observations
!! * PDAFomi_cnt_dim_obs_l_old \n
!!        Set dimension of local obs. vector with isotropic localization
!! * PDAFomi_cnt_dim_obs_l_noniso_old \n
!!        Set dimension of local obs. vector with nonisotropic localization
!! * PDAFomi_init_obsarrays_l_old \n
!!        Initialize arrays for the index of a local observation in 
!!        the full observation vector and its corresponding distance.
!! * PDAFomi_init_obsarrays_l_noniso_old \n
!!        Initialize arrays for the index of a local observation in 
!!        the full observation vector and its corresponding distance
!!        with onoisotrppic localization.
!! * PDAFomi_g2l_obs \n
!!        Initialize local observation vector from full observation vector
!! * PDAFomi_init_obs_l \n
!!        Initialize the local vector of observations
!! * PDAFomi_prodRinvA_l \n
!!        Multiply an intermediate matrix of the local filter analysis
!!        with the inverse of the observation error covariance matrix
!!        and apply observation localization
!! * PDAFomi_prodRinvA_hyb_l \n
!!        Multiply an intermediate matrix of the local filter analysis
!!        with the inverse of the observation error covariance matrix
!!        and apply observation localization. In addition apply the 
!!        hybrid weight
!! * PDAFomi_init_obsvar_l \n
!!        Compute mean observation error variance
!! * PDAFomi_likelihood_l \n
!!        Compute local likelihood for an ensemble member
!! * PDAFomi_likelihood_hyb_l \n
!!        Compute local likelihood for an ensemble member taking into
!!        account a hybrid weight for tempering
!! * PDAFomi_localize_covar_iso \n
!!        Apply covariance isotropic localization in LEnKF
!! * PDAFomi_localize_covar_noniso \n
!!        Apply non-isotropic covariance localization in LEnKF
!! * PDAFomi_g2l_obs_internal \n
!!        Internal routine to initialize local observation vector from full
!!        observation vector (used by PDAFomi_init_obs_l and PDAFomi_g2l_obs)
!! * PDAFomi_comp_dist2 \n
!!        Compute squared distance
!! * PDAFomi_check_dist2 \n
!!        Compute and check distance for isotropic localization
!! * PDAFomi_check_dist2_noniso \n
!!        Compute and check distance for non-isotropic localization
!! * PDAFomi_weights_l \n
!!        Compute a vector of localization weights
!! * PDAFomi_deallocate_obs \n
!!        Deallocate arrays in observation type
!! * PDAFomi_dealloc \n
!!        Deallocate arrays in all observation types
!! * PDAFomi_omit_by_innovation_l \
!!        Exclude observations if innovation is too large (thisobs%inno_exclude)
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE PDAFomi_obs_l

  USE PDAFomi_obs_f, ONLY: obs_f, r_earth, pi, debug, n_obstypes, error
  USE PDAF_mod_filter, ONLY: screen, obs_member
  USE PDAF_mod_filtermpi, ONLY: mype, npes_filter

  IMPLICIT NONE
  SAVE

! *** Module internal variables

  ! Data type to define the local observations by internally shared variables of the module
  TYPE obs_l
     INTEGER :: dim_obs_l                 !< number of local observations
     INTEGER :: off_obs_l                 !< Offset of this observation in overall local obs. vector
     INTEGER, ALLOCATABLE :: id_obs_l(:)  !< Indices of local observations in full obs. vector 
     REAL, ALLOCATABLE :: distance_l(:)   !< Distances of local observations
     REAL, ALLOCATABLE :: cradius_l(:)    !< directional cut-off radii of local observations
     REAL, ALLOCATABLE :: sradius_l(:)    !< directional support radii of local observations
     REAL, ALLOCATABLE :: ivar_obs_l(:)   !< Inverse variance of local observations
     REAL, ALLOCATABLE :: dist_l_v(:)     !< Vertical distances of local observations
     INTEGER :: locweight                 !< Specify localization function
     INTEGER :: locweight_v=0             !< Specify localization function for vertical direction
     INTEGER :: nradii                    !< Length of CRADIUS and SRADIUS
     REAL, ALLOCATABLE :: cradius(:)      !< Localization cut-off radius (single value or vector)
     REAL, ALLOCATABLE :: sradius(:)      !< support radius for localization function (single value or vector)
  END TYPE obs_l

  TYPE obs_arr_l                          ! Type for pointer array over all observation types
     TYPE(obs_l), POINTER :: ptr
  END TYPE obs_arr_l

  TYPE(obs_arr_l), ALLOCATABLE :: obs_l_all(:) ! Declare pointer array

  INTEGER :: firstobs = 0                 ! Flag for very first call to init_dim_obs_l
  INTEGER :: offset_obs_l = 0             ! offset of current observation in overall local obs. vector

!$OMP THREADPRIVATE(obs_l_all, firstobs, offset_obs_l)

  INTERFACE PDAFomi_init_dim_obs_l_old
     MODULE PROCEDURE PDAFomi_init_dim_obs_l_iso_old
     MODULE PROCEDURE PDAFomi_init_dim_obs_l_noniso_old
     MODULE PROCEDURE PDAFomi_init_dim_obs_l_noniso_locweights_old
  END INTERFACE

  INTERFACE PDAFomi_localize_covar
     MODULE PROCEDURE PDAFomi_localize_covar_iso
     MODULE PROCEDURE PDAFomi_localize_covar_noniso
  END INTERFACE


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
!! with coordinates COORD_l and localization cut-off radius CRADIUS.
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
  SUBROUTINE PDAFomi_init_dim_obs_l_iso_old(thisobs_l, thisobs, coords_l, locweight, cradius, &
       sradius, cnt_obs_l)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    TYPE(obs_l), TARGET, INTENT(inout) :: thisobs_l  !< Data type with local observation
    REAL, INTENT(in) :: coords_l(:)          !< Coordinates of current analysis domain
    INTEGER, INTENT(in) :: locweight         !< Type of localization function
    REAL, INTENT(in) :: cradius              !< Localization cut-off radius (single or vector)
    REAL, INTENT(in) :: sradius              !< Support radius of localization function (single or vector)
    INTEGER, INTENT(inout) :: cnt_obs_l      !< Local dimension of current observation vector

! *** Local variables ***
    REAL :: maxcoords_l, mincoords_l         ! Min/Max domain coordinates to check geographic coords
    REAL :: maxocoords_l, minocoords_l       ! Min/Max observation coordinates to check geographic coords


    doassim: IF (thisobs%doassim == 1) THEN


! ***********************************************
! *** Check offset in full observation vector ***
! ***********************************************

       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_init_dim_obs_l -- START'

       IF (thisobs%ncoord/=3 .AND. thisobs%disttype>=10) THEN
          WRITE (*,*) '+++++ ERROR PDAF-OMI: factorized 2+1D localization can only be used for thisobs%ncoord=3'
          error = 14
       END IF

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

       IF (ALLOCATED(thisobs_l%cradius)) DEALLOCATE(thisobs_l%cradius)
       ALLOCATE(thisobs_l%cradius(1))
       IF (ALLOCATED(thisobs_l%sradius)) DEALLOCATE(thisobs_l%sradius)
       ALLOCATE(thisobs_l%sradius(1))

       thisobs_l%nradii = 1
       thisobs_l%cradius(1) = cradius
       thisobs_l%sradius(1) = sradius

       ! Store offset
       thisobs_l%off_obs_l = offset_obs_l


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

          ! For geographic coordinates check whether their range is reasonable
          IF (thisobs%disttype==2 .OR. thisobs%disttype==3 .OR. thisobs%disttype==12 .OR. thisobs%disttype==13) THEN
             maxcoords_l = MAXVAL(coords_l)
             mincoords_l = MINVAL(coords_l)
             maxocoords_l = MAXVAL(thisobs%ocoord_f(1:2, :))
             minocoords_l = MINVAL(thisobs%ocoord_f(1:2, :))

             IF (maxcoords_l>2.0*pi .OR. mincoords_l<-pi .OR. maxocoords_l>2.0*pi .OR. minocoords_l<-pi) THEN
                WRITE (*,*) '++ OMI-debug init_dim_obs_l:', debug, &
                     '  WARNING: The unit for geographic coordinates is radian, thus range (0,2*pi) or (-pi,pi)!'
             END IF
          END IF
          WRITE (*,*) '++ OMI-debug init_dim_obs_l:', debug, &
               '  Note: Please ensure that coords_l and observation coordinates have the same unit'
       END IF

       CALL PDAFomi_cnt_dim_obs_l(thisobs_l, thisobs, coords_l)

       ! Store number of local module-type observations for output
       cnt_obs_l = cnt_obs_l + thisobs_l%dim_obs_l


! ************************************************************
! *** Initialize internal local arrays for local distances ***
! *** and indices of local obs. in full obs. vector        ***
! ************************************************************

       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, &
            '   PDAFomi_init_dim_obs_l -- initialize local observation arrays'

       ! Initialize ID_OBS_L and DISTANCE_L and increment offsets
       CALL PDAFomi_init_obsarrays_l(thisobs_l, thisobs, coords_l, offset_obs_l)

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

  END SUBROUTINE PDAFomi_init_dim_obs_l_iso_old




!-------------------------------------------------------------------------------
!> Set dimension of local obs. vector and local obs. arrays (non-isotropic)
!!
!! This routine sets the number of local observations for the
!! current observation type for the local analysis domain
!! with coordinates COORD_l and a vector of localization cut-off
!! radii CRADIUS.
!! Further the routine initializes arrays for the index of a
!! local observation in the full observation vector and its 
!! corresponding distance.
!! The operation are performed by calling the routines 
!! cnt_dim_obs_l and init_obsarrays_l.
!!
!! __Revision history:__
!! * 2024-02 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_init_dim_obs_l_noniso_old(thisobs_l, thisobs, coords_l, locweight, cradius, &
       sradius, cnt_obs_l)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    TYPE(obs_l), TARGET, INTENT(inout) :: thisobs_l  !< Data type with local observation
    REAL, INTENT(in) :: coords_l(:)          !< Coordinates of current analysis domain
    INTEGER, INTENT(in) :: locweight         !< Type of localization function
    REAL, INTENT(in) :: cradius(:)           !< Vector of localization cut-off radii
    REAL, INTENT(in) :: sradius(:)           !< Vector of support radii of localization function
    INTEGER, INTENT(inout) :: cnt_obs_l      !< Local dimension of current observation vector

! *** Local variables ***
    REAL :: maxcoords_l, mincoords_l         ! Min/Max domain coordinates to check geographic coords
    REAL :: maxocoords_l, minocoords_l       ! Min/Max observation coordinates to check geographic coords


    doassim: IF (thisobs%doassim == 1) THEN

! ***********************************************
! *** Check offset in full observation vector ***
! ***********************************************

       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_init_dim_obs_l_noniso -- START'

       ! Check consistency of dimensions
       IF (SIZE(cradius) /= thisobs%ncoord) THEN
          WRITE (*,*) '+++++ ERROR PDAF-OMI: non-isotropic localization: Size of CRADIUS /= thisobs%ncoord'
          error = 12
       END IF
       IF (SIZE(sradius) /= thisobs%ncoord) THEN
          WRITE (*,*) '+++++ ERROR PDAF-OMI: non-isotropic localization: Size of SRADIUS /= thisobs%ncoord'
          error = 13
       END IF
       IF (thisobs%ncoord/=3 .AND. thisobs%disttype>=10) THEN
          WRITE (*,*) '+++++ ERROR PDAF-OMI: factorized 2+1D localization can only be used for thisobs%ncoord=3'
          error = 14
       END IF

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

       ! Allocate vectors for localization radii and store their values
       IF (ALLOCATED(thisobs_l%cradius)) DEALLOCATE(thisobs_l%cradius)
       ALLOCATE(thisobs_l%cradius(thisobs%ncoord))
       IF (ALLOCATED(thisobs_l%sradius)) DEALLOCATE(thisobs_l%sradius)
       ALLOCATE(thisobs_l%sradius(thisobs%ncoord))

       thisobs_l%nradii = thisobs%ncoord
       thisobs_l%cradius(:) = cradius(:)
       thisobs_l%sradius(:) = sradius(:)

       ! Store offset
       thisobs_l%off_obs_l = offset_obs_l


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


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, &
               '   PDAFomi_init_dim_obs_l_noniso -- count local observations'
          IF (thisobs%obsid == firstobs) THEN
             WRITE (*,*) '++ OMI-debug init_dim_obs_l_noniso:', debug, '  Re-init dim_obs_l=0'
          END IF
          WRITE (*,*) '++ OMI-debug init_dim_obs_l_noniso:', debug, '  coords_l', coords_l

          ! For geographic coordinates check whether their range is reasonable
          IF (thisobs%disttype==2 .OR. thisobs%disttype==3 .OR. thisobs%disttype==12 .OR. thisobs%disttype==13) THEN
             maxcoords_l = MAXVAL(coords_l)
             mincoords_l = MINVAL(coords_l)
             maxocoords_l = MAXVAL(thisobs%ocoord_f(1:2, :))
             minocoords_l = MINVAL(thisobs%ocoord_f(1:2, :))

             IF (maxcoords_l>2.0*pi .OR. mincoords_l<-pi .OR. maxocoords_l>2.0*pi .OR. minocoords_l<-pi) THEN
                WRITE (*,*) '++ OMI-debug init_dim_obs_l_noniso:', debug, &
                     '  WARNING: The unit for geographic coordinates is radian, thus range (0,2*pi) or (-pi,pi)!'
             END IF
          END IF
          WRITE (*,*) '++ OMI-debug init_dim_obs_l_noniso:', debug, &
               '  Note: Please ensure that coords_l and observation coordinates have the same unit'
       END IF

       CALL PDAFomi_cnt_dim_obs_l_noniso(thisobs_l, thisobs, coords_l)

       ! Store number of local module-type observations for output
       cnt_obs_l = cnt_obs_l + thisobs_l%dim_obs_l


! ************************************************************
! *** Initialize internal local arrays for local distances ***
! *** and indices of local obs. in full obs. vector        ***
! ************************************************************

       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, &
            '   PDAFomi_init_dim_obs_l_noniso -- initialize local observation arrays'

       ! Initialize ID_OBS_L and DISTANCE_L and increment offsets
       CALL PDAFomi_init_obsarrays_l_noniso(thisobs_l, thisobs, coords_l, offset_obs_l)

       ! Print debug information
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug init_dim_obs_l_noniso:', debug, '  thisobs_l%dim_obs_l', thisobs_l%dim_obs_l
          IF (thisobs_l%dim_obs_l>0) THEN
             WRITE (*,*) '++ OMI-debug init_dim_obs_l_noniso:', debug, '  thisobs_l%id_obs_l', thisobs_l%id_obs_l
             WRITE (*,*) '++ OMI-debug init_dim_obs_l_noniso:', debug, '  thisobs_l%distance_l', thisobs_l%distance_l
          END IF
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_init_dim_obs_l_noniso -- END'
       END IF

    ELSE doassim

       cnt_obs_l = cnt_obs_l + 0

    END IF doassim

  END SUBROUTINE PDAFomi_init_dim_obs_l_noniso_old




!-------------------------------------------------------------------------------
!> Set dimension of local obs. vector and local obs. arrays
!!
!! This routine is a variant of PDAFomi_init_dim_obs_l_noniso with
!! support for a vector of localization weights. This is used
!! to specify different localization functions for the vertical and 
!! horizontal directions. The routine only stores the value of 
!! locweights(2) for the vertical and calls PDAFomi_init_dim_obs_l_iso.
!!
!! __Revision history:__
!! * 2024-04 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_init_dim_obs_l_noniso_locweights_old(thisobs_l, thisobs, coords_l, locweights, cradius, &
       sradius, cnt_obs_l)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    TYPE(obs_l), TARGET, INTENT(inout) :: thisobs_l  !< Data type with local observation
    REAL, INTENT(in) :: coords_l(:)          !< Coordinates of current analysis domain
    INTEGER, INTENT(in) :: locweights(:)     !< Types of localization function
    REAL, INTENT(in) :: cradius(:)           !< Vector of localization cut-off radii
    REAL, INTENT(in) :: sradius(:)           !< Vector of support radii of localization function
    INTEGER, INTENT(inout) :: cnt_obs_l      !< Local dimension of current observation vector


! *** Store vertical locweight and call standard routine

    IF (debug>0) THEN
       WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_init_dim_obs_l_noniso_locweights -- START'
       WRITE (*,*) '++ OMI-debug init_dim_obs_l_noniso_locweights:', debug, '  locweights', locweights
    END IF

    ! Check consistency of dimensions
    IF (SIZE(locweights) /= 2) THEN
       WRITE (*,*) '+++++ ERROR PDAF-OMI: Input for locweights in horizontal and vertical needs size 2'
       error = 15
    END IF
    IF (thisobs%ncoord /= 3) THEN
       WRITE (*,*) '+++++ WARNING PDAF-OMI: separate locweight for vertical is only utilized if thisobs%ncoord=3'
    END IF

    IF (thisobs%ncoord == 3) THEN
       ! locweight for the vertical is treated separately
       thisobs_l%locweight_v = locweights(2)
    END IF

    ! Call to usual routine that handles a single locweight setting
    CALL PDAFomi_init_dim_obs_l_noniso_old(thisobs_l, thisobs, coords_l, locweights(1), cradius, &
         sradius, cnt_obs_l)

    IF (debug>0) &
         WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_init_dim_obs_l_noniso_locweights -- END'

  END SUBROUTINE PDAFomi_init_dim_obs_l_noniso_locweights_old




!-------------------------------------------------------------------------------
!> Set dimension of local observation vector
!!
!! This routine sets the number of local observations for the
!! current observation type for the local analysis domain
!! with coordinates COORDS_L and localization cut-off radius CRADIUS.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_cnt_dim_obs_l(thisobs_l, thisobs, coords_l)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    REAL, INTENT(in) :: coords_l(:)          !< Coordinates of current analysis domain (thisobs%ncoord)

! *** Local variables ***
    INTEGER :: i, cnt       ! Counters
    REAL :: cradius         ! localization cut-off radius
    REAL :: distance2       ! squared distance
    REAL :: sradius         ! support radius
    LOGICAL :: checkdist    ! Flag whether distance nis not larger than cut-off radius
    

! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    ! Count local observations
    thisobs_l%dim_obs_l = 0
    cnt = 0

    IF (debug>0) THEN
       WRITE (*,*) '++ OMI-debug cnt_dim_obs_l: ', debug, '  thisobs%ncoord', thisobs%ncoord
       WRITE (*,*) '++ OMI-debug cnt_dim_obs_l: ', debug, '  thisobs_l%cradius', thisobs_l%cradius
       WRITE (*,*) '++ OMI-debug cnt_dim_obs_l: ', debug, '  Check for observations within radius'
    END IF

    scancount: DO i = 1, thisobs%dim_obs_f

       CALL PDAFomi_check_dist2(thisobs, thisobs_l, coords_l, thisobs%ocoord_f(1:thisobs%ncoord, i), distance2, &
            checkdist, i-1, cnt)

          ! If distance below limit, add observation to local domain
       IF (checkdist) THEN
          IF (debug>0) THEN
             WRITE (*,*) '++ OMI-debug cnt_dim_obs_l: ', debug, &
                  '  valid observation with coordinates', thisobs%ocoord_f(1:thisobs%ncoord, i)
          END IF
          
          thisobs_l%dim_obs_l = thisobs_l%dim_obs_l + 1
       END IF

    END DO scancount

  END SUBROUTINE PDAFomi_cnt_dim_obs_l




!-------------------------------------------------------------------------------
!> Set dimension of local observation vector for nonisotropic localization
!!
!! This routine sets the number of local observations for the
!! current observation type for the local analysis domain
!! with coordinates COORDS_L and localization cut-off radius CRADIUS.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_cnt_dim_obs_l_noniso(thisobs_l, thisobs, coords_l)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    REAL, INTENT(in) :: coords_l(:)          !< Coordinates of current analysis domain (thisobs%ncoord)

! *** Local variables ***
    INTEGER :: i, cnt       ! Counters
    REAL :: cradius         ! localization cut-off radius
    REAL :: distance2       ! squared distance
    REAL :: sradius         ! support radius
    LOGICAL :: checkdist    ! Flag whether distance nis not larger than cut-off radius
    REAL :: dists(thisobs%ncoord)   ! Distance vector between analysis point and observation


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    ! Count local observations
    thisobs_l%dim_obs_l = 0
    cnt = 0

    IF (debug>0) THEN
       WRITE (*,*) '++ OMI-debug cnt_dim_obs_l_noniso: ', debug, '  thisobs%ncoord', thisobs%ncoord
       WRITE (*,*) '++ OMI-debug cnt_dim_obs_l_noniso: ', debug, '  thisobs_l%cradius', thisobs_l%cradius
       WRITE (*,*) '++ OMI-debug cnt_dim_obs_l_noniso: ', debug, '  Check for observations within radius'
    END IF

    IF (thisobs_l%nradii==1) THEN

       ! 1D but with radius specified as array

       scancount: DO i = 1, thisobs%dim_obs_f

          CALL PDAFomi_check_dist2(thisobs, thisobs_l, coords_l, thisobs%ocoord_f(1:thisobs%ncoord, i), distance2, &
               checkdist, i-1, cnt)

          ! If distance below limit, add observation to local domain
          IF (checkdist) THEN
             IF (debug>0) THEN
                WRITE (*,*) '++ OMI-debug cnt_dim_obs_l_noniso: ', debug, &
                     '  valid observation with coordinates', thisobs%ocoord_f(1:thisobs%ncoord, i)
             END IF

             thisobs_l%dim_obs_l = thisobs_l%dim_obs_l + 1
          END IF

       END DO scancount

    ELSEIF (thisobs_l%nradii==2 .OR. thisobs_l%nradii==3) THEN

       ! Nonisotropic in 2 or 3 dimensions

       scancountB: DO i = 1, thisobs%dim_obs_f

          CALL PDAFomi_check_dist2_noniso(thisobs, thisobs_l, coords_l, thisobs%ocoord_f(1:thisobs%ncoord, i), distance2, &
               dists, cradius, sradius, checkdist, i-1, cnt)

          ! If distance below limit, add observation to local domain
          IF (checkdist .AND. debug>0) THEN
             WRITE (*,*) '++ OMI-debug cnt_dim_obs_l_noniso: ', debug, &
                  '  valid observation with coordinates', thisobs%ocoord_f(1:thisobs%ncoord, i)
             WRITE (*,*) '++ OMI-debug cnt_dim_obs_l_noniso: ', debug, &
                  '  valid observation distance, cradius, sradius', SQRT(distance2), cradius, sradius
          END IF

       END DO scancountB

       thisobs_l%dim_obs_l = cnt

    ELSE
       WRITE (*,*) '+++++ ERROR PDAF-OMI: nonisotropic localization is only possible in 1, 2 or 3 dimensions'
       error = 10
    END IF

  END SUBROUTINE PDAFomi_cnt_dim_obs_l_noniso



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
  SUBROUTINE PDAFomi_init_obsarrays_l(thisobs_l, thisobs, coords_l, off_obs_l_all)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    REAL, INTENT(in) :: coords_l(:)          !< Coordinates of current water column (thisobs%ncoord)
    INTEGER, INTENT(inout) :: off_obs_l_all  !< input: offset of current obs. in local obs. vector
                                             !< output: input + thisobs_l%dim_obs_l

! *** Local variables ***
    INTEGER :: i, off_obs   ! Counters
    REAL :: cradius         ! localization cut-off radius
    REAL :: distance2       ! squared distance
    REAL :: sradius         ! support radius
    LOGICAL :: checkdist    ! Flag whether distance nis not larger than cut-off radius


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    ! Allocate module-internal index array for indices in module-type observation vector
    IF (ALLOCATED(thisobs_l%id_obs_l)) DEALLOCATE(thisobs_l%id_obs_l)
    IF (ALLOCATED(thisobs_l%distance_l)) DEALLOCATE(thisobs_l%distance_l)
    IF (ALLOCATED(thisobs_l%cradius_l)) DEALLOCATE(thisobs_l%cradius_l)
    IF (ALLOCATED(thisobs_l%sradius_l)) DEALLOCATE(thisobs_l%sradius_l)
    IF (thisobs_l%dim_obs_l>0) THEN
       ALLOCATE(thisobs_l%id_obs_l(thisobs_l%dim_obs_l))
       ALLOCATE(thisobs_l%distance_l(thisobs_l%dim_obs_l))
       ALLOCATE(thisobs_l%cradius_l(thisobs_l%dim_obs_l))
       ALLOCATE(thisobs_l%sradius_l(thisobs_l%dim_obs_l))
    ELSE
       ALLOCATE(thisobs_l%id_obs_l(1))
       ALLOCATE(thisobs_l%distance_l(1))
       ALLOCATE(thisobs_l%cradius_l(1))
       ALLOCATE(thisobs_l%sradius_l(1))
    END IF

    off_obs = 0

    ! Count local observations and initialize index and distance arrays

    IF (thisobs_l%dim_obs_l>0) THEN

       scancount: DO i = 1, thisobs%dim_obs_f

          CALL PDAFomi_check_dist2(thisobs, thisobs_l, coords_l, thisobs%ocoord_f(1:thisobs%ncoord, i), distance2, &
               checkdist, i-1, off_obs)

          ! If distance below limit, add observation to local domain
          IF (checkdist) THEN
             ! For internal storage (use in prodRinvA_l)
             thisobs_l%id_obs_l(off_obs) = i                       ! node index
             thisobs_l%distance_l(off_obs) = SQRT(distance2)       ! distance
             thisobs_l%cradius_l(off_obs) = thisobs_l%cradius(1)   ! isotropic cut-off radius
             thisobs_l%sradius_l(off_obs) = thisobs_l%sradius(1)   ! isotropic support radius
          END IF
       END DO scancount

       ! Count overall local observations
       off_obs_l_all = off_obs_l_all + off_obs     ! dimension

    END IF

  END SUBROUTINE PDAFomi_init_obsarrays_l



!-------------------------------------------------------------------------------
!> Initialize local arrays for an observation for nonisotropic localization
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
!! * 2024-04 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_init_obsarrays_l_noniso(thisobs_l, thisobs, coords_l, off_obs_l_all)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    REAL, INTENT(in) :: coords_l(:)          !< Coordinates of current water column (thisobs%ncoord)
    INTEGER, INTENT(inout) :: off_obs_l_all  !< input: offset of current obs. in local obs. vector
                                             !< output: input + thisobs_l%dim_obs_l

! *** Local variables ***
    INTEGER :: i, off_obs   ! Counters
    REAL :: cradius         ! localization cut-off radius
    REAL :: distance2       ! squared distance
    REAL :: sradius         ! support radius
    LOGICAL :: checkdist    ! Flag whether distance nis not larger than cut-off radius
    REAL :: dists(thisobs%ncoord)   ! Distance vector between analysis point and observation


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    ! Allocate module-internal index array for indices in module-type observation vector
    IF (ALLOCATED(thisobs_l%id_obs_l)) DEALLOCATE(thisobs_l%id_obs_l)
    IF (ALLOCATED(thisobs_l%distance_l)) DEALLOCATE(thisobs_l%distance_l)
    IF (ALLOCATED(thisobs_l%cradius_l)) DEALLOCATE(thisobs_l%cradius_l)
    IF (ALLOCATED(thisobs_l%sradius_l)) DEALLOCATE(thisobs_l%sradius_l)
    IF (thisobs_l%dim_obs_l>0) THEN
       ALLOCATE(thisobs_l%id_obs_l(thisobs_l%dim_obs_l))
       ALLOCATE(thisobs_l%distance_l(thisobs_l%dim_obs_l))
       ALLOCATE(thisobs_l%cradius_l(thisobs_l%dim_obs_l))
       ALLOCATE(thisobs_l%sradius_l(thisobs_l%dim_obs_l))
    ELSE
       ALLOCATE(thisobs_l%id_obs_l(1))
       ALLOCATE(thisobs_l%distance_l(1))
       ALLOCATE(thisobs_l%cradius_l(1))
       ALLOCATE(thisobs_l%sradius_l(1))
    END IF

    ! Allocate array for vertical distance in factorized localization (disttype>=10)
    IF (thisobs_l%locweight_v>0) THEN
       IF (ALLOCATED(thisobs_l%dist_l_v)) DEALLOCATE(thisobs_l%dist_l_v)
       IF (thisobs_l%dim_obs_l>0) THEN
          ALLOCATE(thisobs_l%dist_l_v(thisobs_l%dim_obs_l))
       ELSE
          ALLOCATE(thisobs_l%dist_l_v(1))
       END IF
    END IF


    off_obs = 0

    ! Count local observations and initialize index and distance arrays
    IF (thisobs_l%nradii==1) THEN

       ! 1D but with radius specified as array

       IF (thisobs_l%dim_obs_l>0) THEN

          scancount: DO i = 1, thisobs%dim_obs_f

             CALL PDAFomi_check_dist2(thisobs, thisobs_l, coords_l, thisobs%ocoord_f(1:thisobs%ncoord, i), distance2, &
                  checkdist, i-1, off_obs)

             ! If distance below limit, add observation to local domain
             IF (checkdist) THEN
                ! For internal storage (use in prodRinvA_l)
                thisobs_l%id_obs_l(off_obs) = i                       ! node index
                thisobs_l%distance_l(off_obs) = SQRT(distance2)       ! distance
                thisobs_l%cradius_l(off_obs) = thisobs_l%cradius(1)   ! isotropic cut-off radius
                thisobs_l%sradius_l(off_obs) = thisobs_l%sradius(1)   ! isotropic support radius
             END IF
          END DO scancount

          ! Count overall local observations
          off_obs_l_all = off_obs_l_all + off_obs     ! dimension

       END IF

    ELSEIF (thisobs_l%nradii==2 .OR. thisobs_l%nradii==3) THEN

       ! Nonisotropic in 2 or 3 dimensions

       IF (thisobs_l%dim_obs_l>0) THEN

          off_obs = 0
          scancountB: DO i = 1, thisobs%dim_obs_f

             CALL PDAFomi_check_dist2_noniso(thisobs, thisobs_l, coords_l, thisobs%ocoord_f(1:thisobs%ncoord, i), distance2, &
                  dists, cradius, sradius, checkdist, i-1, off_obs)

             ! If distance below limit, add observation to local domain
             IF (checkdist) THEN
                ! For internal storage (use in prodRinvA_l)
                thisobs_l%id_obs_l(off_obs) = i             ! node index
                thisobs_l%distance_l(off_obs) = SQRT(distance2) ! distance
                thisobs_l%cradius_l(off_obs) = cradius          ! directional cut-off radius
                thisobs_l%sradius_l(off_obs) = sradius          ! directional support radius
                IF (thisobs_l%locweight_v>0 .AND. thisobs_l%nradii==3) THEN
                   thisobs_l%dist_l_v(off_obs) = dists(3)       ! distance in vertical direction
                END if
             END IF
          END DO scancountB

          ! Count overall local observations
          off_obs_l_all = off_obs_l_all + off_obs     ! dimension

       END IF

    ELSE

       WRITE (*,*) '+++++ ERROR PDAF-OMI: nonisotropic localization is only possible in 1, 2 or 3 dimensions'
       error = 11

    END IF

  END SUBROUTINE PDAFomi_init_obsarrays_l_noniso



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
          IF (obs_member==0) THEN
             WRITE (*,*) '++ OMI-debug: ', debug, &
                  'PDAFomi_g2l_obs -- START Get local observed ensemble mean'
          ELSE
             WRITE (*,*) '++ OMI-debug: ', debug, &
                  'PDAFomi_g2l_obs -- START Get local observed ensemble member', obs_member
          END IF
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
!> Compute weights for localization
!!
!! The routine is called during the analysis step
!! on each local analysis domain.
!!
  SUBROUTINE PDAFomi_observation_localization_weights(thisobs_l, thisobs, ncols, &
                                           A_l, weight, verbose)

  IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    INTEGER, INTENT(in) :: ncols             !< Rank of initial covariance matrix
    REAL, INTENT(in) :: A_l(:, :)            !< Input matrix (thisobs_l%dim_obs_l, ncols)
    INTEGER, INTENT(in) :: verbose           !< Verbosity flag
    REAL, INTENT(out) :: weight(thisobs_l%dim_obs_l) !> Localization weights


! *** local variables ***
    INTEGER :: i                       ! Index of observation component
    REAL, ALLOCATABLE :: weight_v(:)   ! Localization weights for vertical (for locweight_v>0)

! Screen output
       IF (verbose == 1) THEN
          WRITE (*, '(a, 5x, a, 1x, i3)') &
               'PDAFomi', '--- Domain localization for obs. type ID',thisobs%obsid
          IF (thisobs_l%nradii==1) THEN
             IF (thisobs%disttype<10) THEN
                WRITE (*, '(a, 8x, a)') &
                     'PDAFomi', '--- isotropic localization'
             ELSE
                WRITE (*, '(a, 8x, a)') &
                     'PDAFomi', '--- isotropic localization factorized in 2+1 dimensions'
             END IF
             WRITE (*, '(a, 8x, a, 1x, es11.3)') &
                  'PDAFomi', '--- Localization cut-off radius', thisobs_l%cradius
             WRITE (*, '(a, 8x, a, 1x, es11.3)') &
                  'PDAFomi', '--- Support radius', thisobs_l%sradius
          ELSE IF (thisobs_l%nradii==2) THEN
             WRITE (*, '(a, 8x, a)') &
                  'PDAFomi', '--- non-isotropic localization'
             WRITE (*, '(a, 8x, a, 1x, 2es11.3)') &
                  'PDAFomi', '--- Localization cut-off radii', thisobs_l%cradius
             WRITE (*, '(a, 8x, a, 1x, 2es11.3)') &
                  'PDAFomi', '--- Support radii', thisobs_l%sradius
          ELSE IF (thisobs_l%nradii==3) THEN
             IF (thisobs%disttype<10) THEN
                WRITE (*, '(a, 8x, a)') &
                     'PDAFomi', '--- non-isotropic localization in 3 dimensions'
             ELSE
                WRITE (*, '(a, 8x, a)') &
                     'PDAFomi', '--- non-isotropic localization factorized in 2+1 dimensions'
             END IF
             IF (thisobs_l%locweight_v>0) &
                WRITE (*, '(a, 8x, a)') &
                     'PDAFomi', '--- use separate localization function in vertical direction'
             WRITE (*, '(a, 8x, a, 1x, 3es11.3)') &
                  'PDAFomi', '--- Localization cut-off radii', thisobs_l%cradius
             WRITE (*, '(a, 8x, a, 1x, 3es11.3)') &
                  'PDAFomi', '--- Support radii', thisobs_l%sradius
          END IF
       ENDIF


! ***********************************************
! *** Apply a weight matrix with correlations ***
! *** of compact support to matrix A or the   ***
! *** observation error covariance matrix.    ***
! ***********************************************

       ! *** Initialize weight array

       CALL PDAFomi_weights_l(verbose, thisobs_l%dim_obs_l, ncols, thisobs_l%locweight, &
            thisobs_l%cradius_l, thisobs_l%sradius_l, &
            A_l, thisobs_l%ivar_obs_l, thisobs_l%distance_l, weight)

       ! *** For factorized 2+1D localization use product of horizontal and vertical weights
       IF (thisobs_l%locweight_v>0) then

          IF (verbose == 1) THEN
             WRITE (*, '(a, 8x, a)') &
                  'PDAFomi', '--- initialize also weight function for vertical direction'
          END IF

          ALLOCATE(weight_v(thisobs_l%dim_obs_l))

          CALL PDAFomi_weights_l_sgnl(verbose, thisobs_l%dim_obs_l, ncols, thisobs_l%locweight_v, &
               thisobs_l%cradius(3), thisobs_l%sradius(3), &
               A_l, thisobs_l%ivar_obs_l, thisobs_l%dist_l_v, weight_v)

          DO i = 1, thisobs_l%dim_obs_l
             weight(i) = weight(i) * weight_v(i)
          END DO

          DEALLOCATE(weight_v)
       END IF

       ! *** Handling of special weighting types ***

       lw2: IF (thisobs_l%locweight == 26) THEN
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
  END SUBROUTINE

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
          IF (thisobs_l%locweight_v>0) &
               WRITE (*,*) '++ OMI-debug prodrinva_l:    ', debug, '  thisobs_l%locweight_v', thisobs_l%locweight_v
          WRITE (*,*) '++ OMI-debug prodRinvA_l:    ', debug, '  thisobs%dim_obs_l', thisobs_l%dim_obs_l
          WRITE (*,*) '++ OMI-debug prodRinvA_l:    ', debug, '  thisobs%ivar_obs_l', thisobs_l%ivar_obs_l
          WRITE (*,*) '++ OMI-debug prodRinvA_l:    ', debug, '  Input matrix A_l', A_l
       END IF

       ALLOCATE(weight(thisobs_l%dim_obs_l))
       call PDAFomi_observation_localization_weights(thisobs_l, thisobs, ncols, A_l, &
                                         weight, verbose)

       ! *** Apply weight

       doweighting: IF (thisobs_l%locweight >= 11) THEN

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
!> Compute product of inverse of R with some matrix and hybrid weight
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
!! applied to matrix A. In addition the hybrid
!! weight alpha is applied.
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
!! * 2022-03 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_prodRinvA_hyb_l(thisobs_l, thisobs, nobs_all, ncols, &
       gamma, A_l, C_l, verbose)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    INTEGER, INTENT(in) :: nobs_all          !< Dimension of local obs. vector (all obs. types)
    INTEGER, INTENT(in) :: ncols             !< Rank of initial covariance matrix
    REAL, INTENT(in)    :: gamma             !< Hybrid weight
    REAL, INTENT(inout) :: A_l(:, :)         !< Input matrix (thisobs_l%dim_obs_l, ncols)
    REAL, INTENT(out)   :: C_l(:, :)         !< Output matrix (thisobs_l%dim_obs_l, ncols)
    INTEGER, INTENT(in) :: verbose           !< Verbosity flag


! *** local variables ***
    INTEGER :: i, j                    ! Index of observation component
    REAL, ALLOCATABLE :: weight(:)     ! Localization weights
    REAL, ALLOCATABLE :: weight_v(:)   ! Localization weights for vertical (for locweight_v>0)
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
               'PDAFomi_prodrinva_hyb_l -- START Multiply with inverse R and and apply localization'
          WRITE (*,*) '++ OMI-debug prodrinva_hyb_l:    ', debug, '  thisobs_l%locweight', thisobs_l%locweight
          WRITE (*,*) '++ OMI-debug prodrinva_hyb_l:    ', debug, 'thisobs%dim_obs_f', thisobs_l%dim_obs_l
          WRITE (*,*) '++ OMI-debug prodrinva_hyb_l:    ', debug, 'thisobs%ivar_obs_f', thisobs_l%ivar_obs_l
          WRITE (*,*) '++ OMI-debug prodrinva_hyb_l:    ', debug, 'Input matrix A_l', A_l
       END IF

       IF (verbose == 1) THEN
          WRITE (*,'(a, 5x, a, f12.5)') 'PDAFomi', '--- hybrid gamma=', gamma
          WRITE (*, '(a, 5x, a, 1x)') &
               'PDAFomi', '--- Domain localization'
          WRITE (*, '(a, 8x, a, 1x, es11.3)') &
               'PDAFomi', '--- Localization cut-off radius', thisobs_l%cradius
          WRITE (*, '(a, 8x, a, 1x, es11.3)') &
               'PDAFomi', '--- Support radius', thisobs_l%sradius
       ENDIF


! ***********************************************
! *** Apply a weight matrix with correlations ***
! *** of compact support to matrix A or the   ***
! *** observation error covariance matrix.    ***
! ***********************************************

       ! *** Initialize weight array

       ALLOCATE(weight(thisobs_l%dim_obs_l))

       CALL PDAFomi_weights_l(verbose, thisobs_l%dim_obs_l, ncols, thisobs_l%locweight, &
            thisobs_l%cradius_l, thisobs_l%sradius_l, &
            A_l, thisobs_l%ivar_obs_l, thisobs_l%distance_l, weight)

       ! *** For factorized 2+1D localization use product of horizontal and vertical weights
       IF (thisobs_l%locweight_v>0) then

          IF (verbose == 1) THEN
             WRITE (*, '(a, 8x, a)') &
                  'PDAFomi', '--- initialize also weight function for vertical direction'
          END IF

          ALLOCATE(weight_v(thisobs_l%dim_obs_l))

          CALL PDAFomi_weights_l_sgnl(verbose, thisobs_l%dim_obs_l, ncols, thisobs_l%locweight_v, &
               thisobs_l%cradius(3), thisobs_l%sradius(3), &
               A_l, thisobs_l%ivar_obs_l, thisobs_l%dist_l_v, weight_v)

          DO i = 1, thisobs_l%dim_obs_l
             weight(i) = weight(i) * weight_v(i)
          END DO

          DEALLOCATE(weight_v)
       END IF


       ! *** Handling of special weighting types ***

       lw2: IF (thisobs_l%locweight == 26) THEN
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

       doweighting: IF (thisobs_l%locweight >= 11) THEN

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
                C_l(i+off, j) = gamma * thisobs_l%ivar_obs_l(i) * A_l(i+off, j)
             END DO
          END DO
  
       ELSE doweighting

          ! *** Apply weight to matrix R only
          DO j = 1, ncols
             DO i = 1, thisobs_l%dim_obs_l
                C_l(i+off, j) = gamma * thisobs_l%ivar_obs_l(i) * weight(i) * A_l(i+off, j)
             END DO
          END DO
     
       END IF doweighting

       ! *** Clean up ***

       DEALLOCATE(weight)

       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_prodrinva_hyb_l -- END'

    ENDIF doassim

  END SUBROUTINE PDAFomi_prodRinvA_hyb_l



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
!! defined by the variables 'locweight', 'cradius', 
!! 'cradius2' and 'sradius' in the main program.
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
    REAL, ALLOCATABLE :: weight_v(:)      ! Localization weights for vertical (for locweight_v>0)
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
             WRITE (*, '(a, 5x, a)') &
                  'PDAFomi', '--- Assume Gaussian observation errors'
          ELSE
             WRITE (*, '(a, 5x, a)') &
                  'PDAFomi', '--- Assume double-exponential observation errors'
          END IF
          WRITE (*, '(a, 5x, a, 1x)') &
               'PDAFomi', '--- Domain localization'
          WRITE (*, '(a, 8x, a, 1x, es11.3)') &
               'PDAFomi', '--- Localization cut-off radius', thisobs_l%cradius
          WRITE (*, '(a, 8x, a, 1x, es11.3)') &
               'PDAFomi', '--- Support radius', thisobs_l%sradius
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
            thisobs_l%cradius_l, thisobs_l%sradius_l, &
            resid_obs, thisobs_l%ivar_obs_l, thisobs_l%distance_l, weight)

       ! *** For factorized 2+1D localization use product of horizontal and vertical weights
       IF (thisobs_l%locweight_v>0) then

          IF (verbose == 1) THEN
             WRITE (*, '(a, 8x, a)') &
                  'PDAFomi', '--- initialize also weight function for vertical direction'
          END IF

          ALLOCATE(weight_v(thisobs_l%dim_obs_l))

          CALL PDAFomi_weights_l_sgnl(verbose, thisobs_l%dim_obs_l, 1, thisobs_l%locweight_v, &
               thisobs_l%cradius(3), thisobs_l%sradius(3), &
               resid_obs, thisobs_l%ivar_obs_l, thisobs_l%dist_l_v, weight_v)

          DO i = 1, thisobs_l%dim_obs_l
             weight(i) = weight(i) * weight_v(i)
          END DO

          DEALLOCATE(weight_v)
       END IF

       DEALLOCATE(resid_obs)


       ! *** Handling of special weighting types ***

       lw2: IF (thisobs_l%locweight == 26) THEN
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
          Rinvresid_l(i) = thisobs_l%ivar_obs_l(i) * weight(i) * resid_l(thisobs_l%off_obs_l+i)
       END DO


! ********************************
! *** Compute local likelihood ***
! ********************************

       IF (thisobs%obs_err_type == 0) THEN

          ! Gaussian errors
          ! Calculate exp(-0.5*resid^T*R^-1*resid)

          ! Transform back to log likelihood to increment its values
          IF (lhood_l>0.0) lhood_l = - LOG(lhood_l)

          lhood_one = 0.0
          DO i = 1, thisobs_l%dim_obs_l
             lhood_one = lhood_one + 0.5*resid_l(thisobs_l%off_obs_l+i)*Rinvresid_l(i)
          END DO

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
!> Compute local likelihood for an ensemble member using hybrid weight
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
!! defined by the variables 'locweight', 'cradius', 
!! 'cradius2' and 'sradius' in the main program.
!! A tempering is appply by using the hybrid weight 'gamma'.
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
!! * 2022-03 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_likelihood_hyb_l(thisobs_l, thisobs, resid_l, gamma, lhood_l, verbose)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    REAL, INTENT(inout) :: resid_l(:)        !< Input vector of residuum
    REAL, INTENT(inout) :: lhood_l           !< Output vector - log likelihood
    REAL, INTENT(in)    :: gamma             !< Hybrid weight
    INTEGER, INTENT(in) :: verbose           !< Verbosity flag


! *** local variables ***
    INTEGER :: i                          ! Index of observation component
    REAL, ALLOCATABLE :: weight(:)        ! Localization weights
    REAL, ALLOCATABLE :: weight_v(:)      ! Localization weights for vertical (for locweight_v>0)
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
               'PDAFomi_likelihood_hyb_l -- START localization and likelihood, member', obs_member
          WRITE (*,*) '++ OMI-debug likelihood_hyb_l:  ', debug, '  thisobs_l%locweight', thisobs_l%locweight
       END IF

       ! Screen output
       IF (verbose == 1) THEN
          IF (thisobs%obs_err_type==0) THEN
             WRITE (*, '(a, 5x, a)') &
                  'PDAFomi', '--- Assume Gaussian observation errors'
          ELSE
             WRITE (*, '(a, 5x, a)') &
                  'PDAFomi', '--- Assume double-exponential observation errors'
          END IF
          WRITE (*, '(a, 5x, a, f12.5)') &
               'PDAFomi', '--- Apply tempering with 1.0-gamma= ', 1.0 - gamma
          WRITE (*, '(a, 5x, a, 1x)') &
               'PDAFomi', '--- Domain localization'
          WRITE (*, '(a, 8x, a, 1x, es11.3)') &
               'PDAFomi', '--- Localization cut-off radius', thisobs_l%cradius
          WRITE (*, '(a, 8x, a, 1x, es11.3)') &
               'PDAFomi', '--- Support radius', thisobs_l%sradius
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
            thisobs_l%cradius_l, thisobs_l%sradius_l, &
            resid_obs, thisobs_l%ivar_obs_l, thisobs_l%distance_l, weight)

       ! *** For factorized 2+1D localization use product of horizontal and vertical weights
       IF (thisobs_l%locweight_v>0) then

          IF (verbose == 1) THEN
             WRITE (*, '(a, 8x, a)') &
                  'PDAFomi', '--- initialize also weight function for vertical direction'
          END IF

          ALLOCATE(weight_v(thisobs_l%dim_obs_l))

          CALL PDAFomi_weights_l_sgnl(verbose, thisobs_l%dim_obs_l, 1, thisobs_l%locweight_v, &
               thisobs_l%cradius(3), thisobs_l%sradius(3), &
               resid_obs, thisobs_l%ivar_obs_l, thisobs_l%dist_l_v, weight_v)

          DO i = 1, thisobs_l%dim_obs_l
             weight(i) = weight(i) * weight_v(i)
          END DO

          DEALLOCATE(weight_v)
       END IF

       DEALLOCATE(resid_obs)


       ! *** Handling of special weighting types ***

       lw2: IF (thisobs_l%locweight == 26) THEN
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
          Rinvresid_l(i) = (1.0-gamma) * thisobs_l%ivar_obs_l(i) * weight(i) * resid_l(i)
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
          WRITE (*,*) '++ OMI-debug likelihood_hyb_l:  ', debug, '  accumulated likelihood', lhood_l
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_likelihood_hyb)l -- END'
       END IF

    END IF doassim

  END SUBROUTINE PDAFomi_likelihood_hyb_l



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
  SUBROUTINE PDAFomi_localize_covar_iso(thisobs, dim, locweight, cradius, sradius, &
       coords, HP, HPH)

    USE PDAFomi_obs_f, &
         ONLY: obsdims
    USE PDAF_mod_filtermpi, &
       ONLY: npes

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(in) :: thisobs    !< Data type with full observation
    INTEGER, INTENT(in) :: dim            !< State dimension
    INTEGER, INTENT(in) :: locweight      !< Localization weight type
    REAL, INTENT(in)    :: cradius        !< localization radius
    REAL, INTENT(in)    :: sradius        !< support radius for weight functions
    REAL, INTENT(in)    :: coords(:,:)    !< Coordinates of state vector elements
    REAL, INTENT(inout) :: HP(:, :)       !< Matrix HP, dimension (nobs, dim)
    REAL, INTENT(inout) :: HPH(:, :)      !< Matrix HPH, dimension (nobs, nobs)

! *** local variables ***
    INTEGER :: i, j, pe, cnt ! counters
    INTEGER :: ncoord        ! Number of coordinates
    REAL    :: distance      ! Distance between points in the domain 
    REAL    :: weight        ! Localization weight
    REAL, ALLOCATABLE :: weights(:) ! Localization weights array
    REAL    :: tmp(1,1)= 1.0 ! Temporary, but unused array
    INTEGER :: wtype         ! Type of weight function
    INTEGER :: rtype         ! Type of weight regulation
    REAL, ALLOCATABLE :: co(:), oc(:)   ! Coordinates of single point
    INTEGER, ALLOCATABLE :: id_start(:) ! Start index of obs. type in global averall obs. vector
    INTEGER, ALLOCATABLE :: id_end(:)   ! End index of obs. type in global averall obs. vector
    INTEGER, ALLOCATABLE :: obs_map(:)  ! Mapping indiced for observations 


    doassim: IF (thisobs%doassim == 1) THEN

! **********************
! *** INITIALIZATION ***
! **********************

       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_localize_covar -- START'
          WRITE (*,*) '++ OMI-debug localize_covar:', debug, 'thisobs%off_obs_g', thisobs%off_obs_g
          WRITE (*,*) '++ OMI-debug localize_covar:', debug, 'thisobs%dim_obs_g', thisobs%dim_obs_g
       END IF

       ! Screen output
       IF (screen > 0 .AND. mype==0) THEN
          WRITE (*,'(a, 8x, a, 1x, i3)') &
               'PDAFomi', '--- Apply covariance localization, obs. type ID', thisobs%obsid
          WRITE (*, '(a, 12x, a, 1x, f12.2)') &
               'PDAFomi', '--- Local influence radius', cradius

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


       ! *** Initialize mapping indices

       ALLOCATE(id_start(npes), id_end(npes))

       ! thisobs%ocoord_f is global for all observations of one type
       ! while HP and HPH are ordered obstype-first. Thus the observations
       ! of all types of one sub-domain are combined. This mapping
       ! ensures that the correct indices are used in HP and HPH.
       pe = 1
       id_start(1) = 1
       IF (thisobs%obsid>1) id_start(1) = id_start(1) + SUM(obsdims(1, 1:thisobs%obsid-1))
       id_end(1)   = id_start(1) + obsdims(1,thisobs%obsid) - 1
       DO pe = 2, npes
          id_start(pe) = id_start(pe-1) + SUM(obsdims(pe-1,thisobs%obsid:))
          IF (thisobs%obsid>1) id_start(pe) = id_start(pe) + SUM(obsdims(pe,1:thisobs%obsid-1))
          id_end(pe) = id_start(pe) + obsdims(pe,thisobs%obsid) - 1
       END DO

       ALLOCATE(obs_map(thisobs%dim_obs_g))
       cnt = 1
       DO pe = 1, npes
          DO i = id_start(pe), id_end(pe)
             obs_map(cnt) = i
             cnt = cnt + 1
          END DO
       END DO


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

       ALLOCATE(weights(thisobs%dim_obs_g))

       DO i = 1, dim

          ! Initialize coordinate
          co(1:ncoord) = coords(1:thisobs%ncoord, i)

          DO j = 1, thisobs%dim_obs_g

             ! Initialize coordinate
             oc(1:ncoord) = thisobs%ocoord_f(1:thisobs%ncoord, j)

             ! Compute distance
             CALL PDAFomi_comp_dist2(thisobs, co, oc, distance, (i*j)-1)
             distance = SQRT(distance)

             ! Compute weight
             CALL PDAF_local_weight(wtype, rtype, cradius, sradius, distance, &
                  1, 1, tmp, 1.0, weights(j), 0)
          END DO

          IF (debug==i) THEN
             WRITE (*,*) '++ OMI-debug localize_covar:  ', debug, 'weights for row in HP', weights
          END IF

          DO j = 1, thisobs%dim_obs_g

             ! Apply localization
             HP(obs_map(j), i) = weights(j) * HP(obs_map(j), i)

          END DO
       END DO

       DEALLOCATE(weights)


       ! *** Localize HPH^T ***

       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug localize_covar:', debug, '  localize matrix HPH^T'
       END IF

       DO i = 1, thisobs%dim_obs_g

          ! Initialize coordinate
          co(1:ncoord) = thisobs%ocoord_f(1:thisobs%ncoord, i)

          DO j = 1, thisobs%dim_obs_g

             ! Initialize coordinate
             oc(1:ncoord) = thisobs%ocoord_f(1:thisobs%ncoord, j)

             ! Compute distance
             CALL PDAFomi_comp_dist2(thisobs, co, oc, distance, (i*j)-1)
             distance = SQRT(distance)

             ! Compute weight
             CALL PDAF_local_weight(wtype, rtype, cradius, sradius, distance, &
                  1, 1, tmp, 1.0, weight, 0)

             ! Apply localization
             HPH(obs_map(j), obs_map(i)) = weight * HPH(obs_map(j), obs_map(i))

          END DO
       END DO

       ! clean up
       DEALLOCATE(co, oc)
       DEALLOCATE(id_start, id_end, obs_map)

       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_localize_covar -- END'

    END IF doassim

  END SUBROUTINE PDAFomi_localize_covar_iso



!-------------------------------------------------------------------------------
!> Apply covariance localization: 2+1D factorized with vertical localization weight
!!
!! This routine is a variant of PDAFomi_localize_covar_noniso with
!! support for a vector of localization weights. This is used
!! to specify different localization functions for the vertical and 
!! horizontal directions. The routine only stores the value of 
!! locweights(2) for the vertical and calls PDAFomi_init_dim_obs_l_noniso.
!!
!! __Revision history:__
!! * 2024-04 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_localize_covar_noniso_locweights(thisobs, dim, locweights, cradius, sradius, &
       coords, HP, HPH)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs !< Data type with full observation
    INTEGER, INTENT(in) :: dim            !< State dimension
    INTEGER, INTENT(in) :: locweights(:)  !< Types of localization function
    REAL, INTENT(in) :: cradius(:)        !< Vector of localization cut-off radii
    REAL, INTENT(in) :: sradius(:)        !< Vector of support radii of localization function
    REAL, INTENT(in)    :: coords(:,:)    !< Coordinates of state vector elements
    REAL, INTENT(inout) :: HP(:, :)       !< Matrix HP, dimension (nobs, dim)
    REAL, INTENT(inout) :: HPH(:, :)      !< Matrix HPH, dimension (nobs, nobs)

! *** Store vertical locweight and call standard routine

    IF (debug>0) THEN
       WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_localize_covar_noniso_locweights -- START'
       WRITE (*,*) '++ OMI-debug PDAFomi_localize_covar_noniso_locweights:', debug, '  locweights', locweights
    END IF

    ! Check consistency of dimensions
    IF (SIZE(locweights) /= 2) THEN
       WRITE (*,*) '+++++ ERROR PDAF-OMI: Input for locweights in horizontal and vertical needs size 2'
       error = 15
    END IF
    IF (thisobs%ncoord /= 3) THEN
       WRITE (*,*) '+++++ WARNING PDAF-OMI: separate locweight for vertical is only utilized if thisobs%ncoord=3'
    END IF

    IF (thisobs%ncoord == 3) THEN
       ! locweight for the vertical is treated separately
       thisobs%locweight_v = locweights(2)
    END IF

    CALL PDAFomi_localize_covar_noniso(thisobs, dim, locweights(1), cradius, sradius, &
         coords, HP, HPH)

  END SUBROUTINE PDAFomi_localize_covar_noniso_locweights

!-------------------------------------------------------------------------------
!> Apply covariance localization
!!
!! This routine applies a localization matrix B
!! to the matrices HP and HPH^T of the localized EnKF.
!! This variant is for non-iceotropic localization
!!
!! __Revision history:__
!! * 2020-03 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_localize_covar_noniso(thisobs, dim, locweight, cradius, sradius, &
       coords, HP, HPH)

    USE PDAFomi_obs_f, &
         ONLY: obsdims
    USE PDAF_mod_filtermpi, &
       ONLY: npes

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(in) :: thisobs    !< Data type with full observation
    INTEGER, INTENT(in) :: dim            !< State dimension
    INTEGER, INTENT(in) :: locweight      !< Localization weight type
    REAL, INTENT(in) :: cradius(:)        !< Vector of localization cut-off radii
    REAL, INTENT(in) :: sradius(:)        !< Vector of support radii of localization function
    REAL, INTENT(in)    :: coords(:,:)    !< Coordinates of state vector elements
    REAL, INTENT(inout) :: HP(:, :)       !< Matrix HP, dimension (nobs, dim)
    REAL, INTENT(inout) :: HPH(:, :)      !< Matrix HPH, dimension (nobs, nobs)

! *** local variables ***
    INTEGER :: i, j, pe, cnt ! counters
    INTEGER :: ncoord        ! Number of coordinates
    REAL    :: distance      ! Distance between points in the domain 
    REAL    :: weight        ! Localization weight
    REAL, ALLOCATABLE :: weights(:) ! Localization weights array
    REAL    :: weight_v      ! Weights in vertical direction
    REAL    :: tmp(1,1)= 1.0 ! Temporary, but unused array
    INTEGER :: wtype         ! Type of weight function
    INTEGER :: rtype         ! Type of weight regulation
    REAL :: srad, crad       ! localization cut-off radius
    REAL, ALLOCATABLE :: co(:), oc(:)   ! Coordinates of single point
    INTEGER, ALLOCATABLE :: id_start(:) ! Start index of obs. type in global averall obs. vector
    INTEGER, ALLOCATABLE :: id_end(:)   ! End index of obs. type in global averall obs. vector
    INTEGER, ALLOCATABLE :: obs_map(:)  ! Mapping indiced for observations 
    LOGICAL :: checkdist     ! Flag whether distance nis not larger than cut-off radius
    REAL :: dists(thisobs%ncoord)   ! Distance vector between analysis point and observation
    TYPE(obs_l) :: thisobs_l ! local observation


    doassim: IF (thisobs%doassim == 1) THEN

! **********************
! *** INITIALIZATION ***
! **********************

       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_localize_covar_noniso -- START'
          WRITE (*,*) '++ OMI-debug localize_covar_noniso:', debug, 'thisobs%off_obs_g', thisobs%off_obs_g
          WRITE (*,*) '++ OMI-debug localize_covar_noniso:', debug, 'thisobs%dim_obs_g', thisobs%dim_obs_g
       END IF

       ! Check consistency of dimensions
       IF (SIZE(cradius) /= thisobs%ncoord) THEN
          WRITE (*,*) '+++++ ERROR PDAF-OMI: non-isotropic localization: Size of CRADIUS /= thisobs%ncoord'
          error = 12
       END IF
       IF (SIZE(sradius) /= thisobs%ncoord) THEN
          WRITE (*,*) '+++++ ERROR PDAF-OMI: non-isotropic localization: Size of SRADIUS /= thisobs%ncoord'
          error = 13
       END IF

       ! Screen output
       IF (screen > 0 .AND. mype==0) THEN
          WRITE (*,'(a, 8x, a)') &
               'PDAFomi', '--- Apply non-isotropic covariance localization'
          IF (thisobs%ncoord==1) THEN
             WRITE (*, '(a, 12x, a, 1x, es11.3)') &
                  'PDAFomi', '--- Local influence radius', cradius
          ELSEIF (thisobs%ncoord==2) THEN
             WRITE (*, '(a, 12x, a, 1x, 2es11.3)') &
                  'PDAFomi', '--- Local influence radii', cradius
          ELSE
             WRITE (*, '(a, 12x, a, 1x, 3f11.3)') &
                  'PDAFomi', '--- Local influence radii', cradius
          END IF

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

       ! *** Initialize mapping indices

       ALLOCATE(id_start(npes), id_end(npes))

       ! thisobs%ocoord_f is global for all observations of one type
       ! while HP and HPH are ordered obstype-first. Thus the observations
       ! of all types of one sub-domain are combined. This mapping
       ! ensures that the correct indices are used in HP and HPH.
       pe = 1
       id_start(1) = 1
       IF (thisobs%obsid>1) id_start(1) = id_start(1) + SUM(obsdims(1, 1:thisobs%obsid-1))
       id_end(1)   = id_start(1) + obsdims(1,thisobs%obsid) - 1
       DO pe = 2, npes
          id_start(pe) = id_start(pe-1) + SUM(obsdims(pe-1,thisobs%obsid:))
          IF (thisobs%obsid>1) id_start(pe) = id_start(pe) + SUM(obsdims(pe,1:thisobs%obsid-1))
          id_end(pe) = id_start(pe) + obsdims(pe,thisobs%obsid) - 1
       END DO

       ALLOCATE(obs_map(thisobs%dim_obs_g))
       cnt = 1
       DO pe = 1, npes
          DO i = id_start(pe), id_end(pe)
             obs_map(cnt) = i
             cnt = cnt + 1
          END DO
       END DO

       ! Allocate vectors for localization radii and store their values
       IF (ALLOCATED(thisobs_l%cradius)) DEALLOCATE(thisobs_l%cradius)
       ALLOCATE(thisobs_l%cradius(thisobs%ncoord))
       IF (ALLOCATED(thisobs_l%sradius)) DEALLOCATE(thisobs_l%sradius)
       ALLOCATE(thisobs_l%sradius(thisobs%ncoord))

       ! Set number of localization radii
       thisobs_l%nradii = thisobs%ncoord
       thisobs_l%cradius(:) = cradius(:)
       thisobs_l%sradius(:) = sradius(:)



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

       ALLOCATE(weights(thisobs%dim_obs_g))

       cnt=0
       DO i = 1, dim

          ! Initialize coordinate
          co(1:ncoord) = coords(1:thisobs%ncoord, i)

          DO j = 1, thisobs%dim_obs_g

             ! Initialize coordinate
             oc(1:ncoord) = thisobs%ocoord_f(1:thisobs%ncoord, j)

             ! Compute distance
             CALL PDAFomi_check_dist2_noniso(thisobs, thisobs_l, co, oc, distance, &
                  dists, crad, srad, checkdist, (i*j)-1, cnt)
             distance = SQRT(distance)

             ! Compute weight
             CALL PDAF_local_weight(wtype, rtype, crad, srad, distance, &
                  1, 1, tmp, 1.0, weights(j), 0)

             ! Compute separate weight for vertical direction (for factorized 2+1D localization)
             IF (thisobs%locweight_v>0 .AND. thisobs%ncoord==3) THEN

                CALL PDAF_local_weight(wtype, rtype, cradius(3), sradius(3), dists(3), &
                     1, 1, tmp, 1.0, weight_v, 0)
                weights(j) = weights(j) * weight_v

             END IF
          END DO

          IF (debug==i) THEN
             WRITE (*,*) '++ OMI-debug localize_covar:  ', debug, 'weights for row in HP', weights
          END IF

          DO j = 1, thisobs%dim_obs_g

             ! Apply localization
             HP(obs_map(j), i) = weights(j) * HP(obs_map(j), i)

          END DO
       END DO

       DEALLOCATE(weights)


       ! *** Localize HPH^T ***

       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug localize_covar:', debug, '  localize matrix HPH^T'
       END IF

       cnt=0
       DO i = 1, thisobs%dim_obs_g

          ! Initialize coordinate
          co(1:ncoord) = thisobs%ocoord_f(1:thisobs%ncoord, i)

          DO j = 1, thisobs%dim_obs_g

             ! Initialize coordinate
             oc(1:ncoord) = thisobs%ocoord_f(1:thisobs%ncoord, j)

             ! Compute distance
             CALL PDAFomi_check_dist2_noniso(thisobs, thisobs_l, co, oc, distance, &
                  dists, crad, srad, checkdist, (i*j)-1, cnt)
             distance = SQRT(distance)

             ! Compute weight
             CALL PDAF_local_weight(wtype, rtype, crad, srad, distance, &
                  1, 1, tmp, 1.0, weight, 0)

             ! Compute separate weight for vertical direction (for factorized 2+1D localization)
             IF (thisobs%locweight_v>0 .AND. thisobs%ncoord==3) THEN

                CALL PDAF_local_weight(wtype, rtype, cradius(3), sradius(3), dists(3), &
                     1, 1, tmp, 1.0, weight_v, 0)
                weights(j) = weights(j) * weight_v

             END IF

             ! Apply localization
             HPH(obs_map(j), obs_map(i)) = weight * HPH(obs_map(j), obs_map(i))

          END DO
       END DO

       ! clean up
       DEALLOCATE(co, oc)
       DEALLOCATE(id_start, id_end, obs_map)

       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_localize_covar -- END'

    END IF doassim

  END SUBROUTINE PDAFomi_localize_covar_noniso



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
       IF (thisobs_l%dim_obs_l>0) THEN
          WRITE (*,*) '++ OMI-debug g2l_obs:       ', debug, '  thisobs_l%id_obs_l', thisobs_l%id_obs_l
          WRITE (*,*) '++ OMI-debug g2l_obs:       ', debug, '  obs_l', &
               obs_l_all(1+offset_obs_l_all:offset_obs_l_all+thisobs_l%dim_obs_l)
       ELSE
          WRITE (*,*) '++ OMI-debug g2l_obs:       ', debug, '  no local observations present'
       END IF
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
    REAL, INTENT(in) :: coordsA(:)        !< Coordinates of current analysis domain (ncoord)
    REAL, INTENT(in) :: coordsB(:)        !< Coordinates of observation (ncoord)
    REAL, INTENT(out) :: distance2        !< Squared distance
    INTEGER, INTENT(in) :: verbose        !< Control screen output

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

    norm: IF ((thisobs%disttype==0 .OR. thisobs%disttype==10) .OR. &
         ((thisobs%disttype==1 .OR. thisobs%disttype==11) .AND. domsize==0)) THEN

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

    ELSEIF ((thisobs%disttype==1 .OR. thisobs%disttype==11) .AND. domsize==1) THEN norm

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
       IF (thisobs%disttype<10) THEN
          ! full 3D localization
          DO k = 1, thisobs%ncoord
             distance2 = distance2 + dists(k)*dists(k)
          END DO
       ELSE
          ! factorized 2+1D localization
          DO k = 1, thisobs%ncoord-1
             distance2 = distance2 + dists(k)*dists(k)
          END DO
       END IF

    ELSEIF (thisobs%disttype==2 .OR. thisobs%disttype==12) THEN norm

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

    ELSEIF (thisobs%disttype==3 .OR. thisobs%disttype==13) THEN norm

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
!> Check distance in case of isotropic localization
!!
!! This routine computes the distance between two locations.
!! The computation can be for Cartesian grids with and without
!! periodicity and for geographic coordinates. For Cartesian
!! grids, the coordinates can be in any unit, while geographic
!! coordinates must be provided in radians and the resulting
!! distance will be in meters. Finally, the routine checks
!! whether the distance is not larger than the cut-off radius.
!!
!! Choices for distance computation - disttype:
!! 0: Cartesian distance in ncoord dimensions
!! 1: Cartesian distance in ncoord dimensions with periodicity
!!    (Needs specification of domsize(ncoord))
!! 2: Aproximate geographic distance with horizontal coordinates in radians (-pi/+pi)
!! 3: Geographic distance computation using haversine formula
!! 10-13: Variants of distance types 0-3, but particularly for 3 dimensions in which 
!!    a 2+1 dimensional localization is applied (distance weighting only in the horizontal)
!!
!! __Revision history:__
!! * 2024-04 - Lars Nerger - Initial code based on PDAFomi_comp_dist2
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_check_dist2(thisobs, thisobs_l, coordsA, coordsB, distance2, &
       checkdist, verbose, cnt_obs)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(in) :: thisobs    !< Data type with full observation
    TYPE(obs_l), INTENT(in) :: thisobs_l  !< Data type with local observation
    REAL, INTENT(in) :: coordsA(:)        !< Coordinates of current analysis domain (ncoord)
    REAL, INTENT(in) :: coordsB(:)        !< Coordinates of observation (ncoord)
    REAL, INTENT(out) :: distance2        !< Squared distance
    LOGICAL, INTENT(out) :: checkdist     !< Flag whether distance is within cut-off radius
    INTEGER, INTENT(in) :: verbose        !< Control screen output
    INTEGER, INTENT(inout) :: cnt_obs     !< Count number of local observations

! *** Local variables ***
    INTEGER :: k                    ! Counters
    REAL :: dists(thisobs%ncoord)   ! Distance vector between analysis point and observation
    REAL :: slon, slat              ! sine of distance in longitude or latitude
    REAL :: cradius2                ! squared localization cut-off radius
    INTEGER :: domsize              ! Flag whether domainsize is set
    LOGICAL :: distflag             ! Flag whether distance in a coordinate direction is within cradius


! **********************
! *** Initialization ***
! **********************

    ! Initialize distance flag
    checkdist = .FALSE.
    distflag = .TRUE.


! ************************
! *** Compute distance ***
! ************************

    IF (.NOT.ALLOCATED(thisobs%domainsize)) THEN
       domsize = 0
    ELSE
       domsize = 1
    END IF

    norm: IF ((thisobs%disttype==0 .OR. thisobs%disttype==10) .OR. &
         ((thisobs%disttype==1 .OR. thisobs%disttype==11) .AND. domsize==0)) THEN

       ! *** Compute Cartesian distance ***

       IF (debug>0 .AND. verbose==0) THEN
          WRITE (*,*) '++ OMI-debug check_dist2:    ', debug, '  compute Cartesian distance'
       END IF

       IF (thisobs%ncoord>=3) THEN
          dists(3) = ABS(coordsA(3) - coordsB(3))
          IF (dists(3)>thisobs_l%cradius(1)) THEN
             distflag = .FALSE.
          ELSE
             dists(2) = ABS(coordsA(2) - coordsB(2))
             IF (dists(2)>thisobs_l%cradius(1)) THEN
                distflag = .FALSE.
             ELSE
                dists(1) = ABS(coordsA(1) - coordsB(1))
                IF (dists(1)>thisobs_l%cradius(1)) THEN
                   distflag = .FALSE.
                ELSE
                   ! full squared distance
                   distance2 = 0.0
                   IF (thisobs%disttype<10) THEN
                      ! full 3D localization
                      DO k = 1, thisobs%ncoord
                         distance2 = distance2 + dists(k)*dists(k)
                      END DO
                   ELSE
                      ! factorized 2+1D localization
                      DO k = 1, thisobs%ncoord-1
                         distance2 = distance2 + dists(k)*dists(k)
                      END DO
                   END IF
                END IF
             END IF
          END IF
       ELSEIF (thisobs%ncoord==2) THEN
          dists(2) = ABS(coordsA(2) - coordsB(2))
          IF (dists(2)>thisobs_l%cradius(1)) THEN
             distflag = .FALSE.
          ELSE
             dists(1) = ABS(coordsA(1) - coordsB(1))
             IF (dists(1)>thisobs_l%cradius(1)) THEN
                distflag = .FALSE.
             ELSE
                ! full squared distance
                distance2 = 0.0
                DO k = 1, thisobs%ncoord
                   distance2 = distance2 + dists(k)*dists(k)
                END DO
             END IF
          END IF
       ELSEIF (thisobs%ncoord==1) THEN
          dists(1) = ABS(coordsA(1) - coordsB(1))
          IF (dists(1)>thisobs_l%cradius(1)) THEN
             distflag = .FALSE.
          ELSE
             ! full squared distance
             distance2 = 0.0
             DO k = 1, thisobs%ncoord
                distance2 = distance2 + dists(k)*dists(k)
             END DO
          END IF
       END IF

    ELSEIF ((thisobs%disttype==1 .OR. thisobs%disttype==11) .AND. domsize==1) THEN norm

       ! *** Compute periodic Cartesian distance ***

       IF (debug>0 .AND. verbose==0) THEN
          WRITE (*,*) '++ OMI-debug check_dist2:    ', debug, '  compute periodic Cartesian distance'
       END IF

       IF (thisobs%ncoord>=3) THEN
          IF (thisobs%domainsize(3)<=0.0) THEN 
             dists(3) = ABS(coordsA(3) - coordsB(3))
          ELSE
             dists(3) = MIN(ABS(coordsA(3) - coordsB(3)), &
                  ABS(ABS(coordsA(3) - coordsB(3))-thisobs%domainsize(3)))
          END IF
          IF (dists(3)>thisobs_l%cradius(1)) THEN
             distflag = .FALSE.
          ELSE
             IF (thisobs%domainsize(2)<=0.0) THEN 
                dists(2) = ABS(coordsA(2) - coordsB(2))
             ELSE
                dists(2) = MIN(ABS(coordsA(2) - coordsB(2)), &
                     ABS(ABS(coordsA(2) - coordsB(2))-thisobs%domainsize(2)))
             END IF
             IF (dists(2)>thisobs_l%cradius(1)) THEN
                distflag = .FALSE.
             ELSE
                IF (thisobs%domainsize(1)<=0.0) THEN 
                   dists(1) = ABS(coordsA(1) - coordsB(1))
                ELSE
                   dists(1) = MIN(ABS(coordsA(1) - coordsB(1)), &
                        ABS(ABS(coordsA(1) - coordsB(1))-thisobs%domainsize(1)))
                END IF
                IF (dists(1)>thisobs_l%cradius(1)) THEN
                   distflag = .FALSE.
                ELSE
                   ! full squared distance
                   distance2 = 0.0
                   IF (thisobs%disttype<10) THEN
                      ! full 3D localization
                      DO k = 1, thisobs%ncoord
                         distance2 = distance2 + dists(k)*dists(k)
                      END DO
                   ELSE
                      ! factorized 2+1D localization
                      DO k = 1, thisobs%ncoord-1
                         distance2 = distance2 + dists(k)*dists(k)
                      END DO
                   END IF
                END IF
             END IF
          END IF
       ELSEIF (thisobs%ncoord==2) THEN
          IF (thisobs%domainsize(2)<=0.0) THEN 
             dists(2) = ABS(coordsA(2) - coordsB(2))
          ELSE
             dists(2) = MIN(ABS(coordsA(2) - coordsB(2)), &
                  ABS(ABS(coordsA(2) - coordsB(2))-thisobs%domainsize(2)))
          END IF
          IF (dists(2)>thisobs_l%cradius(1)) THEN
             distflag = .FALSE.
          ELSE
             IF (thisobs%domainsize(1)<=0.0) THEN 
                dists(1) = ABS(coordsA(1) - coordsB(1))
             ELSE
                dists(1) = MIN(ABS(coordsA(1) - coordsB(1)), &
                     ABS(ABS(coordsA(1) - coordsB(1))-thisobs%domainsize(1)))
             END IF
             IF (dists(1)>thisobs_l%cradius(1)) THEN
                distflag = .FALSE.
             ELSE
                ! full squared distance
                distance2 = 0.0
                DO k = 1, thisobs%ncoord
                   distance2 = distance2 + dists(k)*dists(k)
                END DO
             END IF
          END IF
       ELSEIF (thisobs%ncoord==1) THEN
          IF (thisobs%domainsize(1)<=0.0) THEN 
             dists(1) = ABS(coordsA(1) - coordsB(1))
          ELSE
             dists(1) = MIN(ABS(coordsA(1) - coordsB(1)), &
                  ABS(ABS(coordsA(1) - coordsB(1))-thisobs%domainsize(1)))
          END IF
          IF (dists(1)>thisobs_l%cradius(1)) THEN
             distflag = .FALSE.
          ELSE
             ! full squared distance
             distance2 = 0.0
             DO k = 1, thisobs%ncoord
                distance2 = distance2 + dists(k)*dists(k)
             END DO
          END IF
       END IF

    ELSEIF (thisobs%disttype==2 .OR. thisobs%disttype==12) THEN norm

       ! *** Compute distance from geographic coordinates ***

       IF (debug>0 .AND. verbose==0) THEN
          WRITE (*,*) '++ OMI-debug check_dist2:    ', debug, '  compute geographic distance'
       END IF

       IF (thisobs%ncoord==3) THEN
          dists(3) = ABS(coordsA(3) - coordsB(3))
          IF (dists(3)>thisobs_l%cradius(1)) THEN
             distflag = .FALSE.
          ELSE
             dists(2) = r_earth * ABS(coordsA(2) - coordsB(2))
             IF (dists(2)>thisobs_l%cradius(1)) THEN
                distflag = .FALSE.
             ELSE
                dists(1) = r_earth * MIN( ABS(coordsA(1) - coordsB(1))* COS(coordsA(2)), &
                     ABS(ABS(coordsA(1) - coordsB(1)) - 2.0*pi) * COS(coordsA(2)))
                IF (dists(1)>thisobs_l%cradius(1)) THEN
                   distflag = .FALSE.
                ELSE
                   ! full squared distance
                   distance2 = 0.0
                   IF (thisobs%disttype<10) THEN
                      ! full 3D localization
                      DO k = 1, thisobs%ncoord
                         distance2 = distance2 + dists(k)*dists(k)
                      END DO
                   ELSE
                      ! factorized 2+1D localization
                      DO k = 1, thisobs%ncoord-1
                         distance2 = distance2 + dists(k)*dists(k)
                      END DO
                   END IF
                END IF
             END IF
          END IF
       ELSE
          dists(2) = r_earth * ABS(coordsA(2) - coordsB(2))
          IF (dists(2)>thisobs_l%cradius(1)) THEN
             distflag = .FALSE.
          ELSE
             dists(1) = r_earth * MIN( ABS(coordsA(1) - coordsB(1))* COS(coordsA(2)), &
                  ABS(ABS(coordsA(1) - coordsB(1)) - 2.0*pi) * COS(coordsA(2)))
             IF (dists(1)>thisobs_l%cradius(1)) THEN
                distflag = .FALSE.
             ELSE
                ! full squared distance
                distance2 = 0.0
                DO k = 1, thisobs%ncoord
                   distance2 = distance2 + dists(k)*dists(k)
                END DO
             END IF
          END IF
       END IF

    ELSEIF (thisobs%disttype==3 .OR. thisobs%disttype==13) THEN norm

       ! *** Compute distance from geographic coordinates with haversine formula ***

       IF (debug>0 .AND. verbose==0) THEN
          WRITE (*,*) '++ OMI-debug check_dist2:    ', debug, &
               '  compute geographic distance using haversine function'
       END IF

       IF (thisobs%ncoord==3) THEN
          dists(3) = ABS(coordsA(3) - coordsB(3))
          IF (dists(3)>thisobs_l%cradius(1)) THEN
             distflag = .FALSE.
          ELSE
             dists(2) = r_earth * ABS(coordsA(2) - coordsB(2))
             IF (dists(2)>thisobs_l%cradius(1)) THEN
                distflag = .FALSE.
             ELSE
                ! Haversine formula
                slon = SIN((coordsA(1) - coordsB(1))/2)
                slat = SIN((coordsA(2) - coordsB(2))/2)

                dists(2) = SQRT(slat*slat + COS(coordsA(2))*COS(coordsB(2))*slon*slon)
                IF (dists(2)<=1.0) THEN
                   dists(2) = 2.0 * r_earth* ASIN(dists(2))
                ELSE
                   dists(2) = r_earth* pi
                END IF
                IF (dists(2)>thisobs_l%cradius(1)) THEN
                   distflag = .FALSE.
                ELSE
                   ! full squared distance
                   distance2 = 0.0
                   IF (thisobs%disttype<10) THEN
                      ! full 3D localization
                      DO k = 2, thisobs%ncoord
                         distance2 = distance2 + dists(k)*dists(k)
                      END DO
                   ELSE
                      ! factorized 2+1D localization
                      DO k = 2, thisobs%ncoord-1
                         distance2 = distance2 + dists(k)*dists(k)
                      END DO
                   END IF
                END IF
            END IF
          END IF
       ELSE
          dists(2) = r_earth * ABS(coordsA(2) - coordsB(2))
          IF (dists(2)>thisobs_l%cradius(1)) THEN
             distflag = .FALSE.
          ELSE
             ! Haversine formula
             slon = SIN((coordsA(1) - coordsB(1))/2)
             slat = SIN((coordsA(2) - coordsB(2))/2)

             dists(2) = SQRT(slat*slat + COS(coordsA(2))*COS(coordsB(2))*slon*slon)
             IF (dists(2)<=1.0) THEN
                dists(2) = 2.0 * r_earth* ASIN(dists(2))
             ELSE
                dists(2) = r_earth* pi
             END IF
             IF (dists(2)>thisobs_l%cradius(1)) THEN
                distflag = .FALSE.
             ELSE
                ! full squared distance
                distance2 = 0.0
                DO k = 1, thisobs%ncoord
                   distance2 = distance2 + dists(k)*dists(k)
                END DO
             END IF
          END IF
       END IF

    END IF norm

    IF (distflag) THEN
       cradius2 = thisobs_l%cradius(1)*thisobs_l%cradius(1)

       IF (distance2 <= cradius2) THEN
          ! Set flag for valid observation
          checkdist = .TRUE.

          ! Increment counter
          cnt_obs = cnt_obs + 1
       END IF
    END IF

  END SUBROUTINE PDAFomi_check_dist2




!-------------------------------------------------------------------------------
!> Check distance in case of nonisotropic localization
!!
!! This routine computes the distance between the observation and a 
!! model grid point and the cut-off radius of an ellipse (in 2D)
!! or ellipsoid (in 3D) in the direction of the distance. Finally,
!! the routine checks whether the distance is not larger than the
!! cut-off radius.
!!
!! Choices for distance computation - disttype:
!! 0: Cartesian distance in ncoord dimensions
!! 1: Cartesian distance in ncoord dimensions with periodicity
!!    (Needs specification of domsize(ncoord))
!! 2: Aproximate geographic distance with horizontal coordinates in radians (-pi/+pi)
!! 3: Geographic distance computation using haversine formula
!! 10-13: Variants of distance types 0-3, but particularly for 3 dimensions in which 
!!    a 2+1 dimensional localization is applied (distance weighting only in the horizontal)
!!
!! __Revision history:__
!! * 2024-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_check_dist2_noniso(thisobs, thisobs_l, coordsA, coordsB, distance2, &
       dists, cradius, sradius, checkdist, verbose, cnt_obs)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(in) :: thisobs    !< Data type with full observation
    TYPE(obs_l), INTENT(in) :: thisobs_l  !< Data type with local observation
    REAL, INTENT(in) :: coordsA(:)        !< Coordinates of current analysis domain (ncoord)
    REAL, INTENT(in) :: coordsB(:)        !< Coordinates of observation (ncoord)
    REAL, INTENT(out) :: distance2        !< Squared distance
    REAL, INTENT(inout) :: dists(:)       !< Vector of distance in each coordinate direction
    REAL, INTENT(out) :: cradius          !< Directional cut-off radius
    REAL, INTENT(inout) :: sradius        !< Directional support radius
    LOGICAL, INTENT(out) :: checkdist     !< Flag whether distance is within cut-off radius
    INTEGER, INTENT(in) :: verbose        !< Control screen output
    INTEGER, INTENT(inout) :: cnt_obs     !< Count number of local observations

! *** Local variables ***
    INTEGER :: k                    ! Counters
    REAL :: slon, slat              ! sine of distance in longitude or latitude
    INTEGER :: domsize              ! Flag whether domainsize is set
    REAL :: cradius2                ! cut-off radius on ellipse or ellipsoid
    REAL :: phi, theta              ! Angles in ellipse or ellipsoid
    REAL :: dist_xy                 ! Distance in xy-plan in 3D case
    LOGICAL :: distflag             ! Flag whether distance in a coordinate direction is within cradius


! **********************
! *** Initialization ***
! **********************

    ! Initialize distance flag
    checkdist = .FALSE.    ! Whether an observation lies within the local box
    distflag = .TRUE.      ! Whether an observation lies within the local radius (ellipse, ellipsoid)


! ************************
! *** Compute distance ***
! ************************

    IF (.NOT.ALLOCATED(thisobs%domainsize)) THEN
       domsize = 0
    ELSE
       domsize = 1
    END IF

    ! Debug output
    IF (debug>0 .AND. verbose==0) THEN
       WRITE (*,*) '++ OMI-debug check_dist2_noniso: ', debug, '  use non-isotropic localization'
    END IF

    norm: IF ((thisobs%disttype==0 .OR. thisobs%disttype==10) .OR. &
         ((thisobs%disttype==1 .OR. thisobs%disttype==11) .AND. domsize==0)) THEN

       ! *** Compute Cartesian distance ***

       IF (debug>0 .AND. verbose==0) THEN
          WRITE (*,*) '++ OMI-debug check_dist2_noniso: ', debug, '  compute Cartesian distance'
       END IF

       IF (thisobs%ncoord==3) THEN
          dists(3) = ABS(coordsA(3) - coordsB(3))
          IF (dists(3)>thisobs_l%cradius(3)) THEN
             distflag = .FALSE.
          ELSE
             dists(2) = ABS(coordsA(2) - coordsB(2))
             IF (dists(2)>thisobs_l%cradius(2)) THEN
                distflag = .FALSE.
             ELSE
                dists(1) = ABS(coordsA(1) - coordsB(1))
                IF (dists(1)>thisobs_l%cradius(1)) THEN
                   distflag = .FALSE.
                ELSE
                   ! full squared distance
                   distance2 = 0.0
                   IF (thisobs%disttype<10) THEN
                      ! full 3D localization
                      DO k = 1, thisobs%ncoord
                         distance2 = distance2 + dists(k)*dists(k)
                      END DO
                   ELSE
                      ! factorized 2+1D localization
                      DO k = 1, thisobs%ncoord-1
                         distance2 = distance2 + dists(k)*dists(k)
                      END DO
                   END IF
                END IF
             END IF
          END IF
       ELSEIF (thisobs%ncoord==2) THEN
          dists(2) = ABS(coordsA(2) - coordsB(2))
          IF (dists(2)>thisobs_l%cradius(2)) THEN
             distflag = .FALSE.
          ELSE
             dists(1) = ABS(coordsA(1) - coordsB(1))
             IF (dists(1)>thisobs_l%cradius(1)) THEN
                distflag = .FALSE.
             ELSE
                ! full squared distance
                distance2 = 0.0
                DO k = 1, thisobs%ncoord
                   distance2 = distance2 + dists(k)*dists(k)
                END DO
             END IF
          END IF
       ELSEIF (thisobs%ncoord==1) THEN
          dists(1) = ABS(coordsA(1) - coordsB(1))
          IF (dists(1)>thisobs_l%cradius(1)) THEN
             distflag = .FALSE.
          ELSE
             ! full squared distance
             distance2 = 0.0
             DO k = 1, thisobs%ncoord
                distance2 = distance2 + dists(k)*dists(k)
             END DO
          END IF
       END IF

    ELSEIF ((thisobs%disttype==1 .OR. thisobs%disttype==11) .AND. domsize==1) THEN norm

       ! *** Compute periodic Cartesian distance ***

       IF (debug>0 .AND. verbose==0) THEN
          WRITE (*,*) '++ OMI-debug check_dist2_noniso: ', debug, '  compute periodic Cartesian distance'
       END IF

       IF (thisobs%ncoord==3) THEN
          IF (thisobs%domainsize(3)<=0.0) THEN 
             dists(3) = ABS(coordsA(3) - coordsB(3))
          ELSE
             dists(3) = MIN(ABS(coordsA(3) - coordsB(3)), &
                  ABS(ABS(coordsA(3) - coordsB(3))-thisobs%domainsize(3)))
          END IF
          IF (dists(3)>thisobs_l%cradius(3)) THEN
             distflag = .FALSE.
          ELSE
             IF (thisobs%domainsize(2)<=0.0) THEN 
                dists(2) = ABS(coordsA(2) - coordsB(2))
             ELSE
                dists(2) = MIN(ABS(coordsA(2) - coordsB(2)), &
                     ABS(ABS(coordsA(2) - coordsB(2))-thisobs%domainsize(2)))
             END IF
             IF (dists(2)>thisobs_l%cradius(2)) THEN
                distflag = .FALSE.
             ELSE
                IF (thisobs%domainsize(1)<=0.0) THEN 
                   dists(1) = ABS(coordsA(1) - coordsB(1))
                ELSE
                   dists(1) = MIN(ABS(coordsA(1) - coordsB(1)), &
                        ABS(ABS(coordsA(1) - coordsB(1))-thisobs%domainsize(1)))
                END IF
                IF (dists(1)>thisobs_l%cradius(1)) THEN
                   distflag = .FALSE.
                ELSE
                   ! full squared distance
                   distance2 = 0.0
                   IF (thisobs%disttype<10) THEN
                      ! full 3D localization
                      DO k = 1, thisobs%ncoord
                         distance2 = distance2 + dists(k)*dists(k)
                      END DO
                   ELSE
                      ! factorized 2+1D localization
                      DO k = 1, thisobs%ncoord-1
                         distance2 = distance2 + dists(k)*dists(k)
                      END DO
                   END IF
                END IF
             END IF
          END IF
       ELSEIF (thisobs%ncoord==2) THEN
          IF (thisobs%domainsize(2)<=0.0) THEN 
             dists(2) = ABS(coordsA(2) - coordsB(2))
          ELSE
             dists(2) = MIN(ABS(coordsA(2) - coordsB(2)), &
                  ABS(ABS(coordsA(2) - coordsB(2))-thisobs%domainsize(2)))
          END IF
          IF (dists(2)>thisobs_l%cradius(2)) THEN
             distflag = .FALSE.
          ELSE
             IF (thisobs%domainsize(1)<=0.0) THEN 
                dists(1) = ABS(coordsA(1) - coordsB(1))
             ELSE
                dists(1) = MIN(ABS(coordsA(1) - coordsB(1)), &
                     ABS(ABS(coordsA(1) - coordsB(1))-thisobs%domainsize(1)))
             END IF
             IF (dists(1)>thisobs_l%cradius(1)) THEN
                distflag = .FALSE.
             ELSE
                ! full squared distance
                distance2 = 0.0
                DO k = 1, thisobs%ncoord
                   distance2 = distance2 + dists(k)*dists(k)
                END DO
             END IF
          END IF
       ELSEIF (thisobs%ncoord==1) THEN
          IF (thisobs%domainsize(1)<=0.0) THEN 
             dists(1) = ABS(coordsA(1) - coordsB(1))
          ELSE
             dists(1) = MIN(ABS(coordsA(1) - coordsB(1)), &
                  ABS(ABS(coordsA(1) - coordsB(1))-thisobs%domainsize(1)))
          END IF
          IF (dists(1)>thisobs_l%cradius(1)) THEN
             distflag = .FALSE.
          ELSE
             ! full squared distance
             distance2 = 0.0
             DO k = 1, thisobs%ncoord
                distance2 = distance2 + dists(k)*dists(k)
             END DO
          END IF
       END IF

    ELSEIF (thisobs%disttype==2 .OR. thisobs%disttype==12) THEN norm

       ! *** Compute distance from geographic coordinates ***

       IF (debug>0 .AND. verbose==0) THEN
          WRITE (*,*) '++ OMI-debug check_dist2_noniso: ', debug, '  compute geographic distance'
       END IF

       IF (thisobs%ncoord==3) THEN
          dists(3) = ABS(coordsA(3) - coordsB(3))
          IF (dists(3)>thisobs_l%cradius(3)) THEN
             distflag = .FALSE.
          ELSE
             dists(2) = r_earth * ABS(coordsA(2) - coordsB(2))
             IF (dists(2)>thisobs_l%cradius(2)) THEN
                distflag = .FALSE.
             ELSE
                dists(1) = r_earth * MIN( ABS(coordsA(1) - coordsB(1))* COS(coordsA(2)), &
                     ABS(ABS(coordsA(1) - coordsB(1)) - 2.0*pi) * COS(coordsA(2)))
                IF (dists(1)>thisobs_l%cradius(1)) THEN
                   distflag = .FALSE.
                ELSE
                   ! full squared distance
                   distance2 = 0.0
                   IF (thisobs%disttype<10) THEN
                      ! full 3D localization
                      DO k = 1, thisobs%ncoord
                         distance2 = distance2 + dists(k)*dists(k)
                      END DO
                   ELSE
                      ! factorized 2+1D localization
                      DO k = 1, thisobs%ncoord-1
                         distance2 = distance2 + dists(k)*dists(k)
                      END DO
                   END IF
                END IF
             END IF
          END IF
       ELSE
          dists(2) = r_earth * ABS(coordsA(2) - coordsB(2))
          IF (dists(2)>thisobs_l%cradius(2)) THEN
             distflag = .FALSE.
          ELSE
             dists(1) = r_earth * MIN( ABS(coordsA(1) - coordsB(1))* COS(coordsA(2)), &
                  ABS(ABS(coordsA(1) - coordsB(1)) - 2.0*pi) * COS(coordsA(2)))
             IF (dists(1)>thisobs_l%cradius(1)) THEN
                distflag = .FALSE.
             ELSE
                ! full squared distance
                distance2 = 0.0
                DO k = 1, thisobs%ncoord
                   distance2 = distance2 + dists(k)*dists(k)
                END DO
             END IF
          END IF
       END IF

    ELSEIF (thisobs%disttype==3 .OR. thisobs%disttype==13) THEN norm

       ! *** Compute distance from geographic coordinates with haversine formula ***

       IF (debug>0 .AND. verbose==0) THEN
          WRITE (*,*) '++ OMI-debug check_dist2_noniso:    ', debug, &
               '  compute geographic distance using haversine function'
       END IF

       IF (thisobs%ncoord==3) THEN
          dists(3) = ABS(coordsA(3) - coordsB(3))
          IF (dists(3)>thisobs_l%cradius(3)) THEN
             distflag = .FALSE.
          ELSE
             dists(2) = r_earth * ABS(coordsA(2) - coordsB(2))
             IF (dists(2)>thisobs_l%cradius(2)) THEN
                distflag = .FALSE.
             ELSE
                dists(1) = r_earth * MIN( ABS(coordsA(1) - coordsB(1))* COS(coordsA(2)), &
                     ABS(ABS(coordsA(1) - coordsB(1)) - 2.0*pi) * COS(coordsA(2)))

                ! Haversine formula
                slon = SIN((coordsA(1) - coordsB(1))/2)
                slat = SIN((coordsA(2) - coordsB(2))/2)

                dists(2) = SQRT(slat*slat + COS(coordsA(2))*COS(coordsB(2))*slon*slon)
                IF (dists(2)<=1.0) THEN
                   dists(2) = 2.0 * r_earth* ASIN(dists(2))
                ELSE
                   dists(2) = r_earth* pi
                END IF
                IF (dists(2)>thisobs_l%cradius(1)) THEN
                   distflag = .FALSE.
                ELSE
                   ! full squared distance
                   distance2 = 0.0
                   IF (thisobs%disttype<10) THEN
                      ! full 3D localization
                      DO k = 2, thisobs%ncoord
                         distance2 = distance2 + dists(k)*dists(k)
                      END DO
                   ELSE
                      ! factorized 2+1D localization
                      DO k = 2, thisobs%ncoord-1
                         distance2 = distance2 + dists(k)*dists(k)
                      END DO
                   END IF
                END IF
            END IF
          END IF
       ELSE
          dists(2) = r_earth * ABS(coordsA(2) - coordsB(2))
          IF (dists(2)>thisobs_l%cradius(2)) THEN
             distflag = .FALSE.
          ELSE
             dists(1) = r_earth * MIN( ABS(coordsA(1) - coordsB(1))* COS(coordsA(2)), &
                  ABS(ABS(coordsA(1) - coordsB(1)) - 2.0*pi) * COS(coordsA(2)))

             ! Haversine formula
             slon = SIN((coordsA(1) - coordsB(1))/2)
             slat = SIN((coordsA(2) - coordsB(2))/2)

             dists(2) = SQRT(slat*slat + COS(coordsA(2))*COS(coordsB(2))*slon*slon)
             IF (dists(2)<=1.0) THEN
                dists(2) = 2.0 * r_earth* ASIN(dists(2))
             ELSE
                dists(2) = r_earth* pi
             END IF
             IF (dists(2)>thisobs_l%cradius(1)) THEN
                distflag = .FALSE.
             ELSE
                ! full squared distance
                distance2 = 0.0
                DO k = 1, thisobs%ncoord
                   distance2 = distance2 + dists(k)*dists(k)
                END DO
             END IF
          END IF
       END IF

    END IF norm


! ***************************************************************************
! *** Compute directional cut-off and support radii and set distance flag ***
! ***************************************************************************

    dflag: IF (distflag) THEN
       nrad: IF (thisobs_l%nradii == 2 .OR. (thisobs_l%nradii == 3 .AND. thisobs%disttype >= 10)) THEN

          IF ((thisobs_l%cradius(1) == thisobs_l%cradius(2)) .OR. &
               (thisobs_l%sradius(1) == thisobs_l%sradius(2))) THEN
             ! 2D isotropic case

             cradius2 = thisobs_l%cradius(1) * thisobs_l%cradius(1)

             IF (distance2 <= cradius2) THEN
                ! Set flag for valid observation
                checkdist = .TRUE.
                cnt_obs = cnt_obs + 1

                cradius = thisobs_l%cradius(1)
                sradius = thisobs_l%sradius(1)

                IF (debug>0) THEN
                   WRITE (*,*) '++ OMI-debug check_dist2_noniso: ', debug, &
                        '  2D isotropic with separately specified, but equal, radii'
                   WRITE (*,*) '++ OMI-debug check_dist2_noniso: ', debug, '  theta, cradius, sradius', &
                        theta*180/pi, cradius, sradius
                END IF
             END IF
          ELSE

             ! *** 2D anisotropic case: Use polar radius of ellipse in 2 dimensions ***

             ! Compute angle
             IF (dists(1) /= 0.0) THEN
                theta = ATAN(dists(2) / dists(1))
             ELSE
                theta = pi / 2.0
             END IF

             ! Compute radius in direction of theta
             IF (thisobs_l%cradius(1)>0.0 .OR. thisobs_l%cradius(2)>0.0) THEN
                cradius = thisobs_l%cradius(1) * thisobs_l%cradius(2) / &
                     SQRT( (thisobs_l%cradius(2)*COS(theta))**2  &
                     + (thisobs_l%cradius(1)*SIN(theta))**2 )
             ELSE
                cradius = 0.0
             END IF

             cradius2 = cradius * cradius

             IF (distance2 <= cradius2) THEN
                ! Set flag for valid observation
                checkdist = .TRUE.
                cnt_obs = cnt_obs + 1

                ! Compute support radius in direction of theta
                IF (thisobs_l%sradius(1)>0.0 .OR. thisobs_l%sradius(2)>0.0) THEN
                   sradius = thisobs_l%sradius(1) * thisobs_l%sradius(2) / &
                        SQRT( (thisobs_l%sradius(2)*COS(theta))**2 &
                        + (thisobs_l%sradius(1)*SIN(theta))**2 )
                ELSE
                   sradius = 0.0
                END IF

                IF (debug>0) THEN
                   WRITE (*,*) '++ OMI-debug check_dist2_noniso: ', debug, &
                        '  2D nonisotropic localization'
                   WRITE (*,*) '++ OMI-debug check_dist2_noniso: ', debug, '  theta, cradius, sradius', &
                        theta*180/pi, cradius, sradius
                END IF
             END IF

           END IF

       ELSE IF (thisobs_l%nradii == 3  .AND. thisobs%disttype < 10) THEN nrad

          ! To save computing time, we here distinguish whether 
          ! - the horizontal radii are equal and only direction 3 has a different radius
          ! - whether all radii are equal (isotropic but specified with separate radii)
          ! - the anisotropy is in all 3 dimensions (all radii different)

          aniso: IF ((thisobs_l%cradius(1) == thisobs_l%cradius(2)) .AND. &
               (thisobs_l%cradius(1) /= thisobs_l%cradius(3)) .AND. &
               (thisobs_l%sradius(1) == thisobs_l%sradius(2))) THEN

             ! *** Isotropic in horizontal direction, distinct radius in the third direction (vertical) ***

             dist_xy = SQRT(dists(1)*dists(1) + dists(2)*dists(2))

             ! 2D anisotropy: Polar radius of ellipse in 2 dimensions

             ! Compute angle
             IF (dist_xy /= 0.0) THEN
                theta = ATAN(dists(3) / dist_xy)
             ELSE
                theta = pi / 2.0
             END IF

             ! Compute radius in direction of theta
             IF (thisobs_l%cradius(1)>0.0 .OR. thisobs_l%cradius(3)>0.0) THEN
                cradius = thisobs_l%cradius(1) * thisobs_l%cradius(3) / &
                     SQRT( (thisobs_l%cradius(3)*COS(theta))**2  &
                     + (thisobs_l%cradius(1)*SIN(theta))**2 )
             ELSE
                cradius = 0.0
             END IF

             cradius2 = cradius * cradius

             IF (distance2 <= cradius2) THEN
                ! Set flag for valid observation
                checkdist = .TRUE.
                cnt_obs = cnt_obs + 1

                ! Compute support radius in direction of theta
                IF (thisobs_l%sradius(1)>0.0 .OR. thisobs_l%sradius(3)>0.0) THEN
                   sradius = thisobs_l%sradius(1) * thisobs_l%sradius(3) / &
                        SQRT( (thisobs_l%sradius(3)*COS(theta))**2 &
                        + (thisobs_l%sradius(1)*SIN(theta))**2 )
                ELSE
                   sradius = 0.0
                END IF
                   
                IF (debug>0) THEN
                   WRITE (*,*) '++ OMI-debug check_dist2_noniso: ', debug, &
                        '  3D: isotropic in directions 1 and 2, nonisotropic in direction 3'
                   WRITE (*,*) '++ OMI-debug check_dist2_noniso: ', debug, '  theta, cradius, sradius', &
                        theta*180/pi, cradius, sradius
                END IF
             END IF

          ELSEIF ((thisobs_l%cradius(1) == thisobs_l%cradius(2)) .AND. &
               (thisobs_l%cradius(1) == thisobs_l%cradius(3)) .AND. &
               (thisobs_l%sradius(1) == thisobs_l%sradius(2)) .AND. &
               (thisobs_l%sradius(2) == thisobs_l%sradius(3))) THEN aniso

             ! *** 3D isotropic case (all radii equal) ***

             cradius = thisobs_l%cradius(1)
             cradius2 = thisobs_l%cradius(1) * thisobs_l%cradius(1)
             sradius = thisobs_l%sradius(1)
             
             IF (distance2 <= cradius2) THEN
                ! Set flag for valid observation
                checkdist = .TRUE.
                cnt_obs = cnt_obs + 1
             END IF

             IF (debug>0) THEN
                WRITE (*,*) '++ OMI-debug check_dist2_noniso: ', debug, &
                     '  3D isotropic case specified with vector of radii'
                WRITE (*,*) '++ OMI-debug check_dist2_noniso: ', debug, '  theta, cradius, sradius', &
                     theta*180/pi, cradius, sradius
             END IF
          ELSE aniso

             ! *** general 3D anisotropic case ***

             ! Polar radius of ellipsoid in 3 dimensions

             ! Compute angle in x-y direction
             IF (dists(1) /= 0.0) THEN
                theta = ATAN(dists(2) / dists(1))
             ELSE
                theta = pi / 2.0
             END IF

             ! Distance in xy-plane
             dist_xy = SQRT(dists(1)**2 + dists(2)**2)

             ! Compute angle of xy-plane to z direction
             IF (dist_xy /= 0.0) THEN
                phi = ATAN(dists(3) / dist_xy)
             ELSE
                phi = 0.0
             END IF

             ! Compute radius in direction of theta
             IF (thisobs_l%cradius(1)>0.0 .OR. thisobs_l%cradius(2)>0.0 .OR. thisobs_l%cradius(3)>0.0) THEN
                cradius = thisobs_l%cradius(1) * thisobs_l%cradius(2) * thisobs_l%cradius(3) / &
                     SQRT( (thisobs_l%cradius(2)*thisobs_l%cradius(3)*COS(phi)*COS(theta))**2 &
                     + (thisobs_l%cradius(1)*thisobs_l%cradius(3)*COS(phi)*SIN(theta))**2 &
                     + (thisobs_l%cradius(1)*thisobs_l%cradius(2)*SIN(phi))**2 )
             ELSE
                cradius = 0.0
             END IF

             cradius2 = cradius * cradius

             IF (distance2 <= cradius2) THEN
                ! Set flag for valid observation
                checkdist = .TRUE.
                cnt_obs = cnt_obs + 1

                ! Compute support radius in direction of theta
                IF (thisobs_l%sradius(1)>0.0 .OR. thisobs_l%sradius(2)>0.0 .OR. thisobs_l%sradius(3)>0.0) THEN
                   sradius = thisobs_l%sradius(1) * thisobs_l%sradius(2) * thisobs_l%sradius(3) / &
                        SQRT( (thisobs_l%sradius(2)*thisobs_l%sradius(3)*COS(phi)*COS(theta))**2 &
                        + (thisobs_l%sradius(1)*thisobs_l%sradius(3)*COS(phi)*SIN(theta))**2 &
                        + (thisobs_l%sradius(1)*thisobs_l%sradius(2)*SIN(phi))**2 )
                ELSE
                   sradius = 0.0
                END IF

                IF (debug>0) THEN
                   WRITE (*,*) '++ OMI-debug check_dist2_noniso: ', debug, &
                        '  3D nonisotropic localization'
                   WRITE (*,*) '++ OMI-debug check_dist2_noniso: ', debug, '  theta, phi, distance, cradius, sradius', &
                        theta*180/pi, phi*180/pi, SQRT(distance2), cradius, sradius
                END IF
             END IF

          END IF aniso
       ELSEIF (thisobs_l%nradii == 1) THEN nrad
          cradius = thisobs_l%cradius(1)
          cradius2 = thisobs_l%cradius(1) * thisobs_l%cradius(1)
          sradius = thisobs_l%sradius(1)
          
          IF (distance2 <= cradius2) THEN
             ! Set flag for valid observation
             checkdist = .TRUE.
             cnt_obs = cnt_obs + 1
          END IF

       END IF nrad
    END IF dflag

  END SUBROUTINE PDAFomi_check_dist2_noniso




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
  SUBROUTINE PDAFomi_weights_l(verbose, nobs_l, ncols, locweight, cradius, sradius, &
        matA, ivar_obs_l, dist_l, weight_l)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: verbose        !< Verbosity flag
    INTEGER, INTENT(in) :: nobs_l         !< Number of local observations
    INTEGER, INTENT(in) :: ncols          !< 
    INTEGER, INTENT(in) :: locweight      !< Localization weight type
    REAL, INTENT(in)    :: cradius(:)     !< Localization cut-off radius
    REAL, INTENT(in)    :: sradius(:)     !< support radius for weight functions
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
            .OR. locweight == 4 .OR. locweight == 5. .OR. locweight == 15 &
            .OR. locweight == 16) THEN
          WRITE (*, '(a, 8x, a)') &
               'PDAFomi', '--- Use distance-dependent weight for observation errors'

          IF (locweight == 3 .OR. locweight == 15) THEN
             WRITE (*, '(a, 8x, a)') &
                  'PDAFomi', '--- Use regulated weight with mean error variance'
          ELSE IF (locweight == 4 .OR. locweight == 16) THEN
             WRITE (*, '(a, 8x, a)') &
                  'PDAFomi', '--- Use regulated weight with single-point error variance'
          END IF
       ELSE IF (locweight == 11 .OR. locweight == 26 .OR. locweight == 27) THEN
          WRITE (*, '(a, 8x, a)') &
               'PDAFomi', '--- Use distance-dependent weight for observed ensemble'
       END IF
    ENDIF


! *** Initialize weight vector

    IF (locweight == 0) THEN
       ! Uniform (unit) weighting
       wtype = 0
       rtype = 0
    ELSE IF (locweight == 1 .OR. locweight == 11) THEN
       ! Exponential weighting
       wtype = 1
       rtype = 0
    ELSE IF (locweight == 2 .OR. locweight == 3 .OR. locweight == 4 &
         .OR. locweight == 16 .OR. locweight == 17) THEN
       ! 5th-order polynomial (Gaspari&Cohn, 1999)
       wtype = 2

       IF (locweight == 3 .OR. locweight == 4) THEN
          ! Use regulated weight
          rtype = 1
       ELSE   
          ! No regulated weight
          rtype = 0
       END IF

    ELSE IF (locweight == 5 .OR. locweight == 15 .OR. locweight ==16) THEN
       ! 5th-order polynomial (Gaspari&Cohn, 1999)
       wtype = 3

       IF (locweight == 15 .OR. locweight == 16) THEN
          ! Use regulated weight
          rtype = 1
       ELSE   
          ! No regulated weight
          rtype = 0
       END IF

    END IF

    IF (locweight == 4 .OR. locweight == 16) THEN
       ! Allocate array for single observation point
       ALLOCATE(A_obs(1, ncols))
    END IF


    ! Control verbosity of PDAF_local_weight
    IF (verbose==1) THEN
       verbose_w = 1
    ELSE
       verbose_w = 0
    END IF

    IF (locweight /= 4) THEN
       ! All localizations except regulated weight based on variance at 
       ! single observation point
       DO i = 1, nobs_l

          ! set observation variance value
          var_obs_l = 1.0 / ivar_obs_l(i)

          CALL PDAF_local_weight(wtype, rtype, cradius(i), sradius(i), dist_l(i), &
               nobs_l, ncols, matA, var_obs_l, weight_l(i), verbose_w)

          verbose_w = 0

       END DO

    ELSE
       ! Regulated weight using variance at single observation point
       DO i = 1, nobs_l

          ! set observation variance value
          var_obs_l = 1.0 / ivar_obs_l(i)

          A_obs(1,:) = matA(i,:)
          CALL PDAF_local_weight(wtype, rtype, cradius(i), sradius(i), dist_l(i), &
               1, ncols, A_obs, var_obs_l, weight_l(i), verbose_w)

          verbose_w = 0

       END DO
    END IF

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
  SUBROUTINE PDAFomi_weights_l_sgnl(verbose, nobs_l, ncols, locweight, cradius, sradius, &
        matA, ivar_obs_l, dist_l, weight_l)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: verbose        !< Verbosity flag
    INTEGER, INTENT(in) :: nobs_l         !< Number of local observations
    INTEGER, INTENT(in) :: ncols          !< 
    INTEGER, INTENT(in) :: locweight      !< Localization weight type
    REAL, INTENT(in)    :: cradius        !< Localization cut-off radius
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
            .OR. locweight == 4 .OR. locweight == 5. .OR. locweight == 15 &
            .OR. locweight == 16) THEN
          WRITE (*, '(a, 8x, a)') &
               'PDAFomi', '--- Use distance-dependent weight for observation errors'

          IF (locweight == 3 .OR. locweight == 15) THEN
             WRITE (*, '(a, 8x, a)') &
                  'PDAFomi', '--- Use regulated weight with mean error variance'
          ELSE IF (locweight == 4 .OR. locweight == 16) THEN
             WRITE (*, '(a, 8x, a)') &
                  'PDAFomi', '--- Use regulated weight with single-point error variance'
          END IF
       ELSE IF (locweight == 11 .OR. locweight == 26 .OR. locweight == 27) THEN
          WRITE (*, '(a, 8x, a)') &
               'PDAFomi', '--- Use distance-dependent weight for observed ensemble'
       END IF
    ENDIF


! *** Initialize weight vector

    IF (locweight == 0) THEN
       ! Uniform (unit) weighting
       wtype = 0
       rtype = 0
    ELSE IF (locweight == 1 .OR. locweight == 11) THEN
       ! Exponential weighting
       wtype = 1
       rtype = 0
    ELSE IF (locweight == 2 .OR. locweight == 3 .OR. locweight == 4 &
         .OR. locweight == 16 .OR. locweight == 17) THEN
       ! 5th-order polynomial (Gaspari&Cohn, 1999)
       wtype = 2

       IF (locweight == 3 .OR. locweight == 4) THEN
          ! Use regulated weight
          rtype = 1
       ELSE   
          ! No regulated weight
          rtype = 0
       END IF

    ELSE IF (locweight == 5 .OR. locweight == 15 .OR. locweight ==16) THEN
       ! 5th-order polynomial (Gaspari&Cohn, 1999)
       wtype = 3

       IF (locweight == 15 .OR. locweight == 16) THEN
          ! Use regulated weight
          rtype = 1
       ELSE   
          ! No regulated weight
          rtype = 0
       END IF

    END IF

    IF (locweight == 4 .OR. locweight == 16) THEN
       ! Allocate array for single observation point
       ALLOCATE(A_obs(1, ncols))
    END IF


    ! Control verbosity of PDAF_local_weight
    IF (verbose==1) THEN
       verbose_w = 1
    ELSE
       verbose_w = 0
    END IF

    IF (locweight /= 4) THEN
       ! All localizations except regulated weight based on variance at 
       ! single observation point
       DO i = 1, nobs_l

          ! set observation variance value
          var_obs_l = 1.0 / ivar_obs_l(i)

          CALL PDAF_local_weight(wtype, rtype, cradius, sradius, dist_l(i), &
               nobs_l, ncols, matA, var_obs_l, weight_l(i), verbose_w)

          verbose_w = 0

       END DO

    ELSE
       ! Regulated weight using variance at single observation point
       DO i = 1, nobs_l

          ! set observation variance value
          var_obs_l = 1.0 / ivar_obs_l(i)

          A_obs(1,:) = matA(i,:)
          CALL PDAF_local_weight(wtype, rtype, cradius, sradius, dist_l(i), &
               1, ncols, A_obs, var_obs_l, weight_l(i), verbose_w)

          verbose_w = 0

       END DO
    END IF

    ! Print debug information
    IF (debug>0) THEN
       WRITE (*,*) '++ OMI-debug weights_l:     ', debug, 'thisobs_l%dim_obs_l', nobs_l
       WRITE (*,*) '++ OMI-debug weights_l:     ', debug, 'thisobs%distance_l', dist_l(1:nobs_l)
       WRITE (*,*) '++ OMI-debug weights_l:     ', debug, 'weight_l', weight_l(1:nobs_l)
       WRITE (*,*) '++ OMI-debug weights_l:     ', debug, 'thisobs_l%ivar_obs_l', ivar_obs_l(1:nobs_l)
    END IF
  
    IF (locweight == 4) DEALLOCATE(A_obs)

  END SUBROUTINE PDAFomi_weights_l_sgnl



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
         offset_obs_g, obsdims, map_obs_id

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

    ! Reset assim flag
    thisobs%doassim = 0

    IF (ALLOCATED(obs_f_all)) DEALLOCATE(obs_f_all)
    IF (ALLOCATED(obsdims)) DEALLOCATE(obsdims)
    IF (ALLOCATED(map_obs_id)) DEALLOCATE(map_obs_id)

    ! Reset counters over all observation types
    n_obstypes = 0
    obscnt = 0
    offset_obs = 0
    offset_obs_g = 0

    ! Reset flag for first call to local observations
    firstobs = 0

  END SUBROUTINE PDAFomi_deallocate_obs



!-------------------------------------------------------------------------------
!> Exclude observations for too high innovation
!!
!! The routine is called during the analysis step
!! on each local analysis domain. It checks the
!! size of the innovation and sets the observation
!! error to a high value if the squared innovation
!! exceeds a limit relative to the observation error
!! variance.
!!
!! __Revision history:__
!! * 2022-12 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_omit_by_inno_l(thisobs_l, thisobs, inno_l, obs_l_all, obsid, cnt_all, verbose)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    REAL, INTENT(in)    :: inno_l(:)         !< Input vector of observation innovation
    REAL, INTENT(in)    :: obs_l_all(:)      !< Input vector of local observations
    INTEGER, INTENT(in) :: obsid             !< ID of observation type
    INTEGER, INTENT(inout) :: cnt_all        !< Count of omitted observation over all types
    INTEGER, INTENT(in) :: verbose           !< Verbosity flag


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
             WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_omit_by_inno_l -- START   obs-ID', obsid
             WRITE (*,*) '++ OMI-debug omit_by_inno_l:', debug, 'limit for innovation', &
                  thisobs%inno_omit
             WRITE (*,*) '++ OMI-debug omit_by_inno_l:', debug, 'inno_omit_ivar', &
                  thisobs%inno_omit_ivar
             WRITE (*,*) '++ OMI-debug omit_by_inno_l:', debug, 'innovation_l', inno_l
          ENDIF

          ! Squared limit factor
          limit2 = thisobs%inno_omit * thisobs%inno_omit

          ! Check for observations to be excluded
          cnt = 0
          DO i = 1, thisobs_l%dim_obs_l

             ! Squared innovation
             inno2 = inno_l(i + thisobs_l%off_obs_l)* inno_l(i + thisobs_l%off_obs_l)

             IF (inno2 > limit2 * 1.0/thisobs_l%ivar_obs_l(i)) THEN

                IF (debug>0) THEN
                   WRITE (*,*) '++ OMI-debug omit_by_inno_l:', debug, 'omit: innovation:', &
                        inno_l(i + thisobs_l%off_obs_l), 'observation:', obs_l_all(i + thisobs_l%off_obs_l)
                END IF

                ! Exclude observation by increased its observation error
                thisobs_l%ivar_obs_l(i) = thisobs%inno_omit_ivar

                ! Count excluded obs
                cnt = cnt + 1
             END IF
          ENDDO

          IF (debug>0 .AND. cnt>0) THEN
             WRITE (*,*) '++ OMI-debug omit_by_inno_l:', debug, 'count of excluded obs.: ', cnt
             WRITE (*,*) '++ OMI-debug omit_by_inno_l:', debug, 'updated thisobs_l%ivar_obs_l ', &
                  thisobs_l%ivar_obs_l
          ENDIF

          IF (debug>0) &
               WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_omit_by_inno_l -- END   obs-ID', obsid

          cnt_all = cnt_all + cnt

       END IF

    ENDIF doassim

  END SUBROUTINE PDAFomi_omit_by_inno_l



!-------------------------------------------------------------------------------
!> Get statistics on local observations
!!
!! The routine is called in the update routine of
!! local filters and write statistics on locally
!! used and excluded observations.
!!
!! __Revision history:__
!! * 2023-03 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_obsstats_l(screen)

    USE MPI
    USE PDAFomi_obs_f, ONLY: ostats_omit
    USE PDAF_mod_filtermpi, &
         ONLY: COMM_filter, MPIerr

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: screen            !< Verbosity flag

! *** Local variables ***
    INTEGER :: ostats_omit_g(7)

    IF (npes_filter>1) THEN
       CALL MPI_Reduce(ostats_omit, ostats_omit_g, 5, MPI_INTEGER, MPI_SUM, &
            0, COMM_filter, MPIerr)
       CALL MPI_Reduce(ostats_omit(6:7), ostats_omit_g(6:7), 2, MPI_INTEGER, MPI_MAX, &
            0, COMM_filter, MPIerr)
    ELSE
       ! This is a work around for working with nullmpi.F90
       ostats_omit_g = ostats_omit
    END IF

    IF (mype == 0 .AND. screen > 0 .AND. ostats_omit_g(1)>0) THEN
       WRITE (*, '(a, 5x, a)') 'PDAFomi', '--- Global statistics for locally omitted observations:'
       WRITE (*, '(a, 8x, a, i10)') &
            'PDAFomi', 'Local domains with omitted observations:        ', ostats_omit_g(1)
       WRITE (*, '(a, 8x, a, i10)') &
            'PDAFomi', 'Local domains without omitted observations:     ', ostats_omit_g(2)
       WRITE (*, '(a, 8x, a, i10)') &
            'PDAFomi', 'Maximum number of locally omitted observations: ', ostats_omit_g(6)
       WRITE (*, '(a, 8x, a, i10)') &
            'PDAFomi', 'Maximum number of locally used observations:    ', ostats_omit_g(7)
       IF (ostats_omit_g(2) > 0) THEN
          WRITE (*, '(a, 8x, a, f9.1)') &
               'PDAFomi', 'Total avg. locally omitted observations:', &
               REAL(ostats_omit_g(3)) / REAL(ostats_omit_g(1) + ostats_omit_g(2))
          WRITE (*, '(a, 8x, a, f9.1)') &
               'PDAFomi', 'Total avg. locally used observations:   ', &
               REAL(ostats_omit_g(4)) / REAL(ostats_omit_g(1) + ostats_omit_g(2))
          WRITE (*, '(a, 8x, a, f9.1)') &
               'PDAFomi', 'Avg. omitted for domains with omitted observations:    ', &
               REAL(ostats_omit_g(3)) / REAL(ostats_omit_g(1))
          WRITE (*, '(a, 8x, a, f9.1)') &
               'PDAFomi', 'Avg. used for domains with omitted observations:       ', &
            REAL(ostats_omit_g(5)) / REAL(ostats_omit_g(1))
       END IF
    ELSEIF (mype == 0 .AND. screen > 0 .AND. ostats_omit_g(2)>0) THEN
       WRITE (*, '(a, 5x, a)') 'PDAFomi', '--- Global statistics for locally omitted observations:'
       WRITE (*, '(a, 9x, a)') &
            'PDAFomi', 'Zero observations omitted'
    END IF

  END SUBROUTINE PDAFomi_obsstats_l


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
         offset_obs_g, obsdims, map_obs_id

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

          ! Reset assim flag
          obs_f_all(i)%ptr%doassim = 0
       END DO

       IF (ALLOCATED(obs_f_all)) DEALLOCATE(obs_f_all)
       IF (ALLOCATED(obsdims)) DEALLOCATE(obsdims)
       IF (ALLOCATED(map_obs_id)) DEALLOCATE(map_obs_id)

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
