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

!> PDAF-OMI routines for determining local observations
!!
!! This module contains generic routines for several observation-related
!! operations for local filters. The routines are
!!
!! * PDAFomi_init_dim_obs_l \n
!!        Initialize dimension of local obs. vetor and arrays for
!!        local observations
!! * PDAFomi_init_dim_obs_l_iso \n
!!        Initialize dimension of local obs. vetor and arrays for
!!        local observations for isotropic localization
!! * PDAFomi_init_dim_obs_l_noniso \n
!!        Initialize dimension of local obs. vetor and arrays for
!!        local observations for nonisotropic localization
!! * PDAFomi_init_dim_obs_l_noniso_locweights \n
!!        Initialize dimension of local obs. vetor and arrays for
!!        local observations for nonisotropic localization 
!!        and different locweight for horizontal and vertical
!! * PDAFomi_check_dist2_loop \n
!!        Compute and check distance for isotropic localization
!! * PDAFomi_check_dist2_noniso_loop \n
!!        Compute and check distance for non-isotropic localization
!! * PDAFomi_set_localization \n
!!        Store localization parameters in OMI (for isotropic localization)
!! * PDAFomi_set_localization_noniso \n
!!        Store localization parameters in OMI (for non-isotropic localization)
!! * PDAFomi_set_dim_obs_l \n
!!        Register local observation with OMI
!! * PDAFomi_store_obs_l_index \n
!!        Store index, distance, cradius, and sradius of a local observation
!! * PDAFomi_store_obs_l_index_vdist \n
!!        Store index, distance, cradius, sradius, and vertical distance of
!!        a local observation for 2+1D factorized localization
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE PDAFomi_dim_obs_l

  USE PDAFomi_obs_f, ONLY: obs_f, r_earth, pi, debug, n_obstypes, error
  USE PDAFomi_obs_l, ONLY: obs_l, obs_l_all, firstobs, offset_obs_l
  USE PDAF_mod_filtermpi, ONLY: mype, npes_filter

  IMPLICIT NONE
  SAVE

  INTERFACE PDAFomi_init_dim_obs_l
     MODULE PROCEDURE PDAFomi_init_dim_obs_l_iso
     MODULE PROCEDURE PDAFomi_init_dim_obs_l_noniso
     MODULE PROCEDURE PDAFomi_init_dim_obs_l_noniso_locweights
  END INTERFACE

CONTAINS

!-------------------------------------------------------------------------------
!> Set dimension of local obs. vector and local obs. arrays
!!
!! This routine sets the number of local observations for the
!! current observation type for the local analysis domain
!! with coordinates COORD_l and localization cut-off radius CRADIUS.
!! Further the routine initializes arrays for the index of a
!! local observation in the full observation vector and its 
!! corresponding distance.
!! The operations are performed by calling the routine
!! PDAFomi_check_dist2_loop once for counting and a second time
!! for initializing the arrays.  
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_init_dim_obs_l_iso(thisobs_l, thisobs, coords_l, locweight, cradius, &
       sradius, cnt_obs_l_all)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    TYPE(obs_l), TARGET, INTENT(inout) :: thisobs_l  !< Data type with local observation
    REAL, INTENT(in) :: coords_l(:)          !< Coordinates of current analysis domain
    INTEGER, INTENT(in) :: locweight         !< Type of localization function
    REAL, INTENT(in) :: cradius              !< Localization cut-off radius
    REAL, INTENT(in) :: sradius              !< Support radius of localization function
    INTEGER, INTENT(inout) :: cnt_obs_l_all  !< Local dimension of current observation vector

! *** Local variables ***
    REAL :: maxcoords_l, mincoords_l         ! Min/Max domain coordinates to check geographic coords
    REAL :: maxocoords_l, minocoords_l       ! Min/Max observation coordinates to check geographic coords
    INTEGER :: cnt_obs                       ! Counter for valid local observations


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


! **************************************
! *** Store localization information ***
! **************************************

       thisobs_l%locweight = locweight

       ! Allocate vectors for localization radii and store their values
       ! For isotropic localization the size of the arrays is just 1
       IF (ALLOCATED(thisobs_l%cradius)) DEALLOCATE(thisobs_l%cradius)
       ALLOCATE(thisobs_l%cradius(1))
       IF (ALLOCATED(thisobs_l%sradius)) DEALLOCATE(thisobs_l%sradius)
       ALLOCATE(thisobs_l%sradius(1))

       thisobs_l%nradii = 1
       thisobs_l%cradius(1) = cradius
       thisobs_l%sradius(1) = sradius


! **************************************
! *** Count valid local observations ***
! **************************************

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

          WRITE (*,*) '++ OMI-debug init_dim_obs_l: ', debug, '  thisobs%ncoord', thisobs%ncoord
          WRITE (*,*) '++ OMI-debug init_dim_obs_l: ', debug, '  thisobs_l%cradius', thisobs_l%cradius
          WRITE (*,*) '++ OMI-debug init_dim_obs_l: ', debug, '  Check for observations within radius'
       END IF

       cnt_obs = 0
       CALL PDAFomi_check_dist2_loop(thisobs_l, thisobs, coords_l, cnt_obs, 1)


! ************************************************
! *** Initialize local observation for PDAFomi ***
! ************************************************

       CALL PDAFomi_set_dim_obs_l(thisobs_l, thisobs, cnt_obs_l_all, cnt_obs)


! ************************************************************
! *** Initialize internal local arrays for local distances ***
! *** and indices of local obs. in full obs. vector        ***
! ************************************************************

       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, &
            '   PDAFomi_init_dim_obs_l -- initialize local observation arrays'

       ! Count local observations and initialize index and distance arrays
       IF (thisobs_l%dim_obs_l>0) THEN
          cnt_obs = 0
          CALL PDAFomi_check_dist2_loop(thisobs_l, thisobs, coords_l, cnt_obs, 2)
       END IF

       ! Print debug information
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug init_dim_obs_l:', debug, '  thisobs_l%dim_obs_l', thisobs_l%dim_obs_l
          IF (thisobs_l%dim_obs_l>0) THEN
             WRITE (*,*) '++ OMI-debug init_dim_obs_l:', debug, '  thisobs_l%id_obs_l', thisobs_l%id_obs_l
             WRITE (*,*) '++ OMI-debug init_dim_obs_l:', debug, '  thisobs_l%distance_l', thisobs_l%distance_l
          END IF
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_init_dim_obs_l -- END'
       END IF

    END IF doassim

  END SUBROUTINE PDAFomi_init_dim_obs_l_iso


!-------------------------------------------------------------------------------
!> Check distance in case of isotropic localization
!!
!! This routine computes the distance between the location of
!! a local analysis domains and all full observations and checks
!! whether the observations lies within the localization radius.
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
  SUBROUTINE PDAFomi_check_dist2_loop(thisobs_l, thisobs, coordsA, cnt_obs, mode)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    TYPE(obs_f), INTENT(in) :: thisobs       !< Data type with full observation
    REAL, INTENT(in) :: coordsA(:)           !< Coordinates of current analysis domain (ncoord)
    INTEGER, INTENT(inout) :: cnt_obs        !< Count number of local observations
    INTEGER, INTENT(in) :: mode              !< 1: count local observations
                                             !< 2: initialize local arrays

! *** Local variables ***
    INTEGER :: i, k                 ! Counters
    INTEGER :: verbose              ! verbosity flag
    INTEGER :: domsize              ! Flag whether domainsize is set
    LOGICAL :: distflag             ! Flag whether distance in a coordinate direction is within cradius
    REAL :: slon, slat              ! sine of distance in longitude or latitude
    REAL :: distance2               ! square distance
    REAL :: cradius2                ! squared localization cut-off radius
    REAL :: dists(thisobs%ncoord)   ! Distance vector between analysis point and observation
    REAL :: coordsB(thisobs%ncoord) ! Array for coordinates of a single observation


! **********************
! *** Initialization ***
! **********************

    scancount: DO i = 1, thisobs%dim_obs_f

       ! Initialize distance flag
       distflag = .TRUE.

       verbose = i

       coordsB = thisobs%ocoord_f(1:thisobs%ncoord, i)


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

             ! Increment counter
             cnt_obs = cnt_obs + 1

             IF (debug>0) THEN
                WRITE (*,*) '++ OMI-debug cnt_dim_obs_l: ', debug, &
                     '  valid observation with coordinates', thisobs%ocoord_f(1:thisobs%ncoord, i)
             END IF

             IF (mode == 2) THEN
                ! For internal storage (use in prodRinvA_l)
                thisobs_l%id_obs_l(cnt_obs) = i                       ! node index
                thisobs_l%distance_l(cnt_obs) = SQRT(distance2)       ! distance
                thisobs_l%cradius_l(cnt_obs) = thisobs_l%cradius(1)   ! isotropic cut-off radius
                thisobs_l%sradius_l(cnt_obs) = thisobs_l%sradius(1)   ! isotropic support radius
             END IF

          END IF
       END IF
    END DO scancount

  END SUBROUTINE PDAFomi_check_dist2_loop




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
  SUBROUTINE PDAFomi_init_dim_obs_l_noniso(thisobs_l, thisobs, coords_l, locweight, cradius, &
       sradius, cnt_obs_l_all)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    TYPE(obs_l), TARGET, INTENT(inout) :: thisobs_l  !< Data type with local observation
    REAL, INTENT(in) :: coords_l(:)          !< Coordinates of current analysis domain
    INTEGER, INTENT(in) :: locweight         !< Type of localization function
    REAL, INTENT(in) :: cradius(:)           !< Vector of localization cut-off radii
    REAL, INTENT(in) :: sradius(:)           !< Vector of support radii of localization function
    INTEGER, INTENT(inout) :: cnt_obs_l_all  !< Local dimension of current observation vector

! *** Local variables ***
    REAL :: maxcoords_l, mincoords_l         ! Min/Max domain coordinates to check geographic coords
    REAL :: maxocoords_l, minocoords_l       ! Min/Max observation coordinates to check geographic coords
    INTEGER :: cnt_obs                       ! Counter for valid local observations


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


! **************************************
! *** Count valid local observations ***
! **************************************

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

          WRITE (*,*) '++ OMI-debug init_dim_obs_l_noniso: ', debug, '  thisobs%ncoord', thisobs%ncoord
          WRITE (*,*) '++ OMI-debug init_dim_obs_l_noniso: ', debug, '  thisobs_l%cradius', thisobs_l%cradius
          WRITE (*,*) '++ OMI-debug init_dim_obs_l_noniso: ', debug, '  Check for observations within radius'
       END IF

       cnt_obs = 0
       IF (thisobs_l%nradii==1) THEN
          ! 1D but with radius specified as array
          CALL PDAFomi_check_dist2_loop(thisobs_l, thisobs, coords_l, cnt_obs, 1)
       ELSEIF (thisobs_l%nradii==2 .OR. thisobs_l%nradii==3) THEN
          ! Nonisotropic in 2 or 3 dimensions
          CALL PDAFomi_check_dist2_noniso_loop(thisobs_l, thisobs, coords_l, cnt_obs, 1)
       ELSE
          WRITE (*,*) '+++++ ERROR PDAF-OMI: nonisotropic localization is only possible in 1, 2 or 3 dimensions'
          error = 10
       END IF


! ************************************************
! *** Initialize local observation for PDAFomi ***
! ************************************************

       CALL PDAFomi_set_dim_obs_l(thisobs_l, thisobs, cnt_obs_l_all, cnt_obs)


! ************************************************************
! *** Initialize internal local arrays for local distances ***
! *** and indices of local obs. in full obs. vector        ***
! ************************************************************

       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, &
            '   PDAFomi_init_dim_obs_l_noniso -- initialize local observation arrays'

       ! Count local observations and initialize index and distance arrays
       IF (thisobs_l%dim_obs_l>0) THEN

          cnt_obs = 0
          IF (thisobs_l%nradii==1) THEN
             ! 1D but with radius specified as array
             CALL PDAFomi_check_dist2_loop(thisobs_l, thisobs, coords_l, cnt_obs, 2)
          ELSEIF (thisobs_l%nradii==2 .OR. thisobs_l%nradii==3) THEN
             ! Nonisotropic in 2 or 3 dimensions
             CALL PDAFomi_check_dist2_noniso_loop(thisobs_l, thisobs, coords_l, cnt_obs, 2)
          ELSE
             WRITE (*,*) '+++++ ERROR PDAF-OMI: nonisotropic localization is only possible in 1, 2 or 3 dimensions'
             error = 11
          END IF
       END IF

       ! Print debug information
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug init_dim_obs_l_noniso:', debug, '  thisobs_l%dim_obs_l', thisobs_l%dim_obs_l
          IF (thisobs_l%dim_obs_l>0) THEN
             WRITE (*,*) '++ OMI-debug init_dim_obs_l_noniso:', debug, '  thisobs_l%id_obs_l', thisobs_l%id_obs_l
             WRITE (*,*) '++ OMI-debug init_dim_obs_l_noniso:', debug, '  thisobs_l%distance_l', thisobs_l%distance_l
          END IF
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_init_dim_obs_l_noniso -- END'
       END IF

    END IF doassim

  END SUBROUTINE PDAFomi_init_dim_obs_l_noniso




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
  SUBROUTINE PDAFomi_init_dim_obs_l_noniso_locweights(thisobs_l, thisobs, coords_l, locweights, cradius, &
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
       WRITE (*,*) '+++++ ERROR PDAF-OMI: Input for locweight in horizontal and vertical directions needs size 2'
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
    CALL PDAFomi_init_dim_obs_l_noniso(thisobs_l, thisobs, coords_l, locweights(1), cradius, &
         sradius, cnt_obs_l)

    IF (debug>0) &
         WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_init_dim_obs_l_noniso_locweights -- END'

  END SUBROUTINE PDAFomi_init_dim_obs_l_noniso_locweights


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
  SUBROUTINE PDAFomi_check_dist2_noniso_loop(thisobs_l, thisobs, coordsA, cnt_obs, mode)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(in) :: thisobs    !< Data type with full observation
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    REAL, INTENT(in) :: coordsA(:)        !< Coordinates of current analysis domain (ncoord)
    INTEGER, INTENT(inout) :: cnt_obs     !< Count number of local observations
    INTEGER, INTENT(in) :: mode              !< 1: count local observations
                                             !< 2: initialize local arrays

! *** Local variables ***
    INTEGER :: i, k                 ! Counters
    INTEGER :: verbose              ! verbosity flag
    INTEGER :: domsize              ! Flag whether domainsize is set
    LOGICAL :: distflag             ! Flag whether distance in a coordinate direction is within cradius
    REAL :: slon, slat              ! sine of distance in longitude or latitude
    REAL :: distance2               ! square distance
    REAL :: cradius2                ! cut-off radius on ellipse or ellipsoid
    REAL :: phi, theta              ! Angles in ellipse or ellipsoid
    REAL :: dist_xy                 ! Distance in xy-plan in 3D case
    REAL :: dists(thisobs%ncoord)   ! Distance vector between analysis point and observation
    REAL :: coordsB(thisobs%ncoord) ! Array for coordinates of a single observation
    REAL :: cradius                 ! Directional cut-off radius
    REAL :: sradius                 ! Directional support radius
    LOGICAL :: checkdist            ! Flag whether distance is within cut-off radius


! **********************
! *** Initialization ***
! **********************

    scancount: DO i = 1, thisobs%dim_obs_f

       ! Initialize distance flag
       checkdist = .FALSE.    ! Whether an observation lies within the local box
       distflag = .TRUE.      ! Whether an observation lies within the local radius (ellipse, ellipsoid)

       ! Verbosity flag
       verbose = i

       ! Observation coordinates
       coordsB = thisobs%ocoord_f(1:thisobs%ncoord, i)


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

          IF (mode==2 .AND. checkdist) THEN
             ! For internal storage (use in prodRinvA_l and likelihood_l)
             thisobs_l%id_obs_l(cnt_obs) = i                       ! node index
             thisobs_l%distance_l(cnt_obs) = SQRT(distance2)       ! distance
             thisobs_l%cradius_l(cnt_obs) = cradius                ! directional cut-off radius
             thisobs_l%sradius_l(cnt_obs) = sradius                ! directional support radius
             IF (thisobs_l%locweight_v>0 .AND. thisobs_l%nradii==3) THEN
                thisobs_l%dist_l_v(cnt_obs) = dists(3)             ! distance in vertical direction
             END if
          END IF
    END IF dflag

 END DO scancount

  END SUBROUTINE PDAFomi_check_dist2_noniso_loop



!-------------------------------------------------------------------------------
!> Set localization parameters for isotropic localization
!!
!! This routine stores localization information (locweight, cradius, sradius)
!! in OMI and allocates local arrays for cradius and sradius. This variant 
!! is for isotropic localization. The routine is used by user-supplied 
!! implementations of PDAFomi_init_dim_obs_l.
!!
!! The routine is called by all filter processes.
!!
!! __Revision history:__
!! * 2024-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_set_localization(thisobs_l, cradius, sradius, locweight)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    REAL, INTENT(in) :: cradius         !< Localization cut-off radius
    REAL, INTENT(in) :: sradius         !< Support radius of localization function
    INTEGER, INTENT(in) :: locweight    !< Type of localization function


! *** Allocate vectors for localization radii ***

       IF (ALLOCATED(thisobs_l%cradius)) DEALLOCATE(thisobs_l%cradius)
       ALLOCATE(thisobs_l%cradius(1))
       IF (ALLOCATED(thisobs_l%sradius)) DEALLOCATE(thisobs_l%sradius)
       ALLOCATE(thisobs_l%sradius(1))

       thisobs_l%locweight = locweight
       thisobs_l%nradii = 1
       thisobs_l%cradius(:) = cradius
       thisobs_l%sradius(:) = sradius

  END SUBROUTINE PDAFomi_set_localization



!-------------------------------------------------------------------------------
!> Set localization parameters for non-isotropic localization
!!
!! This routine stores localization information (locweight, cradius, sradius)
!! in OMI and allocates local arrays for cradius and sradius. This variant 
!! is for non-isotropic localization. The routine is used by user-supplied 
!! implementations of PDAFomi_init_dim_obs_l.
!!
!! The routine is called by all filter processes.
!!
!! __Revision history:__
!! * 2024-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_set_localization_noniso(thisobs_l, nradii, cradius, sradius, locweight, locweight_v)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    INTEGER, INTENT(in) :: nradii            !< Number of radii to consider for localization
    REAL, INTENT(in) :: cradius(nradii)      !< Localization cut-off radius
    REAL, INTENT(in) :: sradius(nradii)      !< Support radius of localization function
    INTEGER, INTENT(in) :: locweight         !< Type of localization function
    INTEGER, INTENT(in) :: locweight_v       !< Type of localization function in vertical direction (only for nradii=3)



! *** Allocate vectors for localization radii ***

       IF (ALLOCATED(thisobs_l%cradius)) DEALLOCATE(thisobs_l%cradius)
       ALLOCATE(thisobs_l%cradius(nradii))
       IF (ALLOCATED(thisobs_l%sradius)) DEALLOCATE(thisobs_l%sradius)
       ALLOCATE(thisobs_l%sradius(nradii))

       thisobs_l%locweight = locweight
       thisobs_l%nradii = nradii
       thisobs_l%cradius(1:nradii) = cradius(1:nradii)
       thisobs_l%sradius(1:nradii) = sradius(1:nradii)
       IF (nradii==3) thisobs_l%locweight_v = locweight_v

     END SUBROUTINE PDAFomi_set_localization_noniso


!-------------------------------------------------------------------------------
!> Initialization for dim_obs_l
!!
!! This routine initializes information on local observation vectors.
!! It is used by a user-supplied implementations of PDAFomi_init_dim_obs_l.
!!
!! The routine is called by all filter processes.
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_set_dim_obs_l(thisobs_l, thisobs, cnt_obs_l_all, cnt_obs_l)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
    TYPE(obs_l), TARGET, INTENT(inout) :: thisobs_l  !< Data type with local observation
    INTEGER, INTENT(inout) :: cnt_obs_l_all  !< Local dimension of observation vector over all obs. types
    INTEGER, INTENT(inout) :: cnt_obs_l      !< Local dimension of single observation type vector


    ! Store ID of first observation type that calls the routine
    ! This is reset in PDAFomi_deallocate_obs
    IF (firstobs == 0) THEN
       firstobs = thisobs%obsid
    END IF
 
    ! Reset offset of currrent observation in overall local obs. vector
    IF (thisobs%obsid == firstobs) THEN
       offset_obs_l = 0
       cnt_obs_l_all = 0
    END IF

    ! Store offset
    thisobs_l%off_obs_l = offset_obs_l

    ! Initialize pointer array
    IF (thisobs%obsid == firstobs) THEN
       IF (ALLOCATED(obs_l_all)) DEALLOCATE(obs_l_all)
       ALLOCATE(obs_l_all(n_obstypes))
    END IF

    ! Set pointer to current observation
    obs_l_all(thisobs%obsid)%ptr => thisobs_l

    ! Store local observation dimension and increment offset
    thisobs_l%dim_obs_l = cnt_obs_l
    offset_obs_l = offset_obs_l + cnt_obs_l
    cnt_obs_l_all = cnt_obs_l_all + cnt_obs_l

    ! Allocate arrays to store information on local observations
    IF (ALLOCATED(thisobs_l%id_obs_l)) DEALLOCATE(thisobs_l%id_obs_l)
    IF (ALLOCATED(thisobs_l%distance_l)) DEALLOCATE(thisobs_l%distance_l)
    IF (ALLOCATED(thisobs_l%cradius_l)) DEALLOCATE(thisobs_l%cradius_l)
    IF (ALLOCATED(thisobs_l%sradius_l)) DEALLOCATE(thisobs_l%sradius_l)

    haveobs: IF (cnt_obs_l>0) THEN
       ALLOCATE(thisobs_l%id_obs_l(cnt_obs_l))
       ALLOCATE(thisobs_l%distance_l(cnt_obs_l))
       ALLOCATE(thisobs_l%cradius_l(cnt_obs_l))
       ALLOCATE(thisobs_l%sradius_l(cnt_obs_l))
       IF (thisobs_l%locweight_v>0) THEN
          IF (ALLOCATED(thisobs_l%dist_l_v)) DEALLOCATE(thisobs_l%dist_l_v)
          ALLOCATE(thisobs_l%dist_l_v(cnt_obs_l))
       END IF

    ELSE
       ALLOCATE(thisobs_l%id_obs_l(1))
       ALLOCATE(thisobs_l%distance_l(1))
       ALLOCATE(thisobs_l%cradius_l(1))
       ALLOCATE(thisobs_l%sradius_l(1))
       IF (ALLOCATED(thisobs_l%dist_l_v)) DEALLOCATE(thisobs_l%dist_l_v)
       ALLOCATE(thisobs_l%dist_l_v(1))
    END IF haveobs

  END SUBROUTINE PDAFomi_set_dim_obs_l



!-------------------------------------------------------------------------------
!> Store local index, distance and radii
!!
!! This routine stores the mapping index between the global and local
!! observation vectors, the distance and the cradius and sradius
!! for a single observations in OMI. This variant is for non-factorized
!! localization. The routine is used by user-supplied implementations 
!! of PDAFomi_init_dim_obs_l.
!!
!! The routine is called by all filter processes.
!!
!! __Revision history:__
!! * 2024-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_store_obs_l_index(thisobs_l, idx, id_obs_l, distance, &
       cradius_l, sradius_l)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    INTEGER, INTENT(in) :: idx       !< Element of local observation array to be filled
    INTEGER, INTENT(in) :: id_obs_l  !< Index of local observation in full observation array
    REAL, INTENT(in) :: distance     !< Distance between local analysis domain and observation
    REAL, INTENT(in) :: cradius_l    !< cut-off radius for this local observation
                                     !  (directional radius in case of non-isotropic localization)
    REAL, INTENT(in) :: sradius_l    !< support radius for this local observation
                                     !  (directional radius in case of non-isotropic localization)


! *** Store values ***

    thisobs_l%id_obs_l(idx)   = id_obs_l   ! element of local obs. vector in full obs. vector
    thisobs_l%distance_l(idx) = distance   ! distance
    thisobs_l%cradius_l(idx)  = cradius_l  ! cut-off radius
    thisobs_l%sradius_l(idx)  = sradius_l  ! support radius


  END SUBROUTINE PDAFomi_store_obs_l_index



!-------------------------------------------------------------------------------
!> Store local index, dsitance and radii for factorized localization
!!
!! This routine stores the mapping index between the global and local
!! observation vectors, the distance and the cradius and sradius
!! for a single observations in OMI. This variant is for 2+1D factorized
!! localization. The routine is used by user-supplied implementations 
!! of PDAFomi_init_dim_obs_l.
!!
!! The routine is called by all filter processes.
!!
!! __Revision history:__
!! * 2024-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_store_obs_l_index_vdist(thisobs_l, idx, id_obs_l, distance, &
       cradius_l, sradius_l, vdist)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_l), INTENT(inout) :: thisobs_l  !< Data type with local observation
    INTEGER, INTENT(in) :: idx       !< Element of local observation array to be filled
    INTEGER, INTENT(in) :: id_obs_l  !< Index of local observation in full observation array
    REAL, INTENT(in) :: distance     !< Distance between local analysis domain and observation
    REAL, INTENT(in) :: cradius_l    !< cut-off radius for this local observation
                                     !  (directional radius in case of non-isotropic localization)
    REAL, INTENT(in) :: sradius_l    !< support radius for this local observation
                                     !  (directional radius in case of non-isotropic localization)
    REAL, INTENT(in) :: vdist        !< support radius in vertical direction for 2+1D factorized localization


! *** Store values ***

    thisobs_l%id_obs_l(idx)   = id_obs_l   ! element of local obs. vector in full obs. vector
    thisobs_l%distance_l(idx) = distance   ! distance
    thisobs_l%cradius_l(idx)  = cradius_l  ! cut-off radius
    thisobs_l%sradius_l(idx)  = sradius_l  ! support radius
    IF (thisobs_l%locweight_v>0 .AND. thisobs_l%nradii==3) THEN
       thisobs_l%dist_l_v(idx) = vdist    ! distance in vertical direction
    END if


  END SUBROUTINE PDAFomi_store_obs_l_index_vdist

END MODULE PDAFomi_dim_obs_l
