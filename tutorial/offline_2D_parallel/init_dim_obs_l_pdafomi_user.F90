!> Set dimension of local obs. vector and local obs. arrays
!!
!! This routine sets the number of local observations for the
!! current observation type for the local analysis domain
!! with coordinates COORD_l and a vector of localization cut-off
!! radii CRADIUS.
!! Further the routine initializes arrays for the index of a
!! local observation in the full observation vector and its 
!! corresponding distance.
!!
!! This is a user-coded variant of PDAFomi_init_dim_obs_l
!! which yields better performance than the generic routine
!! provided by PDAF-OMI. This example uses an isotropic
!! localization with Cartesian distance calculation and
!! no periodicity of the model domain.
!!
!! __Revision history:__
!! * 2024-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_dim_obs_l_pdafomi_user(thisobs_l, thisobs, coords_l, locweight, cradius, &
     sradius, cnt_obs_l_all)

  USE PDAFomi, &
       ONLY: obs_f, obs_l, PDAFomi_set_dim_obs_l, PDAFomi_set_localization

  IMPLICIT NONE

! *** Arguments ***
  TYPE(obs_f), INTENT(inout) :: thisobs    !< Data type with full observation
  TYPE(obs_l), TARGET, INTENT(inout) :: thisobs_l  !< Data type with local observation
  REAL, INTENT(in) :: coords_l(thisobs%ncoord)     !< Coordinates of current analysis domain
  INTEGER, INTENT(in) :: locweight         !< Type of localization function
  REAL, INTENT(in) :: cradius              !< Vector of localization cut-off radii
  REAL, INTENT(in) :: sradius              !< Vector of support radii of localization function
  INTEGER, INTENT(inout) :: cnt_obs_l_all  !< Local dimension of observation vector over all obs. types

! *** Local variables ***
  INTEGER :: i, cnt_l             ! Counters
  INTEGER :: cnt_obs              ! Counted number of local observations
  REAL :: distance2               ! squared distance
  REAL :: dists(thisobs%ncoord)   ! Distance vector between analysis point and observation
  REAL :: crad2                   ! square cut-off radius


  doassim: IF (thisobs%doassim == 1) THEN

! **************************************
! *** Store localization information ***
! **************************************

     CALL PDAFomi_set_localization(thisobs_l, cradius, sradius, locweight)


! ********************************
! *** Count local observations ***
! ********************************

     crad2 = thisobs_l%cradius(1) * thisobs_l%cradius(1)

     cnt_obs = 0
     countobs: DO i = 1, thisobs%dim_obs_f

        dists(1) = ABS(coords_l(1) - thisobs%ocoord_f(1, i))
        dists(2) = ABS(coords_l(2) - thisobs%ocoord_f(2, i))

        distance2 = dists(1)*dists(1) + dists(2)*dists(2)

        ! Count observations within squared radius
        IF (distance2 <= crad2) cnt_obs = cnt_obs + 1
     END DO countobs


! ************************************************
! *** Initialize local observation for PDAFomi ***
! ************************************************

     CALL PDAFomi_set_dim_obs_l(thisobs_l, thisobs, cnt_obs_l_all, cnt_obs)


! ****************************************************************
! *** Initialize local arrays in thisobs_l for local distances ***
! *** and indices of local obs. in full obs. vector            ***
! ****************************************************************

     haveobs: IF (cnt_obs>0) THEN

        cnt_l = 0
        scanobs: DO i = 1, thisobs%dim_obs_f

           dists(1) = ABS(coords_l(1) - thisobs%ocoord_f(1, i))
           dists(2) = ABS(coords_l(2) - thisobs%ocoord_f(2, i))

           distance2 = dists(1)*dists(1) + dists(2)*dists(2)

           IF (distance2 <= crad2) THEN
                 
              cnt_l = cnt_l + 1

              thisobs_l%id_obs_l(cnt_l) = i                      ! Index of local obs. in full obs. vector
              thisobs_l%distance_l(cnt_l) = SQRT(distance2)      ! distance
              thisobs_l%cradius_l(cnt_l) = thisobs_l%cradius(1)  ! cut-off radius
              thisobs_l%sradius_l(cnt_l) = thisobs_l%sradius(1)  ! support radius
           END IF
        END DO scanobs
     END IF haveobs

  END IF doassim

END SUBROUTINE init_dim_obs_l_pdafomi_user
