!$Id: localize_covar_pdaf.F90 1676 2016-12-10 14:55:45Z lnerger $
!BOP
!
! !ROUTINE: localize_covar_pdaf --- apply localization matrix in LEnKF
!
! !INTERFACE:
SUBROUTINE localize_covar_pdaf(dim_p, dim_obs, HP, HPH)

! !DESCRIPTION:
! User-supplied routine for PDAF (local EnKF)
!
! This routine applies a localization matrix B
! to the matrices HP and HPH^T of the localized EnKF.
!
! !REVISION HISTORY:
! 2016-11 - Lars Nerger
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
!    ONLY: cradius, sradius, locweight, obs_pdaf2nc
! hcp
! we need to store the coordinates of the state vector 
! and obs array in longxy, latixy, and longxy_obs, latixy_obs
! respectively
#if defined CLMSA
    ONLY: cradius, sradius, locweight, obs_pdaf2nc, &
          longxy, latixy, longxy_obs, latixy_obs
!hc  end
#else
    ONLY: cradius, sradius, locweight, obs_pdaf2nc
#endif
!fin hcp

  USE mod_read_obs,&
  ONLY: x_idx_obs_nc, y_idx_obs_nc, z_idx_obs_nc
#if defined CLMSA
   USE shr_kind_mod , only : r8 => shr_kind_r8
   USE mod_read_obs, ONLY: clmobs_lon
   USE mod_read_obs, ONLY: clmobs_lat
   USE enkf_clm_mod, ONLY: init_clm_l_size, clmupdate_T
   USE enkf_clm_mod, ONLY: clm_begc
   USE enkf_clm_mod, ONLY: clm_endc
   USE enkf_clm_mod, ONLY: state_pdaf2clm_c_p
   USE enkf_clm_mod, ONLY: state_pdaf2clm_j_p
   USE mod_parallel_pdaf, ONLY: filterpe
#endif
   USE mod_parallel_pdaf, ONLY: mype_world
   USE mod_parallel_pdaf, ONLY: abort_parallel
!fin hcp

  USE mod_tsmp,&
#if defined CLMSA
    ONLY:  tag_model_parflow, tag_model_clm, &
           enkf_subvecsize, model
#else
    ONLY: tag_model_parflow, tag_model_clm, &
          enkf_subvecsize, &
          nx_glob, ny_glob, nz_glob, &
          xcoord, ycoord, zcoord, &
          xcoord_fortran, ycoord_fortran, zcoord_fortran, &
          model
#endif
  USE mod_tsmp, ONLY: point_obs

#ifdef CLMFIVE
  USE GridcellType, ONLY: grc
  USE ColumnType, ONLY : col
#else
  USE clmtype, ONLY : clm3
#endif

  USE, INTRINSIC :: iso_c_binding, ONLY: C_F_POINTER

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p                  ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_obs                ! number of observations
  REAL, INTENT(inout) :: HP(dim_obs, dim_p)     ! Matrix HP
  REAL, INTENT(inout) :: HPH(dim_obs, dim_obs)  ! Matrix HPH

! *** local variables ***
  INTEGER :: i, j          ! Index of observation component
  REAL    :: dx,dy,distance  ! Distance between points in the domain 
  REAL    :: weight        ! Localization weight
  REAL    :: tmp(1,1)= 1.0 ! Temporary, but unused array
  INTEGER :: wtype         ! Type of weight function
  INTEGER :: rtype         ! Type of weight regulation
!hcp
! dim_l is the number of soil layers
! ncellxy is the number of cell in xy (not counted z) plane
! for each processor
#if defined CLMSA
  ! INTEGER :: dim_l
  ! INTEGER :: ncellxy
  ! INTEGER :: k
  real(r8), pointer :: lon(:)
  real(r8), pointer :: lat(:)
  integer, pointer :: mycgridcell(:) !Pointer for CLM3.5/CLM5.0 col->gridcell index arrays
  REAL :: yhalf
#endif
  INTEGER :: icoord

! **********************
! *** INITIALIZATION ***
! **********************

  ! Screen output
  WRITE (*,'(8x, a)') &
       '--- Apply covariance localization'
  WRITE (*, '(12x, a, 1x, f12.2)') &
       '--- Local influence radius', cradius

  IF (locweight == 1) THEN
     WRITE (*, '(12x, a)') &
          '--- Use exponential distance-dependent weight'
  ELSE IF (locweight == 2) THEN
     WRITE (*, '(12x, a)') &
          '--- Use distance-dependent weight by 5th-order polynomial'
  END IF

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



#ifndef CLMSA
  IF(model==tag_model_parflow)THEN
    call C_F_POINTER(xcoord,xcoord_fortran,[enkf_subvecsize])
    call C_F_POINTER(ycoord,ycoord_fortran,[enkf_subvecsize])
    call C_F_POINTER(zcoord,zcoord_fortran,[enkf_subvecsize])

    ! Check that point observations are used
    if (.not. point_obs .eq. 1) then
      print *, "TSMP-PDAF mype(w)=", mype_world, ": ERROR(2) `point_obs.eq.1` needed for using obs_pdaf2nc."
      call abort_parallel()
    end if

    ! localize HP
    DO j = 1, dim_obs
       DO i = 1, dim_p

         ! Index in coordinate array only spans `enkf_subvecsize`.
         !
         ! Necessary condition: the full state vector consists of
         ! sections of size `enkf_subvecsize`, where each section
         ! corresponds to a single coordinate array.
         icoord = modulo(i,enkf_subvecsize)

         dx = abs(x_idx_obs_nc(obs_pdaf2nc(j)) - int(xcoord_fortran(icoord))-1)
         dy = abs(y_idx_obs_nc(obs_pdaf2nc(j)) - int(ycoord_fortran(icoord))-1)
         distance = sqrt(real(dx)**2 + real(dy)**2)
    
         ! Compute weight
         CALL PDAF_local_weight(wtype, rtype, cradius, sradius, distance, 1, 1, tmp, 1.0, weight, 0)
    
         ! Apply localization
         HP(j,i) = weight * HP(j,i)

       END DO
    END DO
    
    ! localize HPH^T
    DO j = 1, dim_obs
       DO i = 1, dim_obs
    
         ! Compute distance
         dx = abs(x_idx_obs_nc(obs_pdaf2nc(j)) - x_idx_obs_nc(obs_pdaf2nc(i)))
         dy = abs(y_idx_obs_nc(obs_pdaf2nc(j)) - y_idx_obs_nc(obs_pdaf2nc(i)))
         distance = sqrt(real(dx)**2 + real(dy)**2)
    
         ! Compute weight
         CALL PDAF_local_weight(wtype, rtype, cradius, sradius, distance, 1, 1, tmp, 1.0, weight, 0)
    
         ! Apply localization
         HPH(j,i) = weight * HPH(j,i)

       END DO
    END DO
    
  ENDIF ! model==tag_model_parflow
#endif

!by hcp to computer the localized covariance matrix in CLMSA case
#if defined CLMSA

   IF(model==tag_model_clm)THEN

    ! localize HP
    ! ----------- 

    ! Lon/Lat information from CLM
#ifdef CLMFIVE
    ! Obtain CLM lon/lat information
    lon   => grc%londeg
    lat   => grc%latdeg
    ! Obtain CLM column-gridcell information
    mycgridcell => col%gridcell
#else
    lon   => clm3%g%londeg
    lat   => clm3%g%latdeg
    mycgridcell => clm3%g%l%c%gridcell
#endif

    DO j = 1, dim_obs
      DO i = 1, dim_p

         ! Distance: obs - state

         ! Units: Index numbering
         ! dx = abs(longxy_obs(j) - longxy(i)-1)
         ! dy = abs(latixy_obs(j) - latixy(i)-1)
         ! dx = abs(longxy_obs(j) - longxy(state_pdaf2clm_c_p(i)))
         ! dy = abs(latixy_obs(j) - latixy(state_pdaf2clm_c_p(i)))

         ! Units: lat/lon
         dx = abs(clmobs_lon(obs_pdaf2nc(j)) - lon(mycgridcell(state_pdaf2clm_c_p(i))))
         dy = abs(clmobs_lat(obs_pdaf2nc(j)) - lat(mycgridcell(state_pdaf2clm_c_p(i))))

         ! Check for longitude differences that may yield differences
         ! larger than 180deg depending on conventions. Example:
         ! crossing the prime meridian (lon=0deg), when convention is
         ! all-positive lons
         IF (dx > 180.0) THEN
           dx = 360.0 - dx
         END IF

         ! Intermediate latitude
         yhalf = ( clmobs_lat(obs_pdaf2nc(j)) + lat(mycgridcell(state_pdaf2clm_c_p(i))) ) / 2.0

         ! Latitude-dependent factor for longitude difference
         dx = dx * cos(yhalf * 3.14159265358979323846 / 180.0)

         ! Factor ca. 111km comes from R*pi/180, where R is earth
         ! radius and pi/180 is because we input lat/lon in degrees
         distance = 111.19492664455873 * sqrt(real(dx)**2 + real(dy)**2)
    
         ! Compute weight
         CALL PDAF_local_weight(wtype, rtype, cradius, sradius, distance, 1, 1, tmp, 1.0, weight, 0)
    
         ! Apply localization
         HP(j,i) = weight * HP(j,i)

      END DO
    END DO
    
    ! localize HPH^T
    DO j = 1, dim_obs
       DO i = 1, dim_obs
    
         ! Compute distance: obs - obs

         ! Units: Index numbering
         ! dx = abs(longxy_obs(j) - longxy_obs(i))
         ! dy = abs(latixy_obs(j) - latixy_obs(i))

         ! Units: lat/lon
         dx = abs(clmobs_lon(obs_pdaf2nc(j)) - clmobs_lon(obs_pdaf2nc(i)))
         dy = abs(clmobs_lat(obs_pdaf2nc(j)) - clmobs_lat(obs_pdaf2nc(i)))

         ! Check for longitude differences that may yield differences
         ! larger than 180deg depending on conventions. Example:
         ! crossing the prime meridian (lon=0deg), when convention is
         ! all-positive lons
         IF (dx > 180.0) THEN
           dx = 360.0 - dx
         END IF

         ! Intermediate latitude
         yhalf = (clmobs_lat(obs_pdaf2nc(j)) + clmobs_lat(obs_pdaf2nc(i))) / 2.0

         ! Latitude-dependent factor for longitude difference
         dx = dx * cos(yhalf * 3.14159265358979323846 / 180.0)

         ! Factor ca. 111km comes from R*pi/180, where R is earth
         ! radius and pi/180 is because we input lat/lon in degrees
         distance = 111.19492664455873 * sqrt(real(dx)**2 + real(dy)**2)
    
         ! Compute weight
         CALL PDAF_local_weight(wtype, rtype, cradius, sradius, distance, 1, 1, tmp, 1.0, weight, 0)
    
         ! Apply localization
         HPH(j,i) = weight * HPH(j,i)

       END DO
    END DO

    ! Deallocate lon/lat arrays

    ! For other filters these are deallocated in clean_obs_nc, invoked
    ! at the end of init_dim_obs_pdaf. For LEnKF, deallocation is
    ! moved to end of localize_covar_pdaf as these arrays are still
    ! used here.
    if(allocated(clmobs_lon))deallocate(clmobs_lon)
    if(allocated(clmobs_lat))deallocate(clmobs_lat)
    
  ENDIF ! model==tag_model_clm
#endif
!hcp end

END SUBROUTINE localize_covar_pdaf
