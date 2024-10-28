!$Id: localize_covar_pdaf.F90 1676 2016-12-10 14:55:45Z lnerger $
!BOP
!
! !ROUTINE: localize_covar_pdaf --- apply localization matrix in LEnKF
!
! !INTERFACE:
SUBROUTINE localize_covar_pdaf(dim_state, dim_obs, HP, HPH)

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
!    ONLY: cradius, sradius, locweight, obs_nc2pdaf
! hcp
! we need to store the coordinates of the state vector 
! and obs array in longxy, latixy, and longxy_obs, latixy_obs
! respectively
#if defined CLMSA
    ONLY: cradius, sradius, locweight, obs_nc2pdaf, &
          longxy, latixy, longxy_obs, latixy_obs
!hc  end
#else
    ONLY: cradius, sradius, locweight, obs_nc2pdaf
#endif
!fin hcp

  USE mod_read_obs,&
  ONLY: x_idx_obs_nc, y_idx_obs_nc, z_idx_obs_nc
#if defined CLMSA
   USE enkf_clm_mod, ONLY: init_clm_l_size, clmupdate_T
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

  USE, INTRINSIC :: iso_c_binding, ONLY: C_F_POINTER

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_state              ! State dimension
  INTEGER, INTENT(in) :: dim_obs                ! number of observations
  REAL, INTENT(inout) :: HP(dim_obs, dim_state) ! Matrix HP
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
  INTEGER :: dim_l, ncellxy,k
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
      print *, "TSMP-PDAF mype(w)=", mype_world, ": ERROR(2) `point_obs.eq.1` needed for using obs_nc2pdaf."
      call abort_parallel()
    end if

    ! localize HP
    DO j = 1, dim_obs
       DO i = 1, dim_state

         ! Index in coordinate array only spans `enkf_subvecsize`.
         !
         ! Necessary condition: the full state vector consists of
         ! sections of size `enkf_subvecsize`, where each section
         ! corresponds to a single coordinate array.
         icoord = modulo(i,enkf_subvecsize)

         dx = abs(x_idx_obs_nc(obs_nc2pdaf(j)) - int(xcoord_fortran(icoord))-1)
         dy = abs(y_idx_obs_nc(obs_nc2pdaf(j)) - int(ycoord_fortran(icoord))-1)
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
         dx = abs(x_idx_obs_nc(obs_nc2pdaf(j)) - x_idx_obs_nc(obs_nc2pdaf(i)))
         dy = abs(y_idx_obs_nc(obs_nc2pdaf(j)) - y_idx_obs_nc(obs_nc2pdaf(i)))
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
   !hcp
   !get the number of soil layer for evalute the arrangement 
   ! of the element in the state vector.
   call init_clm_l_size(dim_l)

   IF(model==tag_model_clm)THEN

    ! localize HP
    ! ncellxy=dim_state/dim_l
    if (.not. clmupdate_T.EQ.0) then
        dim_l = 2
    end if

    ncellxy=dim_state/dim_l
    DO j = 1, dim_obs
      DO k=1,dim_l
        DO i = 1, ncellxy
         ! Compute distance: obs - gridcell
!         dx = abs(longxy_obs(j) - longxy(i)-1)
!         dy = abs(latixy_obs(j) - latixy(i)-1)
         dx = abs(longxy_obs(j) - longxy(i))
         dy = abs(latixy_obs(j) - latixy(i))
         distance = sqrt(real(dx)**2 + real(dy)**2)
    
         ! Compute weight
         CALL PDAF_local_weight(wtype, rtype, cradius, sradius, distance, 1, 1, tmp, 1.0, weight, 0)
    
         ! Apply localization
         HP(j,i+(k-1)*ncellxy) = weight * HP(j,i+(k-1)*ncellxy)

        END DO
      END DO
    END DO
    
    ! localize HPH^T
    DO j = 1, dim_obs
       DO i = 1, dim_obs
    
         ! Compute distance: obs - obs
         dx = abs(longxy_obs(j) - longxy_obs(i))
         dy = abs(latixy_obs(j) - latixy_obs(i))
!         dy = abs(y_idx_obs_nc(obs_nc2pdaf(j)) - y_idx_obs_nc(obs_nc2pdaf(i)))
         distance = sqrt(real(dx)**2 + real(dy)**2)
    
         ! Compute weight
         CALL PDAF_local_weight(wtype, rtype, cradius, sradius, distance, 1, 1, tmp, 1.0, weight, 0)
    
         ! Apply localization
         HPH(j,i) = weight * HPH(j,i)

       END DO
    END DO
    
  ENDIF ! model==tag_model_parflow
#endif
!hcp end

END SUBROUTINE localize_covar_pdaf
