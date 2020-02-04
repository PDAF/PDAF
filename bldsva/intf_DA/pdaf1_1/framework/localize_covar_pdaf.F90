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
!    ONLY: local_range, srange, locweight, obs_nc2pdaf
! hcp
! we need to store the coordinates of the state vector 
! and obs array in longxy, latixy, and longxy_obs, latixy_obs
! respectively
#if defined CLMSA
    ONLY: local_range, srange, locweight, obs_nc2pdaf, &
          longxy, latixy, longxy_obs, latixy_obs
!hc  end
#else
    ONLY: local_range, srange, locweight, obs_nc2pdaf
#endif
!fin hcp

  USE mod_read_obs,&
  ONLY: x_idx_obs_nc, y_idx_obs_nc, z_idx_obs_nc
#if defined CLMSA
   USE enkf_clm_mod, ONLY: init_clm_l_size
   USE mod_parallel_pdaf, ONLY: filterpe
#endif
!fin hcp
  USE mod_parallel_model, &
    ONLY: model

  USE mod_tsmp,&
#if defined CLMSA
    ONLY:  tag_model_parflow, tag_model_clm, &
           enkf_subvecsize
#else
    ONLY: tag_model_parflow, tag_model_clm, &
          enkf_subvecsize, &
          nx_glob, ny_glob, nz_glob, &
          xcoord, ycoord, zcoord, &
          xcoord_fortran, ycoord_fortran, zcoord_fortran
#endif

  USE, INTRINSIC :: iso_c_binding

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

! **********************
! *** INITIALIZATION ***
! **********************

  ! Screen output
  WRITE (*,'(8x, a)') &
       '--- Apply covariance localization'
  WRITE (*, '(12x, a, 1x, f12.2)') &
       '--- Local influence radius', local_range

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
!hcp enkf_subvecsize is the number of grid cells,
!which is the size of state vector only if soil moisture
!or pressure is the entire state vector, which is not true
!when other parameters are also included in the state vector
!in which the size of x/ycoord is NOT enkf_subvecsize.
    call C_F_POINTER(xcoord,xcoord_fortran,[dim_state]) ![enkf_subvecsize])
    call C_F_POINTER(ycoord,ycoord_fortran,[dim_state]) ![enkf_subvecsize])
    call C_F_POINTER(zcoord,zcoord_fortran,[dim_state]) ![enkf_subvecsize])

    ! localize HP
    DO j = 1, dim_obs
       DO i = 1, dim_state
    
         dx = abs(x_idx_obs_nc(obs_nc2pdaf(j)) - int(xcoord_fortran(i))-1)
         dy = abs(y_idx_obs_nc(obs_nc2pdaf(j)) - int(ycoord_fortran(i))-1)
         distance = sqrt(real(dx)**2 + real(dy)**2)
    
         ! Compute weight
         CALL PDAF_local_weight(wtype, rtype, local_range, srange, distance, 1, 1, tmp, 1.0, weight, 0)
    
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
         CALL PDAF_local_weight(wtype, rtype, local_range, srange, distance, 1, 1, tmp, 1.0, weight, 0)
    
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
!    call C_F_POINTER(xcoord,xcoord_fortran,[enkf_subvecsize])
!    call C_F_POINTER(ycoord,ycoord_fortran,[enkf_subvecsize])
!    call C_F_POINTER(zcoord,zcoord_fortran,[enkf_subvecsize])

    ! localize HP
    ncellxy=dim_state/dim_l
    DO j = 1, dim_obs
      DO k=1,dim_l
        DO i = 1, ncellxy
         dx = abs(longxy_obs(j) - longxy(i)-1)
         dy = abs(latixy_obs(j) - latixy(i)-1)
         distance = sqrt(real(dx)**2 + real(dy)**2)
    
         ! Compute weight
         CALL PDAF_local_weight(wtype, rtype, local_range, srange, distance, 1, 1, tmp, 1.0, weight, 0)
    
         ! Apply localization
         HP(j,i+(k-1)*ncellxy) = weight * HP(j,i+(k-1)*ncellxy)

        END DO
      END DO
    END DO
    
    ! localize HPH^T
    DO j = 1, dim_obs
       DO i = 1, dim_obs
    
         ! Compute distance
         dx = abs(longxy_obs(j) - longxy_obs(i))
         dy = abs(latixy_obs(j) - latixy_obs(i))
!         dy = abs(y_idx_obs_nc(obs_nc2pdaf(j)) - y_idx_obs_nc(obs_nc2pdaf(i)))
         distance = sqrt(real(dx)**2 + real(dy)**2)
    
         ! Compute weight
         CALL PDAF_local_weight(wtype, rtype, local_range, srange, distance, 1, 1, tmp, 1.0, weight, 0)
    
         ! Apply localization
         HPH(j,i) = weight * HPH(j,i)

       END DO
    END DO
    
  ENDIF ! model==tag_model_parflow
#endif
!hcp end

END SUBROUTINE localize_covar_pdaf
