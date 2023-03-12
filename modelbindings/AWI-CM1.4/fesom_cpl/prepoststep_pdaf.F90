!$Id: prepoststep_pdaf.F90 2136 2019-11-22 18:56:35Z lnerger $
!BOP
!
! !ROUTINE: prepoststep_pdaf - Routine controlling ensemble integration for PDAF
!
! !INTERFACE:
SUBROUTINE prepoststep_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! This variant is used with the simplified interface of
! PDAF. In this case, the name of the routine is defined
! within PDAF. This routine just calls the prepoststep
! routine corresponding to the selected filter algorithm.
!
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, npes_filter, COMM_filter, writepe
  USE mod_assim_pdaf, & ! Variables for assimilation
       ONLY: step_null, filtertype, ocoord_n2d, offset, &
       obs_sst, obs_sst_error, ivariance_obs, loc_radius, &
       id_nod2D_ice, mean_ice_p, mean_sst_p
  USE g_parfe, &
       ONLY: MPI_DOUBLE_PRECISION, MPI_SUM, MPIerr, MPI_STATUS_SIZE, &
       MPI_INTEGER, MPI_MAX, MPI_MIN, mydim_nod2d, mydim_nod3d, &
       MPI_COMM_FESOM
  USE g_config, &
       ONLY: r_restart
  USE o_mesh, &
       ONLY: nod2D
  USE output_pdaf, &
       ONLY: write_da, write_netcdf_pdaf, write_netcdf_pdaf_ens, &
       write_pos_da, write_ens, write_pos_da_ens

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
     ! (When the routine is called before the analysis -step is provided.)
  INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   ! PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   ! PE-local dimension of observation vector
  REAL, INTENT(inout) :: state_p(dim_p) ! PE-local forecast/analysis state
  ! The array 'state_p' is not generally not initialized in the case of SEIK.
  ! It can be used freely here.
  REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) ! Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      ! PE-local state ensemble
  INTEGER, INTENT(in) :: flag        ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_get_state      (as U_prepoststep)
! Called by: PDAF_X_update       (as U_prepoststep)
!EOP

! *** Local variables ***
  INTEGER :: i, j, member, field    ! Counters
  INTEGER :: offset_field           ! Offset of a field in the state vector
  INTEGER :: dim_field              ! Dimension of a field
  REAL :: invdim_ens                ! Inverse ensemble size
  REAL :: invdim_ensm1              ! Inverse of ensemble size minus 1 
  REAL :: rmse_p(13)                ! PE-local estimated rms errors
  REAL :: rmse(13)                  ! Global estimated rms errors
  REAL, ALLOCATABLE :: var_p(:)     ! Estimated local model state variances
  CHARACTER(len=1) :: typestr       ! Character indicating call type
  REAL :: min_eff_dim_obs, max_eff_dim_obs       ! Stats on effective observation dimensions
  REAL :: min_eff_dim_obs_g, max_eff_dim_obs_g   ! Stats on effective observation dimensions
  REAL :: sum_eff_dim_obs, avg_eff_dim_obs_g     ! Stats on effective observation dimensions
  INTEGER :: dim_nod2D_ice_p        ! PE-local ice node dimension
  
! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype_filter==0) THEN
     IF (step-step_null==0) THEN
        WRITE (*,'(a, i7,3x,a)') 'FESOM-PDAF', step,'Analyze initial state ensemble'
        WRITE (typestr,'(a1)') 'i'
     ELSE IF (step>0) THEN
        WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', 'Analyze assimilated state ensemble'
        WRITE (typestr,'(a1)') 'a'
     ELSE IF (step<0) THEN
        WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', 'Analyze forecast state ensemble'
        WRITE (typestr,'(a1)') 'f'
     END IF
  END IF

  ! After analysis deallocate observation-related fields
  ! These are allocated in COMPUTE_DIM_OBS_FULL
  IF (ALLOCATED(obs_sst)) DEALLOCATE(obs_sst)
  IF (ALLOCATED(ocoord_n2d)) DEALLOCATE(ocoord_n2d)
  IF (ALLOCATED(obs_sst_error)) DEALLOCATE(obs_sst_error)
  IF (ALLOCATED(ivariance_obs)) DEALLOCATE(ivariance_obs)


  ! Allocate array for variances
  ALLOCATE(var_p(dim_p))

  ! Initialize numbers
  rmse_p = 0.0
  rmse   = 0.0
  invdim_ens = 1.0 / REAL(dim_ens)  
  invdim_ensm1 = 1.0 / REAL(dim_ens-1)  


! ****************************
! *** Perform pre/poststep ***
! ****************************

  ! *** Compute mean state
  IF (mype_filter==0) WRITE (*,'(a, 8x,a)') &
       'FESOM-PDAF', '--- compute ensemble mean'

  ! local 
  state_p = 0.0
  DO member = 1, dim_ens
     DO i = 1, dim_p
        state_p(i) = state_p(i) + ens_p(i,member)
     END DO
  END DO
  state_p(:) = invdim_ens * state_p(:)


! **********************************************************************
! *** store ice concentration and mean SST for observation selection ***
! **********************************************************************

  IF (step < 0) THEN
     if (allocated(mean_ice_p)) deallocate(mean_ice_p)
     ALLOCATE (mean_ice_p(myDim_nod2D))
     mean_ice_p = state_p(1+offset (7):offset(7)+myDim_nod2D)
  END IF

  ! store mean_sst
  IF (step < 0) THEN
     if (allocated(mean_sst_p)) deallocate(mean_sst_p)
     ALLOCATE (mean_sst_p(myDim_nod2D))
     mean_sst_p = state_p(1+offset (5):offset(5)+myDim_nod2D)
  END IF


! ********************************************************
! *** Compute estimate RMS errors for different fields ***
! ********************************************************

  ! *** Compute local sampled variances of state vector ***
  var_p(:) = 0.0
  DO member = 1, dim_ens
     DO j = 1, dim_p
        var_p(j) = var_p(j) + &
             (ens_p(j, member) - state_p(j))* &
             (ens_p(j, member) - state_p(j))
     END DO
  END DO
  var_p(:) = invdim_ensm1 * var_p(:)


  offset_field = 0
  DO field = 1, 12
     ! Specify dimension of field
     ! SSH 
     IF (field == 1) THEN
        dim_field = myDim_nod2D
     ! u,v,w,T,S
     ELSEIF (field >= 2 .AND. field < 7) THEN
        dim_field = myDim_nod3D
     ! ice
     ELSEIF (field >= 7 .AND. field < 12) THEN
        dim_field = myDim_nod2D
     ELSE
     ! field 12 is the total RMS of ocean fields
        !dim_field = dim_p
        dim_field = myDim_nod2D + 5 * myDim_nod3D
        offset_field = 0
     END IF
     
     DO i = 1, dim_field
        rmse_p(field) = rmse_p(field) + var_p(i + offset_field)
     ENDDO
     rmse_p(field) = rmse_p(field) / real(dim_field)


     ! Set offset for next field
     offset_field = offset_field + dim_field
     
  END DO

  ! Calculate the SST rmse
  field = 13
  dim_field = myDim_nod2D
  offset_field = myDim_nod2D + 3 * myDim_nod3D
     DO i = 1, dim_field
        rmse_p(field) = rmse_p(field) + var_p(i + offset_field)
     END DO
  rmse_p(field) = rmse_p(field) / real(dim_field)

  ! Global sum of RMS errors
  CALL MPI_Allreduce (rmse_p, rmse, 13, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)

  rmse = SQRT(rmse)

  ! Display RMS errors
  IF (mype_filter==0) THEN
     WRITE (*,'(a, 10x,a)') &
          'FESOM-PDAF', 'RMS error according to sampled covariance'
     WRITE (*,'(a,7x,a9,1x,a13,a14,a16,a14,/a, 10x,66a)') &
          'FESOM-PDAF', 'ssh','u','v','temp','salt',&
          'FESOM-PDAF', ('-',i=1,66)
     WRITE (*,'(a,10x,es11.4,5es14.4,1x,a5,a1,/a, 10x,66a)') &
          'FESOM-PDAF', rmse(1), rmse(2), rmse(3), rmse(5), rmse(6), rmse(13), 'RMSe-', typestr,&
          'FESOM-PDAF', ('-',i=1,66)
  END IF


! **************************
! *** Write output files ***
! **************************

  write_pos_da_ens = write_pos_da

  output: IF (write_da) THEN
     ! *** Write output to NetCDF files

     IF ((step - step_null)==0) THEN
        ! *** write initial state fields ***
        CALL write_netcdf_pdaf('i', write_pos_da, step, dim_p, state_p, 12, rmse, writepe)
     ELSE IF ((step - step_null) > 0) THEN
        ! *** write assimilated state fields ***
        CALL write_netcdf_pdaf('a', write_pos_da, step, dim_p, state_p, 12, rmse, writepe)

        ! Increment write position
        write_pos_da = write_pos_da + 1

     ELSE IF ((step - step_null) < 0) THEN
        ! *** write forecasted state fields ***
        CALL write_netcdf_pdaf('f', write_pos_da, step, dim_p, state_p, 12, rmse, writepe)
     END IF

    ! *** Write output for individual ensemble member

    IF (write_ens) THEN
     IF ((step - step_null)==0) THEN
        ! *** write initial state fields ***
        CALL write_netcdf_pdaf_ens('i', write_pos_da_ens, step, dim_p, ens_p, 12, rmse, writepe, dim_ens)
     ELSE IF ((step - step_null) > 0) THEN
        ! *** write assimilated state fields ***
        CALL write_netcdf_pdaf_ens('a', write_pos_da_ens, step, dim_p, ens_p, 12, rmse, writepe, dim_ens)

     ELSE IF ((step - step_null) < 0) THEN
        ! *** write forecasted state fields ***
        CALL write_netcdf_pdaf_ens('f', write_pos_da_ens, step, dim_p, ens_p, 12, rmse, writepe, dim_ens)
     END IF
        
    END IF
 ENDIF output


! ********************
! *** finishing up ***
! ********************

  DEALLOCATE(var_p)

END SUBROUTINE prepoststep_pdaf
