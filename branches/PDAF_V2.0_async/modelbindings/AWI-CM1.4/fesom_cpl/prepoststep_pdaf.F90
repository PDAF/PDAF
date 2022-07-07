!$Id: prepoststep_pdaf.F90 2271 2020-04-08 13:04:09Z lnerger $
!>  Routine for pre- and postsstep operations for PDAF
!!
!! User-supplied call-back routine for PDAF.
!!
!! The routine is called for all filters
!! before the analysis and after the ensemble transformation.
!! Also it is called once at the initial time
!! before any forecasts are computed.
!! The routine provides full access to the state 
!! estimate and the state ensemble to the user.
!! Thus, user-controlled pre- and poststep 
!! operations can be performed here. For example 
!! the forecast and the analysis states and ensemble
!! covariance matrix can be analized, e.g. by 
!! computing the estimated variances. In addition, 
!! the estimates can be written to disk. If a user 
!! considers to perform adjustments to the 
!! estimates (e.g. for balances), this routine is 
!! the right place for it.
!!
!! The routine is called by all filter processes.
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
SUBROUTINE prepoststep_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

  USE mod_parallel_pdaf, &
       ONLY: mype_filter_fesom, npes_filter, COMM_filter_fesom, writepe, &
       MPI_DOUBLE_PRECISION, MPI_SUM, MPIerr
  USE mod_assim_pdaf, & ! Variables for assimilation
       ONLY: step_null, filtertype, dim_lag, eff_dim_obs, loctype, &
       off_fields_p, n_fields, dim_fields_p, dim_fields
  USE obs_sst_cmems_pdafomi, &
       ONLY: assim_o_sst, sst_exclude_ice, sst_exclude_diff, mean_ice_p, mean_sst_p
  USE g_parfe, &
       ONLY: mydim_nod2d
  USE o_mesh, &
       ONLY: nod2D
  USE output_pdaf, &
       ONLY: write_da, write_netcdf_pdaf, write_netcdf_pdaf_ens, &
       write_pos_da, write_ens, write_pos_da_ens

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step        !< Current time step
     ! (When the routine is called before the analysis -step is provided.)
  INTEGER, INTENT(in) :: dim_p       !< Process-local state dimension
  INTEGER, INTENT(in) :: dim_ens     !< Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   !< Process-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   !< Process-local dimension of observation vector
  REAL, INTENT(inout) :: state_p(dim_p) !< Process-local forecast/analysis state
  ! The array 'state_p' is not generally not initialized in the case of SEIK.
  ! It can be used freely here.
  REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) !< Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      !< Process-local state ensemble
  INTEGER, INTENT(in) :: flag        !< PDAF status flag


! *** Local variables ***
  INTEGER :: i, j, member, field    ! Counters
  REAL :: invdim_ens                ! Inverse ensemble size
  REAL :: invdim_ensm1              ! Inverse of ensemble size minus 1 
  REAL :: rmse_p(12)                ! Process-local estimated rms errors
  REAL :: rmse(12)                  ! Global estimated rms errors
  REAL, ALLOCATABLE :: var_p(:)     ! Estimated local model state variances
  CHARACTER(len=1) :: typestr       ! Character indicating call type

  
! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype_filter_fesom==0) THEN
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

  ! Allocate array for variances
  ALLOCATE(var_p(dim_p))

  ! Initialize numbers
  rmse_p = 0.0
  rmse   = 0.0
  invdim_ens = 1.0 / REAL(dim_ens)  
  invdim_ensm1 = 1.0 / REAL(dim_ens-1)  

  ! Allocate eff_dim_obs
  IF (.NOT.ALLOCATED(eff_dim_obs)) ALLOCATE(eff_dim_obs(mydim_nod2d))
  eff_dim_obs = 0.0


! ****************************
! *** Perform pre/poststep ***
! ****************************

  ! *** Compute mean state
  IF (mype_filter_fesom==0) WRITE (*,'(a, 8x,a)') &
       'FESOM-PDAF', '--- compute ensemble mean'

  ! local 
  state_p = 0.0
  DO member = 1, dim_ens
     DO i = 1, dim_p
        state_p(i) = state_p(i) + ens_p(i,member)
     END DO
  END DO
  state_p(:) = invdim_ens * state_p(:)


! *********************************************************************
! *** Store ensemble mean values for observation exclusion criteria ***
! *********************************************************************

  IF (assim_o_sst .AND. step<0) THEN

     ! ice concentration
     IF (sst_exclude_ice) THEN 
        IF (ALLOCATED(mean_ice_p)) DEALLOCATE(mean_ice_p)
        ALLOCATE (mean_ice_p(myDim_nod2D))
        mean_ice_p = state_p(1+off_fields_p(7) : myDim_nod2D+off_fields_p(7))
     END IF

     ! SST
     IF (sst_exclude_ice .OR. sst_exclude_diff > 0.0) THEN
        IF (ALLOCATED(mean_sst_p)) DEALLOCATE(mean_sst_p)
        ALLOCATE (mean_sst_p(myDim_nod2D))
        mean_sst_p = state_p(1+off_fields_p(5) : myDim_nod2D+off_fields_p(5))
     END IF

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


  ! *** Compute RMS errors ***

  ! All fields in state vector
  DO field = 1, n_fields

     DO i = 1, dim_fields_p(field)
        rmse_p(field) = rmse_p(field) + var_p(i + off_fields_p(field))
     ENDDO

     rmse_p(field) = rmse_p(field) / REAL(dim_fields(field))

  END DO

  ! RMSS error for SST
  field = 12
  DO i = 1, myDim_nod2D
     rmse_p(field) = rmse_p(field) + var_p(i + off_fields_p(5))
  END DO
  rmse_p(field) = rmse_p(field) / REAL(nod2D)

  ! Global sum of RMS errors
  CALL MPI_Allreduce (rmse_p, rmse, 12, MPI_DOUBLE_PRECISION, MPI_SUM, &
       COMM_filter_fesom, MPIerr)

  rmse = SQRT(rmse)

  ! Display RMS errors
  IF (mype_filter_fesom==0) THEN
     WRITE (*,'(a, 10x,a)') &
          'FESOM-PDAF', 'RMS error according to sampled covariance'
     WRITE (*,'(a,7x,a10,a13,a13,a14,a13,a12,/a, 10x, 77a)') &
          'FESOM-PDAF', 'ssh','u','v','temp','salt','SST', &
          'FESOM-PDAF', ('-',i=1,77)
     WRITE (*,'(a,10x,es11.4,5es13.4,1x,a5,a1,/a, 10x, 77a)') &
          'FESOM-PDAF', rmse(1), rmse(2), rmse(3), rmse(5), rmse(6), rmse(12), 'RMSe-', typestr,&
          'FESOM-PDAF', ('-',i=1,77)
  END IF


! ***************************************************************
! *** Compute statistics for effective observation dimensions ***
! ***************************************************************

  IF (step>0) CALL adaptive_lradius_stats_pdaf()


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


! *********************************************************
! *** Deallocate observation-related arrays             ***
! *** which were allocate in the INIT_DIMOBS_F routines ***
! *********************************************************

  IF (step > 0) THEN
    CALL deallocate_obs_pdafomi()

    IF (ALLOCATED(mean_ice_p)) DEALLOCATE(mean_ice_p)
    IF (ALLOCATED(mean_sst_p)) DEALLOCATE(mean_sst_p)
    IF (ALLOCATED(eff_dim_obs)) DEALLOCATE(eff_dim_obs)
  END IF


! ********************
! *** finishing up ***
! ********************

  DEALLOCATE(var_p)

END SUBROUTINE prepoststep_pdaf
