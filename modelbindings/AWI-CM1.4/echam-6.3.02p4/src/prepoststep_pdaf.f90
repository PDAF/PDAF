!$Id$
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
       ONLY: mype_filter_echam, MPIerr, MPI_SUM, comm_filter_echam, &
       MPI_DOUBLE_PRECISION, writepe
  USE mod_assim_pdaf, &
       ONLY: n_fields, dim_fields_p, dim_fields, off_fields_p
  USE mod_assim_atm_pdaf, ONLY: dp
  USE mo_decomposition, ONLY: dc=>local_decomposition
  USE output_pdaf, &
       ONLY: write_da, write_netcdf_pdaf, write_pos_da


  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step        !< Current time step
     ! (When the routine is called before the analysis -step is provided.)
  INTEGER, INTENT(in) :: dim_p       !< Process-local state dimension
  INTEGER, INTENT(in) :: dim_ens     !< Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   !< Process-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   !< Process-local dimension of observation vector
  REAL(dp), INTENT(inout) :: state_p(dim_p) !< Process-local forecast/analysis state
  ! The array 'state_p' is not generally not initialized in the case of SEIK.
  ! It can be used freely here.
  REAL(dp), INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) !< Inverse of matrix U
  REAL(dp), INTENT(inout) :: ens_p(dim_p, dim_ens)      !< Process-local state ensemble
  INTEGER, INTENT(in) :: flag        !< PDAF status flag


! *** Local variables ***
  INTEGER :: i, j, member, field    ! Counters
  REAL :: invdim_ens                ! Inverse ensemble size
  REAL :: invdim_ensm1              ! Inverse of ensemble size minus 1 
  CHARACTER(len=1) :: typestr       ! Character indicating call type
  REAL :: rmse_p(7)                 ! PE-local estimated rms errors
  REAL :: rmse(7)                   ! Global estimated rms errors
  REAL, ALLOCATABLE :: var_p(:)     ! Estimated local model state variances


! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype_filter_echam==0) THEN
     IF (step==0) THEN
        WRITE (*,'(a, i7,3x,a)') 'ECHAM-PDAF', step,'Analyze initial state ensemble'
        WRITE (typestr,'(a1)') 'i'
     ELSE IF (step>0) THEN
        WRITE (*,'(a, 8x,a)') 'ECHAM-PDAF', 'Analyze assimilated state ensemble'
        WRITE (typestr,'(a1)') 'a'
     ELSE IF (step<0) THEN
        WRITE (*,'(a, 8x,a)') 'ECHAM-PDAF', 'Analyze forecast state ensemble'
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

! ****************************
! *** Perform pre/poststep ***
! ****************************

  ! *** Compute mean state
  IF (mype_filter_echam==0) WRITE (*,'(a, 8x,a)') &
       'ECHAM-PDAF', '--- compute ensemble mean'

  ! local 
  state_p = 0.0
  DO member = 1, dim_ens
     DO i = 1, dim_p
        state_p(i) = state_p(i) + ens_p(i,member)
     END DO
  END DO
  state_p(:) = invdim_ens * state_p(:)


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

  DO field = 1, n_fields

     DO i = 1, dim_fields_p(field)
        rmse_p(field) = rmse_p(field) + var_p(i + off_fields_p(field))
     ENDDO

     rmse_p(field) = rmse_p(field) / REAL(dim_fields(field))

  END DO

  ! Global sum of RMS errors
  CALL MPI_Allreduce (rmse_p, rmse, n_fields, MPI_DOUBLE_PRECISION, MPI_SUM, &
       comm_filter_echam, MPIerr)

  rmse = SQRT(rmse)

  ! Display RMS errors
  IF (mype_filter_echam==0) THEN
     WRITE (*,'(a, 10x,a)') &
          'ECHAM-PDAF', 'RMS error according to sampled covariance'
     WRITE (*,'(a,7x,a10,a13,a13,a14,a13,a12,a13,/a, 10x, 77a)') &
          'ECHAM-PDAF', 'T','lsp','vo','d','q','u','v',&
          'ECHAM-PDAF', ('-',i=1,66)
     WRITE (*,'(a,10x,es11.4,6es13.4,1x,a5,a1,/a, 10x,77a)') &
          'ECHAM-PDAF', rmse(1), rmse(2), rmse(3), rmse(4), rmse(5), rmse(6),rmse(7),'RMSe-', typestr,&
          'ECHAM-PDAF', ('-',i=1,77)
  END IF



! **************************
! *** Write output files ***
! **************************

  ! *** Write output to NetCDF files

  IF (step == 0) THEN
     ! *** write initial state fields ***
     CALL write_netcdf_pdaf('i', write_pos_da, step, dim_p, state_p, n_fields, rmse, writepe)
  ELSE IF (step > 0) THEN
     ! *** write assimilated state fields ***
     CALL write_netcdf_pdaf('a', write_pos_da, step, dim_p, state_p, n_fields, rmse, writepe)

     ! Increment write position
     write_pos_da = write_pos_da + 1

  ELSE IF (step < 0) THEN
     ! *** write forecasted state fields ***
     CALL write_netcdf_pdaf('f', write_pos_da, step, dim_p, state_p, n_fields, rmse, writepe)
  END IF


! *********************************************************
! *** Deallocate observation-related arrays             ***
! *** which were allocate in the INIT_DIMOBS_F routines ***
! *********************************************************

  IF (step > 0) THEN
     CALL deallocate_obs_pdafomi()
  END IF


! ********************
! *** finishing up ***
! ********************

  DEALLOCATE(var_p)

END SUBROUTINE prepoststep_pdaf
