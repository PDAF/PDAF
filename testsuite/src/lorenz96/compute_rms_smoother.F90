!$Id: compute_rms_smoother.F90 1590 2015-06-12 11:59:35Z lnerger $
!BOP
!
! !ROUTINE: compute_rms_smoother --- Compute RMS errors for smoothed states
!
! !INTERFACE:
SUBROUTINE compute_rms_smoother(step, dim_lag, dim, dim_ens, state, &
     variance, rmse_s, trmse_s, mrmse_s_null, mtrmse_s_null, &
     mrmse_s_step, mtrmse_s_step)


! !DESCRIPTION:
! Helper routine for for pre/poststep routine.
! Used in the filters: ETKF/LETKF/ESTKF\LESTKF
! 
! The routine is called by prepoststep_seik. It computes
! the estimated and true RMS errors of the smoothed
! state estimated.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2012-06 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_memcount, &
       ONLY: memcount
!   USE mod_model, &
!        ONLY: dim_state, local_dims, dt, step_null
!   USE mod_parallel, &
!        ONLY: mype_filter, npes_filter, COMM_filter, MPI_DOUBLE_PRECISION, &
!        MPIerr, MPIstatus
  USE mod_assimilation, &
       ONLY: delt_obs, stepnull_means, fileid_state, state_true
  USE PDAF_interfaces_module, &    ! Required to get pointer to smoother state array
       ONLY: PDAF_get_smootherens

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
  INTEGER, INTENT(in) :: dim_lag     ! Size of lag for smoother
  INTEGER, INTENT(in) :: dim         ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
  REAL, INTENT(inout) :: state(dim)            ! Array for state vector
  REAL, INTENT(inout) :: variance(dim)         ! Array for variance
  REAL, INTENT(out)   :: rmse_s(dim_lag)   ! Estimated rms errors
  REAL, INTENT(out)   :: trmse_s(dim_lag)  ! True rms errors
  REAL, INTENT(inout) :: mrmse_s_null(dim_lag)   ! Time-mean estimated rms errors (from step 0)
  REAL, INTENT(inout) :: mtrmse_s_null(dim_lag)  ! Time-mean true rms errors (from step 0 )
  REAL, INTENT(inout) :: mrmse_s_step(dim_lag)   ! Time-mean estimated rms errors
  REAL, INTENT(inout) :: mtrmse_s_step(dim_lag)  ! Time-mean true rms errors


! !CALLING SEQUENCE:
! Called by: prepoststep_seik

! *** local variables ***
  INTEGER :: i, j, s, lag              ! counters
  REAL :: invdim_ens                   ! Inverse ensemble size
  REAL :: invdim_ensm1                 ! Inverse of ensemble size minus 1
  REAL, POINTER :: sens_pointer(:,:,:) ! Pointer to smoother ensemble
  INTEGER :: lastlag                   ! Number initialized lags in sens_pointer
  INTEGER :: status                    ! Status flag for PDAF_get_smootherens
  LOGICAL, SAVE :: firstcall = .true.  ! Flag, whether routine is called for first time
  INTEGER :: id_dim, id_step            ! File dimension IDs
  INTEGER :: id_state                   ! File variable ID
  INTEGER :: stat(50)                   ! Array for status flag
  INTEGER :: nsteps_file                ! Number of time steps in trajectory file
  INTEGER :: pos(2)                     ! Position index for writing
  INTEGER :: cnt(2)                     ! Count index for writing
  INTEGER :: state_step1and2(2)         ! First and second time step index in state file
  INTEGER, SAVE :: statefile_laststep   ! Last time step stored in state file
  INTEGER, SAVE :: delt_state_file      ! Interval between sucessively stored states
  INTEGER :: read_pos                   ! Which time step to read from the state file 
  ! Variables for mean errors from step 0
  REAL, ALLOCATABLE, SAVE :: sum_rmse_null(:)  ! Estimated RMS error accumulated over time
  REAL, ALLOCATABLE, SAVE :: sum_trmse_null(:) ! True RMS error accumulated over time
  INTEGER, ALLOCATABLE, SAVE :: nsum_null(:)   ! Length of sums over time
  ! Variables for sum from step stepnull_means
  REAL, ALLOCATABLE, SAVE :: sum_rmse_step(:)  ! Estimated RMS error accumulated over time
  REAL, ALLOCATABLE, SAVE :: sum_trmse_step(:) ! True RMS error accumulated over time
  INTEGER, ALLOCATABLE, SAVE :: nsum_step(:)   ! Length of sums over time


! **********************
! *** INITIALIZATION ***
! **********************

  ! Get access to smoother array
  IF (dim_lag > 0 .AND. step > 0) THEN
     CALL PDAF_get_smootherens(sens_pointer, lastlag, status)
  END IF

  ! Initialize numbers
  invdim_ens    = 1.0 / REAL(dim_ens)  
  invdim_ensm1  = 1.0 / REAL(dim_ens - 1)

  ! Initialize RMS errors to -1 to show that they are uninitialized
  rmse_s  = -1.0
  trmse_s = -1.0


  initcall: IF (firstcall) THEN

     ! *** Read information from file holding true state.  ***
     ! *** The file was already opened in compute_truerms. ***

     ! Get dimensions
     s = 1
     stat(s) = NF_INQ_DIMID(fileid_state, 'timesteps', id_dim)
     s = s + 1
     stat(s) = NF_INQ_DIMLEN(fileid_state, id_dim, nsteps_file)

     ! Read time step information
     s = s + 1
     stat(s) = NF_INQ_VARID(fileid_state, 'step', id_step)

     pos(1) = 1
     cnt(1) = 2
     s = s + 1
     stat(s) = NF_GET_VARA_INT(fileid_state, id_step, pos, cnt, state_step1and2)
  
     pos(1) = nsteps_file
     cnt(1) = 1
     s = s + 1
     stat(s) = NF_GET_VARA_INT(fileid_state, id_step, pos, cnt, statefile_laststep)
     s = s + 1

     DO i = 1,  s - 1
        IF (stat(i) /= NF_NOERR) &
            WRITE(*, *) 'NetCDF error reading trajectory file, no.', i
     END DO

     ! Initialize observation interval in file
     delt_state_file = state_step1and2(2) - state_step1and2(1)


     ! *** Prepare arrays for computing time-mean RMS errors ***

     ALLOCATE(sum_rmse_null(dim_lag))
     ALLOCATE(sum_trmse_null(dim_lag))
     ALLOCATE(nsum_null(dim_lag))
     ALLOCATE(sum_rmse_step(dim_lag))
     ALLOCATE(sum_trmse_step(dim_lag))
     ALLOCATE(nsum_step(dim_lag))

     sum_rmse_null = 0.0
     sum_trmse_null = 0.0
     nsum_null = 0
     sum_rmse_step = 0.0
     sum_trmse_step = 0.0
     nsum_step = 0

     ! *** Set firstcall flag ***
     firstcall = .false.

  END IF initcall


! **********************************************
! *** Compute RMS errors for smoothed states ***
! **********************************************

  ! Initialize screen output
  IF (lastlag > 0) THEN
     WRITE (*, '(8x, a)') 'Analize smoothed states'
     WRITE (*, '(8x, a10,7x,a12,14x,a11,13x,a10,i6)') &
          'RMS errors:','current step','mean from 0','mean from ',stepnull_means
     write (*, '(12x,a6,1x,a10,6x,a4,4x,a1,3x,a9,3x,a9,1x,a1,3x,a9,3x,a9)') &
          'lag ','estimate','true','|','mean est.','mean true','|','mean est.','mean true'
  END IF


  lagloopA: DO lag = 1, lastlag

     ! ************************************
     ! *** Compute estimated RMS errors ***
     ! ************************************

     ! Compute mean state 
     state = 0.0
     DO i = 1, dim_ens
        DO j = 1, dim
           state(j) = state(j) + sens_pointer(j, i, lag)
        END DO
     END DO
     state(:) = invdim_ens * state(:)

     ! *** Compute local sampled variances ***
     variance(:) = 0.0
     DO i = 1, dim_ens
        DO j = 1, dim
           variance(j) = variance(j) &
             + (sens_pointer(j, i, lag) - state(j)) &
             * (sens_pointer(j, i, lag) - state(j))
        END DO
     END DO
     variance(:) = invdim_ensm1 * variance(:)

     ! *** Compute RMS error according to sampled variances ***
     rmse_s(lag)  = 0.0
     DO i = 1, dim
        rmse_s(lag) = rmse_s(lag) + variance(i)
     ENDDO
     rmse_s(lag) = SQRT(rmse_s(lag) / dim)


     ! *******************************
     ! *** Compute true RMS errors ***
     ! *******************************

     comprms: IF (abs(step) < statefile_laststep) THEN

        ! Initialize file position corresponding to current time step
        read_pos =  (step - lag * delt_obs) / delt_state_file + 1

        ! Read true state
        s = 1
        stat(s) = NF_INQ_VARID(fileid_state, 'state', id_state)

        pos(2) = read_pos
        cnt(2) = 1
        pos(1) = 1
        cnt(1) = dim
        s = s + 1
        stat(s) = NF_GET_VARA_DOUBLE(fileid_state, id_state, pos, cnt, state_true)

        DO i = 1,  s
           IF (stat(i) /= NF_NOERR) &
                WRITE(*, *) 'NetCDF error in reading true state from file, no.', i
        END DO

        ! Compute RMS
        trmse_s(lag) = 0.0
        DO j = 1, dim
           trmse_s(lag) = trmse_s(lag) + (state_true(j) - state(j))**2
        END DO
        trmse_s(lag) = SQRT(trmse_s(lag) / dim)

     END IF comprms

  END DO lagloopA


! ********************************************************
! *** Compute time-mean RMS errors for smoothed states ***
! ********************************************************

  lagloopB: DO lag = 1, lastlag

     ! Sums from step 0
     sum_rmse_null(lag) = sum_rmse_null(lag) + rmse_s(lag)
     sum_trmse_null(lag) = sum_trmse_null(lag) + trmse_s(lag)
     nsum_null(lag) = nsum_null(lag) + 1
     mrmse_s_null(lag) = sum_rmse_null(lag) / REAL(nsum_null(lag))
     mtrmse_s_null(lag) = sum_trmse_null(lag) / REAL(nsum_null(lag))

     ! Sums from step stepnull_means
     IF (ABS(step) >= stepnull_means) THEN
        sum_rmse_step(lag) = sum_rmse_step(lag) + rmse_s(lag)
        sum_trmse_step(lag) = sum_trmse_step(lag) + trmse_s(lag)
        nsum_step(lag) = nsum_step(lag) + 1
        mrmse_s_step(lag) = sum_rmse_step(lag) / REAL(nsum_step(lag))
        mtrmse_s_step(lag) = sum_trmse_step(lag) / REAL(nsum_step(lag))
     END IF

     ! Screen output
     write (*,'(12x,i6,2es12.4,1x,a1,2es12.4,1x,a1,2es12.4)') &
          lag, rmse_s(lag), trmse_s(lag), '|',mrmse_s_null(lag), &
          mtrmse_s_null(lag), '|',mrmse_s_step(lag), mtrmse_s_step(lag)
  END DO lagloopB


! ********************
! *** finishing up ***
! ********************

END SUBROUTINE compute_rms_smoother
