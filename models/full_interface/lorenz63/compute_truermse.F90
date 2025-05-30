!$Id$
!BOP
!
! !ROUTINE: compute_truermse --- Compute true rms errors for assimilation
!
! !INTERFACE:
SUBROUTINE compute_truermse(calltype, step, time, dim, state_est, &
       trmse, dim_ens, ens, hist_true, hist_mean, &
       skewness, kurtosis, crps_stats)

! !DESCRIPTION:
! Helper routine for the pre/poststep routines. The routine
! reads the true state from the trajectory file and computes
! the true rms error for the currrent state estimate. In 
! addition, rank histograms about the ensemble mean and
! the true state as well as higher-order statistical moments
! (skewness, kurtosis) are computed.
!
! !REVISION HISTORY:
! 2010-01 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE netcdf
  USE PDAF, &
       ONLY: PDAF_diag_histogram, PDAF_diag_ensstats, PDAF_diag_crps
  USE mod_memcount, &
       ONLY: memcount
  USE mod_assimilation, &
       ONLY: fileid_state, state_true, stepnull_means, &
       observation_g, use_obs_mask, obsindx
  USE output_netcdf, &
       ONLY: file_state

  IMPLICIT NONE

! !ARGUMENTS:
  CHARACTER(len=3), INTENT(in) :: calltype ! Whether routine is called at 
        !(ini) initial time, (for) after forecast, or (ana) after analysis
  INTEGER, INTENT(in) :: step           ! Current time step
        ! (When the routine is called before the analysis -step is provided.)
  REAL, INTENT(in)    :: time           ! Currrent model time
  INTEGER, INTENT(in) :: dim            ! PE-local state dimension
  REAL, INTENT(in)    :: state_est(dim) ! Forecast/analysis state
  REAL, INTENT(out)   :: trmse          ! True rms error
  INTEGER, INTENT(in) :: dim_ens        ! Ensemble size
  REAL, INTENT(in)    :: ens(dim, dim_ens)     ! State ensemble
  ! Arrays for rank histograms
  INTEGER, INTENT(inout) :: hist_true(dim_ens+1, 2) ! Histogram about true state
  INTEGER, INTENT(inout) :: hist_mean(dim_ens+1, 2) ! Histogram about ensemble mean
  REAL, INTENT(out) :: skewness         ! Skewnes of ensemble
  REAL, INTENT(out) :: kurtosis         ! Kurtosis of ensemble
  REAL ,INTENT(out) :: crps_stats(4)    ! CRPS, reli, resol, uncert

! !CALLING SEQUENCE:
! Called by: prepoststep routine
!EOP

! *** local variables ***
  INTEGER :: i, j, s                    ! Counters
  INTEGER :: id_dim, id_step            ! File dimension IDs
  INTEGER :: id_state                   ! File variable ID
  INTEGER :: stat(50)                   ! Array for status flag
  INTEGER :: nsteps_file                ! Number of time steps in trajectory file
  INTEGER :: pos(2)                     ! Position index for writing
  INTEGER :: cnt(2)                     ! Count index for writing
  INTEGER :: state_step1and2(2)         ! First and second time step index in state file
  INTEGER, SAVE :: statefile_laststep(1)   ! Last time step stored in state file
  INTEGER, SAVE :: delt_state_file      ! Interval between sucessively stored states
  INTEGER :: read_pos                   ! Which time step to read from the state file 
  INTEGER :: status                     ! Status output for PDAF_computehistogram
  REAL :: delta_hist                    ! Delta value for histogram
  REAL :: crps, reli, resol, uncert     ! CRPS = reli + resol


! **********************
! *** INITIALIZATION ***
! **********************

  initcall: IF (calltype == 'ini') THEN
     
     ! Open Netcdf file holding true trajectory
     s = 1
     stat(s) = NF90_OPEN(TRIM(file_state), NF90_NOWRITE, fileid_state)

     ! Get dimensions
     s = s + 1
     stat(s) = NF90_INQ_DIMID(fileid_state, 'timesteps', id_dim)
     s = s + 1
     stat(s) = NF90_Inquire_dimension(fileid_state, id_dim, len=nsteps_file)

     ! Read time step information
     s = s + 1
     stat(s) = NF90_INQ_VARID(fileid_state, 'step', id_step)

     pos(1) = 1
     cnt(1) = 2
     s = s + 1
     stat(s) = NF90_GET_VAR(fileid_state, id_step, state_step1and2, start=pos(1:1), count=cnt(1:1))
  
     pos(1) = nsteps_file
     cnt(1) = 1
     s = s + 1
     stat(s) = NF90_GET_VAR(fileid_state, id_step, statefile_laststep, start=pos(1:1), count=cnt(1:1))
     s = s + 1

     DO i = 1,  s - 1
        IF (stat(i) /= NF90_NOERR) &
            WRITE(*, *) 'NetCDF error reading trajectory file, no.', i
     END DO

     ! Initialize observation interval in file
     delt_state_file = state_step1and2(2) - state_step1and2(1)

     ! Allocate array for true state
     ALLOCATE(state_true(dim))
     ! count allocated memory
     CALL memcount(3, 'r', dim)

  END IF initcall

  ! Initialize file position corresponding to current time step
  IF (calltype == 'for') THEN
     read_pos = -step / delt_state_file + 1
  ELSE
     read_pos =  step / delt_state_file + 1
  END IF
  write (*,'(8x,a,i6)') &
       '--- Read true state at file position', read_pos


! ******************************
! *** Compute true RMS error ***
! ******************************
  comprms: IF (abs(step) < statefile_laststep(1)) THEN

     ! Read true state
     s = 1
     stat(s) = NF90_INQ_VARID(fileid_state, 'state', id_state)

     pos(2) = read_pos
     cnt(2) = 1
     pos(1) = 1
     cnt(1) = dim
     s = s + 1
     stat(s) = NF90_GET_VAR(fileid_state, id_state, state_true, start=pos(1:2), count=cnt(1:2))

     DO i = 1,  s
        IF (stat(i) /= NF90_NOERR) &
             WRITE(*, *) 'NetCDF error in reading true state from file, no.', i
     END DO

     ! Compute RMS
     trmse = 0.0
     DO j = 1, dim
        trmse = trmse + (state_true(j) - state_est(j))**2
     END DO
     trmse = SQRT(trmse / dim)

     
     ! *** Compute histogram information
     IF(calltype=='ana') THEN
        ! Histograms from beginning
        CALL PDAF_diag_histogram(1, dim, dim_ens, 1, &
             state_true, ens, hist_true(:,1), delta_hist, status)
        CALL PDAF_diag_histogram(1, dim, dim_ens, 1, &
             state_est, ens, hist_mean(:,1), delta_hist, status)
        IF (ABS(step) >= stepnull_means) THEN
           ! Histograms for after stepnull_means
           CALL PDAF_diag_histogram(1, dim, dim_ens, 1, &
                state_true, ens, hist_true(:,2), delta_hist, status)
           CALL PDAF_diag_histogram(1, dim, dim_ens, 1, &
                state_est, ens, hist_mean(:,2), delta_hist, status)
        END IF
     END IF

     ! *** Compute statistics of ensemble (skewness, kurtosis)
     CALL PDAF_diag_ensstats(dim, dim_ens, 1, &
          state_est, ens, skewness, kurtosis, status)

  ELSE
     ! We don't have true state information any more - set rmse to zero
     trmse = 0.0

  END IF comprms


! ********************
! *** Compute CRPS ***
! ********************

  IF (calltype=='ana') THEN
     CALL PDAF_diag_CRPS(dim, dim_ens, 0, ens, state_true, &
          crps, reli, resol, uncert, status)
     crps_stats (1) = crps
     crps_stats (2) = reli
     crps_stats (3) = resol
     crps_stats (4) = uncert

     WRITE (*,'(a,4es13.5)') 'CRPS:', crps, reli, resol, uncert
  END IF

END SUBROUTINE compute_truermse
!BOP
!
! !ROUTINE: close_netcdf_state  --- close netcdf file for true state
!
! !INTERFACE:
SUBROUTINE close_netcdf_state()

! !DESCRIPTION:
! This routine closes the netcdf file holding the true trajectory.
! It is called in the routine next_observation.

! !USES:
  USE netcdf
  USE mod_assimilation, &
       ONLY: fileid_state

  IMPLICIT NONE
!EOP

! Local variables
  INTEGER :: stat(50)                ! Array for status flag

! Close file
  stat(1) = NF90_CLOSE(fileid_state)
  IF (stat(1) /= NF90_NOERR) &
       WRITE(*, *) 'NetCDF error in closing file with true states, no. 1'

END SUBROUTINE close_netcdf_state
