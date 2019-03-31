!$Id: output_netcdf_asml.F90 1590 2015-06-12 11:59:35Z lnerger $
!BOP
!
! !MODULE:
MODULE output_netcdf_asml

! !DESCRIPTION: 
! This module provides routines to initialize
! NetCDF output files for the assimilation in the
! Lorenz96 model and to write output into the files.
!
! !REVISION HISTORY:
! 2010-01 - Lars Nerger - Initial code
! Later revisions - see SVN log
!
! !USES:
  IMPLICIT NONE
  SAVE
  PUBLIC

! !PUBLIC DATA MEMBERS:
  CHARACTER(len=100) :: file_asml = 'asml.nc' ! Name of the NetCDF output file
  INTEGER :: delt_write_asml = 1              ! Output interval in assimilation intervals
  LOGICAL :: write_states    = .TRUE.         ! Whether to write estimated states
  LOGICAL :: write_stats     = .FALSE.        ! Whether to write ensemble statistics
  LOGICAL :: write_ens       = .FALSE.        ! Whether to write full ensemble 

!EOP

! Private variables
  INTEGER, PRIVATE :: file_pos   ! File position to write to
  INTEGER, PRIVATE :: fileid     ! Id of netcdf file
  INTEGER, PRIVATE :: cnt_steps  ! Count time step for delt_write_asml

CONTAINS
!BOP
!
! !ROUTINE: init_netcdf_asml  --- initialize netcdf output
!
! !INTERFACE:
  SUBROUTINE init_netcdf_asml(step, dt, dim, filtertype, subtype, &
       dim_ens, forget, type_ensinit, local_range, local_range2, &
       locweight, srange, rms_obs, delt_obs, total_steps, &
       seedset, stepnull_means, dim_lag)

! !DESCRIPTION:
! This routine initializes the netcdf file 

! !USES:
    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

! !ARGUMENTS:
    INTEGER, INTENT(IN) :: step       ! Initial time step
    REAL, INTENT(IN)    :: dt         ! Size of time step
    INTEGER, INTENT(IN) :: dim        ! Dimension of model state
    INTEGER, INTENT(IN) :: filtertype ! Type of filter
    INTEGER, INTENT(IN) :: subtype    ! Sub-type of filter
    INTEGER, INTENT(IN) :: dim_ens    ! ensemble_size
    REAL, INTENT(IN)    :: forget     ! forgetting factor
    CHARACTER(len=3), INTENT(IN) :: type_ensinit ! Type of ensemble init
    REAL, INTENT(IN)    :: local_range  ! Localization radius
    REAL, INTENT(IN)    :: local_range2  ! Localization radius - right side
    INTEGER, INTENT(IN) :: locweight    ! Type of localization
    REAL, INTENT(IN)    :: rms_obs      ! RMS error of observations
    REAL, INTENT(IN)    :: srange       ! Support range for 5th order polynomial
                                        !   and range for 1/e for exponential weighting
    INTEGER, INTENT(IN) :: delt_obs     ! Number of time steps between two analysis steps
    INTEGER, INTENT(IN) :: total_steps  ! Total number of time steps in experiment
    INTEGER, INTENT(IN) :: seedset      ! Set id of seeds for random numbers in initialization
    INTEGER, INTENT(IN) :: stepnull_means ! Step at which computation of time mean errors is started
    INTEGER, INTENT(IN) :: dim_lag    ! Lag for smoothing
!EOP

! Local variables
    INTEGER :: i, s                   ! Counters
    INTEGER :: type_ensinit_id        ! ID for type of ensemble initialization
    INTEGER :: dimid_state, dimid_1   ! Dimension IDs
    INTEGER :: dimid_step             ! Dimension ID
    INTEGER :: dimid_ens, dimid_ensp1 ! Dimension IDs
    INTEGER :: dimid_lag              ! Dimension ID
    INTEGER :: ID_tmp                 ! Variable IDs
    INTEGER :: ID_state               ! Variable IDs
    INTEGER :: stat(250)              ! Array for status flag
    INTEGER :: dimarray(2)            ! Array for dimensions
    INTEGER :: dimarray3(3)           ! Array for dimensions
    CHARACTER(len=100) :: attstr      ! String to write attributes

! *** Initialize file ***    

! Print screen information
    WRITE (*, '(/1x, a)') 'Initialize NetCDF output'
    WRITE (*,'(5x,a,i6)') 'Writing interval (steps) ',delt_write_asml
    IF (write_states) WRITE (*,'(5x,a)') '--> Write model states'
    IF (write_ens) WRITE (*,'(5x,a)') '--> Write ensemble states'
    IF (write_stats) WRITE (*,'(5x,a)') '--> Write higher-order ensemble statistics'

! Initialize file position
    file_pos = 1

! Initialize counter for output interval
    cnt_steps = 1

! Initialize file and write global atributes

    s = 1

    stat(s) = NF_CREATE(TRIM(file_asml), 0, fileid) 
    s = s + 1

    attstr  = 'Assimilation into Lorenz96 model'
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'title', LEN_TRIM(attstr), &
         TRIM(attstr)) 
    s = s + 1

! Define Dimensions

    stat(s) = NF_DEF_DIM(fileid, 'dim_state', dim, dimid_state)             
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'iteration', NF_UNLIMITED, dimid_step)
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'one', 1, dimid_1)
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'dim_ens', dim_ens, dimid_ens)
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'dim_ensp1', dim_ens+1, dimid_ensp1)
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'dim_lag', dim_lag, dimid_lag)
    s = s + 1

! Define variables characterizing the experiment
    stat(s) = NF_DEF_VAR(fileid, 'filtertype', NF_INT, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'subtype', NF_INT, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'dim_ens', NF_INT, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'forget', NF_DOUBLE, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'step_null', NF_INT, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'total_steps', NF_INT, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'type_ensinit_id', NF_INT, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'local_range', NF_DOUBLE, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'local_range2', NF_DOUBLE, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'locweight', NF_INT, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'srange', NF_DOUBLE, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'rms_obs', NF_DOUBLE, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'delt_obs', NF_INT, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'seedset', NF_INT, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'stepnull_means', NF_INT, 1, DimId_1, Id_tmp) 
    s = s + 1


! Define variables
    
    s = 1
    stat(s) = NF_DEF_VAR(fileid, 'mrmse_for_null', NF_DOUBLE, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'mrmse_ana_null', NF_DOUBLE, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'mtrmse_for_null', NF_DOUBLE, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'mtrmse_ana_null', NF_DOUBLE, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'mrmse_for_step', NF_DOUBLE, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'mrmse_ana_step', NF_DOUBLE, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'mtrmse_for_step', NF_DOUBLE, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'mtrmse_ana_step', NF_DOUBLE, 1, DimId_1, Id_tmp) 
    s = s + 1

    smootherA: IF (dim_lag > 0) THEN
       stat(s) = NF_DEF_VAR(fileid, 'mrmse_smoother_null', NF_DOUBLE, 1, DimId_lag, Id_tmp) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'mtrmse_smoother_null', NF_DOUBLE, 1, DimId_lag, Id_tmp) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'mrmse_smoother_step', NF_DOUBLE, 1, DimId_lag, Id_tmp) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'mtrmse_smoother_step', NF_DOUBLE, 1, DimId_lag, Id_tmp) 
       s = s + 1
    END IF smootherA

    stat(s) = NF_DEF_VAR(fileid, 'step', NF_INT, 1, DimId_step, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'step_ini', NF_INT, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'time', NF_DOUBLE, 1, DimId_step, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'time_ini', NF_DOUBLE, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'rmse_ini', NF_DOUBLE, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'rmse_for', NF_DOUBLE, 1, DimId_step, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'rmse_ana', NF_DOUBLE, 1, DimId_step, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'trmse_ini', NF_DOUBLE, 1, DimId_1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'trmse_for', NF_DOUBLE, 1, DimId_step, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'trmse_ana', NF_DOUBLE, 1, DimId_step, Id_tmp) 
    s = s + 1

    smootherB: IF (dim_lag > 0) THEN
       dimarray(1) = dimid_lag
       dimarray(2) = dimid_step
       stat(s) = NF_DEF_VAR(fileid, 'rmse_smoother', NF_DOUBLE, 2, dimarray, Id_tmp) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'trmse_smoother', NF_DOUBLE, 2, dimarray, Id_tmp) 
       s = s + 1
    END IF smootherB

    stat(s) = NF_DEF_VAR(fileid, 'hist_true_null', NF_INT, 1, Dimid_ensp1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'hist_mean_null', NF_INT, 1, Dimid_ensp1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'hist_true_step', NF_INT, 1, Dimid_ensp1, Id_tmp) 
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'hist_mean_step', NF_INT, 1, Dimid_ensp1, Id_tmp) 
    s = s + 1

    writestats: IF (write_stats) THEN
       stat(s) = NF_DEF_VAR(fileid, 'skewness_ini', NF_DOUBLE, 1, DimId_1, Id_tmp) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'kurtosis_ini', NF_DOUBLE, 1, DimId_1, Id_tmp) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'skewness_for', NF_DOUBLE, 1, DimId_step, Id_tmp) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'kurtosis_for', NF_DOUBLE, 1, DimId_step, Id_tmp) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'skewness_ana', NF_DOUBLE, 1, DimId_step, Id_tmp) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'kurtosis_ana', NF_DOUBLE, 1, DimId_step, Id_tmp) 
       s = s + 1
    END IF writestats

    writestates: IF (write_states) THEN
       dimarray(1) = dimid_state
       dimarray(2) = dimid_1
       stat(s) = NF_DEF_VAR(fileid, 'state_ini', NF_DOUBLE, 2, dimarray, Id_state)
       s = s + 1
       dimarray(1) = dimid_state
       dimarray(2) = dimid_step
       stat(s) = NF_DEF_VAR(fileid, 'state_for', NF_DOUBLE, 2, dimarray, Id_state)
       s = s + 1
       dimarray(1) = dimid_state
       dimarray(2) = dimid_step
       stat(s) = NF_DEF_VAR(fileid, 'state_ana', NF_DOUBLE, 2, dimarray, Id_state)
       s = s + 1
    END IF writestates

    writeens: IF (write_ens) THEN
       dimarray3(1) = dimid_state
       dimarray3(2) = dimid_ens
       dimarray3(3) = dimid_1
       stat(s) = NF_DEF_VAR(fileid, 'ens_ini', NF_DOUBLE, 3, dimarray3, Id_state)
       s = s + 1
       dimarray3(1) = dimid_state
       dimarray3(2) = dimid_ens
       dimarray3(3) = dimid_step
       stat(s) = NF_DEF_VAR(fileid, 'ens_for', NF_DOUBLE, 3, dimarray3, Id_state)
       s = s + 1
       dimarray3(1) = dimid_state
       dimarray3(2) = dimid_ens
       dimarray3(3) = dimid_step
       stat(s) = NF_DEF_VAR(fileid, 'ens_ana', NF_DOUBLE, 3, dimarray3, Id_state)
       s = s + 1
    END IF writeens

    stat(s) = NF_ENDDEF(fileid) 
    s = s + 1

! Write variables characterizing the experiment

    stat(s) = NF_INQ_VARID(fileid, 'filtertype', Id_tmp) 
    s = s + 1
    stat(s) = NF_PUT_VAR_INT(fileid, Id_tmp, filtertype)
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, 'subtype', Id_tmp) 
    s = s + 1
    stat(s) = NF_PUT_VAR_INT(fileid, Id_tmp, subtype)
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, 'dim_ens', Id_tmp) 
    s = s + 1
    stat(s) = NF_PUT_VAR_INT(fileid, Id_tmp, dim_ens)
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, 'forget', Id_tmp) 
    s = s + 1
    stat(s) = NF_PUT_VAR_DOUBLE(fileid, Id_tmp, forget)
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, 'step_null', Id_tmp) 
    s = s + 1
    stat(s) = NF_PUT_VAR_INT(fileid, Id_tmp, step)
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, 'total_steps', Id_tmp) 
    s = s + 1
    stat(s) = NF_PUT_VAR_INT(fileid, Id_tmp, total_steps)
    s = s + 1
    IF (type_ensinit == 'eof') THEN
       type_ensinit_id = 1
    ELSE IF (type_ensinit == 'rnd') THEN
       type_ensinit_id = 2
    END IF
    stat(s) = NF_INQ_VARID(fileid, 'type_ensinit_id', Id_tmp) 
    s = s + 1
    stat(s) = NF_PUT_VAR_INT(fileid, Id_tmp, type_ensinit_id)
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, 'local_range', Id_tmp) 
    s = s + 1
    stat(s) = NF_PUT_VAR_DOUBLE(fileid, Id_tmp, local_range)
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, 'local_range2', Id_tmp) 
    s = s + 1
    stat(s) = NF_PUT_VAR_DOUBLE(fileid, Id_tmp, local_range2)
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, 'locweight', Id_tmp) 
    s = s + 1
    stat(s) = NF_PUT_VAR_INT(fileid, Id_tmp, locweight)
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, 'srange', Id_tmp) 
    s = s + 1
    stat(s) = NF_PUT_VAR_DOUBLE(fileid, Id_tmp, srange)
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, 'rms_obs', Id_tmp) 
    s = s + 1
    stat(s) = NF_PUT_VAR_DOUBLE(fileid, Id_tmp, rms_obs)
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, 'delt_obs', Id_tmp) 
    s = s + 1
    stat(s) = NF_PUT_VAR_INT(fileid, Id_tmp, delt_obs)
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, 'seedset', Id_tmp) 
    s = s + 1
    stat(s) = NF_PUT_VAR_INT(fileid, Id_tmp, seedset)
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, 'stepnull_means', Id_tmp) 
    s = s + 1
    stat(s) = NF_PUT_VAR_INT(fileid, Id_tmp, stepnull_means)
    s = s + 1

    DO i = 1,  s - 1
       IF (stat(i) /= NF_NOERR) &
            WRITE(*, *) 'NetCDF error in file initialization, no.', i
    END DO

  END SUBROUTINE init_netcdf_asml
!BOP
!
! !ROUTINE: write_netcdf_asml  --- write netcdf output during assimilation
!
! !INTERFACE:
  SUBROUTINE write_netcdf_asml(calltype, step, time, dim, state, &
       rmse, trmse, mrmse_null, mtrmse_null, mrmse_step, &
       mtrmse_step, dim_ens, ens, hist_true, hist_mean, &
       skewness, kurtosis, dim_lag, rmse_s, trmse_s, mrmse_s_null, &
       mtrmse_s_null, mrmse_s_step, mtrmse_s_step)

! !DESCRIPTION:
! This routine initializes the netcdf file 

! !USES:
    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

! !ARGUMENTS:
    CHARACTER(len=3)    :: calltype    ! Type of output call
    INTEGER, INTENT(IN) :: step        ! Current time step
    REAL, INTENT(IN)    :: time        ! Current model time
    INTEGER, INTENT(IN) :: dim         ! Dimension of model state
    REAL, INTENT(IN)    :: state(dim)  ! Model state
    REAL, INTENT(IN)    :: rmse        ! Estimated RMS error
    REAL, INTENT(IN)    :: trmse       ! True RMS error
    REAL, INTENT(IN)    :: mrmse_null  ! Time-mean estimated RMS error from step 0
    REAL, INTENT(IN)    :: mtrmse_null ! Time-mean true RMS error from step 0
    REAL, INTENT(IN)    :: mrmse_step  ! Time-mean estimated RMS error from stepnull_means
    REAL, INTENT(IN)    :: mtrmse_step ! Time-mean true RMS error from stepnull_means
    INTEGER, INTENT(IN) :: dim_ens     ! Ensemble size
    REAL, INTENT(IN)    :: ens(dim, dim_ens)       ! Ensemble
    INTEGER, INTENT(IN) :: hist_true(dim_ens+1, 2) ! Rank histogram about true state
    INTEGER, INTENT(IN) :: hist_mean(dim_ens+1, 2) ! Rank histogram about ensemble mean
    REAL, INTENT(IN)    :: skewness                ! Skewness of ensemble
    REAL, INTENT(IN)    :: kurtosis                ! Kurtosis of ensemble
    ! RMS errors for smoother
    INTEGER, INTENT(IN) :: dim_lag             ! Size of lag for smoothing
    REAL, INTENT(IN) :: rmse_s(dim_lag)        ! Estimated RMS error
    REAL, INTENT(IN) :: trmse_s(dim_lag)       ! True RMS error
    REAL, INTENT(IN) :: mrmse_s_null(dim_lag)  ! Time-mean estimated RMS error from step 0
    REAL, INTENT(IN) :: mtrmse_s_null(dim_lag) ! Time-mean true RMS error from step 0
    REAL, INTENT(IN) :: mrmse_s_step(dim_lag)  ! Time-mean estimated RMS error from stepnull_means
    REAL, INTENT(IN) :: mtrmse_s_step(dim_lag) ! Time-mean true RMS error from stepnull_means
!EOP

! Local variables
    INTEGER :: i, s                    ! Counters
    INTEGER :: ID_time, ID_step        ! Variable IDs
    INTEGER :: ID_state                ! Variable IDs
    INTEGER :: ID_ens                  ! Variable IDs
    INTEGER :: ID_rmse, ID_trmse       ! Variable IDs
    INTEGER :: ID_mrmseN, ID_mtrmseN   ! Variable IDs
    INTEGER :: ID_mrmseS, ID_mtrmseS   ! Variable IDs
    INTEGER :: ID_hist_true_null, ID_hist_mean_null ! Variable IDs
    INTEGER :: ID_hist_true_step, ID_hist_mean_step ! Variable IDs
    INTEGER :: ID_skew, ID_kurt        ! Variable IDs
    INTEGER :: Id_rmse_s, Id_trmse_s, Id_mrmseN_s      ! Variable IDs for smoother output
    INTEGER :: Id_mtrmseN_s, Id_mrmseS_s, Id_mtrmseS_s ! Variable IDs for smoother output
    INTEGER :: stat(100)               ! Array for status flag
    INTEGER :: pos(2)                  ! Position index for writing
    INTEGER :: cnt(2)                  ! Count index for writing
    INTEGER :: pos3(3)                 ! Position index for writing
    INTEGER :: cnt3(3)                 ! Count index for writing
    LOGICAL :: dowrite                 ! Flag whether to write at the current call


! Check, if we have to write states at this time step
    IF (cnt_steps==delt_write_asml) THEN
       dowrite = .TRUE.
       IF (calltype == 'ana') THEN
          cnt_steps = 1
       END IF
    ELSE
       dowrite = .FALSE.
       IF (calltype == 'ana') THEN
          cnt_steps = cnt_steps + 1
       END IF
    END IF

! Inquire variable Ids

    s = 1
    IF (calltype == 'ini') THEN

       IF (write_states) THEN
          stat(s) = NF_INQ_VARID(fileid, "state_ini", Id_state)
          s = s + 1
       END IF
       IF (write_ens) THEN
          stat(s) = NF_INQ_VARID(fileid, "ens_ini", Id_ens)
          s = s + 1
       END IF
       stat(s) = NF_INQ_VARID(fileid, "time_ini", Id_time)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "step_ini", Id_step)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "rmse_ini", Id_rmse)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "trmse_ini", Id_trmse)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "mrmse_ana_null", Id_mrmseN)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "mtrmse_ana_null", Id_mtrmseN)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "mrmse_ana_step", Id_mrmseS)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "mtrmse_ana_step", Id_mtrmseS)
       s = s + 1
       writestats: IF (write_stats) THEN
          stat(s) = NF_INQ_VARID(fileid, "skewness_ini", Id_skew)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "kurtosis_ini", Id_kurt)
          s = s + 1
       END IF writestats

    ELSE IF (calltype == 'for') THEN

       IF (write_states) THEN
          stat(s) = NF_INQ_VARID(fileid, "state_for", Id_state)
          s = s + 1
       END IF
       IF (write_ens) THEN
          stat(s) = NF_INQ_VARID(fileid, "ens_for", Id_ens)
          s = s + 1
       END IF
       stat(s) = NF_INQ_VARID(fileid, "time", Id_time) 
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "step", Id_step)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "rmse_for", Id_rmse)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "trmse_for", Id_trmse)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "mrmse_for_null", Id_mrmseN)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "mtrmse_for_null", Id_mtrmseN)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "mrmse_for_step", Id_mrmseS)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "mtrmse_for_step", Id_mtrmseS)
       s = s + 1
       writestatsB: IF (write_stats) THEN
          stat(s) = NF_INQ_VARID(fileid, "skewness_for", Id_skew)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "kurtosis_for", Id_kurt)
          s = s + 1
       END IF writestatsB

    ELSE

       IF (write_states) THEN
          stat(s) = NF_INQ_VARID(fileid, "state_ana", Id_state)
          s = s + 1
       END IF
       IF (write_ens) THEN
          stat(s) = NF_INQ_VARID(fileid, "ens_ana", Id_ens)
          s = s + 1
       END IF
       stat(s) = NF_INQ_VARID(fileid, "time", Id_time) 
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "step", Id_step)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "rmse_ana", Id_rmse)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "trmse_ana", Id_trmse)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "mrmse_ana_null", Id_mrmseN)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "mtrmse_ana_null", Id_mtrmseN)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "mrmse_ana_step", Id_mrmseS)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "mtrmse_ana_step", Id_mtrmseS)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "hist_true_null", Id_hist_true_null)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "hist_mean_null", Id_hist_mean_null)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "hist_true_step", Id_hist_true_step)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "hist_mean_step", Id_hist_mean_step)
       s = s + 1
       writestatsC: IF (write_stats) THEN
          stat(s) = NF_INQ_VARID(fileid, "skewness_ana", Id_skew)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "kurtosis_ana", Id_kurt)
          s = s + 1
       END IF writestatsC

    END IF

    DO i = 1,  s - 1
       IF (stat(i) /= NF_NOERR) &
            WRITE(*, *) 'NetCDF error in preparing output, no.', i
    END DO
    s = 1

    smootherA: IF (dim_lag > 0 .AND. calltype=='ana') THEN
       stat(s) = NF_INQ_VARID(fileid, "rmse_smoother", Id_rmse_s)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "trmse_smoother", Id_trmse_s)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "mrmse_smoother_null", Id_mrmseN_s)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "mtrmse_smoother_null", Id_mtrmseN_s)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "mrmse_smoother_step", Id_mrmseS_s)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "mtrmse_smoother_step", Id_mtrmseS_s)
       s = s + 1
    END IF smootherA

    DO i = 1,  s - 1
       IF (stat(i) /= NF_NOERR) &
            WRITE(*, *) 'NetCDF error in preparing smoother output, no.', i
    END DO
    s = 1

! Write variables
    IF (calltype == 'ini') THEN
       pos(1) = file_pos
       cnt(1) = 1
       stat(s) = NF_PUT_VARA_DOUBLE(fileid, Id_time, pos(1), cnt(1), time)
       s = s + 1
       stat(s) = NF_PUT_VARA_INT(fileid, Id_step, pos(1), cnt(1), step)
       s = s + 1
    END IF

    IF (calltype == 'for' .AND. dowrite) THEN
       pos(1) = file_pos
       cnt(1) = 1
       stat(s) = NF_PUT_VARA_DOUBLE(fileid, Id_time, pos(1), cnt(1), time)
       s = s + 1
       stat(s) = NF_PUT_VARA_INT(fileid, Id_step, pos(1), cnt(1), -step)
       s = s + 1
    END IF

    ! Write state information and RMS errors only at specified intervals
    ! and at initialization
    writetimedep: IF (dowrite .OR. calltype=='ini') THEN
       IF (write_states) THEN
          pos(1) = 1
          pos(2) = file_pos
          cnt(1) = dim
          cnt(2) = 1
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, Id_state, pos, cnt, state)
          s = s + 1
       END IF
       IF (write_ens) THEN
          pos3(1) = 1
          pos3(2) = 1
          pos3(3) = file_pos
          cnt3(1) = dim
          cnt3(2) = dim_ens
          cnt3(3) = 1
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, Id_ens, pos3, cnt3, ens)
          s = s + 1
       END IF

       pos(1) = file_pos
       cnt(1) = 1
       stat(s) = NF_PUT_VARA_DOUBLE(fileid, Id_rmse, pos(1), cnt(1), rmse)
       s = s + 1
    
       pos(1) = file_pos
       cnt(1) = 1
       stat(s) = NF_PUT_VARA_DOUBLE(fileid, Id_trmse, pos(1), cnt(1), trmse)
       s = s + 1

       ! Write RMS errors from smoothing
       smootherB: IF (dim_lag > 0 .AND. calltype=='ana') THEN
          pos(1) = 1
          pos(2) = file_pos
          cnt(1) = dim_lag
          cnt(2) = 1
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, Id_rmse_s, pos, cnt, rmse_s)
          s = s + 1
    
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, Id_trmse_s, pos, cnt, trmse_s)
          s = s + 1
       END IF smootherB

       IF (write_stats) THEN
          ! Ensemble statistics
          pos(1) = file_pos
          cnt(1) = 1
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, Id_skew, pos(1), cnt(1), skewness)
          s = s + 1
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, Id_kurt, pos(1), cnt(1), kurtosis)
          s = s + 1
       END IF

    END IF writetimedep

    IF (calltype=='ana') THEN
       ! Histogram information - written at each time
       stat(s) = NF_PUT_VAR_INT(fileid, Id_hist_true_null, hist_true(:, 1))
       s = s + 1
       stat(s) = NF_PUT_VAR_INT(fileid, Id_hist_mean_null, hist_mean(:, 1))
       s = s + 1
       stat(s) = NF_PUT_VAR_INT(fileid, Id_hist_true_step, hist_true(:, 2))
       s = s + 1
       stat(s) = NF_PUT_VAR_INT(fileid, Id_hist_mean_step, hist_mean(:, 2))
       s = s + 1
    END IF

    ! Mean errors are written at each time
    stat(s) = NF_PUT_VAR_DOUBLE(fileid, Id_mrmseN, mrmse_null)
    s = s + 1
    stat(s) = NF_PUT_VAR_DOUBLE(fileid, Id_mtrmseN, mtrmse_null)
    s = s + 1
    stat(s) = NF_PUT_VAR_DOUBLE(fileid, Id_mrmseS, mrmse_step)
    s = s + 1
    stat(s) = NF_PUT_VAR_DOUBLE(fileid, Id_mtrmseS, mtrmse_step)
    s = s + 1

    ! Write mean smoother errors at each time
    smootherC: IF (dim_lag > 0 .AND. calltype=='ana') THEN
       stat(s) = NF_PUT_VAR_DOUBLE(fileid, Id_mrmseN_s, mrmse_s_null)
       s = s + 1
       stat(s) = NF_PUT_VAR_DOUBLE(fileid, Id_mtrmseN_s, mtrmse_s_null)
       s = s + 1
       stat(s) = NF_PUT_VAR_DOUBLE(fileid, Id_mrmseS_s, mrmse_s_step)
       s = s + 1
       stat(s) = NF_PUT_VAR_DOUBLE(fileid, Id_mtrmseS_s, mtrmse_s_step)
       s = s + 1
    END IF smootherC

    DO i = 1,  s - 1
       IF (stat(i) /= NF_NOERR) &
            WRITE(*, *) 'NetCDF error in writing output, no.', i
    END DO

    ! Increment file position counter
    IF (dowrite .AND. calltype == 'ana') THEN
       file_pos = file_pos + 1
    END IF

  END SUBROUTINE write_netcdf_asml
!BOP
!
! !ROUTINE: close_netcdf_asml  --- close netcdf file
!
! !INTERFACE:
  SUBROUTINE close_netcdf_asml()

! !DESCRIPTION:
! This routine closes the netcdf file.
! It is called in the routine next_observation.

! !USES:
    IMPLICIT NONE

    INCLUDE 'netcdf.inc'
!EOP

! Local variables
    INTEGER :: stat(50)                ! Array for status flag

! Close file

    stat(1) = NF_CLOSE(fileid)
    IF (stat(1) /= NF_NOERR) &
         WRITE(*, *) 'NetCDF error in closing output file, no. 1'

  END SUBROUTINE close_netcdf_asml

END MODULE output_netcdf_asml
