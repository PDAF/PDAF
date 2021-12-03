!$Id$
!> Module holding routines to output assimilation results into NetCDF
!!
!! This module provides routines to initialize NetCDF
!! output files for ECHAM and to write output into the files.
!!
!! __Revision history:__
!! 2020-04 - Qi Tang - Initial code based on output_netcdf module of FESOM
!! * Later revisions - see repository log
!!
MODULE output_pdaf

  IMPLICIT NONE
  SAVE
  PUBLIC

! *** Public variables ***
  CHARACTER(len=100) :: str_daspec='DA'        ! String to identify assimilation experiment
  CHARACTER(len=1)   :: prec_nc = 's'          ! Precision of NetCDF output
                                               ! (s) single, (d) double
  INTEGER :: write_pos_da = 1                  ! Counter for next time slice to be written
  INTEGER :: write_pos_da_ens
  LOGICAL :: write_da = .true.                 ! Whether to write output file from assimilation
  LOGICAL :: write_ens = .true.                ! Whether to write output file for each individual ensemble member


! *** Private variables ***
  LOGICAL, PRIVATE :: debugoutput=.FALSE.   ! Write output for debugging
                                            ! (file contains only last writing)
  INTEGER, PRIVATE :: nf_prec               ! Precision of NetCDF output

CONTAINS

!----------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_output_pdaf - Initialize NetCDf output file
!
! !INTERFACE: 
  SUBROUTINE init_output_pdaf(dim_lag, writepe)

! !USES:
    USE mod_assim_pdaf, &
         ONLY: dim_ens

    IMPLICIT NONE

! !ARGUMENTS:    
    INTEGER, INTENT(in) :: dim_lag              ! Smoother lag
    LOGICAL, INTENT(in) :: writepe
!EOP

    IF (writepe .AND. write_da) THEN
       ! Initialize ehcam file
       CALL init_ncfile_echam_pdaf(dim_lag, writepe)
    END IF

  END SUBROUTINE init_output_pdaf
!----------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_ncfile_echam_pdaf - Initialize NetCDf output file for
! atmosphere fields
!
! !INTERFACE: 

  SUBROUTINE init_ncfile_echam_pdaf(dim_lag, writepe)

! !USES:

    USE mo_filename, &
        ONLY: out_datapath
    USE mo_decomposition, &
        ONLY: gc => global_decomposition
    USE mo_control, &
        ONLY: nvclev

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

! !ARGUMENTS:    
    INTEGER, INTENT(in) :: dim_lag              ! Smoother lag
    LOGICAL, INTENT(in) :: writepe
!EOP 

! Local variables
    CHARACTER(200)     :: filename          ! Full name of output file
    CHARACTER(len=100) :: attstr            ! String to write attributes
    INTEGER :: stat(500)                    ! auxiliary: status array
    INTEGER :: fileid                       ! ID of netCDF file
    INTEGER :: s                            ! auxiliary: status counter
    INTEGER :: dimid_lat, dimid_lon, dimid_lev, dimid_time      ! IDs for dimension
    INTEGER :: dimarray(4)                  ! auxiliary: array dimension
    INTEGER :: varid_tf,varid_ta,varid_ti,varid_lon,varid_lat,varid_lev,varid_time                      ! IDs for variables
    INTEGER :: varid_vof, varid_voa, varid_voi
    INTEGER :: varid_df, varid_da, varid_di
    INTEGER :: varid_qf, varid_qa, varid_qi
    INTEGER :: varid_uf, varid_ua, varid_ui
    INTEGER :: varid_vf, varid_va, varid_vi
    INTEGER :: i                            ! Counters
    INTEGER :: dimid_nhym, dimid_nhyi
    INTEGER :: dimid_one
    INTEGER :: varid_hyai, varid_hybi, varid_hyam, varid_hybm
    INTEGER :: varid_lspf, varid_lspa, varid_lspi
    INTEGER :: varid_rmstini, varid_rmslspini, varid_rmsvoini, varid_rmsdini
    INTEGER :: varid_rmsqini, varid_rmsuini, varid_rmsvini
    INTEGER :: varid_rmstf, varid_rmslspf, varid_rmsvof, varid_rmsdf
    INTEGER :: varid_rmsqf, varid_rmsuf, varid_rmsvf
    INTEGER :: varid_rmsta, varid_rmslspa, varid_rmsvoa, varid_rmsda
    INTEGER :: varid_rmsqa, varid_rmsua, varid_rmsva




    pe0: IF (writepe) THEN

! Print screen information
          WRITE (*, '(/a, 1x, a)') 'FESOM-PDAF', 'Initialize assimilation NetCDF atmosphere file - single precision'
          nf_prec = NF_FLOAT

! ----- open file and write global attributes
       filename=TRIM(out_datapath)//'echam_'//TRIM(str_daspec)//'.nc'

       s = 1

       stat(s) = NF_CREATE(filename, 0, fileid)
       s = s + 1

       attstr  = 'ECHAM Assimilation'
       stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'title', LEN_TRIM(attstr), &
            TRIM(attstr))
       s = s + 1


! ----- DEFINE DIMENSIONS ------------------------------

       stat(s) = NF_DEF_DIM(fileid, 'lon', gc%nlon, dimid_lon)
       s = s + 1
       stat(s) = NF_DEF_DIM(fileid, 'lat', gc%nlat, dimid_lat)
       s = s + 1
       stat(s) = NF_DEF_DIM(fileid, 'lev', gc%nlev, dimid_lev)
       s = s + 1
       stat(s) = NF_DEF_DIM(fileid, 'nhym', gc%nlev, dimid_nhym)
       s = s + 1
       stat(s) = NF_DEF_DIM(fileid, 'nhyi', nvclev, dimid_nhyi)
       s = s + 1
       stat(s) = NF_DEF_DIM(fileid, 'time', NF_UNLIMITED, dimid_time)
       s = s + 1
       stat(s) = NF_DEF_DIM(fileid, 'one', 1, dimid_one)

       DO i = 1,  s
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in ATM dimension definitions, no.', i
       END DO


! ----- DEFINE VARIABLES ---------------------------------

       s = 1

       !- numbers

       stat(s) = NF_DEF_VAR(fileid, 'lon', NF_DOUBLE, 1, dimid_lon, varid_lon)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'lat', NF_DOUBLE, 1, dimid_lat, varid_lat)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'lev', NF_DOUBLE, 1, dimid_lev, varid_lev)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'hyai', NF_DOUBLE, 1, dimid_nhyi, varid_hyai)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'hybi', NF_DOUBLE, 1, dimid_nhyi, varid_hybi)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'hyam', NF_DOUBLE, 1, dimid_nhym, varid_hyam)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'hybm', NF_DOUBLE, 1, dimid_nhym, varid_hybm)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'time', NF_DOUBLE, 1, dimid_time, varid_time)
       s = s + 1

! ----- F I E L D S 

       ! 3D temperature 
       dimarray(4) = dimid_time
       dimarray(3) = dimid_lev
       dimarray(2) = dimid_lat
       dimarray(1) = dimid_lon

       stat(s) = NF_DEF_VAR(fileid, 't_i', NF_FLOAT, 3, dimarray(1:3),varid_ti)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 't_f', NF_FLOAT, 4, dimarray,varid_tf)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 't_a', NF_FLOAT, 4,dimarray,varid_ta)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'vo_i', NF_FLOAT, 3, dimarray(1:3),varid_voi)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'vo_f', NF_FLOAT, 4, dimarray,varid_vof)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'vo_a', NF_FLOAT, 4,dimarray,varid_voa)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'd_i', NF_FLOAT, 3, dimarray(1:3),varid_di)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'd_f', NF_FLOAT, 4, dimarray,varid_df)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'd_a', NF_FLOAT, 4,dimarray,varid_da)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'q_i', NF_FLOAT, 3, dimarray(1:3),varid_qi)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'q_f', NF_FLOAT, 4, dimarray,varid_qf)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'q_a', NF_FLOAT, 4,dimarray,varid_qa)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'u_i', NF_FLOAT, 3, dimarray(1:3),varid_ui)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'u_f', NF_FLOAT, 4, dimarray,varid_uf)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'u_a', NF_FLOAT, 4,dimarray,varid_ua)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'v_i', NF_FLOAT, 3, dimarray(1:3),varid_vi)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'v_f', NF_FLOAT, 4, dimarray,varid_vf)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'v_a', NF_FLOAT, 4,dimarray,varid_va)

       dimarray(3) = dimid_time

       ! 2D log surface temperature
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'lsp_i', NF_FLOAT, 2, dimarray(1:2),varid_lspi)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'lsp_f', NF_FLOAT, 3, dimarray(1:3),varid_lspf)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'lsp_a', NF_FLOAT, 3,dimarray(1:3),varid_lspa)

       DO i = 1,  s
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in ATM variable definition, no.', i
       END DO

       !- RMS errors

       s = 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_t_ini', NF_FLOAT, 1, dimid_one, varid_rmstini)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_vo_ini', NF_FLOAT, 1, dimid_one, varid_rmsvoini)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_d_ini', NF_FLOAT, 1, dimid_one, varid_rmsdini)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_q_ini', NF_FLOAT, 1, dimid_one, varid_rmsqini)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_u_ini', NF_FLOAT, 1, dimid_one, varid_rmsuini)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_v_ini', NF_FLOAT, 1, dimid_one, varid_rmsvini)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_lsp_ini', NF_FLOAT, 1, dimid_one, varid_rmslspini)

       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_t_f', NF_FLOAT, 1, dimid_time, varid_rmstf)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_vo_f', NF_FLOAT, 1, dimid_time, varid_rmsvof)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_d_f', NF_FLOAT, 1, dimid_time, varid_rmsdf)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_q_f', NF_FLOAT, 1, dimid_time, varid_rmsqf)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_u_f', NF_FLOAT, 1, dimid_time, varid_rmsuf)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_v_f', NF_FLOAT, 1, dimid_time, varid_rmsvf)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_lsp_f', NF_FLOAT, 1, dimid_time, varid_rmslspf)

       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_t_a', NF_FLOAT, 1, dimid_time, varid_rmsta)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_vo_a', NF_FLOAT, 1, dimid_time, varid_rmsvoa)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_d_a', NF_FLOAT, 1, dimid_time, varid_rmsda)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_q_a', NF_FLOAT, 1, dimid_time, varid_rmsqa)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_u_a', NF_FLOAT, 1, dimid_time, varid_rmsua)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_v_a', NF_FLOAT, 1, dimid_time, varid_rmsva)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_lsp_a', NF_FLOAT, 1, dimid_time, varid_rmslspa)

! ----- DEFINE ATTRIBUTES -----------------------------------------

       s = 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_lon, 'standard_name', 9,'longitude');
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_lon, 'long_name', 9,'longitude');
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_lon, 'units', 12,'degrees_east');
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_lon, 'axis', 1,'X');
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_lat, 'standard_name', 8,'latitude');
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_lat, 'long_name', 8,'latitude');
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_lat, 'units', 13,'degrees_north');
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_lat, 'axis', 1,'Y');
       s = s + 1
       attstr  = 'hybrid_sigma_pressure'
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_lev, 'standard_name',LEN_TRIM(attstr),attstr);
       s = s + 1
       attstr  = 'hybrid level at layer midpoints'
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_lev, 'long_name',LEN_TRIM(attstr),attstr);
       s = s + 1
       attstr  = 'hyam hybm (mlev=hyam+hybm*aps)'
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_lev, 'formula',LEN_TRIM(attstr),attstr);
       s = s + 1
       attstr  = 'ap: hyam b: hybm ps: aps'
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_lev, 'formula_terms',LEN_TRIM(attstr),attstr);
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_lev, 'units', 5, 'level');
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_lev, 'positive', 4, 'down');
       s = s + 1
       attstr  = 'hybrid A coefficient at layer interfaces'
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_hyai, 'long_name',LEN_TRIM(attstr),attstr);
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_hyai, 'units', 2, 'Pa');
       s = s + 1
       attstr  = 'hybrid B coefficient at layer interfaces'
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_hybi,'long_name',LEN_TRIM(attstr),attstr);
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_hybi, 'units', 1, '1');
       s = s + 1
       attstr  = 'hybrid A coefficient at layer midpoints'
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_hyam,'long_name',LEN_TRIM(attstr),attstr);
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_hyam, 'units', 2, 'Pa');
       s = s + 1
       attstr  = 'hybrid B coefficient at layer midpoints'
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_hybm,'long_name',LEN_TRIM(attstr),attstr);
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_hybm, 'units', 1, '1');
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_time, 'standard_name', 4, 'time');
       s = s + 1
       attstr  = 'days since the starting day'
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_time, 'units',LEN_TRIM(attstr),attstr);
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_time, 'axis', 1, 'T');

       ! Fields
       s = s + 1
       attstr  = 'initial temperature'
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_ti,'long_name', LEN_TRIM(attstr),attstr);
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_ti,'units', 1 , 'K');
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_ti,'grid_type', 8, 'gaussian');
       s = s + 1
       attstr  = 'forecast temperature'
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_tf,'long_name', LEN_TRIM(attstr),attstr);
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_tf,'units', 1 , 'K');
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_tf,'grid_type', 8, 'gaussian');
       s = s + 1
       attstr  = 'analysed temperature'
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_ta,'long_name', LEN_TRIM(attstr),attstr);
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_ta,'units', 1 , 'K');
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_ta,'grid_type', 8, 'gaussian');
       s = s + 1
       attstr = 'forecast log surface pressure'
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_lspf, 'long_name',LEN_TRIM(attstr),attstr);
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_lspf, 'grid_type', 8, 'gaussian');
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_lspf, 'axis', 3, 'T--');
       s = s + 1
       attstr = 'analysed log surface pressure'
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_lspa, 'long_name',LEN_TRIM(attstr),attstr);
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_lspa, 'grid_type', 8, 'gaussian');
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, varid_lspa, 'axis', 3, 'T--');

       s = s + 1
       stat(s) = NF_ENDDEF(fileid)
       s = s + 1
       stat(1) = NF_CLOSE(fileid)

       IF (stat(1) /= NF_NOERR) THEN
          WRITE(*, *) 'NetCDF error in closing ATM NetCDF file'
       END IF
    END IF pe0

  END SUBROUTINE init_ncfile_echam_pdaf
!------------------------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: write_netcdf_pdaf - Write global fields into NetCDF files
!
! !INTERFACE: 
  SUBROUTINE write_netcdf_pdaf(writetype, write_pos_da, iteration, dim, state_l, nfields, rmse, writepe)

! ! USES:
    IMPLICIT NONE

! !ARGUMENTS:    
    CHARACTER(len=1), INTENT(in) :: writetype     ! Write (i) initial, (a) assimilated, (f) forecast fields
    INTEGER, INTENT(in) :: write_pos_da           ! Write position
    INTEGER, INTENT(in) :: iteration              ! Current model time step
    INTEGER, INTENT(in) :: dim                    ! Size of state vector
    REAL, INTENT(in) :: state_l(dim)              ! State vector
    INTEGER, INTENT(in) :: nfields                ! number of fields in state vector
    REAL, INTENT(in) :: rmse(nfields)             ! Array of RMS errors
    LOGICAL, INTENT(in) :: writepe
!EOP

    ! Write temperature fields
    CALL write_nc_airt_pdaf(writetype, write_pos_da, iteration, dim, state_l, &
         nfields, rmse, writepe)

    ! Write other fields

  END SUBROUTINE write_netcdf_pdaf
!-----------------------------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: write_nc_airt_pdaf - Write global atm temperature fields into NetCDF file
!
! !INTERFACE: 
  SUBROUTINE write_nc_airt_pdaf(writetype, write_pos_da, iteration, dim, state_l, nfields, rmse, writepe)

! !USES:

    USE mo_filename, &
        ONLY: out_datapath
    USE mo_decomposition, &
        ONLY: ldc => local_decomposition, gdc => global_decomposition
    USE mo_tr_gather, &
        ONLY: gather_field
    USE mo_time_control,&
        ONLY: get_time_step, get_date_components, current_date
    USE mo_geoloc, &
        ONLY: philat_2d, philon_2d
    USE mo_transpose, &
        ONLY: gather_gp
    USE mo_control, &
        ONLY: nvclev, vct, nlev
    USE mod_assim_pdaf, &
        ONLY: off_fields_p
    USE mod_assim_atm_pdaf, ONLY: dp, wp


    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

! !ARGUMENTS:    
    CHARACTER(len=1), INTENT(in) :: writetype     ! Write (i) initial, (a) assimilated, (f) forecast fields
    INTEGER, INTENT(in) :: write_pos_da           ! Write position
    INTEGER, INTENT(in) :: iteration              ! Current model time step
    INTEGER, INTENT(in) :: dim                    ! Size of state vector
    REAL, INTENT(in) :: state_l(dim)              ! State vector
    INTEGER, INTENT(in) :: nfields                ! number of fields in state vector
    REAL, INTENT(in) :: rmse(nfields)             ! Array of RMS errors
    LOGICAL, INTENT(in) :: writepe
!EOP

! Local variables
    CHARACTER(200)            :: filename
    INTEGER :: stat(50)                  ! Status array
    INTEGER :: FileId                    ! Id of Netcdf file
    INTEGER :: VarId_iter, VarId_time    ! Id numbers
    INTEGER :: varid_lon,varid_lat,varid_lev   !Ids
    INTEGER :: varid_t, varid_lsp, varid_vo, varid_d, varid_q, varid_u, varid_v
    INTEGER :: varid_hyai, varid_hybi, varid_hyam, varid_hybm
    INTEGER, ALLOCATABLE ::tmp_i(:)
    INTEGER :: s                         ! status counter
    INTEGER :: i                         ! Counters
    INTEGER :: nlon, nlat
    REAL(dp), POINTER :: gl_t(:,:,:),g_t(:,:,:)
    REAL(wp), POINTER :: gl_lon(:,:)
    REAL(wp), POINTER :: gl_lat(:,:)
    REAL(dp), POINTER :: gl_lsp(:,:)
    REAL(dp), POINTER :: gl_vo(:,:,:),g_vo(:,:,:)
    REAL(dp), POINTER :: gl_d(:,:,:),g_d(:,:,:)
    REAL(dp), POINTER :: gl_q(:,:,:),g_q(:,:,:)
    REAL(dp), POINTER :: gl_u(:,:,:),g_u(:,:,:)
    REAL(dp), POINTER :: gl_v(:,:,:),g_v(:,:,:)
    REAL(dp), POINTER :: local_3d(:,:,:), local_2d(:,:)

    INTEGER :: global_dim(3)
    INTEGER :: istep
    INTEGER :: pos1, nmb                 ! Position index for writing
    INTEGER :: current_day,current_year,current_month
    INTEGER :: pos1vec4(4), nmbvec4(4)   ! Position index arrays for writing
    INTEGER :: pos1vec3(3), nmbvec3(3)   ! Position index arrays for writing
    REAL(dp), DIMENSION(nvclev) :: hyai, hybi
    REAL(wp) :: hyam(nlev), hybm(nlev)
    INTEGER :: k, jk, jrow, jl          ! Counter
    INTEGER :: nproma, ngpblks
    INTEGER :: varid_rmst, varid_rmslsp, varid_rmsvo, varid_rmsd, varid_rmsq
    INTEGER :: varid_rmsu, varid_rmsv

  filename=TRIM(out_datapath)//'echam_'//TRIM(str_daspec)//'.nc'
! Get current time step
  istep = get_time_step()

! Print screen information
    IF (writepe) THEN
       IF (writetype == 'i') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'ECHAM-PDAF', 'Write initial atmosphere state to NetCDF at step ', &
               istep, ' position ', write_pos_da
       ELSE IF (writetype== 'f') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'ECHAM-PDAF', 'Write atmosphere forecast to NetCDF at step ', &
               istep-1, ' position ', write_pos_da
       ELSE IF (writetype== 'a') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'ECHAM-PDAF', 'Write atmosphere analysis to NetCDF at step ', &
               istep-1, ' position ', write_pos_da
       END IF
    END IF


  nlon=ldc%nlon
  nlat=ldc%nlat
  ngpblks=ldc%ngpblks

! Gather full fields and write into files
  ! Gather coordinates
  ALLOCATE(gl_lon(nlon,nlat))
  ALLOCATE(gl_lat(nlon,nlat))
  CALL gather_field(gl_lon,philon_2d)
  CALL gather_field(gl_lat,philat_2d) 

  hyai = vct(1:nvclev)
  hybi = vct(nvclev+1:2*nvclev)
  DO i=1,nlev
    hyam(i)=0.5_wp*(hyai(i)+hyai(i+1))
    hybm(i)=0.5_wp*(hybi(i)+hybi(i+1))
  END DO

  ! Gather global temperature field
  ALLOCATE (gl_t(nlon,nlev,nlat))
  global_dim(1) = nlon
  global_dim(3) = nlev
  global_dim(2) = nlat

  ALLOCATE(local_3d(ldc%nproma,nlev,ngpblks))
! 3D temperature
  k = 1
  DO jk = nlev,1,-1
    DO jrow = 1, ngpblks

      IF ( jrow == ngpblks ) THEN
        nproma = ldc%npromz
      ELSE
        nproma = ldc%nproma
      END IF

      DO jl = 1, nproma
        local_3d(jl,jk,jrow) = state_l(k)
        k = k + 1
      END DO

    END DO
  END DO

  CALL gather_field(gl_t,global_dim,local_3d)
  DEALLOCATE(local_3d)

  ! Gather global log surface pressure field
  ALLOCATE(gl_lsp(nlon,nlat))
  ALLOCATE(local_2d(ldc%nproma,ngpblks))

! 2D log surface pressure
  k = 1
  DO jrow = 1, ngpblks
    IF ( jrow == ngpblks ) THEN
       nproma = ldc%npromz
    ELSE
       nproma = ldc%nproma
    END IF

    DO jl = 1, nproma
      local_2d(jl,jrow) = state_l(k+off_fields_p(2))
      k = k + 1
    END DO

  END DO

  CALL gather_field(gl_lsp,local_2d)
  DEALLOCATE(local_2d)  

  ! Gather global vorticity
  ALLOCATE (gl_vo(nlon,nlev,nlat))
  ALLOCATE(local_3d(ldc%nproma,nlev,ngpblks))

! 3D vorticity
  k = 1
  DO jk = nlev,1,-1
    DO jrow = 1, ngpblks

      IF ( jrow == ngpblks ) THEN
        nproma = ldc%npromz
      ELSE
        nproma = ldc%nproma
      END IF

      DO jl = 1, nproma
        local_3d(jl,jk,jrow) = state_l(k+off_fields_p(3))
        k = k + 1
      END DO

    END DO
  END DO

  CALL gather_field(gl_vo,global_dim,local_3d)
  DEALLOCATE(local_3d)

  ! Gather global divergence
  ALLOCATE (gl_d(nlon,nlev,nlat))
  ALLOCATE(local_3d(ldc%nproma,nlev,ngpblks))

! 3D divergence
  k = 1
  DO jk = nlev,1,-1
    DO jrow = 1, ngpblks

      IF ( jrow == ngpblks ) THEN
        nproma = ldc%npromz
      ELSE
        nproma = ldc%nproma
      END IF

      DO jl = 1, nproma
        local_3d(jl,jk,jrow) = state_l(k+off_fields_p(4))
        k = k + 1
      END DO

    END DO
  END DO

  CALL gather_field(gl_d,global_dim,local_3d)
  DEALLOCATE(local_3d)

  ! Gather global specific humidity
  ALLOCATE (gl_q(nlon,nlev,nlat))
  ALLOCATE(local_3d(ldc%nproma,nlev,ngpblks))

! 3D specific humidity
  k = 1
  DO jk = nlev,1,-1
    DO jrow = 1, ngpblks

      IF ( jrow == ngpblks ) THEN
        nproma = ldc%npromz
      ELSE
        nproma = ldc%nproma
      END IF

      DO jl = 1, nproma
        local_3d(jl,jk,jrow) = state_l(k+off_fields_p(5))
        k = k + 1
      END DO

    END DO
  END DO

  CALL gather_field(gl_q,global_dim,local_3d)
  DEALLOCATE(local_3d)

  ! Gather global u
  ALLOCATE (gl_u(nlon,nlev,nlat))
  ALLOCATE(local_3d(ldc%nproma,nlev,ngpblks))

! 3D u
  k = 1
  DO jk = nlev,1,-1
    DO jrow = 1, ngpblks

      IF ( jrow == ngpblks ) THEN
        nproma = ldc%npromz
      ELSE
        nproma = ldc%nproma
      END IF

      DO jl = 1, nproma
        local_3d(jl,jk,jrow) = state_l(k+off_fields_p(6))
        k = k + 1
      END DO

    END DO
  END DO

  CALL gather_field(gl_u,global_dim,local_3d)
  DEALLOCATE(local_3d)

  ! Gather global v
  ALLOCATE (gl_v(nlon,nlev,nlat))
  ALLOCATE(local_3d(ldc%nproma,nlev,ngpblks))

! 3D v
  k = 1
  DO jk = nlev,1,-1
    DO jrow = 1, ngpblks

      IF ( jrow == ngpblks ) THEN
        nproma = ldc%npromz
      ELSE
        nproma = ldc%nproma
      END IF

      DO jl = 1, nproma
        local_3d(jl,jk,jrow) = state_l(k+off_fields_p(7))
        k = k + 1
      END DO

    END DO
  END DO

  CALL gather_field(gl_v,global_dim,local_3d)
  DEALLOCATE(local_3d)



! Reshape gl_t 
  ALLOCATE(g_t(nlon,nlat,nlev))
  DO i = 1, nlev
     g_t(:,:,i)=gl_t(:,i,:)
  END DO
  DEALLOCATE(gl_t)
!Reshape gl_vo
  ALLOCATE(g_vo(nlon,nlat,nlev))
  DO i = 1, nlev
     g_vo(:,:,i)=gl_vo(:,i,:)
  END DO
  DEALLOCATE(gl_vo)
!Reshape gl_d
  ALLOCATE(g_d(nlon,nlat,nlev))
  DO i = 1, nlev
     g_d(:,:,i)=gl_d(:,i,:)
  END DO
  DEALLOCATE(gl_d)
!Reshape gl_q
  ALLOCATE(g_q(nlon,nlat,nlev))
  DO i = 1, nlev
     g_q(:,:,i)=gl_q(:,i,:)
  END DO
  DEALLOCATE(gl_q)
!Reshape gl_u
  ALLOCATE(g_u(nlon,nlat,nlev))
  DO i = 1, nlev
     g_u(:,:,i)=gl_u(:,i,:)
  END DO
  DEALLOCATE(gl_u)
!Reshape gl_v
  ALLOCATE(g_v(nlon,nlat,nlev))
  DO i = 1, nlev
     g_v(:,:,i)=gl_v(:,i,:)
  END DO
  DEALLOCATE(gl_v)

  if_write:IF (writepe) THEN

! ----- Open Netcdf File
       stat(1) = NF_OPEN(filename, NF_WRITE, FileId)

       IF (stat(1) /= NF_NOERR) STOP 'nc-file error'

! ----- INQUIRE VARIABLE IDs
       s = 1
       stat(s) = NF_INQ_VARID(fileid, "time", varid_time)
       s = s + 1

       IF (writetype=='i') THEN
          stat(s) = NF_INQ_VARID(fileid, "lon", varid_lon)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "lat", varid_lat)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "lev", varid_lev)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "hyai", varid_hyai)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "hybi", varid_hybi)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "hyam", varid_hyam)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "hybm", varid_hybm)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "t_i",varid_t)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "vo_i",varid_vo)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "d_i",varid_d)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "q_i",varid_q)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "u_i",varid_u)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "v_i",varid_v)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "lsp_i",varid_lsp)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_t_ini",varid_rmst)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_vo_ini",varid_rmsvo)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_d_ini",varid_rmsd)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_q_ini",varid_rmsq)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_u_ini",varid_rmsu)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_v_ini",varid_rmsv)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_lsp_ini",varid_rmslsp)

       ELSE IF (writetype=='f') THEN
          stat(s) = NF_INQ_VARID(fileid, "t_f",varid_t)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "lsp_f",varid_lsp) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "vo_f",varid_vo)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "d_f",varid_d)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "q_f",varid_q)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "u_f",varid_u)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "v_f",varid_v)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_t_f",varid_rmst)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_vo_f",varid_rmsvo)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_d_f",varid_rmsd)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_q_f",varid_rmsq)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_u_f",varid_rmsu)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_v_f",varid_rmsv)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_lsp_f",varid_rmslsp)

       ELSE IF (writetype=='a') THEN
          stat(s) = NF_INQ_VARID(fileid, "t_a",varid_t)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "lsp_a",varid_lsp)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "vo_a",varid_vo)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "d_a",varid_d)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "q_a",varid_q)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "u_a",varid_u)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "v_a",varid_v)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_t_a",varid_rmst)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_vo_a",varid_rmsvo)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_d_a",varid_rmsd)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_q_a",varid_rmsq)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_u_a",varid_rmsu)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_v_a",varid_rmsv)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_lsp_a",varid_rmslsp)

       END IF
       DO i = 1, s
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error inquiring variable IDs, no.', i
       END DO

! ----- WRITE VARIABLES
       s = 1
       pos1 = write_pos_da
       nmb  = 1

       CALL get_date_components (current_date,year=current_year,month=current_month, day=current_day)
       stat(s) = NF_PUT_VARA_INT(fileid,varid_time, pos1, nmb, current_day)
       s = s + 1

       IF (writetype=='i') THEN
          ! WRITE lon & lat & lev
          stat(s) = NF_PUT_VAR_DOUBLE(fileid,varid_lon,gl_lon(:,1))
          s = s + 1
          stat(s) = NF_PUT_VAR_DOUBLE(fileid,varid_lat,gl_lat(1,:))
          s = s + 1
          ALLOCATE(tmp_i(nlev))
          DO i = 1, nlev
             tmp_i(i)=i
          END DO
          stat(s) = NF_PUT_VAR_INT(fileid,varid_lev,tmp_i)
          s = s + 1
          DEALLOCATE(tmp_i)
          stat(s) = NF_PUT_VAR_DOUBLE(fileid,varid_hyai,hyai)
          s = s + 1
          stat(s) = NF_PUT_VAR_DOUBLE(fileid,varid_hybi,hybi)
          s = s + 1
          stat(s) = NF_PUT_VAR_DOUBLE(fileid,varid_hyam,hyam)
          s = s + 1
          stat(s) = NF_PUT_VAR_DOUBLE(fileid,varid_hybm,hybm)
          s = s + 1
          stat(s) = NF_PUT_VAR_REAL(fileid, varid_t, REAL(g_t,4))
          s = s + 1
          stat(s) = NF_PUT_VAR_REAL(fileid, varid_vo, REAL(g_vo,4))
          s = s + 1
          stat(s) = NF_PUT_VAR_REAL(fileid, varid_d, REAL(g_d,4))
          s = s + 1
          stat(s) = NF_PUT_VAR_REAL(fileid, varid_q, REAL(g_q,4))
          s = s + 1
          stat(s) = NF_PUT_VAR_REAL(fileid, varid_u, REAL(g_u,4))
          s = s + 1
          stat(s) = NF_PUT_VAR_REAL(fileid, varid_v, REAL(g_v,4))
          s = s + 1
          stat(s) = NF_PUT_VAR_REAL(fileid, varid_rmst, REAL(rmse(1),4))
          s = s + 1
          stat(s) = NF_PUT_VAR_REAL(fileid, varid_rmslsp, REAL(rmse(2),4))
          s = s + 1
          stat(s) = NF_PUT_VAR_REAL(fileid, varid_rmsvo, REAL(rmse(3),4))
          s = s + 1
          stat(s) = NF_PUT_VAR_REAL(fileid, varid_rmsd, REAL(rmse(4),4))
          s = s + 1
          stat(s) = NF_PUT_VAR_REAL(fileid, varid_rmsq, REAL(rmse(5),4))
          s = s + 1
          stat(s) = NF_PUT_VAR_REAL(fileid, varid_rmsu, REAL(rmse(6),4))
          s = s + 1
          stat(s) = NF_PUT_VAR_REAL(fileid, varid_rmsv, REAL(rmse(7),4))
       END IF
   
       ! Write 3D fields
       pos1vec4 = (/1,1,1, write_pos_da/)
       nmbvec4  = (/nlon,nlat,nlev, 1 /)

       IF (writetype=='f' .OR. writetype=='a') THEN 
       
          ! write global temperature field
          stat(s) = NF_PUT_VARA_REAL(fileid, varid_t, pos1vec4, nmbvec4, REAL(g_t, 4))
          s = s + 1
          ! write global vorticity field
          stat(s) = NF_PUT_VARA_REAL(fileid, varid_vo, pos1vec4, nmbvec4,REAL(g_vo, 4))
          s = s + 1
          ! write global divergence field
          stat(s) = NF_PUT_VARA_REAL(fileid, varid_d, pos1vec4, nmbvec4,REAL(g_d, 4))
          s = s + 1
          ! write global specific humidity field
          stat(s) = NF_PUT_VARA_REAL(fileid, varid_q, pos1vec4, nmbvec4,REAL(g_q, 4))
          s = s + 1
          ! write global u field
          stat(s) = NF_PUT_VARA_REAL(fileid, varid_u, pos1vec4, nmbvec4,REAL(g_u, 4))
          s = s + 1
          ! write global v field
          stat(s) = NF_PUT_VARA_REAL(fileid, varid_v, pos1vec4, nmbvec4,REAL(g_v, 4))
          s = s + 1
       END IF

       ! Write 2D fields
       pos1vec3 = (/1,1,write_pos_da/)
       nmbvec3  = (/nlon,nlat,1/)

       IF (writetype=='f' .OR. writetype=='a') THEN

          ! write global logical surface pressure field
          stat(s) = NF_PUT_VARA_REAL(fileid, varid_lsp, pos1vec3, nmbvec3,REAL(gl_lsp, 4))

       END IF

       IF (writetype=='f' .OR. writetype=='a') THEN
          ! Write rms
          s = s + 1
          stat(s) = NF_PUT_VARA_REAL(fileid, varid_rmst, pos1, nmb, REAL(rmse(1),4))
          s = s + 1
          stat(s) = NF_PUT_VARA_REAL(fileid, varid_rmslsp, pos1, nmb, REAL(rmse(2),4))
          s = s + 1
          stat(s) = NF_PUT_VARA_REAL(fileid, varid_rmsvo, pos1, nmb, REAL(rmse(3),4))
          s = s + 1
          stat(s) = NF_PUT_VARA_REAL(fileid, varid_rmsd, pos1, nmb, REAL(rmse(4),4))
          s = s + 1
          stat(s) = NF_PUT_VARA_REAL(fileid, varid_rmsq, pos1, nmb, REAL(rmse(5),4))
          s = s + 1
          stat(s) = NF_PUT_VARA_REAL(fileid, varid_rmsu, pos1, nmb, REAL(rmse(6),4))
          s = s + 1
          stat(s) = NF_PUT_VARA_REAL(fileid, varid_rmsv, pos1, nmb, REAL(rmse(7),4))
       END IF
        
       DO i = 1, s
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in writing variables, no.', i
       END DO
! ----- CLOSE THE FILE
!       DEALLOCATE(gl_t)
       stat(1) = NF_CLOSE(fileid)

       IF (stat(1) /= NF_NOERR) THEN
          WRITE(*, *) 'NetCDF error in closing NetCDF file'
       END IF
    

  END IF if_write

  DEALLOCATE(g_t)
  DEALLOCATE(gl_lsp)
  DEALLOCATE(g_vo)
  DEALLOCATE(g_d)
  DEALLOCATE(g_q)
  DEALLOCATE(g_u)
  DEALLOCATE(g_v)

  END SUBROUTINE write_nc_airt_pdaf

END MODULE output_pdaf

