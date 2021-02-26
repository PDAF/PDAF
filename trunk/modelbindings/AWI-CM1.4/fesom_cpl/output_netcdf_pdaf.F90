!$Id: output_netcdf_pdaf.F90 2300 2020-05-17 05:40:30Z lnerger $
!> Module holding routines to output assimilation results into NetCDF
!!
!! This module provides routines to initialize NetCDF
!! output files for FESOM and to write output into the files.
!!
!! __Revision history:__
!! 2012-03 - Lars Nerger - Initial code based on output_netcdf module of pFEOM
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
  LOGICAL :: write_da = .TRUE.                 ! Whether to write output file from assimilation
  LOGICAL :: write_ens = .TRUE.                ! Whether to write output file for each individual ensemble member

! *** Private variables ***
  LOGICAL, PRIVATE :: debugoutput=.FALSE.   ! Write output for debugging
                                            ! (file contains only last writing)
  INTEGER, PRIVATE :: nf_prec               ! Precision of NetCDF output

CONTAINS

!----------------------------------------------------------------------------
!>  Routine to initialize NetCDF output file
!!
  SUBROUTINE init_output_pdaf(dim_lag, writepe)

    USE g_clock, &
         ONLY: yearnew, yearold
    USE mod_assim_pdaf, &
         ONLY: dim_ens
 
    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_lag           !< Smoother lag
    LOGICAL, INTENT(in) :: writepe           !< Whether this process writes


! *** Call initialization routines ***

    IF (yearnew /= yearold .AND. writepe .AND. write_da) THEN
       ! Initialize ocean file
       CALL init_ncfile_oce_pdaf(dim_lag, writepe)
#ifdef use_ice
       ! Initialize ice file
       CALL init_ncfile_ice_pdaf(dim_lag, writepe)
#endif
        IF (write_ens) THEN
           CALL init_ncfile_oce_pdaf_ens(dim_lag, writepe, dim_ens)
        END IF
    END IF

  END SUBROUTINE init_output_pdaf


!----------------------------------------------------------------------------
!>  Initialize NetCDF output file for ocean fields
!!
  SUBROUTINE init_ncfile_oce_pdaf(dim_lag, writepe)

    USE g_config, &
         ONLY: runid, ResultPath
    USE g_clock, &
         ONLY: cyearnew
    USE g_parfe, &
         ONLY: mydim_nod2d, mydim_nod3d, edim_nod2d, edim_nod3d
    USE o_mesh, &
         ONLY: nod2d, nod3d, coord_nod2d, coord_nod3d

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_lag           !< Smoother lag
    LOGICAL, INTENT(in) :: writepe           !< Whether this process writes
    
! *** Local variables ***
    CHARACTER(len=100) :: attstr            ! String to write attributes
    INTEGER :: i, n                         ! Counters
    INTEGER :: fileid                       ! ID of netCDF file
    INTEGER :: DimId_iter                   ! dimension: iteration
    INTEGER :: DimId_n2D, DimID_n3D         ! dimension: nodes
    INTEGER :: Dim1                         ! dimensions
    INTEGER :: VarId_time, VarId_iter       ! variables: time, iteration 
    INTEGER :: VarId_asmlstep               ! Number of assimilation step
    INTEGER :: VarId_effdimobs, VarId_locradius      ! Variables: observation dim and loc. radius
    INTEGER :: VarId_ssha, VarId_sshf, VarId_sshi    ! variable: sea surface height
    INTEGER :: VarId_tempa, VarId_tempf, VarId_tempi ! variable: temperature
    INTEGER :: VarId_salta, VarId_saltf, VarId_salti ! variable: salt
    INTEGER :: VarId_ua, VarId_uf, VarId_ui          ! variable: u-velocity in nodes
    INTEGER :: VarId_va, VarId_vf, VarId_vi          ! variable: v-velocity in nodes
    INTEGER :: VarId_wa, VarId_wf, VarId_wi          ! variable: w-velocity in nodes
    INTEGER :: VarID_sshs, VarID_temps, VarID_salts  ! smoothed variables
    INTEGER :: VarID_us, VarID_vs, VarID_ws          ! smoother variables
    INTEGER :: VarID_rmsssha, VARID_rmssshf          ! variable: RMS errors SSH
    INTEGER :: VarId_rmsua, VarId_rmsuf              ! variable: RMS errors u
    INTEGER :: VarId_rmsva, VarId_rmsvf              ! variable: RMS errors v
    INTEGER :: VarId_rmswa, VarId_rmswf              ! variable: RMS errors w
    INTEGER :: VarId_rmstempa, VarId_rmstempf        ! variable: RMS errors temperature
    INTEGER :: VarId_rmssalta, VarId_rmssaltf        ! variable: RMS errors salinity
    INTEGER :: VarId_rmssshini, VarId_rmswini        ! RMS for initial fields
    INTEGER :: VarId_rmsuini, VarId_rmsvini          ! RMS for initial fields
    INTEGER :: VarId_rmstempini, VarId_rmssaltini    ! RMS for initial fields
    INTEGER :: VarId_rmssshs, VarId_rmsws            ! RMS for smoothed fields
    INTEGER :: VarId_rmsus, VarId_rmsvs              ! RMS for smoothed fields
    INTEGER :: VarId_rmstemps, VarId_rmssalts        ! RMS for smoothed fields
    INTEGER :: s       			    ! auxiliary: status counter
    INTEGER :: dimarray(3)                  ! auxiliary: array dimension
    INTEGER :: stat(500)                    ! auxiliary: status array
    CHARACTER(200)            :: filename   ! Full name of output file

    
    pe0: IF (writepe) THEN

! Print screen information
       IF (prec_nc == 's') THEN
          WRITE (*, '(/a, 1x, a)') 'FESOM-PDAF', 'Initialize assimilation NetCDF ocean file - single precision'
          nf_prec = NF_FLOAT
       ELSE
          WRITE (*, '(/a, 1x, a)') 'FESOM-PDAF', 'Initialize assimilation NetCDF ocean file - double precision'
          nf_prec = NF_DOUBLE
       END IF

! ----- open file and write global attributes

       filename=TRIM(ResultPath)//runid//'.'//cyearnew//'.oce.'//TRIM(str_daspec)//'.nc'

       s = 1
    
!       stat(s) = NF_CREATE(filename, 0, fileid) 
       stat(s) = NF_CREATE(filename, NF_NETCDF4, fileid) 
       s = s + 1

       attstr  = 'FESOM Assimilation'
       stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'title', LEN_TRIM(attstr), &
            TRIM(attstr)) 
       s = s + 1


! ----- DEFINE DIMENSIONS ------------------------------

       stat(s) = NF_DEF_DIM(fileid, 'nodes_2D',    nod2D, DimId_n2D)             
       s = s + 1
       stat(s) = NF_DEF_DIM(fileid, 'nodes_3D',    nod3D, DimId_n3D)             
       s = s + 1
       stat(s) = NF_DEF_DIM(fileid, 'one',         1, dim1)
       s = s + 1
       stat(s) = NF_DEF_DIM(fileid, 'iteration',   NF_UNLIMITED, DimId_iter)
       s = s + 1
       
       
       DO i = 1,  s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in ocean dimension definitions, no.', i
       END DO

! ----- DEFINE VARIABLES ---------------------------------

       s = 1

       !- numbers

       stat(s) = NF_DEF_VAR(fileid, 'iter', NF_INT, 1, DimId_iter, VarId_iter) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'time',   nf_prec, 1, DimId_iter, VarId_time) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'asmlstep', NF_INT, 1, DimId_iter, VarId_asmlstep) 
       s = s + 1


! ----- F I E L D S 

       !- scalar variables

       dimarray(1) = DimId_n2D
       dimarray(2) = DimId_iter

       stat(s) = NF_DEF_VAR(fileid, 'ssh_a', nf_prec, 2, dimarray(1:2), VarId_ssha); 
       s = s + 1
       stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_ssha, 0, 1, 5) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'ssh_f', nf_prec, 2, dimarray(1:2), VarId_sshf); 
       s = s + 1
       stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_sshf, 0, 1, 5) 
       s = s + 1
       
       dimarray(2) = dim1
       stat(s) = NF_DEF_VAR(fileid, 'ssh_ini', nf_prec, 2, dimarray(1:2), VarId_sshi); 
       s = s + 1
       stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_sshi, 0, 1, 5) 
       s = s + 1

       dimarray(1) = DimId_n3D
       dimarray(2) = DimId_iter

       stat(s) = NF_DEF_VAR(fileid, 'temp_a', nf_prec, 2, dimarray(1:2), VarId_tempa)
       s = s + 1
       stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_tempa, 0, 1, 5) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'temp_f', nf_prec, 2, dimarray(1:2), VarId_tempf)
       s = s + 1
       stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_tempf, 0, 1, 5) 
       s = s + 1

       dimarray(2) = dim1
       stat(s) = NF_DEF_VAR(fileid, 'temp_ini', nf_prec, 2, dimarray(1:2), VarId_tempi)
       s = s + 1
       stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_tempi, 0, 1, 5) 
       s = s + 1

       dimarray(2) = DimId_iter
       stat(s) = NF_DEF_VAR(fileid, 'salt_a', nf_prec, 2, dimarray(1:2), VarId_salta)
       s = s + 1
       stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_salta, 0, 1, 5) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'salt_f', nf_prec, 2, dimarray(1:2), VarId_saltf)
       s = s + 1
       stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_saltf, 0, 1, 5) 
       s = s + 1

       dimarray(2) = dim1
       stat(s) = NF_DEF_VAR(fileid, 'salt_ini', nf_prec, 2, dimarray(1:2), VarId_salti)
       s = s + 1
       stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_salti, 0, 1, 5) 
       s = s + 1

       !- vector variables

       !- velocity in nodes

       dimarray(1) = DimId_n3D
       dimarray(2) = DimId_iter

       stat(s) = NF_DEF_VAR(fileid, 'u_a', nf_prec, 2, dimarray(1:2), VarId_ua)
       s = s + 1
       stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_ua, 0, 1, 5) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'u_f', nf_prec, 2, dimarray(1:2), VarId_uf)
       s = s + 1
       stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_uf, 0, 1, 5) 
       s = s + 1

       dimarray(2) = dim1
       stat(s) = NF_DEF_VAR(fileid, 'u_ini', nf_prec, 2, dimarray(1:2), VarId_ui)
       s = s + 1
       stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_ui, 0, 1, 5) 
       s = s + 1

       dimarray(2) = DimId_iter
       stat(s) = NF_DEF_VAR(fileid, 'v_a', nf_prec, 2, dimarray(1:2), VarId_va)
       s = s + 1
       stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_va, 0, 1, 5) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'v_f', nf_prec, 2, dimarray(1:2), VarId_vf)
       s = s + 1
       stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_vf, 0, 1, 5) 
       s = s + 1

       dimarray(2) = dim1
       stat(s) = NF_DEF_VAR(fileid, 'v_ini', nf_prec, 2, dimarray(1:2), VarId_vi)
       s = s + 1
       stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_vi, 0, 1, 5) 
       s = s + 1

       dimarray(2) = DimId_iter
       stat(s) = NF_DEF_VAR(fileid, 'wpot_a', nf_prec, 2, dimarray(1:2), VarId_wa)
       s = s + 1
       stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_wa, 0, 1, 5) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'wpot_f', nf_prec, 2, dimarray(1:2), VarId_wf)
       s = s + 1
       stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_wf, 0, 1, 5) 
       s = s + 1

       dimarray(2) = dim1
       stat(s) = NF_DEF_VAR(fileid, 'wpot_ini', nf_prec, 2, dimarray(1:2), VarId_wi)
       s = s + 1
       stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_wi, 0, 1, 5) 
       s = s + 1


! ----- Smoother F I E L D S 

       smootherA: IF (dim_lag > 0) THEN
          !- scalar variables

          dimarray(1) = DimId_n2D
          dimarray(2) = DimId_iter
          stat(s) = NF_DEF_VAR(fileid, 'ssh_s', nf_prec, 2, dimarray(1:2), VarId_sshs); 
          s = s + 1
          stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_sshs, 0, 1, 5) 
          s = s + 1

          dimarray(1) = DimId_n3D
          dimarray(2) = DimId_iter
          stat(s) = NF_DEF_VAR(fileid, 'temp_s', nf_prec, 2, dimarray(1:2), VarId_temps)
          s = s + 1
          stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_temps, 0, 1, 5) 
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'salt_s', nf_prec, 2, dimarray(1:2), VarId_salts)
          s = s + 1
          stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_salts, 0, 1, 5) 
          s = s + 1

          !- vector variables

          dimarray(1) = DimId_n3D
          dimarray(2) = DimId_iter

          stat(s) = NF_DEF_VAR(fileid, 'u_s', nf_prec, 2, dimarray(1:2), VarId_us)
          s = s + 1
          stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_us, 0, 1, 5) 
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'v_s', nf_prec, 2, dimarray(1:2), VarId_vs)
          s = s + 1
          stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_vs, 0, 1, 5) 
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'wpot_s', nf_prec, 2, dimarray(1:2), VarId_ws)
          s = s + 1
          stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_ws, 0, 1, 5) 
          s = s + 1
       END IF smootherA

! ----- Effective observation dimensions and localization radii

       dimarray(1) = DimId_n2D
       dimarray(2) = DimId_iter
       stat(s) = NF_DEF_VAR(fileid, 'eff_dim_obs', nf_prec, 2, dimarray(1:2), VarId_effdimobs); 
       s = s + 1
       stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_effdimobs, 0, 1, 5) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'loc_radius', nf_prec, 2, dimarray(1:2), VarId_locradius); 
       s = s + 1
       stat(s) = NF_DEF_VAR_DEFLATE(fileid, VarId_locradius, 0, 1, 5) 
       s = s + 1

       DO i = 1,  s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in ocean variable definition, no.', i
       END DO


       !- RMS errors

       s = 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_ssh_ini', nf_prec, 1, dim1, VarId_rmssshini) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_u_ini', nf_prec, 1, dim1, VarId_rmsuini) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_v_ini', nf_prec, 1, dim1, VarId_rmsvini) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_wpot_ini', nf_prec, 1, dim1, VarId_rmswini) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_temp_ini', nf_prec, 1, dim1, VarId_rmstempini) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_salt_ini', nf_prec, 1, dim1, VarId_rmssaltini) 
       s = s + 1

       stat(s) = NF_DEF_VAR(fileid, 'rms_ssh_a', nf_prec, 1, DimId_iter, VarId_rmsssha) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_u_a', nf_prec, 1, DimId_iter, VarId_rmsua) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_v_a', nf_prec, 1, DimId_iter, VarId_rmsva) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_wpot_a', nf_prec, 1, DimId_iter, VarId_rmswa) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_temp_a', nf_prec, 1, DimId_iter, VarId_rmstempa) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_salt_a', nf_prec, 1, DimId_iter, VarId_rmssalta) 
       s = s + 1

       stat(s) = NF_DEF_VAR(fileid, 'rms_ssh_f', nf_prec, 1, DimId_iter, VarId_rmssshf) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_u_f', nf_prec, 1, DimId_iter, VarId_rmsuf) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_v_f', nf_prec, 1, DimId_iter, VarId_rmsvf) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_wpot_f', nf_prec, 1, DimId_iter, VarId_rmswf) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_temp_f', nf_prec, 1, DimId_iter, VarId_rmstempf) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_salt_f', nf_prec, 1, DimId_iter, VarId_rmssaltf) 
       s = s + 1

       smootherB: IF (dim_lag > 0) THEN
          stat(s) = NF_DEF_VAR(fileid, 'rms_ssh_s', nf_prec, 1, DimId_iter, VarId_rmssshs)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_u_s', nf_prec, 1, DimId_iter, VarId_rmsus)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_v_s', nf_prec, 1, DimId_iter, VarId_rmsvs)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_wpot_s', nf_prec, 1, DimId_iter, VarId_rmsws)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_temp_s', nf_prec, 1, DimId_iter, VarId_rmstemps)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_salt_s', nf_prec, 1, DimId_iter, VarId_rmssalts)
          s = s + 1
       END IF smootherB

       DO i = 1,  s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in rms error variable definition, no.', i
       END DO

! ----- DEFINE ATTRIBUTES -----------------------------------------

       s = 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_time,     'long_name',    4, 'time'); 
       s = s + 1  
       attstr  = 'seconds'
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_time,     'units',        LEN_TRIM(attstr), attstr) 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_iter,     'long_name',     9, 'iteration') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_asmlstep, 'long_name',    17, 'Assimilation step') 
       s = s + 1

       ! Fields
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sshf,      'long_name',   32, 'sea surface elevation - forecast') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sshf,      'units',        5, 'meter') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sshf,      'field',       19, 'ssh, scalar, series') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sshf,      'connections', 20, 'triangles, triangles') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sshf,      'positions',   17, 'surface_locations') 
       s = s + 1

       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ssha,      'long_name',   32, 'sea surface elevation - analysis') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ssha,      'units',        5, 'meter') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ssha,      'field',       19, 'ssh, scalar, series') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ssha,      'connections', 20, 'triangles, triangles') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ssha,      'positions',   17, 'surface_locations') 
       s = s + 1

       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sshi,      'long_name',   31, 'sea surface elevation - initial') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sshi,      'units',        5, 'meter') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sshi,      'field',       19, 'ssh, scalar, series') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sshi,      'connections', 20, 'triangles, triangles') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sshi,      'positions',   17, 'surface_locations') 
       s = s + 1
       
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tempf,      'long_name',   22, 'temperature - forecast') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tempf,      'units',       15, 'degrees Celcius') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tempf,      'field',       27, 'temperature, scalar, series') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tempf,      'connections', 20, 'triangles, triangles') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tempf,      'positions',   17, 'surface_locations') 
       s = s + 1
       
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tempa,      'long_name',   22, 'temperature - analysis') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tempa,      'units',       15, 'degrees Celcius') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tempa,      'field',       27, 'temperature, scalar, series') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tempa,      'connections', 20, 'triangles, triangles') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tempa,      'positions',   17, 'surface_locations') 
       s = s + 1
       
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tempi,      'long_name',   21, 'temperature - initial') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tempi,      'units',       15, 'degrees Celcius') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tempi,      'field',       27, 'temperature, scalar, series') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tempi,      'connections', 20, 'triangles, triangles') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tempi,      'positions',   17, 'surface_locations') 
       s = s + 1
       
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_saltf,      'long_name',   19, 'salinity - forecast') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_saltf,      'units',        3, 'psu') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_saltf,      'field',       24, 'salinity, scalar, series') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_saltf,      'connections', 20, 'triangles, triangles') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_saltf,      'positions',   17, 'surface_locations') 
       s = s + 1
       
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_salta,      'long_name',   19, 'salinity - analysis') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_salta,      'units',        3, 'psu') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_salta,      'field',       24, 'salinity, scalar, series') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_salta,      'connections', 20, 'triangles, triangles') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_salta,      'positions',   17, 'surface_locations') 
       s = s + 1
       
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_salti,      'long_name',   18, 'salinity - initial') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_salti,      'units',        3, 'psu') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_salti,      'field',       24, 'salinity, scalar, series') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_salti,      'connections', 20, 'triangles, triangles') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_salti,      'positions',   17, 'surface_locations') 
       s = s + 1       
       
       DO i = 1, s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in ocean file attribute assignments, no.', i
       END DO

       s = 1
       
       stat(s) = NF_ENDDEF(fileid) 
       s = s + 1
       stat(1) = NF_CLOSE(fileid)

       IF (stat(1) /= NF_NOERR) THEN
          WRITE(*, *) 'NetCDF error in closing ocean NetCDF file'
       END IF
    END IF pe0

  END SUBROUTINE init_ncfile_oce_pdaf


!----------------------------------------------------------------------------
!> Initialize NetCDF output file for ocean fields for each ensemble member
!!
  SUBROUTINE init_ncfile_oce_pdaf_ens(dim_lag, writepe, dim_ens)

    USE g_config, &
         ONLY: runid, ResultPath
    USE g_clock, &
         ONLY: cyearnew
    USE g_parfe, &
         ONLY: mydim_nod2d, mydim_nod3d, edim_nod2d, edim_nod3d
    USE o_mesh, &
         ONLY: nod2d, nod3d, coord_nod2d, coord_nod3d

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_lag          !< Smoother lag
    LOGICAL, INTENT(in) :: writepe          !< Whether this process writes
    INTEGER, INTENT(IN) :: dim_ens          !< Ensemble size

! *** Local variables ***
    CHARACTER(len=100) :: attstr            ! String to write attributes
    INTEGER :: i, n, member                 ! Counters
    CHARACTER(len=10) :: member_string
    INTEGER :: fileid                       ! ID of netCDF file
    INTEGER :: DimId_iter                   ! dimension: iteration
    INTEGER :: DimId_n2D, DimID_n3D         ! dimension: nodes
    INTEGER :: Dim1                         ! dimensions
    INTEGER :: VarId_time, VarId_iter       ! variables: time, iteration 
    INTEGER :: VarId_asmlstep               ! Number of assimilation step
    INTEGER :: VarId_effdimobs, VarId_locradius      ! Variables: observation dim and loc. radius
    INTEGER :: VarId_ssha, VarId_sshf, VarId_sshi    ! variable: sea surface height
    INTEGER :: VarId_tempa, VarId_tempf, VarId_tempi ! variable: temperature
    INTEGER :: VarId_salta, VarId_saltf, VarId_salti ! variable: salt
    INTEGER :: VarId_ua, VarId_uf, VarId_ui          ! variable: u-velocity in nodes
    INTEGER :: VarId_va, VarId_vf, VarId_vi          ! variable: v-velocity in nodes
    INTEGER :: VarId_wa, VarId_wf, VarId_wi          ! variable: w-velocity in nodes
    INTEGER :: VarID_sshs, VarID_temps, VarID_salts  ! smoothed variables
    INTEGER :: VarID_us, VarID_vs, VarID_ws          ! smoother variables
    INTEGER :: VarID_rmsssha, VARID_rmssshf          ! variable: RMS errors SSH
    INTEGER :: VarId_rmsua, VarId_rmsuf              ! variable: RMS errors u
    INTEGER :: VarId_rmsva, VarId_rmsvf              ! variable: RMS errors v
    INTEGER :: VarId_rmswa, VarId_rmswf              ! variable: RMS errors w
    INTEGER :: VarId_rmstempa, VarId_rmstempf        ! variable: RMS errors temperature
    INTEGER :: VarId_rmssalta, VarId_rmssaltf        ! variable: RMS errors salinity
    INTEGER :: VarId_rmssshini, VarId_rmswini        ! RMS for initial fields
    INTEGER :: VarId_rmsuini, VarId_rmsvini          ! RMS for initial fields
    INTEGER :: VarId_rmstempini, VarId_rmssaltini    ! RMS for initial fields
    INTEGER :: VarId_rmssshs, VarId_rmsws            ! RMS for smoothed fields
    INTEGER :: VarId_rmsus, VarId_rmsvs              ! RMS for smoothed fields
    INTEGER :: VarId_rmstemps, VarId_rmssalts        ! RMS for smoothed fields
    INTEGER :: s                            ! auxiliary: status counter
    INTEGER :: dimarray(3)                  ! auxiliary: array dimension
    INTEGER :: stat(500)                    ! auxiliary: status array
    CHARACTER(200)            :: filename   ! Full name of output file


    pe0: IF (writepe) THEN
! Print screen information
       IF (prec_nc == 's') THEN
          WRITE (*, '(/a, 1x, a)') 'FESOM-PDAF', 'Initialize assimilation NetCDF ocean file for each ensemble member - single precision'
          nf_prec = NF_FLOAT
       ELSE
          WRITE (*, '(/a, 1x, a)') 'FESOM-PDAF', 'Initialize assimilation NetCDF ocean file for each ensemble member - double precision'
          nf_prec = NF_DOUBLE
       END IF

! ----- open file and write global attributes
       DO member = 1, dim_ens
          WRITE(member_string,'(i4.4)') member
          filename=TRIM(ResultPath)//runid//'.'//cyearnew//'.oce.'//TRIM(str_daspec)//'_'//TRIM(member_string)//'.nc'

          s = 1

          stat(s) = NF_CREATE(filename, 0, fileid)
          s = s + 1

          attstr  = 'FESOM Assimilation'
          stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'title', LEN_TRIM(attstr), &
            TRIM(attstr))
          s = s + 1


! ----- DEFINE DIMENSIONS ------------------------------

          stat(s) = NF_DEF_DIM(fileid, 'nodes_2D',    nod2D, DimId_n2D)
          s = s + 1
          stat(s) = NF_DEF_DIM(fileid, 'nodes_3D',    nod3D, DimId_n3D)
          s = s + 1
          stat(s) = NF_DEF_DIM(fileid, 'one',         1, dim1)
          s = s + 1
          stat(s) = NF_DEF_DIM(fileid, 'iteration',   NF_UNLIMITED, DimId_iter)
          s = s + 1


          DO i = 1,  s - 1
             IF (stat(i) /= NF_NOERR) &
                  WRITE(*, *) 'NetCDF error in ocean dimension definitions, no.', i
          END DO

! ----- DEFINE VARIABLES ---------------------------------

          s = 1

          !- numbers

          stat(s) = NF_DEF_VAR(fileid, 'iter', NF_INT, 1, DimId_iter, VarId_iter)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'time',   nf_prec, 1, DimId_iter, VarId_time)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'asmlstep', NF_INT, 1, DimId_iter, VarId_asmlstep)
          s = s + 1

! ----- F I E L D S 

          !- scalar variables

          dimarray(1) = DimId_n2D
          dimarray(2) = DimId_iter

          stat(s) = NF_DEF_VAR(fileid, 'ssh_a', nf_prec, 2, dimarray(1:2), VarId_ssha);
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'ssh_f', nf_prec, 2, dimarray(1:2), VarId_sshf);
          s = s + 1

          dimarray(2) = dim1
          stat(s) = NF_DEF_VAR(fileid, 'ssh_ini', nf_prec, 2, dimarray(1:2), VarId_sshi);
          s = s + 1

          dimarray(1) = DimId_n3D
          dimarray(2) = DimId_iter

          stat(s) = NF_DEF_VAR(fileid, 'temp_a', nf_prec, 2, dimarray(1:2), VarId_tempa)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'temp_f', nf_prec, 2, dimarray(1:2), VarId_tempf)
          s = s + 1

          dimarray(2) = dim1
          stat(s) = NF_DEF_VAR(fileid, 'temp_ini', nf_prec, 2, dimarray(1:2), VarId_tempi)
          s = s + 1

          dimarray(2) = DimId_iter
          stat(s) = NF_DEF_VAR(fileid, 'salt_a', nf_prec, 2, dimarray(1:2), VarId_salta)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'salt_f', nf_prec, 2, dimarray(1:2), VarId_saltf)
          s = s + 1

          dimarray(2) = dim1
          stat(s) = NF_DEF_VAR(fileid, 'salt_ini', nf_prec, 2, dimarray(1:2), VarId_salti)
          s = s + 1

          !- vector variables

          !- velocity in nodes

          dimarray(1) = DimId_n3D
          dimarray(2) = DimId_iter

          stat(s) = NF_DEF_VAR(fileid, 'u_a', nf_prec, 2, dimarray(1:2), VarId_ua)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'u_f', nf_prec, 2, dimarray(1:2), VarId_uf)
          s = s + 1

          dimarray(2) = dim1
          stat(s) = NF_DEF_VAR(fileid, 'u_ini', nf_prec, 2, dimarray(1:2), VarId_ui)
          s = s + 1

          dimarray(2) = DimId_iter
          stat(s) = NF_DEF_VAR(fileid, 'v_a', nf_prec, 2, dimarray(1:2), VarId_va)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'v_f', nf_prec, 2, dimarray(1:2), VarId_vf)
          s = s + 1

          dimarray(2) = dim1
          stat(s) = NF_DEF_VAR(fileid, 'v_ini', nf_prec, 2, dimarray(1:2), VarId_vi)
          s = s + 1

          dimarray(2) = DimId_iter
          stat(s) = NF_DEF_VAR(fileid, 'wpot_a', nf_prec, 2, dimarray(1:2), VarId_wa)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'wpot_f', nf_prec, 2, dimarray(1:2), VarId_wf)
          s = s + 1

          dimarray(2) = dim1
          stat(s) = NF_DEF_VAR(fileid, 'wpot_ini', nf_prec, 2, dimarray(1:2), VarId_wi)
          s = s + 1


! ----- Smoother F I E L D S 

          smootherA: IF (dim_lag > 0) THEN
             !- scalar variables

             dimarray(1) = DimId_n2D
             dimarray(2) = DimId_iter
             stat(s) = NF_DEF_VAR(fileid, 'ssh_s', nf_prec, 2, dimarray(1:2), VarId_sshs);
             s = s + 1

             dimarray(1) = DimId_n3D
             dimarray(2) = DimId_iter
             stat(s) = NF_DEF_VAR(fileid, 'temp_s', nf_prec, 2, dimarray(1:2), VarId_temps)
             s = s + 1
             stat(s) = NF_DEF_VAR(fileid, 'salt_s', nf_prec, 2, dimarray(1:2), VarId_salts)
             s = s + 1

             !- vector variables

             dimarray(1) = DimId_n3D
             dimarray(2) = DimId_iter

             stat(s) = NF_DEF_VAR(fileid, 'u_s', nf_prec, 2, dimarray(1:2), VarId_us)
             s = s + 1
             stat(s) = NF_DEF_VAR(fileid, 'v_s', nf_prec, 2, dimarray(1:2), VarId_vs)
             s = s + 1
             stat(s) = NF_DEF_VAR(fileid, 'wpot_s', nf_prec, 2, dimarray(1:2), VarId_ws)
             s = s + 1
          END IF smootherA

! ----- Effective observation dimensions and localization radii

          dimarray(1) = DimId_n2D
          dimarray(2) = DimId_iter
          stat(s) = NF_DEF_VAR(fileid, 'eff_dim_obs', nf_prec, 2, dimarray(1:2), VarId_effdimobs);
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'loc_radius', nf_prec, 2, dimarray(1:2), VarId_locradius);
          s = s + 1

          DO i = 1,  s - 1
             IF (stat(i) /= NF_NOERR) &
                  WRITE(*, *) 'NetCDF error in ocean variable definition, no.', i
          END DO


          !- RMS errors

          s = 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_ssh_ini', nf_prec, 1, dim1, VarId_rmssshini)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_u_ini', nf_prec, 1, dim1, VarId_rmsuini)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_v_ini', nf_prec, 1, dim1, VarId_rmsvini)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_wpot_ini', nf_prec, 1, dim1, VarId_rmswini)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_temp_ini', nf_prec, 1, dim1, VarId_rmstempini)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_salt_ini', nf_prec, 1, dim1, VarId_rmssaltini)
          s = s + 1

          stat(s) = NF_DEF_VAR(fileid, 'rms_ssh_a', nf_prec, 1, DimId_iter, VarId_rmsssha)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_u_a', nf_prec, 1, DimId_iter, VarId_rmsua)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_v_a', nf_prec, 1, DimId_iter, VarId_rmsva)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_wpot_a', nf_prec, 1, DimId_iter, VarId_rmswa)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_temp_a', nf_prec, 1, DimId_iter, VarId_rmstempa)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_salt_a', nf_prec, 1, DimId_iter, VarId_rmssalta)
          s = s + 1

          stat(s) = NF_DEF_VAR(fileid, 'rms_ssh_f', nf_prec, 1, DimId_iter, VarId_rmssshf)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_u_f', nf_prec, 1, DimId_iter, VarId_rmsuf)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_v_f', nf_prec, 1, DimId_iter, VarId_rmsvf)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_wpot_f', nf_prec, 1, DimId_iter, VarId_rmswf)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_temp_f', nf_prec, 1, DimId_iter, VarId_rmstempf)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_salt_f', nf_prec, 1, DimId_iter, VarId_rmssaltf)
          s = s + 1

          smootherB: IF (dim_lag > 0) THEN
             stat(s) = NF_DEF_VAR(fileid, 'rms_ssh_s', nf_prec, 1, DimId_iter, VarId_rmssshs)
             s = s + 1
             stat(s) = NF_DEF_VAR(fileid, 'rms_u_s', nf_prec, 1, DimId_iter, VarId_rmsus)
             s = s + 1
             stat(s) = NF_DEF_VAR(fileid, 'rms_v_s', nf_prec, 1, DimId_iter, VarId_rmsvs)
             s = s + 1
             stat(s) = NF_DEF_VAR(fileid, 'rms_wpot_s', nf_prec, 1, DimId_iter, VarId_rmsws)
             s = s + 1
             stat(s) = NF_DEF_VAR(fileid, 'rms_temp_s', nf_prec, 1, DimId_iter, VarId_rmstemps)
             s = s + 1
             stat(s) = NF_DEF_VAR(fileid, 'rms_salt_s', nf_prec, 1, DimId_iter, VarId_rmssalts)
             s = s + 1
          END IF smootherB

          DO i = 1,  s - 1
             IF (stat(i) /= NF_NOERR) &
                  WRITE(*, *) 'NetCDF error in rms error variable definition, no.', i
          END DO

! ----- DEFINE ATTRIBUTES -----------------------------------------

          s = 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_time,     'long_name',    4, 'time');
          s = s + 1
          attstr  = 'seconds'
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_time,     'units',        LEN_TRIM(attstr), attstr)
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_iter,     'long_name',     9, 'iteration')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_asmlstep, 'long_name',    17, 'Assimilation step')
          s = s + 1

          ! Fields
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sshf,      'long_name',   32, 'sea surface elevation - forecast')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sshf,      'units',        5, 'meter')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sshf,      'field',       19, 'ssh, scalar, series')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sshf,      'connections', 20, 'triangles, triangles')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sshf,      'positions',   17, 'surface_locations')
          s = s + 1

          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ssha,      'long_name',   32, 'sea surface elevation - analysis')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ssha,      'units',        5, 'meter')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ssha,      'field',       19, 'ssh, scalar, series')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ssha,      'connections', 20, 'triangles, triangles')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ssha,      'positions',   17, 'surface_locations')
          s = s + 1

          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sshi,      'long_name',   31, 'sea surface elevation - initial')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sshi,      'units',        5, 'meter')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sshi,      'field',       19, 'ssh, scalar, series')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sshi,      'connections', 20, 'triangles, triangles')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sshi,      'positions',   17, 'surface_locations')
          s = s + 1

          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tempf,      'long_name',   22, 'temperature - forecast')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tempf,      'units',       15, 'degrees Celcius')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tempf,      'field',       27, 'temperature, scalar, series')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tempf,      'connections', 20, 'triangles, triangles')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tempf,      'positions',   17, 'surface_locations')
          s = s + 1

          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_saltf,      'long_name',   19, 'salinity - forecast')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_saltf,      'units',        3, 'psu')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_saltf,      'field',       24, 'salinity, scalar, series')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_saltf,      'connections', 20, 'triangles, triangles')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_saltf,      'positions',   17, 'surface_locations')
          s = s + 1

          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_salta,      'long_name',   19, 'salinity - analysis')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_salta,      'units',        3, 'psu')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_salta,      'field',       24, 'salinity, scalar, series')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_salta,      'connections', 20, 'triangles, triangles')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_salta,      'positions',   17, 'surface_locations')
          s = s + 1

          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_salti,      'long_name',   18, 'salinity - initial')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_salti,      'units',        3, 'psu')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_salti,      'field',       24, 'salinity, scalar, series')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_salti,      'connections', 20, 'triangles, triangles')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_salti,      'positions',   17, 'surface_locations')
          s = s + 1

          DO i = 1, s - 1
             IF (stat(i) /= NF_NOERR) &
                  WRITE(*, *) 'NetCDF error in ocean file attribute assignments, no.', i
          END DO

          s = 1
   
          stat(s) = NF_ENDDEF(fileid)
          s = s + 1
          stat(1) = NF_CLOSE(fileid)

          IF (stat(1) /= NF_NOERR) THEN
             WRITE(*, *) 'NetCDF error in closing ocean NetCDF file'
          END IF
       END DO
    END IF pe0

  END SUBROUTINE init_ncfile_oce_pdaf_ens


!----------------------------------------------------------------------------
!>  Initialize NetCDF output file for ice fields
!!
  SUBROUTINE init_ncfile_ice_pdaf(dim_lag, writepe)

    USE g_config, &
         ONLY: runid, ResultPath
    USE g_clock, &
         ONLY: cyearnew
    USE g_parfe, &
         ONLY: mydim_nod2d, mydim_nod3d, edim_nod2d, edim_nod3d
    USE o_mesh, &
         ONLY: nod2d, nod3d, coord_nod2d, coord_nod3d

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_lag           !< Smoother lag
    LOGICAL, INTENT(in) :: writepe           !< Whether this process writes
    
! *** Local variables ***
    CHARACTER(len=100) :: attstr            ! String to write attributes
    INTEGER :: i, n                         ! Counters
    INTEGER :: fileid                       ! ID of netCDF file
    INTEGER :: DimId_iter                   ! dimension: iteration
    INTEGER :: DimId_n2D                    ! dimension: nodes
    INTEGER :: Dim1                         ! dimensions
    INTEGER :: VarId_time, VarId_iter       ! variables: time, iteration 
    INTEGER :: VarId_asmlstep               ! Number of assimilation step
    INTEGER :: VarId_area_a, VarId_area_f, VarId_area_i    ! variable: ice concentration
    INTEGER :: VarId_hice_a, VarId_hice_f, VarId_hice_i    ! variable: ice thickness
    INTEGER :: VarId_hsnow_a, VarId_hsnow_f, VarId_hsnow_i ! variable: snow thickness
    INTEGER :: VarId_uice_a, VarId_uice_f, VarId_uice_i    ! variable: u-velocity in nodes
    INTEGER :: VarId_vice_a, VarId_vice_f, VarId_vice_i    ! variable: v-velocity in nodes
    INTEGER :: VarID_area_s, VarID_hice_s, VarID_hsnow_s   ! Smoothed ice fields
    INTEGER :: VarID_uice_s, VarID_vice_s                  ! Smoothed ice fields
    INTEGER :: VarId_rmsareaa, VarId_rmsareaf              ! variable: rms error ice conc.
    INTEGER :: VarId_rmshicea, VarId_rmshicef              ! variable: rms error ice thickness
    INTEGER :: VarId_rmshsnowa, VarId_rmshsnowf            ! variable: rms error snow thickness
    INTEGER :: VarId_rmsuicea, VarId_rmsuicef              ! variable: rms error ice u
    INTEGER :: VarId_rmsvicea, VarId_rmsvicef              ! variable: rms error ice v
    INTEGER :: VarId_rmsareaini                            ! RMS of initial fields
    INTEGER :: VarId_rmshiceini, VarId_rmshsnowini         ! RMS of initial fields
    INTEGER :: VarId_rmsuiceini, VarId_rmsviceini          ! RMS of initial fields
    INTEGER :: VarId_rmsareas                              ! RMS of smoothed fields
    INTEGER :: VarId_rmshices, VarId_rmshsnows             ! RMS of smoothed fields
    INTEGER :: VarId_rmsuices, VarId_rmsvices              ! RMS of smoothed fields
    INTEGER :: s       			    ! auxiliary: status counter
    INTEGER :: dimarray(3)                  ! auxiliary: array dimension
    INTEGER :: stat(500)                    ! auxiliary: status array
    CHARACTER(200)            :: filename   ! Full file name

    
    pe0: IF (writepe) THEN

! Print screen information
       IF (prec_nc == 's') THEN
          WRITE (*, '(a, 1x, a)') 'FESOM-PDAF', 'Initialize assimilation NetCDF ice file - single precision'
          nf_prec = NF_FLOAT
       ELSE
          WRITE (*, '(a, 1x, a)') 'FESOM-PDAF', 'Initialize assimilation NetCDF ice file - double precision'
          nf_prec = NF_DOUBLE
       END IF

! ----- open file and write global attributes

       filename=TRIM(ResultPath)//runid//'.'//cyearnew//'.ice.'//TRIM(str_daspec)//'.nc'

       s = 1
    
       stat(s) = NF_CREATE(filename, 0, fileid) 
       s = s + 1

       attstr  = 'FESOM Assimilation - ice fields'
       stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'title', LEN_TRIM(attstr), &
            TRIM(attstr)) 
       s = s + 1


! ----- DEFINE DIMENSIONS ------------------------------

       stat(s) = NF_DEF_DIM(fileid, 'nodes_2D',    nod2D, DimId_n2D)             
       s = s + 1
       stat(s) = NF_DEF_DIM(fileid, 'one',         1, dim1)
       s = s + 1
       stat(s) = NF_DEF_DIM(fileid, 'iter',   NF_UNLIMITED, DimId_iter)
       s = s + 1
       
       
       DO i = 1,  s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in ice dimension definitions, no.', i
       END DO

! ----- DEFINE VARIABLES ---------------------------------

       s = 1

       !- numbers

       stat(s) = NF_DEF_VAR(fileid, 'iter', NF_INT, 1, DimId_iter, VarId_iter) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'time',   nf_prec, 1, DimId_iter, VarId_time) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'asmlstep', NF_INT, 1, DimId_iter, VarId_asmlstep) 
       s = s + 1

! ----- F I E L D S 

       dimarray(1) = DimId_n2D
       dimarray(2) = DimId_iter

       ! Ice concentration
       stat(s) = NF_DEF_VAR(fileid, 'area_a', nf_prec, 2, dimarray(1:2), VarId_area_a);
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'area_f', nf_prec, 2, dimarray(1:2), VarId_area_f); 
       s = s + 1
       
       dimarray(2) = dim1
       stat(s) = NF_DEF_VAR(fileid, 'area_ini', nf_prec, 2, dimarray(1:2), VarId_area_i); 
       s = s + 1

       ! Ice thickness
       dimarray(2) = DimId_iter
       stat(s) = NF_DEF_VAR(fileid, 'hice_a', nf_prec, 2, dimarray(1:2), VarId_hice_a)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'hice_f', nf_prec, 2, dimarray(1:2), VarId_hice_f)
       s = s + 1

       dimarray(2) = dim1
       stat(s) = NF_DEF_VAR(fileid, 'hice_ini', nf_prec, 2, dimarray(1:2), VarId_hice_i)
       s = s + 1

       ! Snow height
       dimarray(2) = DimId_iter
       stat(s) = NF_DEF_VAR(fileid, 'hsnow_a', nf_prec, 2, dimarray(1:2), VarId_hsnow_a)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'hsnow_f', nf_prec, 2, dimarray(1:2), VarId_hsnow_f)
       s = s + 1

       dimarray(2) = dim1
       stat(s) = NF_DEF_VAR(fileid, 'hsnow_ini', nf_prec, 2, dimarray(1:2), VarId_hsnow_i)
       s = s + 1

       !- velocity in nodes

       dimarray(2) = DimId_iter
       stat(s) = NF_DEF_VAR(fileid, 'uice_a', nf_prec, 2, dimarray(1:2), VarId_uice_a)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'uice_f', nf_prec, 2, dimarray(1:2), VarId_uice_f)
       s = s + 1

       dimarray(2) = dim1
       stat(s) = NF_DEF_VAR(fileid, 'uice_ini', nf_prec, 2, dimarray(1:2), VarId_uice_i)
       s = s + 1

       dimarray(2) = DimId_iter
       stat(s) = NF_DEF_VAR(fileid, 'vice_a', nf_prec, 2, dimarray(1:2), VarId_vice_a)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'vice_f', nf_prec, 2, dimarray(1:2), VarId_vice_f)
       s = s + 1

       dimarray(2) = dim1
       stat(s) = NF_DEF_VAR(fileid, 'vice_ini', nf_prec, 2, dimarray(1:2), VarId_vice_i)
       s = s + 1


! ----- SMOOTHER F I E L D S 

       smootherA: IF (dim_lag > 0) THEN
          dimarray(1) = DimId_n2D
          dimarray(2) = DimId_iter

          ! Ice concentration
          stat(s) = NF_DEF_VAR(fileid, 'area_s', nf_prec, 2, dimarray(1:2), VarId_area_s)
          s = s + 1

          ! Ice thickness
          stat(s) = NF_DEF_VAR(fileid, 'hice_s', nf_prec, 2, dimarray(1:2), VarId_hice_s)
          s = s + 1

          ! Snow height
          stat(s) = NF_DEF_VAR(fileid, 'hsnow_s', nf_prec, 2, dimarray(1:2), VarId_hsnow_s)
          s = s + 1

          !- velocity in nodes
          stat(s) = NF_DEF_VAR(fileid, 'uice_s', nf_prec, 2, dimarray(1:2), VarId_uice_s)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'vice_s', nf_prec, 2, dimarray(1:2), VarId_vice_s)
          s = s + 1

       END IF smootherA

       DO i = 1,  s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in ice variable definitions, no.', i
       END DO

       !- RMS errors

       s = 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_area_ini', nf_prec, 1, Dim1, VarId_rmsareaini) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_hice_ini', nf_prec, 1, Dim1, VarId_rmshiceini) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_hsnow_ini', nf_prec, 1, Dim1, VarId_rmshsnowini) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_uice_ini', nf_prec, 1, Dim1, VarId_rmsuiceini) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_vice_ini', nf_prec, 1, Dim1, VarId_rmsviceini) 
       s = s + 1

       stat(s) = NF_DEF_VAR(fileid, 'rms_area_a', nf_prec, 1, DimId_iter, VarId_rmsareaa) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_hice_a', nf_prec, 1, DimId_iter, VarId_rmshicea) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_hsnow_a', nf_prec, 1, DimId_iter, VarId_rmshsnowa) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_uice_a', nf_prec, 1, DimId_iter, VarId_rmsuicea) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_vice_a', nf_prec, 1, DimId_iter, VarId_rmsvicea) 
       s = s + 1

       stat(s) = NF_DEF_VAR(fileid, 'rms_area_f', nf_prec, 1, DimId_iter, VarId_rmsareaf) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_hice_f', nf_prec, 1, DimId_iter, VarId_rmshicef) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_hsnow_f', nf_prec, 1, DimId_iter, VarId_rmshsnowf) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_uice_f', nf_prec, 1, DimId_iter, VarId_rmsuicef) 
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'rms_vice_f', nf_prec, 1, DimId_iter, VarId_rmsvicef) 
       s = s + 1

       smootherB: IF (dim_lag > 0) THEN
          stat(s) = NF_DEF_VAR(fileid, 'rms_area_s', nf_prec, 1, DimId_iter, VarId_rmsareas) 
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_hice_s', nf_prec, 1, DimId_iter, VarId_rmshices) 
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_hsnow_s', nf_prec, 1, DimId_iter, VarId_rmshsnows) 
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_uice_s', nf_prec, 1, DimId_iter, VarId_rmsuices)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'rms_vice_s', nf_prec, 1, DimId_iter, VarId_rmsvices)
          s = s + 1
       END IF smootherB

       DO i = 1,  s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in rms error ice variable definition, no.', i
       END DO

! ----- DEFINE ATTRIBUTES -----------------------------------------

       s = 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_time,     'long_name',    4, 'time'); 
       s = s + 1  
       attstr  = 'seconds'
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_time,     'units',        LEN_TRIM(attstr), attstr) 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_iter,     'long_name',     10, 'iterations') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_asmlstep, 'long_name',     17, 'Assimilation step') 
       s = s + 1

       ! Fields
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_area_f,      'long_name',   37, 'ice concentration [0 to 1] - forecast') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_area_a,      'long_name',   37, 'ice concentration [0 to 1] - analysis') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_area_i,      'long_name',   36, 'ice concentration [0 to 1] - initial') 
       s = s + 1

       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hice_f,      'long_name',   34, 'effective ice thickness - forecast') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hice_f,      'units',       1, 'm') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hice_a,      'long_name',   34, 'effective ice thickness - analysis') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hice_a,      'units',       1, 'm') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hice_i,      'long_name',   33, 'effective ice thickness - initial') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hice_i,      'units',       1, 'm') 
       s = s + 1

       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hsnow_f,      'long_name',   35, 'effective snow thickness - forecast') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hsnow_f,      'units',       1, 'm') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hsnow_a,      'long_name',   35, 'effective snow thickness - analysis') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hsnow_a,      'units',       1, 'm') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hsnow_i,      'long_name',   34, 'effective snow thickness - initial') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hsnow_i,      'units',       1, 'm') 
       s = s + 1

       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_uice_f,      'long_name',   25, 'zonal velocity - forecast') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_uice_f,      'units',       3, 'm/s') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_uice_a,      'long_name',   25, 'zonal velocity - analysis') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_uice_a,      'units',       3, 'm/s') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_uice_i,      'long_name',   24, 'zonal velocity - initial') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_uice_i,      'units',       3, 'm/s') 
       s = s + 1

       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vice_f,      'long_name',   30, 'meridional velocity - forecast') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vice_f,      'units',       3, 'm/s') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vice_a,      'long_name',   30, 'meridional velocity - analysis') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vice_a,      'units',       3, 'm/s') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vice_i,      'long_name',   29, 'meridional velocity - initial') 
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vice_i,      'units',       3, 'm/s') 
       s = s + 1

       DO i = 1, s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in ice file attribute assignments, no.', i
       END DO

       s = 1
       
       stat(s) = NF_ENDDEF(fileid) 
       s = s + 1
       stat(1) = NF_CLOSE(fileid)

       IF (stat(1) /= NF_NOERR) THEN
          WRITE(*, *) 'NetCDF error in closing NetCDF ice file'
       END IF
    END IF pe0

  END SUBROUTINE init_ncfile_ice_pdaf


!----------------------------------------------------------------------------
!>  Write global fields into NetCDF files
!!
  SUBROUTINE write_netcdf_pdaf(writetype, write_pos_da, iteration, dim, state_l, &
       nfields, rmse, writepe)

    IMPLICIT NONE

! *** Arguments ***
    CHARACTER(len=1), INTENT(in) :: writetype   !< Write (i) initial, (a) assimilated, (f) forecast fields
    INTEGER, INTENT(in) :: write_pos_da         !< Write position
    INTEGER, INTENT(in) :: iteration            !< Current model time step
    INTEGER, INTENT(in) :: dim                  !< Size of state vector
    REAL, INTENT(in) :: state_l(dim)            !< State vector
    INTEGER, INTENT(in) :: nfields              !< number of fields in state vector
    REAL, INTENT(in) :: rmse(nfields)           !< Array of RMS errors
    LOGICAL, INTENT(in) :: writepe              !< Whether this process writes

    ! Write ocean fields
    CALL write_nc_oce_pdaf(writetype, write_pos_da, iteration, dim, state_l, &
         nfields, rmse, writepe)

#ifdef use_ice
    ! Write oce fields
    CALL write_nc_ice_pdaf(writetype, write_pos_da, iteration, dim, state_l, &
         nfields, rmse, writepe)
#endif

  END SUBROUTINE write_netcdf_pdaf


!----------------------------------------------------------------------------
!>  Write global ensemble fields into NetCDF files
!!
  SUBROUTINE write_netcdf_pdaf_ens(writetype, write_pos_da, iteration, dim, state_l, &
       nfields, rmse, writepe, dim_ens)

    IMPLICIT NONE

! *** Arguments ***
    CHARACTER(len=1), INTENT(in) :: writetype   !< Write (i) initial, (a) assimilated, (f) forecast fields
    INTEGER, INTENT(in) :: write_pos_da         !< Write position
    INTEGER, INTENT(in) :: iteration            !< Current model time step
    INTEGER, INTENT(in) :: dim                  !< Size of state vector
    REAL, INTENT(in) :: state_l(dim)            !< State vector
    INTEGER, INTENT(in) :: nfields              !< number of fields in state vector
    REAL, INTENT(in) :: rmse(nfields)           !< Array of RMS errors
    LOGICAL, INTENT(in) :: writepe              !< Whether this process writes
    INTEGER, INTENT(in) :: dim_ens              !< Ensemble size


    ! Write ocean fields
    CALL write_nc_oce_pdaf_ens(writetype, write_pos_da, iteration, dim, state_l, &
         nfields, rmse, writepe, dim_ens)

#ifdef use_ice
    ! Write oce fields
    CALL write_nc_ice_pdaf(writetype, write_pos_da, iteration, dim, state_l, &
         nfields, rmse, writepe)
#endif

  END SUBROUTINE write_netcdf_pdaf_ens


!----------------------------------------------------------------------------
!>  Write global ocean fields into NetCDF file
!!
  SUBROUTINE write_nc_oce_pdaf(writetype, write_pos_da, iteration, dim, state_l, &
       nfields, rmse, writepe)

    USE mod_assim_pdaf, &
         ONLY: off_fields_p, istep_asml, eff_dim_obs
    USE obs_SST_CMEMS_pdafomi, &
         ONLY: loc_radius_sst
    USE g_config, &
         ONLY: runid, ResultPath
    USE o_mesh
    USE g_clock
    USE g_parfe, &
         ONLY: myDim_nod2D, myDim_nod3D, eDim_nod3D
    USE o_array, &
         ONLY: uf, ssh, tracer, w

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

! *** Arguments ***
    CHARACTER(len=1), INTENT(in) :: writetype     !< Write (i) initial, (a) assimilated, (f) forecast fields
    INTEGER, INTENT(in) :: write_pos_da           !< Write position
    INTEGER, INTENT(in) :: iteration              !< Current model time step
    INTEGER, INTENT(in) :: dim                    !< Size of state vector
    REAL, INTENT(in) :: state_l(dim)              !< State vector
    INTEGER, INTENT(in) :: nfields                !< number of fields in state vector
    REAL, INTENT(in) :: rmse(nfields)             !< Array of RMS errors
    LOGICAL, INTENT(in) :: writepe                !< Whether this process writes

! *** Local variables ***
    INTEGER :: FileId                    ! Id of Netcdf file
    INTEGER :: i, n                      ! Counters
    INTEGER :: pos1, nmb                 ! Position index for writing
    INTEGER :: pos1vec(2), nmbvec(2)     ! Position index arrays for writing
    INTEGER :: pos1vec3(3), nmbvec3(3)   ! Position index arrays for writing
    INTEGER :: VarId_iter, VarId_time    ! Id numbers
    INTEGER :: VarId_asmlstep            ! Number of assimilation step
    INTEGER :: VarId_ssh                 ! Ids for fields: ssh
    INTEGER :: VarId_effdimobs, VarId_locradius      ! Variables: observation dim and loc. radius
    INTEGER :: VarId_temp, VarId_salt    ! Ids for fields: temperature, salt
    INTEGER :: VarId_u, VarId_v, VarId_w ! Ids for fields: velocity components
    INTEGER :: VarId_tsurf, VarId_ssurf  ! Ids for forcing fields: T and S
    INTEGER :: VarID_taux, VarId_tauy    ! Ids for forcing fields: wind stress
    INTEGER :: VarId_heatf, VarId_fwaterf   ! Ids for forcing fields: fluxes
    INTEGER :: VarID_rmsssh                        ! variable: RMS errors SSH
    INTEGER :: VarId_rmsu, VarId_rmsv, VarId_rmsw  ! variable: RMS errors u
    INTEGER :: VarId_rmstemp, VarId_rmssalt        ! variable: RMS errors temperature
    INTEGER :: stat(50)                  ! Status array
    INTEGER :: s                         ! status counter
    REAL, ALLOCATABLE :: temp_arr3d(:) ! Temporary array for 3D fields
    REAL, ALLOCATABLE :: temp_arr2d(:) ! Temporary array for 2D fields
   
  REAL(kind=8)              :: sec_in_year
  CHARACTER(200)            :: filename


  filename=TRIM(ResultPath)//runid//'.'//cyearnew//'.oce.'//TRIM(str_daspec)//'.nc'

! Print screen information
    IF (writepe) THEN
       IF (writetype == 'i') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'FESOM-PDAF', 'Write initial ocean state to NetCDF at step ', &
               istep, ' position ', write_pos_da
       ELSE IF (writetype== 'f') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'FESOM-PDAF', 'Write ocean forecast to NetCDF at step ', &
               istep-1, ' position ', write_pos_da
       ELSE IF (writetype== 'a') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'FESOM-PDAF', 'Write ocean analysis to NetCDF at step ', &
               istep-1, ' position ', write_pos_da
       ELSE IF (writetype== 's') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'FESOM-PDAF', 'Write smoothed ocean to NetCDF at step ', &
               istep-1, ' position ', write_pos_da
       END IF
    END IF

    sec_in_year = dt * REAL(istep_asml)


! Gather full fields on PE 0 and write into files
    ALLOCATE(temp_arr3d(nod3D), temp_arr2d(nod2D))
    temp_arr3d = 0.0
    temp_arr2d = 0.0

    pe0: IF (writepe) THEN

! ----- Open Netcdf File
       stat(1) = NF_OPEN(filename, NF_WRITE, FileId)

       IF (stat(1) /= NF_NOERR) STOP 'nc-file error'

! ----- INQUIRE VARIABLE IDs

       s = 1

       stat(s) = NF_INQ_VARID(fileid, "iter", VarId_iter) 
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "time",      VarId_time) 
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "asmlstep", VarId_asmlstep) 
       s = s + 1
       IF (writetype=='i') THEN
          stat(s) = NF_INQ_VARID(fileid, "ssh_ini",       VarId_ssh) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "temp_ini",      VarId_temp) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "salt_ini",      VarId_salt) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "u_ini",         VarId_u) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "v_ini",         VarId_v) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "wpot_ini",      VarId_w) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_ssh_ini",       VarId_rmsssh) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_u_ini",       VarId_rmsu) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_v_ini",       VarId_rmsv) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_wpot_ini",       VarId_rmsw) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_temp_ini",       VarId_rmstemp) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_salt_ini",       VarId_rmssalt) 
          s = s + 1
       ELSE IF (writetype=='a') THEN
          stat(s) = NF_INQ_VARID(fileid, "ssh_a",       VarId_ssh) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "temp_a",      VarId_temp) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "salt_a",      VarId_salt) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "u_a",         VarId_u) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "v_a",         VarId_v) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "wpot_a",      VarId_w) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_ssh_a",       VarId_rmsssh) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_u_a",       VarId_rmsu) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_v_a",       VarId_rmsv) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_wpot_a",       VarId_rmsw) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_temp_a",       VarId_rmstemp) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_salt_a",       VarId_rmssalt) 
          s = s + 1
       ELSE IF (writetype=='f') THEN
          stat(s) = NF_INQ_VARID(fileid, "ssh_f",       VarId_ssh) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "temp_f",      VarId_temp) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "salt_f",      VarId_salt) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "u_f",         VarId_u) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "v_f",         VarId_v) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "wpot_f",        VarId_w) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_ssh_f",       VarId_rmsssh) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_u_f",       VarId_rmsu) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_v_f",       VarId_rmsv) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_wpot_f",       VarId_rmsw) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_temp_f",       VarId_rmstemp) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_salt_f",       VarId_rmssalt) 
          s = s + 1
       ELSE IF (writetype=='s') THEN
          stat(s) = NF_INQ_VARID(fileid, "ssh_s",       VarId_ssh) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "temp_s",      VarId_temp) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "salt_s",      VarId_salt) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "u_s",         VarId_u) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "v_s",         VarId_v) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "wpot_s",      VarId_w) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_ssh_s",       VarId_rmsssh) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_u_s",       VarId_rmsu) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_v_s",       VarId_rmsv) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_wpot_s",       VarId_rmsw) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_temp_s",       VarId_rmstemp) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_salt_s",       VarId_rmssalt) 
          s = s + 1
       END IF
       stat(s) = NF_INQ_VARID(fileid, "eff_dim_obs",      VarId_effdimobs) 
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "loc_radius",       VarId_locradius) 
       s = s + 1

       DO i = 1, s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error inquiring variable IDs, no.', i
       END DO

! ----- WRITE VARIABLES
       s = 1
    
       IF (.NOT. DEBUGOUTPUT) THEN
          ! Normal output
          pos1 = write_pos_da
       ELSE
          ! Write keeping only a single time instance
          pos1 = 1
       END IF
       nmb  = 1
    
       IF (writetype/='s') THEN
          stat(s) = NF_PUT_VARA_INT(fileid,  VarId_asmlstep, pos1, nmb, write_pos_da) 
          s = s + 1
          stat(s) = NF_PUT_VARA_INT(fileid,  VarId_iter, pos1, nmb, istep_asml) 
          s = s + 1
          IF (prec_nc == 's') THEN
             ! Write single precision
             stat(s) = NF_PUT_VARA_REAL(fileid, VarId_time, pos1, nmb, REAL(sec_in_year, 4))
          ELSE
             ! Write double precision
             stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_time, pos1, nmb, sec_in_year)
          END IF
          s = s + 1
       END IF
    END IF pe0

! ----- 2D FIELDS

    !----- topography
       
    !----- sea surface height
    DO i = 1, myDim_nod2D
       ssh(i) = state_l(i + off_fields_p(1))
    END DO
    CALL broadcast2d_pdaf(ssh, temp_arr2D)
    pe0a: IF (writepe) THEN
       IF (.NOT. DEBUGOUTPUT) THEN
          ! Normal output
          pos1vec = (/ 1, write_pos_da  /)
       ELSE
          ! Write keeping only a single time instance
          pos1vec = (/ 1, 1  /)
       END IF
       nmbvec  = (/ nod2D, 1 /)

       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_ssh, pos1vec, nmbvec, &
               REAL(temp_arr2d(:),4)) 
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_ssh, pos1vec, nmbvec, &
               temp_arr2d(:)) 
       END IF
       s = s + 1
    END IF pe0a


    IF (writetype=='a') THEN
    !----- Effective observation dimension

       CALL broadcast2d_pdaf(eff_dim_obs, temp_arr2D)
       pe0z: IF (writepe) THEN
          IF (.NOT. DEBUGOUTPUT) THEN
             ! Normal output
             pos1vec = (/ 1, write_pos_da  /)
          ELSE
             ! Write keeping only a single time instance
             pos1vec = (/ 1, 1  /)
          END IF
          nmbvec  = (/ nod2D, 1 /)

          IF (prec_nc == 's') THEN
             ! Write single precision
             stat(s) = NF_PUT_VARA_REAL(fileid, VarId_effdimobs, pos1vec, nmbvec, &
                  REAL(temp_arr2d(:),4)) 
          ELSE
             ! Write double precision
             stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_effdimobs, pos1vec, nmbvec, &
                  temp_arr2d(:)) 
          END IF
          s = s + 1
       END IF pe0z


    !----- Localization radius

       CALL broadcast2d_pdaf(loc_radius_sst, temp_arr2D)
       pe0z1: IF (writepe) THEN
          IF (.NOT. DEBUGOUTPUT) THEN
             ! Normal output
             pos1vec = (/ 1, write_pos_da  /)
          ELSE
             ! Write keeping only a single time instance
             pos1vec = (/ 1, 1  /)
          END IF
          nmbvec  = (/ nod2D, 1 /)

          IF (prec_nc == 's') THEN
             ! Write single precision
             stat(s) = NF_PUT_VARA_REAL(fileid, VarId_locradius, pos1vec, nmbvec, &
                  REAL(temp_arr2d(:),4)) 
          ELSE
             ! Write double precision
             stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_locradius, pos1vec, nmbvec, &
                  temp_arr2d(:)) 
          END IF
          s = s + 1
       END IF pe0z1
    END IF


    !----- Temperature
    DO i = 1, myDim_nod3D
       tracer(i, 1) = state_l(i + off_fields_p(5))
    END DO
    CALL broadcast3d_pdaf(tracer(:,1), temp_arr3d)
    pe0b: IF (writepe) THEN

       nmbvec  = (/ nod3D, 1 /)

       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_temp, pos1vec, nmbvec, &
               REAL(temp_arr3d(:), 4)) 
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_temp, pos1vec, nmbvec, &
               temp_arr3d(:)) 
       END IF
       s = s + 1
    END IF pe0b

    !----- Salinity
    DO i = 1, myDim_nod3D
       tracer(i, 2) = state_l(i + off_fields_p(6))
    END DO
    CALL broadcast3d_pdaf(tracer(:,2), temp_arr3d)
    pe0c: IF (writepe) THEN
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_salt, pos1vec, nmbvec, &
               REAL(temp_arr3d(:), 4)) 
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_salt, pos1vec, nmbvec, &
               temp_arr3d(:)) 
       END IF
       s = s + 1
    END IF pe0c

! ----- VECTOR FIELDS

    !----- nodal velocity
    DO i = 1, myDim_nod3D
       uf(i) = state_l(i + off_fields_p(2))
    END DO
    CALL broadcast3d_pdaf(uf(1 : myDim_nod3D), temp_arr3d)
    pe0d: IF (writepe) THEN
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_u, pos1vec, nmbvec, &
               REAL(temp_arr3d(:), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_u, pos1vec, nmbvec, &
               temp_arr3d(:))
       END IF
       s = s + 1
    END IF pe0d

    DO i = 1, myDim_nod3D
       uf(i + myDim_nod3D + eDim_nod3D) = state_l(i + off_fields_p(3))
    END DO
    CALL broadcast3d_pdaf(uf(1 + myDim_nod3D + eDim_nod3D : 2*myDim_nod3D + eDim_nod3D), temp_arr3d)
    pe0e: IF (writepe) THEN
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_v, pos1vec, nmbvec, &
               REAL(temp_arr3d(:), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_v, pos1vec, nmbvec, &
               temp_arr3d(:))
       END IF
       s = s + 1
    END IF pe0e

    DO i = 1, myDim_nod3D
       w(i) = state_l(i + off_fields_p(4))
    END DO
    CALL broadcast3d_pdaf(w(1 : myDim_nod3D), temp_arr3d)
    pe0f: IF (writepe) THEN
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_w, pos1vec, nmbvec, &
               REAL(temp_arr3d(:), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_w, pos1vec, nmbvec, &
               temp_arr3d(:))
       END IF
       s = s + 1
       DO i = 1, s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in writing variables, no.', i
       END DO

! ----- RMS errors

       s = 1

       IF (.NOT. DEBUGOUTPUT) THEN
          ! Normal output
          pos1 = write_pos_da
       ELSE
          ! Write keeping only a single time instance
          pos1 = 1
       END IF
       nmb  = 1

       ! RMS of SSH
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_rmsssh, pos1, nmb, REAL(rmse(1), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_rmsssh, pos1, nmb, rmse(1))
       END IF
       s = s + 1

       ! RMS of U
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_rmsu, pos1, nmb, REAL(rmse(2), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_rmsu, pos1, nmb, rmse(2))
       END IF
       s = s + 1

       ! RMS of V
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_rmsv, pos1, nmb, REAL(rmse(3), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_rmsv, pos1, nmb, rmse(3))
       END IF
       s = s + 1

       ! RMS of Wpot
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_rmsw, pos1, nmb, REAL(rmse(4), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_rmsw, pos1, nmb, rmse(4))
       END IF
       s = s + 1

       ! RMS of temperature
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_rmstemp, pos1, nmb, REAL(rmse(5), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_rmstemp, pos1, nmb, rmse(5))
       END IF
       s = s + 1

       ! RMS of salinity
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_rmssalt, pos1, nmb, REAL(rmse(6), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_rmssalt, pos1, nmb, rmse(6))
       END IF
       s = s + 1

       DO i = 1, s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in writing RMS errors, no.', i
       END DO

! ----- CLOSE THE FILE

       stat(1) = NF_CLOSE(fileid)

       IF (stat(1) /= NF_NOERR) THEN
          WRITE(*, *) 'NetCDF error in closing NetCDF file'
       END IF
    END IF pe0f

! ----- Clean up

    DEALLOCATE(temp_arr3d, temp_arr2d)

  END SUBROUTINE write_nc_oce_pdaf


!----------------------------------------------------------------------------
!> Write ensemble of global ocean fields into NetCDF file
!!
  SUBROUTINE write_nc_oce_pdaf_ens(writetype, write_pos_da, iteration, dim, state_l, &
       nfields, rmse, writepe, dim_ens)

    USE g_config, &
         ONLY: runid, ResultPath
    USE o_mesh
    USE g_clock
    USE g_parfe, &
         ONLY: myDim_nod2D, myDim_nod3D, eDim_nod3D
    USE o_array, &
         ONLY: uf, ssh, tracer, w
    USE mod_assim_pdaf, &
         ONLY: off_fields_p, istep_asml, eff_dim_obs
    USE obs_SST_CMEMS_pdafomi, &
         ONLY: loc_radius_sst

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

! *** Arguments ***
    CHARACTER(len=1), INTENT(in) :: writetype     !< Write (i) initial, (a) assimilated, (f) forecast fields
    INTEGER, INTENT(in) :: write_pos_da           !< Write position
    INTEGER, INTENT(in) :: iteration              !< Current model time step
    INTEGER, INTENT(in) :: dim                    !< Size of state vector
    INTEGER, INTENT(in) :: dim_ens                !< Ensemble size
    REAL, INTENT(in) :: state_l(dim, dim_ens)     !< State vector
    INTEGER, INTENT(in) :: nfields                !< number of fields in state vector
    REAL, INTENT(in) :: rmse(nfields)             !< Array of RMS errors
    LOGICAL, INTENT(in) :: writepe                !< Whether this process writes

! *** Local variables ***
    INTEGER :: FileId                    ! Id of Netcdf file
    INTEGER :: i, n, member              ! Counters
    CHARACTER(len=10) :: member_string
    INTEGER :: pos1, nmb                 ! Position index for writing
    INTEGER :: pos1vec(2), nmbvec(2)     ! Position index arrays for writing
    INTEGER :: pos1vec3(3), nmbvec3(3)   ! Position index arrays for writing
    INTEGER :: VarId_iter, VarId_time    ! Id numbers
    INTEGER :: VarId_asmlstep            ! Number of assimilation step
    INTEGER :: VarId_ssh                 ! Ids for fields: ssh
    INTEGER :: VarId_effdimobs, VarId_locradius      ! Variables: observation dim and loc. radius
    INTEGER :: VarId_temp, VarId_salt    ! Ids for fields: temperature, salt
    INTEGER :: VarId_u, VarId_v, VarId_w ! Ids for fields: velocity components
    INTEGER :: VarId_tsurf, VarId_ssurf  ! Ids for forcing fields: T and S
    INTEGER :: VarID_taux, VarId_tauy    ! Ids for forcing fields: wind stress
    INTEGER :: VarId_heatf, VarId_fwaterf   ! Ids for forcing fields: fluxes
    INTEGER :: VarID_rmsssh                        ! variable: RMS errors SSH
    INTEGER :: VarId_rmsu, VarId_rmsv, VarId_rmsw  ! variable: RMS errors u
    INTEGER :: VarId_rmstemp, VarId_rmssalt        ! variable: RMS errors temperature
    INTEGER :: stat(50)                  ! Status array
    INTEGER :: s                         ! status counter
    REAL, ALLOCATABLE :: temp_arr3d(:) ! Temporary array for 3D fields
    REAL, ALLOCATABLE :: temp_arr2d(:) ! Temporary array for 2D fields
    REAL, ALLOCATABLE :: ssh_temp(:)
    REAL, ALLOCATABLE :: w_temp(:)
    REAL, ALLOCATABLE :: uf_temp(:)
    REAL, ALLOCATABLE :: tracer_temp(:,:)

  REAL(kind=8)              :: sec_in_year
  CHARACTER(200)            :: filename


! Allocate local variables
    ALLOCATE(ssh_temp(myDim_nod2D))
    ALLOCATE(w_temp(myDim_nod3D))
    ALLOCATE(uf_temp(2*myDim_nod3D + eDim_nod3D))
    ALLOCATE(tracer_temp(myDim_nod3D,2))
    
DO member = 1, dim_ens
  WRITE(member_string,'(i4.4)') member
  filename=TRIM(ResultPath)//runid//'.'//cyearnew//'.oce.'//TRIM(str_daspec)//'_'//TRIM(member_string)//'.nc'

! Print screen information
    IF (writepe) THEN
       IF (writetype == 'i') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'FESOM-PDAF', 'Write initial ocean state to NetCDF at step ', &
               istep, ' position ', write_pos_da
       ELSE IF (writetype== 'f') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'FESOM-PDAF', 'Write ocean forecast to NetCDF at step ', &
               istep-1, ' position ', write_pos_da
       ELSE IF (writetype== 'a') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'FESOM-PDAF', 'Write ocean analysis to NetCDF at step ', &
               istep-1, ' position ', write_pos_da
       ELSE IF (writetype== 's') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'FESOM-PDAF', 'Write smoothed ocean to NetCDF at step ', &
               istep-1, ' position ', write_pos_da
       END IF
    END IF

    sec_in_year = dt * REAL(istep_asml)


! Gather full fields on PE 0 and write into files
    ALLOCATE(temp_arr3d(nod3D), temp_arr2d(nod2D))
    temp_arr3d = 0.0
    temp_arr2d = 0.0

    pe0: IF (writepe) THEN

! ----- Open Netcdf File
       stat(1) = NF_OPEN(filename, NF_WRITE, FileId)

       IF (stat(1) /= NF_NOERR) STOP 'nc-file error'

! ----- INQUIRE VARIABLE IDs

       s = 1

       stat(s) = NF_INQ_VARID(fileid, "iter", VarId_iter)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "time",      VarId_time)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "asmlstep", VarId_asmlstep)
       s = s + 1
       IF (writetype=='i') THEN
          stat(s) = NF_INQ_VARID(fileid, "ssh_ini",       VarId_ssh)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "temp_ini",      VarId_temp)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "salt_ini",      VarId_salt)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "u_ini",         VarId_u)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "v_ini",         VarId_v)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "wpot_ini",      VarId_w)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_ssh_ini",       VarId_rmsssh)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_u_ini",       VarId_rmsu)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_v_ini",       VarId_rmsv)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_wpot_ini",       VarId_rmsw)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_temp_ini",       VarId_rmstemp)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_salt_ini",       VarId_rmssalt)
          s = s + 1
       ELSE IF (writetype=='a') THEN
          stat(s) = NF_INQ_VARID(fileid, "ssh_a",       VarId_ssh)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "temp_a",      VarId_temp)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "salt_a",      VarId_salt)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "u_a",         VarId_u)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "v_a",         VarId_v)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "wpot_a",      VarId_w)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_ssh_a",       VarId_rmsssh)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_u_a",       VarId_rmsu)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_v_a",       VarId_rmsv)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_wpot_a",       VarId_rmsw)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_temp_a",       VarId_rmstemp)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_salt_a",       VarId_rmssalt)
          s = s + 1
       ELSE IF (writetype=='f') THEN
          stat(s) = NF_INQ_VARID(fileid, "ssh_f",       VarId_ssh)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "temp_f",      VarId_temp)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "salt_f",      VarId_salt)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "u_f",         VarId_u)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "v_f",         VarId_v)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "wpot_f",        VarId_w)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_ssh_f",       VarId_rmsssh)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_u_f",       VarId_rmsu)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_v_f",       VarId_rmsv)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_wpot_f",       VarId_rmsw)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_temp_f",       VarId_rmstemp)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_salt_f",       VarId_rmssalt)
          s = s + 1
       ELSE IF (writetype=='s') THEN
          stat(s) = NF_INQ_VARID(fileid, "ssh_s",       VarId_ssh)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "temp_s",      VarId_temp)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "salt_s",      VarId_salt)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "u_s",         VarId_u)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "v_s",         VarId_v)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "wpot_s",      VarId_w)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_ssh_s",       VarId_rmsssh)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_u_s",       VarId_rmsu)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_v_s",       VarId_rmsv)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_wpot_s",       VarId_rmsw)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_temp_s",       VarId_rmstemp)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_salt_s",       VarId_rmssalt)
          s = s + 1
       END IF
       stat(s) = NF_INQ_VARID(fileid, "eff_dim_obs",      VarId_effdimobs)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "loc_radius",       VarId_locradius)
       s = s + 1

       DO i = 1, s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error inquiring variable IDs, no.', i
       END DO

! ----- WRITE VARIABLES
       s = 1

       IF (.NOT. DEBUGOUTPUT) THEN
          ! Normal output
          pos1 = write_pos_da
       ELSE
          ! Write keeping only a single time instance
          pos1 = 1
       END IF
       nmb  = 1

       IF (writetype/='s') THEN
          stat(s) = NF_PUT_VARA_INT(fileid,  VarId_asmlstep, pos1, nmb, write_pos_da)
          s = s + 1
          stat(s) = NF_PUT_VARA_INT(fileid,  VarId_iter, pos1, nmb, istep_asml)
          s = s + 1
          IF (prec_nc == 's') THEN
             ! Write single precision
             stat(s) = NF_PUT_VARA_REAL(fileid, VarId_time, pos1, nmb, REAL(sec_in_year, 4))
          ELSE
             ! Write double precision
             stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_time, pos1, nmb, sec_in_year)
          END IF
          s = s + 1
       END IF
    END IF pe0

! ----- 2D FIELDS

    !----- topography

    !----- sea surface height
    DO i = 1, myDim_nod2D
       ssh_temp(i) = state_l(i + off_fields_p(1), member)
    END DO
    CALL broadcast2d_pdaf(ssh_temp, temp_arr2D)
    pe0a: IF (writepe) THEN
       IF (.NOT. DEBUGOUTPUT) THEN
          ! Normal output
          pos1vec = (/ 1, write_pos_da  /)
       ELSE
          ! Write keeping only a single time instance
          pos1vec = (/ 1, 1  /)
       END IF
       nmbvec  = (/ nod2D, 1 /)

       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_ssh, pos1vec, nmbvec, &
               REAL(temp_arr2d(:),4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_ssh, pos1vec, nmbvec, &
               temp_arr2d(:))
       END IF
       s = s + 1
    END IF pe0a


    IF (writetype=='a') THEN
    !----- Effective observation dimension

       CALL broadcast2d_pdaf(eff_dim_obs, temp_arr2D)
       pe0z: IF (writepe) THEN
          IF (.NOT. DEBUGOUTPUT) THEN
             ! Normal output
             pos1vec = (/ 1, write_pos_da  /)
          ELSE
             ! Write keeping only a single time instance
             pos1vec = (/ 1, 1  /)
          END IF
          nmbvec  = (/ nod2D, 1 /)

          IF (prec_nc == 's') THEN
             ! Write single precision
             stat(s) = NF_PUT_VARA_REAL(fileid, VarId_effdimobs, pos1vec, nmbvec, &
                  REAL(temp_arr2d(:),4))
          ELSE
             ! Write double precision
             stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_effdimobs, pos1vec, nmbvec, &
                  temp_arr2d(:))
          END IF
          s = s + 1
       END IF pe0z


    !----- Localization radius

       CALL broadcast2d_pdaf(loc_radius_sst, temp_arr2D)
       pe0z1: IF (writepe) THEN
          IF (.NOT. DEBUGOUTPUT) THEN
             ! Normal output
             pos1vec = (/ 1, write_pos_da  /)
          ELSE
             ! Write keeping only a single time instance
             pos1vec = (/ 1, 1  /)
          END IF
          nmbvec  = (/ nod2D, 1 /)

          IF (prec_nc == 's') THEN
             ! Write single precision
             stat(s) = NF_PUT_VARA_REAL(fileid, VarId_locradius, pos1vec, nmbvec, &
                  REAL(temp_arr2d(:),4))
          ELSE
             ! Write double precision
             stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_locradius, pos1vec, nmbvec, &
                  temp_arr2d(:))
          END IF
          s = s + 1
       END IF pe0z1
    END IF


    !----- Temperature
    DO i = 1, myDim_nod3D
       tracer_temp(i, 1) = state_l(i + off_fields_p(5), member)
    END DO
    CALL broadcast3d_pdaf(tracer_temp(:,1), temp_arr3d)
    pe0b: IF (writepe) THEN

       nmbvec  = (/ nod3D, 1 /)

       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_temp, pos1vec, nmbvec, &
               REAL(temp_arr3d(:), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_temp, pos1vec, nmbvec, &
               temp_arr3d(:))
       END IF
       s = s + 1
    END IF pe0b

    !----- Salinity
    DO i = 1, myDim_nod3D
       tracer_temp(i, 2) = state_l(i + off_fields_p(6), member)
    END DO
    CALL broadcast3d_pdaf(tracer_temp(:,2), temp_arr3d)
    pe0c: IF (writepe) THEN
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_salt, pos1vec, nmbvec, &
               REAL(temp_arr3d(:), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_salt, pos1vec, nmbvec, &
               temp_arr3d(:))
       END IF
       s = s + 1
    END IF pe0c

! ----- VECTOR FIELDS

    !----- nodal velocity
    DO i = 1, myDim_nod3D
       uf_temp(i) = state_l(i + off_fields_p(2), member)
    END DO
    CALL broadcast3d_pdaf(uf_temp(1 : myDim_nod3D), temp_arr3d)
    pe0d: IF (writepe) THEN
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_u, pos1vec, nmbvec, &
               REAL(temp_arr3d(:), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_u, pos1vec, nmbvec, &
               temp_arr3d(:))
       END IF
       s = s + 1
    END IF pe0d

    DO i = 1, myDim_nod3D
       uf_temp(i + myDim_nod3D + eDim_nod3D) = state_l(i + off_fields_p(3), member)
    END DO
    CALL broadcast3d_pdaf(uf_temp(1 + myDim_nod3D + eDim_nod3D : 2*myDim_nod3D + eDim_nod3D), temp_arr3d)
    pe0e: IF (writepe) THEN
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_v, pos1vec, nmbvec, &
               REAL(temp_arr3d(:), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_v, pos1vec, nmbvec, &
               temp_arr3d(:))
       END IF
       s = s + 1
    END IF pe0e

    DO i = 1, myDim_nod3D
       w_temp(i) = state_l(i + off_fields_p(4), member)
    END DO
    CALL broadcast3d_pdaf(w_temp(1 : myDim_nod3D), temp_arr3d)
    pe0f: IF (writepe) THEN
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_w, pos1vec, nmbvec, &
               REAL(temp_arr3d(:), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_w, pos1vec, nmbvec, &
               temp_arr3d(:))
       END IF
       s = s + 1
       DO i = 1, s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in writing variables, no.', i
       END DO

! ----- RMS errors

       s = 1

       IF (.NOT. DEBUGOUTPUT) THEN
          ! Normal output
          pos1 = write_pos_da
       ELSE
          ! Write keeping only a single time instance
          pos1 = 1
       END IF
       nmb  = 1

       ! RMS of SSH
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_rmsssh, pos1, nmb, REAL(rmse(1), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_rmsssh, pos1, nmb, rmse(1))
       END IF
       s = s + 1

       ! RMS of U
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_rmsu, pos1, nmb, REAL(rmse(2), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_rmsu, pos1, nmb, rmse(2))
       END IF
       s = s + 1

       ! RMS of V
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_rmsv, pos1, nmb, REAL(rmse(3), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_rmsv, pos1, nmb, rmse(3))
       END IF
       s = s + 1

       ! RMS of Wpot
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_rmsw, pos1, nmb, REAL(rmse(4), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_rmsw, pos1, nmb, rmse(4))
       END IF
       s = s + 1

       ! RMS of temperature
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_rmstemp, pos1, nmb, REAL(rmse(5), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_rmstemp, pos1, nmb, rmse(5))
       END IF
       s = s + 1

       ! RMS of salinity
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_rmssalt, pos1, nmb, REAL(rmse(6), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_rmssalt, pos1, nmb, rmse(6))
       END IF
       s = s + 1

       DO i = 1, s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in writing RMS errors, no.', i
       END DO

! ----- CLOSE THE FILE

       stat(1) = NF_CLOSE(fileid)

       IF (stat(1) /= NF_NOERR) THEN
          WRITE(*, *) 'NetCDF error in closing NetCDF file'
       END IF
    END IF pe0f

! ----- Clean up

    DEALLOCATE(temp_arr3d, temp_arr2d)
   END DO
DEALLOCATE(ssh_temp, uf_temp, tracer_temp, w_temp)

  END SUBROUTINE write_nc_oce_pdaf_ens


!----------------------------------------------------------------------------
!>  Write global ice fields into NetCDF file
!!
  SUBROUTINE write_nc_ice_pdaf(writetype, write_pos_da, iteration, dim, state_l, &
       nfields, rmse, writepe)

    USE mod_assim_pdaf, &
         ONLY: off_fields_p, istep_asml
    USE g_config, &
         ONLY: runid, ResultPath
    USE o_mesh
    USE g_clock
    USE g_parfe, &
         ONLY: myDim_nod2D
    USE i_array, &
         ONLY: a_ice, m_ice, m_snow, u_ice, v_ice

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

! *** Arguments ***
    CHARACTER(len=1), INTENT(in) :: writetype     !< Write (i) initial, (a) assimilated, (f) forecast fields
    INTEGER, INTENT(in) :: write_pos_da           !< Write position
    INTEGER, INTENT(in) :: iteration              !< Current model time step
    INTEGER, INTENT(in) :: dim                    !< Size of state vector
    REAL, INTENT(in) :: state_l(dim)              !< State vector
    INTEGER, INTENT(in) :: nfields                !< number of fields in state vector
    REAL, INTENT(in) :: rmse(nfields)             !< Array of RMS errors
    LOGICAL, INTENT(in) :: writepe                !< Whether this process writes

! *** Local variables ***
    INTEGER :: FileId                    ! Id of Netcdf file
    INTEGER :: i, n                      ! Counters
    INTEGER :: pos1, nmb                 ! Position index for writing
    INTEGER :: pos1vec(2), nmbvec(2)     ! Position index arrays for writing
    INTEGER :: VarId_iter, VarId_time    ! Id numbers
    INTEGER :: VarId_asmlstep               ! Number of assimilation step
    INTEGER :: VarId_aice                ! Ids for fields: ice concentration
    INTEGER :: VarId_hice, VarId_hsnow   ! Ids for fields: ice and snow thickness
    INTEGER :: VarId_uice, VarId_vice    ! Ids for fields: velocity components
    INTEGER :: VarId_rmsaice
    INTEGER :: VarId_rmshice, VarID_rmshsnow
    INTEGER :: VarId_rmsuice, VarID_rmsvice
    INTEGER :: stat(50)                  ! Status array
    INTEGER :: s                         ! status counter
    REAL, ALLOCATABLE :: temp_arr2d(:) ! Temporary array for 2D fields
   
    REAL(kind=8)              :: sec_in_year
    CHARACTER(200)            :: filename


  filename=TRIM(ResultPath)//runid//'.'//cyearnew//'.ice.'//TRIM(str_daspec)//'.nc'

! Print screen information
    IF (writepe) THEN
       IF (writetype == 'i') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'FESOM-PDAF', 'Write initial ice state to NetCDF at step ', &
               istep, ' position ', write_pos_da
       ELSE IF (writetype== 'f') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'FESOM-PDAF', 'Write ice forecast to NetCDF at step ', &
               istep-1, ' position ', write_pos_da
       ELSE IF (writetype== 'a') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'FESOM-PDAF', 'Write ice analysis to NetCDF at step ', &
               istep-1, ' position ', write_pos_da
       ELSE IF (writetype== 's') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'FESOM-PDAF', 'Write smoothed ice to NetCDF at step ', &
               istep-1, ' position ', write_pos_da
       END IF
    END IF

    sec_in_year = dt * REAL(istep_asml)


! Gather full fields on PE 0 and write into files
    ALLOCATE(temp_arr2d(nod2D))
    temp_arr2d = 0.0

    pe0: IF (writepe) THEN

! ----- Open Netcdf File
       stat(1) = NF_OPEN(filename, NF_WRITE, FileId)

       IF (stat(1) /= NF_NOERR) STOP 'nc-file error'

! ----- INQUIRE VARIABLE IDs

       s = 1

       stat(s) = NF_INQ_VARID(fileid, "iter", VarId_iter) 
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "asmlstep", VarId_asmlstep) 
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "time", VarId_time) 
       s = s + 1
       IF (writetype=='i') THEN
          stat(s) = NF_INQ_VARID(fileid, "area_ini",  VarId_aice) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "hice_ini",  VarId_hice) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "hsnow_ini", VarId_hsnow) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "uice_ini",  VarId_uice)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "vice_ini",  VarId_vice)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_area_ini",  VarId_rmsaice)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_hice_ini",  VarId_rmshice)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_hsnow_ini",  VarId_rmshsnow)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_uice_ini",  VarId_rmsuice)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_vice_ini",  VarId_rmsvice)
          s = s + 1
       ELSE IF (writetype=='a') THEN
          stat(s) = NF_INQ_VARID(fileid, "area_a",  VarId_aice) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "hice_a",  VarId_hice) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "hsnow_a", VarId_hsnow) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "uice_a",  VarId_uice)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "vice_a",  VarId_vice)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_area_a",  VarId_rmsaice)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_hice_a",  VarId_rmshice)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_hsnow_a",  VarId_rmshsnow)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_uice_a",  VarId_rmsuice)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_vice_a",  VarId_rmsvice)
          s = s + 1
       ELSE IF (writetype=='f') THEN
          stat(s) = NF_INQ_VARID(fileid, "area_f",  VarId_aice) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "hice_f",  VarId_hice) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "hsnow_f", VarId_hsnow) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "uice_f",  VarId_uice)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "vice_f",  VarId_vice)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_area_f",  VarId_rmsaice)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_hice_f",  VarId_rmshice)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_hsnow_f",  VarId_rmshsnow)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_uice_f",  VarId_rmsuice)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_vice_f",  VarId_rmsvice)
          s = s + 1
       ELSE IF (writetype=='s') THEN
          stat(s) = NF_INQ_VARID(fileid, "area_s",  VarId_aice) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "hice_s",  VarId_hice) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "hsnow_s", VarId_hsnow) 
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "uice_s",  VarId_uice)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "vice_s",  VarId_vice)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_area_s",  VarId_rmsaice)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_hice_s",  VarId_rmshice)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_hsnow_s",  VarId_rmshsnow)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_uice_s",  VarId_rmsuice)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "rms_vice_s",  VarId_rmsvice)
          s = s + 1
       END IF

       DO i = 1, s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error inquiring variable IDs, no.', i
       END DO

! ----- WRITE VARIABLES
       s = 1
    
       IF (.NOT. DEBUGOUTPUT) THEN
          ! Normal output
          pos1 = write_pos_da
       ELSE
          ! Write keeping only a single time instance
          pos1 = 1
       END IF
       nmb  = 1
    
       IF (writetype/='s') THEN
          stat(s) = NF_PUT_VARA_INT(fileid,  VarId_asmlstep, pos1, nmb, write_pos_da) 
          s = s + 1
          stat(s) = NF_PUT_VARA_INT(fileid,  VarId_iter, pos1, nmb, istep_asml) 
          s = s + 1
          IF (prec_nc == 's') THEN
             ! Write single precision
             stat(s) = NF_PUT_VARA_REAL(fileid, VarId_time, pos1, nmb, REAL(sec_in_year, 4))
          ELSE
             ! Write double precision
             stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_time, pos1, nmb, sec_in_year)
          END IF
          s = s + 1
       END IF
    END IF pe0

! ----- FIELDS

    !----- ice concentration
    DO i = 1, myDim_nod2D
       a_ice(i) = state_l(i + off_fields_p(7))
    END DO
    CALL broadcast2d_pdaf(a_ice, temp_arr2D)
    pe0a: IF (writepe) THEN
       IF (.NOT. DEBUGOUTPUT) THEN
          ! Normal output
          pos1vec = (/ 1, write_pos_da  /)
       ELSE
          ! Write keeping only a single time instance
          pos1vec = (/ 1, 1  /)
       END IF
       nmbvec  = (/ nod2D, 1 /)

       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_aice, pos1vec, nmbvec, &
               REAL(temp_arr2d(:),4)) 
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_aice, pos1vec, nmbvec, &
               temp_arr2d(:)) 
       END IF
       s = s + 1
    END IF pe0a

    !----- ice thickness
    DO i = 1, myDim_nod2D
       m_ice(i) = state_l(i + off_fields_p(8))
    END DO
    CALL broadcast2d_pdaf(m_ice, temp_arr2D)
    pe0b: IF (writepe) THEN
       IF (.NOT. DEBUGOUTPUT) THEN
          ! Normal output
          pos1vec = (/ 1, write_pos_da  /)
       ELSE
          ! Write keeping only a single time instance
          pos1vec = (/ 1, 1  /)
       END IF
       nmbvec  = (/ nod2D, 1 /)

       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_hice, pos1vec, nmbvec, &
               REAL(temp_arr2d(:),4)) 
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_hice, pos1vec, nmbvec, &
               temp_arr2d(:)) 
       END IF
       s = s + 1
    END IF pe0b

    !----- snow thickness
    DO i = 1, myDim_nod2D
       m_snow(i) = state_l(i + off_fields_p(9))
    END DO
    CALL broadcast2d_pdaf(m_snow, temp_arr2D)
    pe0c: IF (writepe) THEN
       IF (.NOT. DEBUGOUTPUT) THEN
          ! Normal output
          pos1vec = (/ 1, write_pos_da  /)
       ELSE
          ! Write keeping only a single time instance
          pos1vec = (/ 1, 1  /)
       END IF
       nmbvec  = (/ nod2D, 1 /)

       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_hsnow, pos1vec, nmbvec, &
               REAL(temp_arr2d(:),4)) 
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_hsnow, pos1vec, nmbvec, &
               temp_arr2d(:)) 
       END IF
       s = s + 1
    END IF pe0c

    !----- ice zonal velocity
    DO i = 1, myDim_nod2D
       u_ice(i) = state_l(i + off_fields_p(10))
    END DO
    CALL broadcast2d_pdaf(u_ice, temp_arr2D)
    pe0d: IF (writepe) THEN
       IF (.NOT. DEBUGOUTPUT) THEN
          ! Normal output
          pos1vec = (/ 1, write_pos_da  /)
       ELSE
          ! Write keeping only a single time instance
          pos1vec = (/ 1, 1  /)
       END IF
       nmbvec  = (/ nod2D, 1 /)

       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_uice, pos1vec, nmbvec, &
               REAL(temp_arr2d(:),4)) 
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_uice, pos1vec, nmbvec, &
               temp_arr2d(:)) 
       END IF
       s = s + 1
    END IF pe0d

    !----- ice meridional velocity
    DO i = 1, myDim_nod2D
       v_ice(i) = state_l(i + off_fields_p(11))
    END DO
    CALL broadcast2d_pdaf(v_ice, temp_arr2D)
    pe0e: IF (writepe) THEN
       IF (.NOT. DEBUGOUTPUT) THEN
          ! Normal output
          pos1vec = (/ 1, write_pos_da  /)
       ELSE
          ! Write keeping only a single time instance
          pos1vec = (/ 1, 1  /)
       END IF
       nmbvec  = (/ nod2D, 1 /)

       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_vice, pos1vec, nmbvec, &
               REAL(temp_arr2d(:),4)) 
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_vice, pos1vec, nmbvec, &
               temp_arr2d(:)) 
       END IF
       s = s + 1

       DO i = 1, s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in writing variables, no.', i
       END DO

! ----- RMS errors

       s = 1
    
       IF (.NOT. DEBUGOUTPUT) THEN
          ! Normal output
          pos1 = write_pos_da
       ELSE
          ! Write keeping only a single time instance
          pos1 = 1
       END IF
       nmb  = 1
    
       ! RMS of ice concentration
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_rmsaice, pos1, nmb, REAL(rmse(7), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_rmsaice, pos1, nmb, rmse(7))
       END IF
       s = s + 1
    
       ! RMS of ice thickness
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_rmshice, pos1, nmb, REAL(rmse(8), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_rmshice, pos1, nmb, rmse(8))
       END IF
       s = s + 1

       ! RMS of snow thickness
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_rmshsnow, pos1, nmb, REAL(rmse(9), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_rmshsnow, pos1, nmb, rmse(9))
       END IF
       s = s + 1

       ! RMS of ice thickness
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_rmsuice, pos1, nmb, REAL(rmse(10), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_rmsuice, pos1, nmb, rmse(10))
       END IF
       s = s + 1

       ! RMS of ice thickness
       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_rmsvice, pos1, nmb, REAL(rmse(11), 4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_rmsvice, pos1, nmb, rmse(11))
       END IF
       s = s + 1

       DO i = 1, s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in writing RMS errors, no.', i
       END DO

! ----- CLOSE THE FILE

       stat(1) = NF_CLOSE(fileid)

       IF (stat(1) /= NF_NOERR) THEN
          WRITE(*, *) 'NetCDF error in closing NetCDF file'
       END IF
    END IF pe0e

! ----- Clean up

    DEALLOCATE(temp_arr2d)

  END SUBROUTINE write_nc_ice_pdaf

!===================================================================
  SUBROUTINE broadcast3D_pdaf(arr3D, arr3Dglobal)
  ! Makes nodal information available to all PE
  ! arr3d is any array like TF or SF of local size.
  ! arr3Dglobal is an array of nod3D size which 
  ! should be allocated before calling this routine.
  ! It will be filled with information on other PE in 
  ! natural numbering. The routine can be used to organize
  ! output in the same way as in global memory setup  
  !
  ! Coded by Sergey Danilov
  ! Reviewed by ??
  !===================================================================

    USE g_PARFE
    USE mod_parallel_pdaf, ONLY: comm_filter_fesom, mype_filter_fesom
    USE o_MESH
    USE o_ELEMENTS

    IMPLICIT NONE

    REAL(kind=8), INTENT(inout) ::  arr3D(myDim_nod3D+eDim_nod3D)
    REAL(kind=8), INTENT(out) ::  arr3Dglobal(nod3D)

    INTEGER :: ireals
    INTEGER      ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
    INTEGER, ALLOCATABLE, DIMENSION(:) ::  isendbuf, irecvbuf

    REAL(kind=8), ALLOCATABLE, DIMENSION(:) ::  sendbuf, recvbuf

    CALL MPI_Barrier(COMM_FILTER_FESOM, MPIERR)

    IF ( mype_filter_fesom == 0 ) THEN
       IF (npes>1) THEN
          arr3Dglobal(myList_nod3D(1:myDim_nod3D))=arr3D(1:myDim_nod3D)
       END IF
       DO  n = 1, npes-1

          CALL MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
               0, COMM_FILTER_FESOM, status, MPIerr )
          sender = status(MPI_SOURCE)
          ALLOCATE( recvbuf(1:nTS), irecvbuf(1:nTS) )
          CALL MPI_RECV( irecvbuf(1), nTS, MPI_INTEGER, sender, &
               1, COMM_FILTER_FESOM, status, MPIerr )
          CALL MPI_RECV( recvbuf(1), nTS, MPI_DOUBLE_PRECISION, sender, &
               2, COMM_FILTER_FESOM, status, MPIerr )
          
          DO i = 1, nTS
             arr3Dglobal(irecvbuf(i)) = recvbuf(i)
          ENDDO
          DEALLOCATE( recvbuf, irecvbuf )
          
       ENDDO
       
    ELSE

       ALLOCATE( sendbuf(1:myDim_nod3D), isendbuf(1:myDim_nod3D) )
       DO n = 1, myDim_nod3D
          isendbuf(n) = myList_nod3D(n)
          sendbuf(n)  = arr3D(n)
       ENDDO
       CALL MPI_SEND( myDim_nod3D, 1, MPI_INTEGER, 0, 0, COMM_FILTER_FESOM, MPIerr )
       CALL MPI_SEND( isendbuf(1), myDim_nod3D, MPI_INTEGER, 0, 1, &
            COMM_FILTER_FESOM, MPIerr )
       CALL MPI_SEND( sendbuf(1), myDim_nod3D, MPI_DOUBLE_PRECISION, &
            0, 2, COMM_FILTER_FESOM, MPIerr )
       DEALLOCATE( sendbuf, isendbuf )

    ENDIF

    CALL MPI_BCAST( arr3Dglobal, nod3d, MPI_DOUBLE_PRECISION, 0, &
         COMM_FILTER_FESOM, MPIerr)

  END SUBROUTINE broadcast3D_pdaf
!
!===================================================================
!
  SUBROUTINE broadcast2D_pdaf(arr2D, arr2Dglobal)
  ! Makes nodal information available to all PE 
  ! As the preceeding routine, but for 2D arrays
  !
  ! Coded by Sergey Danilov
  ! Reviewed by ??
  !===================================================================

    USE g_PARFE
    USE mod_parallel_pdaf, ONLY: comm_filter_fesom, mype_filter_fesom
    USE o_MESH
    USE o_ELEMENTS

    IMPLICIT NONE

    INTEGER :: ireals
    INTEGER      ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
    INTEGER, ALLOCATABLE, DIMENSION(:) ::  isendbuf, irecvbuf

    REAL(kind=8) ::  arr2D(myDim_nod2D+eDim_nod2D)
    REAL(kind=8) ::  arr2Dglobal(nod2D)
    REAL(kind=8), ALLOCATABLE, DIMENSION(:) ::  sendbuf, recvbuf

    CALL MPI_Barrier(COMM_FILTER_FESOM, MPIERR)
    IF ( mype_filter_fesom == 0 ) THEN
       IF (npes>1) THEN
          arr2Dglobal(myList_nod2D(1:myDim_nod2D))=arr2D(1:myDim_nod2D)
       END IF
       DO  n = 1, npes-1

          CALL MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
               0, COMM_FILTER_FESOM, status, MPIerr )
          sender = status(MPI_SOURCE)
          ALLOCATE( recvbuf(1:nTS), irecvbuf(1:nTS) )
          CALL MPI_RECV( irecvbuf(1), nTS, MPI_INTEGER, sender, &
               1, COMM_FILTER_FESOM, status, MPIerr )
          CALL MPI_RECV( recvbuf(1), nTS, MPI_DOUBLE_PRECISION, sender, &
               2, COMM_FILTER_FESOM, status, MPIerr )

          DO i = 1, nTS
             arr2Dglobal(irecvbuf(i)) = recvbuf(i)
          ENDDO
          DEALLOCATE( recvbuf, irecvbuf )

       ENDDO

    ELSE

       ALLOCATE( sendbuf(1:myDim_nod2D), isendbuf(1:myDim_nod2D) )
       DO n = 1, myDim_nod2D
          isendbuf(n) = myList_nod2D(n)
          sendbuf(n)  = arr2D(n)
       ENDDO
       CALL MPI_SEND( myDim_nod2D, 1, MPI_INTEGER, 0, 0, COMM_FILTER_FESOM, MPIerr )
       CALL MPI_SEND( isendbuf(1), myDim_nod2D, MPI_INTEGER, 0, 1, &
            COMM_FILTER_FESOM, MPIerr )
       CALL MPI_SEND( sendbuf(1), myDim_nod2D, MPI_DOUBLE_PRECISION, &
            0, 2, COMM_FILTER_FESOM, MPIerr )
       DEALLOCATE( sendbuf, isendbuf )

    ENDIF

    CALL MPI_BCAST( arr2Dglobal, nod2d, MPI_DOUBLE_PRECISION, 0, &
         COMM_FILTER_FESOM, MPIerr)

  END SUBROUTINE broadcast2D_pdaf

END MODULE output_pdaf

