!$Id: output_netcdf_pdaf.F90 2096 2019-08-08 15:31:33Z lnerger $
!BOP
!
! !MODULE:
MODULE output_pdaf

! !DESCRIPTION: 
! This modules provides routines to initialize
! NetCDF output files for FEOM and to write
! output into the files.
!
! !REVISION HISTORY:
! 2012-03 - Lars Nerger - Initial code based on output_netcdf module of pFEOM
! Later revisions - see SVN log
!
! !USES:
  IMPLICIT NONE
  SAVE
  PUBLIC

! !PUBLIC DATA MEMBERS:
  CHARACTER(len=100) :: str_daspec='DA'        ! String to identify assimilation experiment
  CHARACTER(len=1)   :: prec_nc = 's'          ! Precision of NetCDF output
                                               ! (s) single, (d) double
  INTEGER :: write_pos_da = 1                  ! Counter for next time slice to be written
  INTEGER :: write_pos_da_ens
  LOGICAL :: write_da = .true.                 ! Whether to write output file from assimilation
  LOGICAL :: write_ens = .true.                ! Whether to write output file for each individual ensemble member
!EOP

! Private variables
  LOGICAL, PRIVATE :: debugoutput=.FALSE.   ! Write output for debugging
                                            ! (file contains only last writing)
  INTEGER :: nf_prec                        ! Precision of NetCDF output

CONTAINS
!BOP
!
! !ROUTINE: init_output_pdaf - Initialize NetCDf output file
!
! !INTERFACE: 
  SUBROUTINE init_output_pdaf(dim_lag, writepe)

! !USES:
    USE g_clock, &
         ONLY: yearnew, yearold
    USE mod_assim_pdaf, &
         ONLY: dim_ens
 
    IMPLICIT NONE

! !ARGUMENTS:    
    INTEGER, INTENT(in) :: dim_lag              ! Smoother lag
    LOGICAL, INTENT(in) :: writepe
!EOP


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
!BOP
!
! !ROUTINE: init_ncfile_oce_pdaf - Initialize NetCDf output file for ocean fields
!
! !INTERFACE: 

  SUBROUTINE init_ncfile_oce_pdaf(dim_lag, writepe)

! !USES:
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

! !ARGUMENTS:    
    INTEGER, INTENT(in) :: dim_lag              ! Smoother lag
    LOGICAL, INTENT(in) :: writepe
!EOP 
    
! Local variables
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
!BOP
!
! !ROUTINE: init_ncfile_oce_pdaf_ens - Initialize NetCDf output file for ocean fields for each ensemble member
!
! !INTERFACE: 

  SUBROUTINE init_ncfile_oce_pdaf_ens(dim_lag, writepe, dim_ens)

! !USES:
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

! !ARGUMENTS:    
    INTEGER, INTENT(in) :: dim_lag              ! Smoother lag
    LOGICAL, INTENT(in) :: writepe
    INTEGER, INTENT(IN) :: dim_ens
!EOP 

! Local variables
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
!BOP
!
! !ROUTINE: init_ncfile_ice_pdaf - Initialize NetCDf output file for ice fields
!
! !INTERFACE: 

  SUBROUTINE init_ncfile_ice_pdaf(dim_lag, writepe)

! !USES:
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

! !ARGUMENTS:    
    INTEGER, INTENT(in) :: dim_lag              ! Smoother lag
    LOGICAL, INTENT(in) :: writepe

!EOP 
    
! Local variables
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
!BOP
!
! !ROUTINE: write_netcdf_pdaf - Write global fields into NetCDF files
!
! !INTERFACE: 
  SUBROUTINE write_netcdf_pdaf(writetype, write_pos_da, iteration, dim, state_l, &
       nfields, rmse, writepe)

! ! USES:
    IMPLICIT NONE

! !ARGUMENTS:    
    CHARACTER(len=1), INTENT(in) :: writetype     ! Write (i) initial, (a) assimilated, (f) forecast fields
    INTEGER, INTENT(in) :: write_pos_da         ! Write position
    INTEGER, INTENT(in) :: iteration              ! Current model time step
    INTEGER, INTENT(in) :: dim                    ! Size of state vector
    REAL, INTENT(in) :: state_l(dim)              ! State vector
    INTEGER, INTENT(in) :: nfields                ! number of fields in state vector
    REAL, INTENT(in) :: rmse(nfields)             ! Array of RMS errors
    LOGICAL, INTENT(in) :: writepe
!EOP

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
!BOP
!
! !ROUTINE: write_netcdf_pdaf_ens - Write global fields into NetCDF files
!
! !INTERFACE: 
  SUBROUTINE write_netcdf_pdaf_ens(writetype, write_pos_da, iteration, dim, state_l, &
       nfields, rmse, writepe, dim_ens)

! ! USES:
    IMPLICIT NONE

! !ARGUMENTS:    
    CHARACTER(len=1), INTENT(in) :: writetype     ! Write (i) initial, (a) assimilated, (f) forecast fields
    INTEGER, INTENT(in) :: write_pos_da         ! Write position
    INTEGER, INTENT(in) :: iteration              ! Current model time step
    INTEGER, INTENT(in) :: dim                    ! Size of state vector
    REAL, INTENT(in) :: state_l(dim)              ! State vector
    INTEGER, INTENT(in) :: nfields                ! number of fields in state vector
    REAL, INTENT(in) :: rmse(nfields)             ! Array of RMS errors
    LOGICAL, INTENT(in) :: writepe
    INTEGER, INTENT(in) :: dim_ens          
!EOP


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
!BOP
!
! !ROUTINE: write_nc_oce_pdaf - Write global ocean fields into NetCDF file
!
! !INTERFACE: 
  SUBROUTINE write_nc_oce_pdaf(writetype, write_pos_da, iteration, dim, state_l, &
       nfields, rmse, writepe)

! !USES:
    USE g_config, &
         ONLY: runid, ResultPath
    use o_mesh
    use g_clock
    USE g_parfe, &
         ONLY: myDim_nod2D, myDim_nod3D, eDim_nod3D
    USE o_array, &
         ONLY: uf, ssh, tracer, w
    USE mod_assim_pdaf, &
         ONLY: offset, istep_asml, loc_radius

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
   
  real(kind=8)              :: sec_in_year
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
       stat(s) = NF_INQ_VARID(fileid, "loc_radius",       VarId_locradius) 
       s = s + 1

       DO i = 1, s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error inquiring variable IDs, no.', i
       END DO

! ----- WRITE VARIABLES
       s = 1
    
       IF (.not. DEBUGOUTPUT) THEN
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
       ssh(i) = state_l(i + offset(1))
    END DO
    CALL broadcast2d_pdaf(ssh, temp_arr2D)
    pe0a: IF (writepe) THEN
       IF (.not. DEBUGOUTPUT) THEN
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


    !----- Localization radius

       CALL broadcast2d_pdaf(loc_radius, temp_arr2D)
       pe0z1: IF (writepe) THEN
          IF (.not. DEBUGOUTPUT) THEN
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
       tracer(i, 1) = state_l(i + offset(5))
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
       tracer(i, 2) = state_l(i + offset(6))
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
       uf(i) = state_l(i + offset(2))
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
       uf(i + myDim_nod3D + eDim_nod3D) = state_l(i + offset(3))
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
       w(i) = state_l(i + offset(4))
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

       IF (.not. DEBUGOUTPUT) THEN
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
!BOP
!
! !ROUTINE: write_nc_oce_pdaf_ens - Write global ocean fields into NetCDF file
!
! !INTERFACE: 
  SUBROUTINE write_nc_oce_pdaf_ens(writetype, write_pos_da, iteration, dim, state_l, &
       nfields, rmse, writepe, dim_ens)

! !USES:
    USE g_config, &
         ONLY: runid, ResultPath
    use o_mesh
    use g_clock
    USE g_parfe, &
         ONLY: myDim_nod2D, myDim_nod3D, eDim_nod3D
    USE o_array, &
         ONLY: uf, ssh, tracer, w
    USE mod_assim_pdaf, &
         ONLY: offset, istep_asml, loc_radius

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

! !ARGUMENTS:    
    CHARACTER(len=1), INTENT(in) :: writetype     ! Write (i) initial, (a) assimilated, (f) forecast fields
    INTEGER, INTENT(in) :: write_pos_da           ! Write position
    INTEGER, INTENT(in) :: iteration              ! Current model time step
    INTEGER, INTENT(in) :: dim                    ! Size of state vector
    INTEGER, INTENT(in) :: dim_ens
    REAL, INTENT(in) :: state_l(dim, dim_ens)              ! State vector
    INTEGER, INTENT(in) :: nfields                ! number of fields in state vector
    REAL, INTENT(in) :: rmse(nfields)             ! Array of RMS errors
    LOGICAL, INTENT(in) :: writepe
!EOP

! Local variables
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

  real(kind=8)              :: sec_in_year
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
       stat(s) = NF_INQ_VARID(fileid, "loc_radius",       VarId_locradius)
       s = s + 1

       DO i = 1, s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error inquiring variable IDs, no.', i
       END DO

! ----- WRITE VARIABLES
       s = 1

       IF (.not. DEBUGOUTPUT) THEN
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
       ssh_temp(i) = state_l(i + offset(1), member)
    END DO
    CALL broadcast2d_pdaf(ssh_temp, temp_arr2D)
    pe0a: IF (writepe) THEN
       IF (.not. DEBUGOUTPUT) THEN
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

    !----- Localization radius

       CALL broadcast2d_pdaf(loc_radius, temp_arr2D)
       pe0z1: IF (writepe) THEN
          IF (.not. DEBUGOUTPUT) THEN
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
       tracer_temp(i, 1) = state_l(i + offset(5), member)
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
       tracer_temp(i, 2) = state_l(i + offset(6), member)
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
       uf_temp(i) = state_l(i + offset(2), member)
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
       uf_temp(i + myDim_nod3D + eDim_nod3D) = state_l(i + offset(3), member)
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
       w_temp(i) = state_l(i + offset(4), member)
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

       IF (.not. DEBUGOUTPUT) THEN
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
!BOP
!
! !ROUTINE: write_nc_ice_pdaf - Write global icean fields into NetCDF file
!
! !INTERFACE: 
  SUBROUTINE write_nc_ice_pdaf(writetype, write_pos_da, iteration, dim, state_l, &
       nfields, rmse, writepe)

! !USES:
    USE g_config, &
         ONLY: runid, ResultPath
    use o_mesh
    use g_clock
    USE g_parfe, &
         ONLY: myDim_nod2D
    USE i_array, &
         ONLY: a_ice, m_ice, m_snow, u_ice, v_ice
    USE mod_assim_pdaf, &
         ONLY: offset, istep_asml

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
   
    real(kind=8)              :: sec_in_year
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
    
       IF (.not. DEBUGOUTPUT) THEN
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
       a_ice(i) = state_l(i + offset(7))
    END DO
    CALL broadcast2d_pdaf(a_ice, temp_arr2D)
    pe0a: IF (writepe) THEN
       IF (.not. DEBUGOUTPUT) THEN
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
       m_ice(i) = state_l(i + offset(8))
    END DO
    CALL broadcast2d_pdaf(m_ice, temp_arr2D)
    pe0b: IF (writepe) THEN
       IF (.not. DEBUGOUTPUT) THEN
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
       m_snow(i) = state_l(i + offset(9))
    END DO
    CALL broadcast2d_pdaf(m_snow, temp_arr2D)
    pe0c: IF (writepe) THEN
       IF (.not. DEBUGOUTPUT) THEN
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
       u_ice(i) = state_l(i + offset(10))
    END DO
    CALL broadcast2d_pdaf(u_ice, temp_arr2D)
    pe0d: IF (writepe) THEN
       IF (.not. DEBUGOUTPUT) THEN
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
       v_ice(i) = state_l(i + offset(11))
    END DO
    CALL broadcast2d_pdaf(v_ice, temp_arr2D)
    pe0e: IF (writepe) THEN
       IF (.not. DEBUGOUTPUT) THEN
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
    
       IF (.not. DEBUGOUTPUT) THEN
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
  subroutine broadcast3D_pdaf(arr3D, arr3Dglobal)
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

    use g_PARFE
    use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
    use o_MESH
    use o_ELEMENTS

    implicit none

    real(kind=8), intent(inout) ::  arr3D(myDim_nod3D+eDim_nod3D)
    real(kind=8), intent(out) ::  arr3Dglobal(nod3D)

    integer :: ireals
    integer      ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
    integer, allocatable, dimension(:) ::  isendbuf, irecvbuf

    real(kind=8), allocatable, dimension(:) ::  sendbuf, recvbuf

    call MPI_Barrier(COMM_FILTER, MPIERR)

    if ( mype_filter == 0 ) then
       if (npes>1) then
          arr3Dglobal(myList_nod3D(1:myDim_nod3D))=arr3D(1:myDim_nod3D)
       end if
       do  n = 1, npes-1

          call MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
               0, COMM_FILTER, status, MPIerr )
          sender = status(MPI_SOURCE)
          allocate( recvbuf(1:nTS), irecvbuf(1:nTS) )
          call MPI_RECV( irecvbuf(1), nTS, MPI_INTEGER, sender, &
               1, COMM_FILTER, status, MPIerr )
          call MPI_RECV( recvbuf(1), nTS, MPI_DOUBLE_PRECISION, sender, &
               2, COMM_FILTER, status, MPIerr )
          
          do i = 1, nTS
             arr3Dglobal(irecvbuf(i)) = recvbuf(i)
          enddo
          deallocate( recvbuf, irecvbuf )
          
       enddo
       
    else

       allocate( sendbuf(1:myDim_nod3D), isendbuf(1:myDim_nod3D) )
       do n = 1, myDim_nod3D
          isendbuf(n) = myList_nod3D(n)
          sendbuf(n)  = arr3D(n)
       enddo
       call MPI_SEND( myDim_nod3D, 1, MPI_INTEGER, 0, 0, COMM_FILTER, MPIerr )
       call MPI_SEND( isendbuf(1), myDim_nod3D, MPI_INTEGER, 0, 1, &
            COMM_FILTER, MPIerr )
       call MPI_SEND( sendbuf(1), myDim_nod3D, MPI_DOUBLE_PRECISION, &
            0, 2, COMM_FILTER, MPIerr )
       deallocate( sendbuf, isendbuf )

    endif

    call MPI_BCAST( arr3Dglobal, nod3d, MPI_DOUBLE_PRECISION, 0, &
         COMM_FILTER, MPIerr)

  end subroutine broadcast3D_pdaf
!
!===================================================================
!
  subroutine broadcast2D_pdaf(arr2D, arr2Dglobal)
  ! Makes nodal information available to all PE 
  ! As the preceeding routine, but for 2D arrays
  !
  ! Coded by Sergey Danilov
  ! Reviewed by ??
  !===================================================================

    use g_PARFE
    use mod_parallel_pdaf, ONLY: comm_filter, mype_filter
    use o_MESH
    use o_ELEMENTS

    implicit none

    integer :: ireals
    integer      ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
    integer, allocatable, dimension(:) ::  isendbuf, irecvbuf

    real(kind=8) ::  arr2D(myDim_nod2D+eDim_nod2D)
    real(kind=8) ::  arr2Dglobal(nod2D)
    real(kind=8), allocatable, dimension(:) ::  sendbuf, recvbuf

    call MPI_Barrier(COMM_FILTER, MPIERR)
    if ( mype_filter == 0 ) then
       if (npes>1) then
          arr2Dglobal(myList_nod2D(1:myDim_nod2D))=arr2D(1:myDim_nod2D)
       end if
       do  n = 1, npes-1

          call MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
               0, COMM_FILTER, status, MPIerr )
          sender = status(MPI_SOURCE)
          allocate( recvbuf(1:nTS), irecvbuf(1:nTS) )
          call MPI_RECV( irecvbuf(1), nTS, MPI_INTEGER, sender, &
               1, COMM_FILTER, status, MPIerr )
          call MPI_RECV( recvbuf(1), nTS, MPI_DOUBLE_PRECISION, sender, &
               2, COMM_FILTER, status, MPIerr )

          do i = 1, nTS
             arr2Dglobal(irecvbuf(i)) = recvbuf(i)
          enddo
          deallocate( recvbuf, irecvbuf )

       enddo

    else

       allocate( sendbuf(1:myDim_nod2D), isendbuf(1:myDim_nod2D) )
       do n = 1, myDim_nod2D
          isendbuf(n) = myList_nod2D(n)
          sendbuf(n)  = arr2D(n)
       enddo
       call MPI_SEND( myDim_nod2D, 1, MPI_INTEGER, 0, 0, COMM_FILTER, MPIerr )
       call MPI_SEND( isendbuf(1), myDim_nod2D, MPI_INTEGER, 0, 1, &
            COMM_FILTER, MPIerr )
       call MPI_SEND( sendbuf(1), myDim_nod2D, MPI_DOUBLE_PRECISION, &
            0, 2, COMM_FILTER, MPIerr )
       deallocate( sendbuf, isendbuf )

    endif

    call MPI_BCAST( arr2Dglobal, nod2d, MPI_DOUBLE_PRECISION, 0, &
         COMM_FILTER, MPIerr)

  end subroutine broadcast2D_pdaf

  SUBROUTINE init_write_profile(write_var)

    USE g_config, &
         ONLY: runid, ResultPath
    USE g_clock, &
         ONLY: cyearnew

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

    ! !ARGUMENTS:    
      CHARACTER(len=1), INTENT(in) :: write_var     ! Write (t) temperature, (s) salinity

    ! Local variables
      CHARACTER(200) :: filename                    ! Full name of output file
      INTEGER :: stat(100)                          ! auxiliary: status array
      CHARACTER(len=100) :: attstr                  ! String to write attributes
      INTEGER :: fileid
      INTEGER :: s, i                               ! Counters
      INTEGER :: dimid_nprof, dimid_time, dimid_nobs, dimid_ndepth   ! IDs of dimensions
      INTEGER :: nlayer
      INTEGER :: id_time, id_lon, id_lat, id_temp, id_sal, id_nprof, id_nobs, id_lonobs, id_latobs,id_tempobs, id_salobs            ! IDs of variables
      INTEGER :: id_node1, id_node2, id_node3, id_node1depth, id_node2depth, id_node3depth                                          ! IDs of nodes
      INTEGER :: dimids(2),dim3ids(3)
      REAL, PARAMETER :: pi=3.14159265358979
      INTEGER :: startv_out(2), countv_out(2)
      

      nlayer = 46

     ! create file and write global attributes
     IF (write_var=='t') THEN
        filename=TRIM(ResultPath)//runid//'_pro_temp.'//cyearnew//'.nc'
        WRITE(*,*) 'Create profile_temp file:',filename
     ELSEIF (write_var=='s') THEN
        filename=TRIM(ResultPath)//runid//'_pro_sal.'//cyearnew//'.nc'
        WRITE(*,*) 'Create profile_sal file:',filename
     END IF
     
     s = 1
     stat(s) = NF_CREATE(filename, NF_netcdf4, fileid)
     s = s + 1
     IF (write_var=='t') THEN
        attstr  = 'EN4 profile_temperature'
     ELSEIF (write_var=='s') THEN
        attstr  = 'EN4 profile_salinity'
     END IF
     stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'title', LEN_TRIM(attstr), TRIM(attstr))
     s = s + 1

     ! define dimensions
     stat(s) = NF_DEF_DIM(fileid, 'n_prof', NF_UNLIMITED, dimid_nprof)
     s = s + 1
     stat(s) = NF_DEF_DIM(fileid, 'n_obs', NF_UNLIMITED, dimid_nobs)
     s = s + 1
     stat(s) = NF_DEF_DIM(fileid, 'time', NF_UNLIMITED, dimid_time)
     s = s + 1
     stat(s) = NF_DEF_DIM(fileid, 'n_depth', nlayer, dimid_ndepth)

     DO i = 1,  s - 1
        IF (stat(i) /= NF_NOERR) &
           WRITE(*, *) 'NetCDF error in profile dimension definitions, no.', i
     END DO

     ! define variables
     s = 1
     stat(s) = NF_DEF_VAR(fileid, 'time', NF_INT, 1, dimid_time, id_time)
     s = s + 1
     stat(s) = NF_DEF_VAR(fileid, 'nprof', NF_INT, 1, dimid_time, id_nprof)
     s = s + 1 
     stat(s) = NF_DEF_VAR(fileid, 'nobs', NF_INT, 1, dimid_time, id_nobs)


     ! lat & lon per profile
     dimids(1) = dimid_nprof
     dimids(2) = dimid_time
     s = s + 1
     stat(s) = NF_DEF_VAR(fileid, 'lon', NF_REAL, 2, dimids, id_lon)
     s = s + 1
     stat(s) = NF_DEF_VAR(fileid, 'lat', NF_REAL, 2, dimids, id_lat)
     s = s + 1
     
     ! lat & lon per obs node
     dimids(1) = dimid_nobs
     dimids(2) = dimid_time
     stat(s) = NF_DEF_VAR(fileid, 'lon_obs', NF_REAL, 2, dimids, id_lonobs)     
     s = s + 1
     stat(s) = NF_DEF_VAR(fileid, 'lat_obs', NF_REAL, 2, dimids, id_latobs)
     s = s + 1

     ! node ids
     stat(s) = NF_DEF_VAR(fileid, 'node1', NF_INT, 2, dimids, id_node1)     
     s = s + 1
     stat(s) = NF_DEF_VAR(fileid, 'node2', NF_INT, 2, dimids, id_node2)
     s = s + 1
     stat(s) = NF_DEF_VAR(fileid, 'node3', NF_INT, 2, dimids, id_node3)
     s = s + 1

     ! obs per nodes 
     IF (write_var=='t') THEN
        stat(s) = NF_DEF_VAR(fileid, 'temp_obs', NF_REAL, 2, dimids, id_tempobs)
     ELSEIF (write_var=='s') THEN
        stat(s) = NF_DEF_VAR(fileid, 'sal_obs', NF_REAL, 2, dimids, id_salobs)
     END IF

     dim3ids(1) = dimid_ndepth
     dim3ids(2) = dimid_nprof
     dim3ids(3) = dimid_time
     s = s + 1  
     IF (write_var=='t') THEN
        stat(s) = NF_DEF_VAR(fileid, 'temp', NF_REAL, 3, dim3ids, id_temp) 
     ELSEIF (write_var=='s') THEN
        stat(s) = NF_DEF_VAR(fileid, 'sal', NF_REAL, 3, dim3ids, id_sal)
     END IF

     s = s + 1
     stat(s) = NF_DEF_VAR(fileid, 'node1_depth', NF_INT, 3, dim3ids, id_node1depth)
     s = s + 1
     stat(s) = NF_DEF_VAR(fileid, 'node2_depth', NF_INT, 3, dim3ids, id_node2depth)
     s = s + 1
     stat(s) = NF_DEF_VAR(fileid, 'node3_depth', NF_INT, 3, dim3ids, id_node3depth)



     DO i = 1, s
       IF (stat(i) /= NF_NOERR) &
           WRITE(*, *) 'NetCDF error in profile variable definitions, no.', i
     END DO
     
     s = 1
     stat(s) = NF_ENDDEF(fileid)
     s = s + 1
     stat(s) = NF_CLOSE(fileid)

     DO i = 1,  s
         IF (stat(i) /= NF_NOERR) &
           WRITE(*, *) 'NetCDF error in closing profile file, no.', i
     END DO

  END SUBROUTINE init_write_profile


  SUBROUTINE write_profile(write_var,n_prof,n_obs,proocoord_n2d,obs_pro,obs_depth,ocoord_n2d,nod3d_g,iteration)

  USE g_config, &
         ONLY: runid, ResultPath
  USE g_clock, &
         ONLY: cyearnew
  USE mod_assim_pdaf, &
         ONLY: istep_asml    

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'



! !ARGUMENTS:    
    CHARACTER(len=1), INTENT(in) :: write_var     ! Write (t) temperature, (s) salinity
    INTEGER, INTENT(in) :: iteration              ! Current model time step
    INTEGER, INTENT(in) :: n_prof                 ! Size of profile
    INTEGER, INTENT(in) :: n_obs                  ! Size of observation valuse 
    REAL, INTENT(in) :: proocoord_n2d(2,n_prof)   ! coordinates array (one variable per profile)
    REAL, INTENT(in) :: obs_pro(n_obs)            ! observation values
    REAL, INTENT(in) :: obs_depth(n_obs)          ! Observation depth
    REAL, INTENT(in) :: ocoord_n2d(2,n_obs)       ! coordiantes array (one value per observation value)
    INTEGER, INTENT(in) :: nod3d_g(3,n_obs)       ! global node index of model nodes (three node index per each observation)
    
! Local variables
    CHARACTER(200) :: filename                    ! Full name of output file
    INTEGER :: stat(100)                          ! auxiliary: status array
    CHARACTER(len=100) :: attstr                  ! String to write attributes
    INTEGER :: fileid
    INTEGER :: s, i, j
    INTEGER :: dimid_nprof, dimid_time, dimid_nobs, dimid_ndepth   ! IDs of dimensions
    INTEGER :: nlayer, n_depth
    INTEGER :: id_time, id_lon, id_lat, id_temp, id_sal, id_nprof, id_lonobs, id_latobs  ! Ids of variables
    INTEGER :: id_tempobs, id_salobs, id_nobs, id_node1, id_node2, id_node3, id_node1depth, id_node2depth, id_node3depth
    INTEGER :: dimids(2), dim3ids(3)
    REAL, PARAMETER :: pi=3.14159265358979
    INTEGER :: startv_out(2), countv_out(2)
    REAL, ALLOCATABLE :: proocoord_n2d_r(:,:),obs_pro_depth(:,:),ocoord_n2d_r(:,:)
    INTEGER :: pos1, nmb                           ! Position index for writing
    INTEGER :: pos1vec(2), nmbvec(2), pos2vec(3), nmb2vec(3)     ! Position index arrays for writing
    REAL, ALLOCATABLE :: depth_layer(:), occord_pro(:,:)
    INTEGER :: curr_prof_idx, curr_layer_idx
    REAL :: curr_lon, curr_lat, pre_location, curr_location
    INTEGER, ALLOCATABLE :: obs_layer(:), node1_depth(:,:),node2_depth(:,:),node3_depth(:,:)

    nlayer = 46
    n_depth = nlayer

   ! Open the file
     IF (write_var=='t') THEN
        filename=TRIM(ResultPath)//runid//'_pro_temp.'//cyearnew//'.nc'
        WRITE(*,*) 'Write profile_temp into netcdf file:',filename
     ELSEIF (write_var=='s') THEN
        filename=TRIM(ResultPath)//runid//'_pro_sal.'//cyearnew//'.nc'
        WRITE(*,*) 'Write profile_sal into netcdf file:',filename
     END IF

       stat(1) = NF_OPEN(filename, NF_WRITE, fileid)

       IF (stat(1) /= NF_NOERR) THEN
        STOP 'nc-file error'
        print *, trim(nf_strerror(stat(1)))
       END IF


! ----- INQUIRE VARIABLE IDs

       s = 1
       stat(s) = NF_INQ_VARID(fileid, "time", id_time)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "lon", id_lon)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "lat", id_lat)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "lon_obs", id_lonobs)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "lat_obs", id_latobs)
       IF (write_var == 't') THEN
          stat(s) = NF_INQ_VARID(fileid, "temp", id_temp)
          s = s + 1
          stat(s) = NF_INQ_VARID(fileid, "temp_obs", id_tempobs)
       ELSEIF (write_var == 's') THEN
          stat(s) = NF_INQ_VARID(fileid, "sal", id_sal)
          s = s + 1 
          stat(s) = NF_INQ_VARID(fileid, "sal_obs", id_salobs)
       END IF
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "nprof", id_nprof)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "nobs", id_nobs)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "node1", id_node1)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "node2", id_node2)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "node3", id_node3)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "node1_depth", id_node1depth)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "node2_depth", id_node2depth)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "node3_depth", id_node3depth)


 ! ---- WRITE VARIABLES ---------------------------------
   IF (.not. DEBUGOUTPUT) THEN
          ! Normal output
          pos1 = write_pos_da
       ELSE
          ! Write keeping only a single time instance
          pos1 = 1
       END IF
       nmb  = 1

  ! Write time
       s = s + 1
       stat(s) = NF_PUT_VARA_INT(fileid,id_time, pos1, nmb, istep_asml)
   
  ! Write nprof
       s = s + 1
       stat(s) = NF_PUT_VARA_INT(fileid, id_nprof, pos1, nmb, n_prof)

  ! Write nobs
       s = s + 1
       stat(s) = NF_PUT_VARA_INT(fileid, id_nobs, pos1, nmb, n_obs)

  ! Write coordinates for each profile
    ALLOCATE(proocoord_n2d_r(2,n_prof))
    DO i = 1, n_prof
      proocoord_n2d_r(1, i) = proocoord_n2d(1, i)*180.0/pi
      proocoord_n2d_r(2, i) = proocoord_n2d(2, i)*180.0/pi
    ENDDO

    IF (.not. DEBUGOUTPUT) THEN
          ! Normal output
          pos1vec = (/ 1, write_pos_da  /)
    ELSE
          ! Write keeping only a single time instance
          pos1vec = (/ 1, 1 /)
    END IF
    nmbvec  = (/ n_prof, 1 /)

    s = s + 1
    stat(s) = NF_PUT_VARA_REAL(fileid, id_lon, pos1vec, nmbvec, REAL(proocoord_n2d_r(1,:), 4))
    s = s + 1
    stat(s) = NF_PUT_VARA_REAL(fileid, id_lat, pos1vec, nmbvec, REAL(proocoord_n2d_r(2,:), 4))
    s = s + 1

  ! Write coordinates for each obs node
    ALLOCATE(ocoord_n2d_r(2,n_obs))
    DO i = 1, n_obs
      ocoord_n2d_r(1, i) = ocoord_n2d(1, i)*180.0/pi
      ocoord_n2d_r(2, i) = ocoord_n2d(2, i)*180.0/pi
    ENDDO

    IF (.not. DEBUGOUTPUT) THEN
          ! Normal output
          pos1vec = (/ 1, write_pos_da  /)
    ELSE
          ! Write keeping only a single time instance
          pos1vec = (/ 1, 1 /)
    END IF
    nmbvec  = (/ n_obs, 1 /)

    s = s + 1
    stat(s) = NF_PUT_VARA_REAL(fileid, id_lonobs, pos1vec, nmbvec, REAL(ocoord_n2d_r(1,:), 4))
    s = s + 1
    stat(s) = NF_PUT_VARA_REAL(fileid, id_latobs, pos1vec, nmbvec, REAL(ocoord_n2d_r(2,:), 4))
    s = s + 1

    ! Write node ids
    stat(s) = NF_PUT_VARA_INT(fileid, id_node1, pos1vec, nmbvec, nod3d_g(1,:))
    s = s + 1
    stat(s) = NF_PUT_VARA_INT(fileid, id_node2, pos1vec, nmbvec, nod3d_g(2,:))
    s = s + 1
    stat(s) = NF_PUT_VARA_INT(fileid, id_node3, pos1vec, nmbvec, nod3d_g(3,:))
    s = s + 1

  ! Write temperature/salinity per obs node
    IF (write_var=='t') THEN
        stat(s) = NF_PUT_VARA_REAL(fileid, id_tempobs, pos1vec, nmbvec, REAL(obs_pro, 4))
    ELSEIF (write_var == 's') THEN
        stat(s) = NF_PUT_VARA_REAL(fileid, id_salobs, pos1vec, nmbvec, REAL(obs_pro, 4))
    END IF
      
   ! Write temperature/salinity per profile
 
    IF (.not. DEBUGOUTPUT) THEN
          ! Normal output
          pos2vec = (/ 1, 1, write_pos_da  /)
    ELSE
          ! Write keeping only a single time instance
          pos2vec = (/ 1, 1, 1 /)
    END IF
    nmb2vec  = (/ n_depth, n_prof, 1 /)
   
      ! Convert depth into layer
      ALLOCATE(depth_layer(nlayer))
      depth_layer(1) = 0.0
      DO i = 2, 11
        depth_layer(i) = depth_layer(i-1) + 10.0
      END DO
      depth_layer(12) = 115.0
      depth_layer(13) = 135.0
      depth_layer(14) = 160.0
      depth_layer(15) = 190.0
      depth_layer(16) = 230.0
      depth_layer(17) = 280.0
      depth_layer(18) = 340.0
      depth_layer(19) = 410.0
      depth_layer(20) = 490.0
      depth_layer(21) = 580.0
      depth_layer(22) = 680.0
      depth_layer(23) = 790.0
      depth_layer(24) = 910.0
      depth_layer(25) = 1040.0
      depth_layer(26) = 1180.0  
      depth_layer(27) = 1330.0
      depth_layer(28) = 1500.0
      depth_layer(29) = 1700.0
      depth_layer(30) = 1920.0
      depth_layer(31) = 2150.0
      depth_layer(32) = 2400.0
      depth_layer(33) = 2650.0
      depth_layer(34) = 2900.0
      depth_layer(35) = 3150.0
      depth_layer(36) = 3400.0
      depth_layer(37) = 3650.0
      depth_layer(38) = 3900.0
      depth_layer(39) = 4150.0
      depth_layer(40) = 4400.0
      depth_layer(41) = 4650.0
      depth_layer(42) = 4900.0
      depth_layer(43) = 5150.0
      depth_layer(44) = 5400.0
      depth_layer(45) = 5650.0
      depth_layer(46) = 5900.0

      ALLOCATE(obs_layer(n_obs))
      DO i = 1, n_obs
        DO j = 1, nlayer
           IF (obs_depth (i) == depth_layer(j)) THEN
                obs_layer(i) = j
           END IF
        END DO
      END DO  

      ! Store the variables for each profile (with 46 depths per profile)
!      ALLOCATE(obs_pro_depth(n_prof,nlayer))
      ALLOCATE(obs_pro_depth(nlayer,n_prof))
      ALLOCATE(node1_depth(n_prof,nlayer))
      ALLOCATE(node2_depth(n_prof,nlayer))
      ALLOCATE(node3_depth(n_prof,nlayer))
      ! Initialize obs_pro_depth
      obs_pro_depth = 99999.0
      node1_depth = 0
      node2_depth = 0
      node3_depth = 0
      ! Loop over dim_obs_f to fill in the obs_pro_depth matrix
      ! For testing purpose, also store the coordinates
      ALLOCATE(occord_pro(2,n_prof))

      ! Define location = lon + lat * 1000 
      pre_location = 0.0
      curr_prof_idx = 0
      DO i = 1, n_obs
        curr_location = ocoord_n2d(1,i) + 1000 * ocoord_n2d(2,i)
        IF (curr_location /= pre_location) THEN
                curr_prof_idx = curr_prof_idx + 1
                pre_location = curr_location
                occord_pro(:,curr_prof_idx) = ocoord_n2d(:,i)
        END IF
        !obs_pro_depth(curr_prof_idx,obs_layer(i)) = obs_pro(i)
        obs_pro_depth(obs_layer(i),curr_prof_idx) = obs_pro(i)
        node1_depth(curr_prof_idx,obs_layer(i)) = nod3d_g(1,i) 
        node2_depth(curr_prof_idx,obs_layer(i)) = nod3d_g(2,i)
        node3_depth(curr_prof_idx,obs_layer(i)) = nod3d_g(3,i)
      END DO

 
    IF (write_var=='t') THEN
        stat(s) = NF_PUT_VARA_REAL(fileid, id_temp, pos2vec, nmb2vec, REAL(obs_pro_depth, 4))
    ELSEIF (write_var == 's') THEN
        stat(s) = NF_PUT_VARA_REAL(fileid, id_sal, pos2vec, nmb2vec, REAL(obs_pro_depth, 4))
    END IF       

    s = s + 1
    stat(s) = NF_PUT_VARA_INT(fileid, id_node1depth, pos2vec, nmb2vec, node1_depth)
    s = s + 1
    stat(s) = NF_PUT_VARA_INT(fileid, id_node2depth, pos2vec, nmb2vec, node2_depth)
    s = s + 1
    stat(s) = NF_PUT_VARA_INT(fileid, id_node3depth, pos2vec, nmb2vec, node3_depth)

    DEALLOCATE(depth_layer, obs_layer, occord_pro)
    DEALLOCATE(proocoord_n2d_r,obs_pro_depth,ocoord_n2d_r,node1_depth,node2_depth,node3_depth)
    
 ! ---- Close file -------------------------------------
    s = 1
    stat(s) = NF_CLOSE(fileid)

    DO i = 1,  s
        IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in writing of temperature profile file, no.', i
    END DO
    END SUBROUTINE write_profile

END MODULE output_pdaf

