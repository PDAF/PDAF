!$Id$
!BOP
!
! !MODULE:
MODULE output_netcdf

! !DESCRIPTION: 
! This module provides routines to initialize
! NetCDF output files for the Lorenz63 model 
! and to write output into the files.
!
! !REVISION HISTORY:
! 2019-07 - Lars Nerger - Initial code adopted from Lorenz96
! Later revisions - see SVN log
!
! !USES:
  IMPLICIT NONE
  SAVE
  PUBLIC

! !PUBLIC DATA MEMBERS:
  CHARACTER(len=100) :: file_state = 'state_L63.nc' ! Name of the NetCDF output file
  INTEGER :: delt_write = 1                         ! Output interval in time steps

!EOP

! Private variables
  INTEGER, PRIVATE :: file_pos     ! File position to write to
  INTEGER, PRIVATE :: cnt_steps    ! Count time step for delt_write
  INTEGER, PRIVATE :: fileid       ! Id of netcdf file

CONTAINS
!BOP
!
! !ROUTINE: init_netcdf  --- initialize netcdf output
!
! !INTERFACE:
  SUBROUTINE init_netcdf(step, time, dt, gamma, rho, beta, x0, y0, z0, dim, state)

! !DESCRIPTION:
! This routine initializes the netcdf file 

! !USES:
    USE netcdf

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(IN) :: step        ! Initial time step
    REAL, INTENT(IN)    :: time        ! Initial model time
    REAL, INTENT(IN)    :: dt          ! Size of time step
    REAL, INTENT(IN)    :: gamma       ! Forcing variable of Lorenz63 model
    REAL, INTENT(IN)    :: rho         ! Forcing variable of Lorenz63 model
    REAL, INTENT(IN)    :: beta        ! Forcing variable of Lorenz63 model
    REAL, INTENT(IN)    :: x0          ! Initial state value for variable x
    REAL, INTENT(IN)    :: y0          ! Initial state value for variable y
    REAL, INTENT(IN)    :: z0          ! Initial state value for variable z
    INTEGER, INTENT(IN) :: dim         ! Dimension of model state
    REAL, INTENT(IN)    :: state(dim)  ! Model state
!EOP

! Local variables
    INTEGER :: i, s                      ! Counters
    INTEGER :: dimid_state, dimid_1      ! Dimension IDs
    INTEGER :: dimid_step                ! Dimension ID
    INTEGER :: ID_time, ID_dt, ID_step   ! Variable Ids
    INTEGER :: Id_gamma, Id_rho, Id_beta ! Variable Ids
    INTEGER :: Id_x0, Id_y0, Id_z0       ! Variable Ids
    INTEGER :: Id_state                  ! Variable Ids
    INTEGER :: stat(50)                  ! Array for status flag
    INTEGER :: dimarray(2)               ! Array for dimensions
    INTEGER :: pos(2)                    ! Position index for writing
    INTEGER :: cnt(2)                    ! Count index for writing
    CHARACTER(len=100) :: attstr         ! String to write attributes

! *** Initialize file ***    

! Print screen information
    WRITE (*, '(/1x, a)') 'Initialize NetCDF output'

! Initialize file position
    file_pos = 1

! Initialize counter for output interval
    cnt_steps = 1

! Initialize file and write global atributes

    dooutput: IF (delt_write>0) THEN
       s = 1

       stat(s) = NF90_CREATE(TRIM(file_state), 0, fileid) 
       s = s + 1

       attstr  = 'Lorenz63 model'
       stat(s) = NF90_PUT_ATT(fileid, NF90_GLOBAL, 'title', TRIM(attstr)) 
       s = s + 1

! Define Dimensions

       stat(s) = NF90_DEF_DIM(fileid, 'dim_state', dim, dimid_state)             
       s = s + 1
       stat(s) = NF90_DEF_DIM(fileid, 'one', 1, dimid_1)
       s = s + 1
       stat(s) = NF90_DEF_DIM(fileid, 'timesteps', NF90_UNLIMITED, dimid_step)
       s = s + 1

! Define variables
    
       stat(s) = NF90_DEF_VAR(fileid, 'gamma', NF90_DOUBLE, DimId_1, Id_gamma) 
       s = s + 1
       stat(s) = NF90_DEF_VAR(fileid, 'rho', NF90_DOUBLE, DimId_1, Id_rho) 
       s = s + 1
       stat(s) = NF90_DEF_VAR(fileid, 'beta', NF90_DOUBLE, DimId_1, Id_beta) 
       s = s + 1
       stat(s) = NF90_DEF_VAR(fileid, 'x0', NF90_DOUBLE, DimId_1, Id_x0) 
       s = s + 1
       stat(s) = NF90_DEF_VAR(fileid, 'y0', NF90_DOUBLE, DimId_1, Id_y0) 
       s = s + 1
       stat(s) = NF90_DEF_VAR(fileid, 'z0', NF90_DOUBLE, DimId_1, Id_z0) 
       s = s + 1
       stat(s) = NF90_DEF_VAR(fileid, 'dt', NF90_DOUBLE, DimId_1, Id_dt) 
       s = s + 1
       stat(s) = NF90_DEF_VAR(fileid, 'step', NF90_INT, DimId_step, Id_step) 
       s = s + 1
       stat(s) = NF90_DEF_VAR(fileid, 'time', NF90_DOUBLE, DimId_step, Id_time) 
       s = s + 1

       dimarray(1) = dimid_state
       dimarray(2) = dimid_step
       stat(s) = NF90_DEF_VAR(fileid, 'state', NF90_DOUBLE, dimarray, Id_state) 
       s = s + 1

       stat(s) = NF90_ENDDEF(fileid) 
       s = s + 1
       
! Write initial and constant variables

       stat(s) = NF90_PUT_VAR(fileid, Id_gamma, gamma) 
       s = s + 1
       stat(s) = NF90_PUT_VAR(fileid, Id_rho, rho) 
       s = s + 1
       stat(s) = NF90_PUT_VAR(fileid, Id_beta, beta) 
       s = s + 1
       stat(s) = NF90_PUT_VAR(fileid, Id_x0, x0) 
       s = s + 1
       stat(s) = NF90_PUT_VAR(fileid, Id_y0, y0) 
       s = s + 1
       stat(s) = NF90_PUT_VAR(fileid, Id_z0, z0) 
       s = s + 1
       stat(s) = NF90_PUT_VAR(fileid, Id_dt, dt) 
       s = s + 1

       pos(1) = 1
       cnt(1) = 1
       stat(s) = NF90_PUT_VAR(fileid, Id_time, time, start=pos(1:1))
       s = s + 1

       stat(s) = NF90_PUT_VAR(fileid, Id_step, step, start=pos(1:1))
       s = s + 1

       pos(1) = 1
       pos(2) = file_pos
       cnt(1) = dim
       cnt(2) = 1
       stat(s) = NF90_PUT_VAR(fileid, Id_state, state, pos)
       s = s + 1

       DO i = 1,  s - 1
          IF (stat(i) /= NF90_NOERR) &
               WRITE(*, *) 'NetCDF error in file initialization, no.', i
       END DO

    END IF dooutput

  END SUBROUTINE init_netcdf
!BOP
!
! !ROUTINE: write_netcdf  --- write netcdf output
!
! !INTERFACE:
  SUBROUTINE write_netcdf(step, time, dim, state)

! !DESCRIPTION:
! This routine initializes the netcdf file 

! !USES:
    USE netcdf

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(IN) :: step       ! Current time step
    REAL, INTENT(IN)    :: time       ! Current model time
    INTEGER, INTENT(IN) :: dim        ! Dimension of model state
    REAL, INTENT(IN)    :: state(dim) ! Model state
!EOP

! Local variables
    INTEGER :: i, s                    ! Counters
    INTEGER :: ID_time, Id_step        ! Variable IDs
    INTEGER :: ID_state                ! Variable ID
    INTEGER :: stat(50)                ! Array for status flag
    INTEGER :: pos(2)                  ! Position index for writing
    INTEGER :: cnt(2)                  ! Count index for writing
    LOGICAL :: dowrite                 ! Flag whether to write at the current call

! Check, if we have to write at this time step
    IF (cnt_steps==delt_write) THEN
       dowrite = .TRUE.
       cnt_steps = 1
       file_pos = file_pos + 1
    ELSE
       dowrite = .FALSE.
       cnt_steps = cnt_steps + 1
    END IF

    dooutput: IF (dowrite) THEN

! Inquire variable Ids

       s = 1
       stat(s) = NF90_INQ_VARID(fileid, "time", Id_time) 
       s = s + 1
       stat(s) = NF90_INQ_VARID(fileid, "step", Id_step) 
       s = s + 1
       stat(s) = NF90_INQ_VARID(fileid, "state", Id_state) 
       s = s + 1

! Write variables

       pos(1) = file_pos
       cnt(1) = 1
       stat(s) = NF90_PUT_VAR(fileid, Id_time, time, start=pos(1:1))
       s = s + 1

       pos(1) = file_pos
       cnt(1) = 1
       stat(s) = NF90_PUT_VAR(fileid, Id_step, step, start=pos(1:1))
       s = s + 1

       pos(1) = 1
       pos(2) = file_pos
       cnt(1) = dim
       cnt(2) = 1
       stat(s) = NF90_PUT_VAR(fileid, Id_state, state, pos)
       s = s + 1

       DO i = 1,  s - 1
          IF (stat(i) /= NF90_NOERR) &
               WRITE(*, *) 'NetCDF error in writing output, no.', i
       END DO

    END IF dooutput

  END SUBROUTINE write_netcdf
!BOP
!
! !ROUTINE: close_netcdf  --- close netcdf file
!
! !INTERFACE:
  SUBROUTINE close_netcdf()

! !DESCRIPTION:
! This routine closes the netcdf file 

! !USES:
    USE netcdf

    IMPLICIT NONE
!EOP

! Local variables
    INTEGER :: stat(50)                ! Array for status flag

! Close file

    stat(1) = NF90_CLOSE(fileid)
    IF (stat(1) /= NF90_NOERR) &
         WRITE(*, *) 'NetCDF error in closing output file, no. 1'

  END SUBROUTINE close_netcdf

END MODULE output_netcdf
