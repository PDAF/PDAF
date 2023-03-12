!$Id: init_dim_obs_full.F90 1631 2016-08-14 07:41:45Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_full --- Set number of observations for full domain
!
! !INTERFACE:
SUBROUTINE init_dim_obs_full(step, dim_obs)

! !DESCRIPTION:
! User supplied routine for PDAF (LSEIK):
!
! The routine is called in PDAF\_lseik\_update 
! at the beginning of the analysis step before 
! the loop through all local analysis domains. 
! It has to determine the dimension of the 
! observation vector according to the current 
! time step for all observations required for 
! the analyses in the loop over all local 
! analysis domains on the PE-local state domain.
!
! This variant is for the Lorenz96 model without
! parallelization. For simplicity the full 
! observation vector holds a global state.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: file_obs, delt_obs_file, obsfile_laststep, have_obs, &
       observation_g, use_obs_mask, obs_mask, obsindx, obsindx_l
  USE mod_model, &
       ONLY: dim_state

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step      ! Current time step
  INTEGER, INTENT(out) :: dim_obs ! Dimension of full observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_obs)
!EOP

! *** Local variables ***
  INTEGER :: i, s               ! Counters
  INTEGER :: stat(50)           ! Array for status flag
  INTEGER :: fileid             ! Id of netcdf file
  INTEGER :: id_step            ! ID for step information
  INTEGER :: id_obs             ! File ID for observation
  INTEGER :: nsteps_file        ! Number of time steps in file
  INTEGER :: obs_step1and2(2)   ! Array holding first and second time step index
  INTEGER :: pos(2)             ! Position index for writing
  INTEGER :: cnt(2)             ! Count index for writing
  INTEGER, save :: allocflag = 1  ! Whether allocation has already been performed
  INTEGER, save :: allocflag2 = 1 ! Flag for allocation of obsindx


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Retrieve information on observations from file
  s = 1
  stat(s) = NF_OPEN(TRIM(file_obs), NF_NOWRITE, fileid)

  ! Read number of time steps in file
  s = s + 1
  stat(s) = NF_INQ_DIMID(fileid, 'timesteps', id_step)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(fileid, id_step, nsteps_file)
  
  ! Read time step information
  s = s + 1
  stat(s) = NF_INQ_VARID(fileid, 'step', id_step)

  pos(1) = 1
  cnt(1) = 2
  s = s + 1
  stat(s) = NF_GET_VARA_INT(fileid, id_step, pos, cnt, obs_step1and2)
  
  pos(1) = nsteps_file
  cnt(1) = 1
  s = s + 1
  stat(s) = NF_GET_VARA_INT(fileid, id_step, pos, cnt, obsfile_laststep)

  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in reading observation step information, no.', i
  END DO

  ! Initialize observation interval in file
  delt_obs_file = obs_step1and2(2) - obs_step1and2(1)

  ! observation dimension 
  IF (step <= obsfile_laststep) THEN
     dim_obs = dim_state
  ELSE
     dim_obs = 0
  END IF

  ! Set dim_obs to 0, if we are at the end of a forecast phase without obs.
  IF (.NOT.have_obs) dim_obs = 0

  ! Read global observation
  readobs: IF (dim_obs > 0) THEN

     IF (allocflag == 1) THEN
        ALLOCATE(observation_g(dim_state))
        allocflag = 0
     END IF

     s = 1
     stat(s) = NF_INQ_VARID(fileid, 'obs', id_obs)

     write (*,'(8x,a,i6)') &
          '--- Read observation at file position', step / delt_obs_file + 1

     pos(2) = step/delt_obs_file + 1
     cnt(2) = 1
     pos(1) = 1
     cnt(1) = dim_obs
     s = s + 1
     stat(s) = NF_GET_VARA_DOUBLE(fileid, id_obs, pos, cnt, observation_g(1:dim_obs))

     s = s + 1
     stat(s) = nf_close(fileid)

     DO i = 1,  s
        IF (stat(i) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in reading global observation, no.', i
     END DO

  END IF readobs

  ! For gappy observations initialize index array
  ! and reorder global observation array
  obsgaps: IF (use_obs_mask) THEN

     IF (allocflag2 == 1) THEN
        ALLOCATE(obsindx(dim_state))
        ALLOCATE(obsindx_l(dim_state))
        allocflag2 = 0
     END IF

     obsindx = 0

     s = 1
     DO i=1, dim_state
        IF (obs_mask(i) == 1) THEN
           obsindx(s) = i
           observation_g(s) = observation_g(i)
           s = s + 1
        END IF
     END DO
     dim_obs = s - 1

  END IF obsgaps

END SUBROUTINE init_dim_obs_full

