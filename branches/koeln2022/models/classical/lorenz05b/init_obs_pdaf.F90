!$Id: init_obs_pdaf.F90 831 2021-11-06 16:16:30Z lnerger $
!BOP
!
! !ROUTINE: init_obs_pdaf --- Initialize observation vector
!
! !INTERFACE:
SUBROUTINE init_obs_pdaf(step, dim_obs, observation)

! !DESCRIPTION:
! User-supplied routine for PDAF (global filters):
!
! The routine is called during the analysis step. 
! It has to provide the PE-local observation vector 
! for the current time step.
!
! Version for the Lorenz05b model without parallelization.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE netcdf
  USE mod_assimilation, &
       ONLY: file_obs, delt_obs_file, observation_g, use_obs_mask, obsindx
  USE mod_model, &
       ONLY: dim_state

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                 ! Current time step
  INTEGER, INTENT(in) :: dim_obs              ! Dimension of obs. vector
  REAL, INTENT(out)   :: observation(dim_obs) ! Observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_seek_analysis    (as U_init_obs)
! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT
! Called by: PDAF_enkf_obs_ensemble
!EOP

! *** local variables ***
  INTEGER :: i, s               ! Counters
  INTEGER :: stat(50)           ! Array for status flag
  INTEGER :: fileid             ! Id of netcdf file
  INTEGER :: id_obs             ! ID for observation
  INTEGER :: pos(2)             ! Position index for writing
  INTEGER :: cnt(2)             ! Count index for writing
  INTEGER, SAVE :: allocflag = 1 ! Whether allocation has already been performed


! ******************************
! *** Initialize observation ***
! ******************************

  IF (allocflag == 1) THEN
     ALLOCATE(observation_g(dim_state))
     allocflag = 0
  END IF

  ! Read observation information from file
  s = 1
  stat(s) = NF90_OPEN(TRIM(file_obs), NF90_NOWRITE, fileid)

  s = s + 1
  stat(s) = NF90_INQ_VARID(fileid, 'obs', id_obs)

  write (*,'(8x,a,i6)') &
       '--- Read observation at file position', step / delt_obs_file

  pos(2) = step/delt_obs_file
  cnt(2) = 1
  pos(1) = 1
  cnt(1) = dim_state
  s = s + 1
  stat(s) = NF90_GET_VAR(fileid, id_obs, observation_g, start=pos, count=cnt)

  s = s + 1
  stat(s) = NF90_CLOSE(fileid)

  DO i = 1,  s
     IF (stat(i) /= NF90_NOERR) &
          WRITE(*, *) 'NetCDF error in reading observation from file, no.', i
  END DO


! ************************************
! *** Initialize observation array ***
! ************************************
  
  IF (.NOT. use_obs_mask) THEN
     ! Full state is observed
     observation(:) = observation_g(:)
  ELSE
     ! Use gappy observations
     DO i = 1, dim_obs
        observation(i) = observation_g(obsindx(i))
     END DO
  END IF

END SUBROUTINE init_obs_pdaf

