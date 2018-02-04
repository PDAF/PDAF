!$Id: init_obs_full.F90 1631 2016-08-14 07:41:45Z lnerger $
!BOP
!
! !ROUTINE: init_obs_full --- Initialize observation vector
!
! !INTERFACE:
SUBROUTINE init_obs_full(step, dim_obs, observation)

! !DESCRIPTION:
! User-supplied routine for PDAF (LSEIK):
!
! The routine is called in PDAF\_lseik\_update
! before the loop over all local analysis domains
! is entered. It has to provide the full observation 
! vector according to current time step (where 'full' 
! means 'all observations required for the localized 
! analysis on the PE-local domain).  This routine 
! is only used for LSEIK if a globally adaptive 
! forgetting factor is requested, rather than an 
! individual forgetting factor for each analysis 
! domain. This routine has to be implemented 
! consistently with the routines for the full 
! observation dimension and the full observation 
! operator. The forgetting factor will only be 
! globally adaptive, if the full observation vector 
! is the global observation vector.
!
! Version for the Lorenz96 model without parallelization.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: file_obs, delt_obs_file, observation_g, use_obs_mask, obsindx
  USE mod_model, &
       ONLY: dim_state

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

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
  INTEGER :: i, s                ! Counters
  INTEGER :: stat(50)            ! Array for status flag
  INTEGER :: fileid              ! Id of netcdf file
  INTEGER :: id_obs              ! ID for observation
  INTEGER :: pos(2)              ! Position index for writing
  INTEGER :: cnt(2)              ! Count index for writing
  INTEGER, save :: allocflag = 1 ! Whether allocation has already been performed


! ******************************
! *** Initialize observation ***
! ******************************

  IF (allocflag == 1) THEN
     ALLOCATE(observation_g(dim_state))
     allocflag = 0
  END IF

  ! Read observation information from file
  s = 1
  stat(s) = NF_OPEN(TRIM(file_obs), NF_NOWRITE, fileid)

  s = s + 1
  stat(s) = NF_INQ_VARID(fileid, 'obs', id_obs)

  write (*,'(8x,a,i6)') &
       '--- Read full observation at file position', step / delt_obs_file + 1

  pos(2) = step/delt_obs_file + 1
  cnt(2) = 1
  pos(1) = 1
  cnt(1) = dim_state
  s = s + 1
  stat(s) = NF_GET_VARA_DOUBLE(fileid, id_obs, pos, cnt, observation_g)

  s = s + 1
  stat(s) = nf_close(fileid)

  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
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

END SUBROUTINE init_obs_full

