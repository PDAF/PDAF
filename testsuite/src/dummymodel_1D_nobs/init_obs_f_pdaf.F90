!$Id: init_obs_f_pdaf.F90 1249 2012-01-27 13:44:54Z lnerger $
!BOP
!
! !ROUTINE: init_obs_f_pdaf --- Initialize observation vector
!
! !INTERFACE:
SUBROUTINE init_obs_f_pdaf(step, dim_obs_f, observation_f)

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
! The routine is called by all filter processes.
!
! Version for the dummy model with domain 
! decomposition. Here, the state is fully observed.
!
! !REVISION HISTORY:
! 2007-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel, &
       ONLY: mype_filter
  USE mod_assimilation, &
       ONLY: rms_obs, dim_obs
  USE mod_model, &
       ONLY: dim_state, dt

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Currrent time step
  INTEGER, INTENT(in) :: dim_obs_f   ! Dimension of full observation vector
  REAL, INTENT(out)   :: observation_f(dim_obs_f) ! Full observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_obs)
! Calls: dlarnv (LAPACK)
!EOP


! *** local variables ***
  INTEGER :: i  ! counter
  REAL, ALLOCATABLE :: observation(:) ! global observation vector
  REAL, ALLOCATABLE:: obs_errors(:)   ! global array holding obs. errors
  INTEGER, SAVE :: first = 1          ! flag for init of random number seed
  INTEGER, SAVE :: iseed(4)           ! seed array for random number generator
  INTEGER, SAVE :: iseed_save(4)      ! stored seed array
  ! variables and arrays for domain decomposition
  INTEGER :: offset   ! Row-offset according to domain decomposition


! ******************************
! *** Initialize observation ***
! ******************************

  ! For the dummy model the observation is just the true state
  ! plus some Gaussian error.
  ! We generate the global observation vector which is then
  ! distributed. This is motivated by the consistency with the 
  ! mode-decomposition variant with regard to the random-number
  ! generation.
  IF (mype_filter == 0) THEN
     WRITE (*, '(8x, a, i7)') &
          ' initialize full observation at time step', step
  ENDIF

  ! Allocate global observation vector
  ALLOCATE(observation(dim_obs_f))

  ! *** Compute global true state
  observation(:) = 1.0
  DO i = 1, step
     observation(:) = observation(:) + 1.0 * dt
  END DO

  ! *** generate array of observation_f errors
  ALLOCATE(obs_errors(dim_obs_f))

  ! Initialized seed for random number routine
  IF (first == 1) THEN
     iseed(1) = 1000
     iseed(2) = 2034
     iseed(3) = 0
     iseed(4) = 3
     first = 2
  END IF

  ! generate random number with Gaussian distribution variance 1
  CALL dlarnv(3, iseed, dim_obs_f, obs_errors(1 : dim_obs_f))

  ! disturb true state
  DO i = 1, dim_obs_f
     observation(i) = observation(i) + rms_obs * obs_errors(i)
  END DO


! *****************************************************
! *** Initialize substate for local analysis domain ***
! *****************************************************

  ! Here we use simply the full global observation vector
  observation_f(1 : dim_obs_f) = observation(1 : dim_obs_f)


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(obs_errors, observation)

END SUBROUTINE init_obs_f_pdaf

