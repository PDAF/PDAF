!$Id: init_dim_obs_f_pdaf.F90 1711 2016-12-20 10:17:58Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_f_pdaf --- Set full dimension of observation
!
! !INTERFACE:
SUBROUTINE init_dim_obs_f_pdaf(step, dim_obs_f)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
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
! The routine is called by each filter process.
!
! Implementation for the dummy model with domain 
! decomposition. For simplicity the full 
! observation vector holds a global state.
!
! !REVISION HISTORY:
! 2005-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &
       ONLY: dim_state, dt, observation
  USE mod_assimilation, &
       ONLY: rms_obs

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step      ! Current time step
  INTEGER, INTENT(out) :: dim_obs_f ! Dimension of full observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_obs)
! Calls: slarnv (LAPACK)
!EOP

! *** local variables ***
  INTEGER :: i                        ! counter
  INTEGER, SAVE :: first = 1          ! flag for allocations
  INTEGER, SAVE :: iseed(4)           ! seed array for random number generator
  REAL, ALLOCATABLE :: obs_errors(:)  ! global array holding obs. errors


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! dimension for PE-local domain
  dim_obs_f = dim_state


! ******************************
! *** Initialize observation ***
! ******************************
  
  ! For the dummy model the observation is just the true state
  ! plus some Gaussian error.
  ! We generate first the global observation vector which is then
  ! distributed. This is motivated by the consistency with the 
  ! mode-decomposition variant with regard to the random-number
  ! generation.

  ! Allocate global observation vector
  IF (first==1) ALLOCATE(observation(dim_state))

  ! *** Compute global true state
  observation(:) = 1.0
  DO i = 1, step
     observation(:) = observation(:) + 1.0 * dt
  END DO

  ! *** generate array of observation errors
  ALLOCATE(obs_errors(dim_state))

  ! Initialized seed for random number routine
  IF (first == 1) THEN
     iseed(1) = 1000
     iseed(2) = 2034
     iseed(3) = 0
     iseed(4) = 3
     first = 2
  END IF

  ! generate random number with Gaussian distribution variance 1
  CALL slarnv(3, iseed, dim_state, obs_errors(1 : dim_state))

  ! disturb true state
  DO i = 1, dim_state
     observation(i) = observation(i) + rms_obs * obs_errors(i)
  END DO

  DEALLOCATE(obs_errors)

END SUBROUTINE init_dim_obs_f_pdaf

