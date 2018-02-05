!$Id: init_obs_l_pdaf.F90 1590 2015-06-12 11:59:35Z lnerger $
!BOP
!
! !ROUTINE: init_obs_l_pdaf --- Initialize local observation vector
!
! !INTERFACE:
SUBROUTINE init_obs_l_pdaf(domain_p, step, dim_obs_l, observation_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the analysis step
! on each local analysis domain in 
! PDAF\_X\_analysis.  It has to initialize 
! the local vector of observations for the 
! current local analysis domain.
!
! The routine is called by all filter processes.
!
! Implementation for the dummy model with domain 
! decomposition. In this variant a local observation 
! domain is used that is defined by the cut-off 
! distance lseik\_range around the current grid
! point that is updated. (See also the variant  
! using a global observation domain.)
!
! !REVISION HISTORY:
! 2005-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel, &
       ONLY: mype_filter
  USE mod_assimilation, &
       ONLY: rms_obs, local_range
  USE mod_model, &
       ONLY: dim_state, dt, local_dims

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain_p      ! Current local analysis domain
  INTEGER, INTENT(in) :: step        ! Current time step
  INTEGER, INTENT(in) :: dim_obs_l   ! Local dimension of observation vector
  REAL, INTENT(out)   :: observation_l(dim_obs_l) ! Local observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_analysis   (as U_init_obs_l)
! Calls: dlarnv (LAPACK)
!EOP


! *** local variables ***
  INTEGER :: i          ! counter
  INTEGER :: ilow, iup  ! Index for domain range for observations
  INTEGER :: domain_g   ! Index of local domain in global grid
  REAL, SAVE, ALLOCATABLE :: observation(:)  ! global observation vector
  REAL, ALLOCATABLE :: obs_errors(:)         ! global array holding obs. errors
  INTEGER, SAVE :: first = 1      ! flag for init of random number seed
  INTEGER, SAVE :: iseed(4)       ! seed array for random number generator


! ******************************
! *** Initialize observation ***
! ******************************

  ! For the dummy model the observation is just the true state
  ! plus some Gaussian error.
  ! We generate the global observation vector which is then
  ! distributed. This is motivated by the consistency with the 
  ! mode-decomposition variant with regard to the random-number
  ! generation.
  IF (domain_p == 1) THEN
     WRITE (*, '(8x, a, i3, a, I7)') &
          '--- PE-Domain:', mype_filter, &
          ' initialize observation at time step', step
  END IF

  IF (domain_p == 1) THEN
     ! Generate a global vector of observations when called
     ! for the first domain

     ! Allocate global observation vector
     IF (first == 1) ALLOCATE(observation(dim_state))

     ! *** Compute global true state
     observation(:) = 1.0
     DO i = 1, step
        observation(:) = observation(:) + 1.0 * dt
     END DO

     ! *** generate array of observation_l errors
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
     CALL dlarnv(3, iseed, dim_state, obs_errors(1 : dim_state))

     ! disturb true state
     DO i = 1, dim_state
        observation(i) = observation(i) + rms_obs * obs_errors(i)
     END DO
     
     DEALLOCATE(obs_errors)

  ENDIF

! *****************************************************
! *** Initialize substate for local analysis domain ***
! *****************************************************

  ! Get domain index in global grid
  domain_g = domain_p
  DO i = 1, mype_filter
     domain_g = domain_g + local_dims(i)
  ENDDO

  ! Get grid index range for local observations
  IF (domain_g > INT(local_range)) THEN
     ilow = domain_g - INT(local_range)
  ELSE
     ilow = 1
  ENDIF
  IF (domain_g + INT(local_range) <= dim_state) THEN
     iup = domain_g + INT(local_range)
  ELSE
     iup = dim_state
  ENDIF

  ! Perform localization
  DO i = ilow, iup
     observation_l(i - ilow + 1) = observation(i)
  END DO


! ********************
! *** Finishing up ***
! ********************
  
END SUBROUTINE init_obs_l_pdaf

