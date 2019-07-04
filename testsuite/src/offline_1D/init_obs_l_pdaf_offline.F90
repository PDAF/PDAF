!$Id: init_obs_l_pdaf_offline.F90 1678 2016-12-11 12:32:53Z lnerger $
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
! PDAF\_lseik\_analysis.  It has to initialize 
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
! The variant for the offline mode is generally identical
! to the online-variant. In this implementation, both
! routines are different, because the online variant
! explicitly simulates on observation by using the time
! stepping of the dummy model.
!
! !REVISION HISTORY:
! 2008-07 - Lars Nerger - Initial code based on online implementation
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel, &
       ONLY: mype_filter
  USE mod_assimilation, &
       ONLY: rms_obs, local_range
  USE mod_model, &
       ONLY: dim_state, local_dims, observation

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
  REAL, ALLOCATABLE :: obs_errors(:)         ! global array holding obs. errors
  INTEGER, SAVE :: first = 1      ! flag for init of random number seed
  INTEGER, SAVE :: iseed(4)       ! seed array for random number generator


! ******************************
! *** Initialize observation ***
! ******************************

  ! For the dummy model the observation is just the true state
  ! plus some Gaussian error.


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

