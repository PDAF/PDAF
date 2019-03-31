!$Id: init_obs_local.F90 1633 2016-08-15 08:02:27Z lnerger $
!BOP
!
! !ROUTINE: init_obs_local --- Initialize local observation vector
!
! !INTERFACE:
SUBROUTINE init_obs_local(domain, step, dim_obs_l, observation_l)

! !DESCRIPTION:
! User-supplied routine for PDAF (LSEIK):
!
! The routine is called during the analysis step
! on each local analysis domain in 
! PDAF\_lseik\_analysis.  It has to initialize 
! the local vector of observations for the 
! current local analysis domain.
!
! The routine is called by all filter processes.
!
! This variant is for the Lorenz96 model without
! parallelization. A local observation 
! domain is used that is defined by the cut-off 
! distance lseik\_range around the current grid
! point that is updated. (See also the variant  
! using a global observation domain.)
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel, &
       ONLY: mype_filter
  USE mod_assimilation, &
       ONLY: local_range, local_range2, observation_g, use_obs_mask, &
       obsindx_l
  USE mod_model, &
       ONLY: dim_state

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain      ! Current local analysis domain
  INTEGER, INTENT(in) :: step        ! Current time step
  INTEGER, INTENT(in) :: dim_obs_l   ! Local dimension of observation vector
  REAL, INTENT(out)   :: observation_l(dim_obs_l) ! Local observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_analysis   (as U_init_obs_l)
!EOP


! *** local variables ***
  INTEGER :: i          ! counter
  INTEGER :: ilow, iup  ! Index for domain range for observations


! ******************************
! *** Initialize observation ***
! ******************************

  IF (domain == 1) THEN
     WRITE (*, '(8x, a, i3, a, I7)') &
          '--- PE-Domain:', mype_filter, &
          ' initialize observation at time step', step
  END IF

! *****************************************************
! *** Initialize substate for local analysis domain ***
! *****************************************************

  obsgaps: IF (.NOT. use_obs_mask) THEN
     ! Full state is observed

     ! Get grid index range for local observations
     ! and consider periodic boundary conditions
     ilow = domain - INT(local_range)
     iup = domain + INT(local_range2)

     ! Perform localization
     IF (ilow >= 1 .AND. iup <= dim_state) THEN
        ! Observed region completely within observed region
        DO i = ilow, iup
           observation_l(i - ilow + 1) = observation_g(i)
        END DO
     ELSE IF (ilow < 1) THEN
        ! Use lower periodic BC
        DO i = ilow + dim_state, dim_state
           observation_l(i-ilow-dim_state+1) = observation_g(i)
        END DO
        DO i = 1, iup
           observation_l(i-ilow+1) = observation_g(i)
        END DO
     ELSE IF (iup > dim_state) THEN
        ! Use upper periodic BC
        DO i = ilow, dim_state
           observation_l(i-ilow+1) = observation_g(i)
        END DO
        DO i = 1, iup - dim_state
           observation_l(i+dim_state-ilow+1) = observation_g(i)
        END DO
     END IF

  ELSE
     ! Gappy observations

     DO i = 1, dim_obs_l
        observation_l(i) = observation_g(obsindx_l(i))
     END DO
  END IF obsgaps

END SUBROUTINE init_obs_local

