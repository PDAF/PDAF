!$Id: init_obs_l_pdaf.F90 1711 2016-12-20 10:17:58Z lnerger $
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
! distance local\_range around the current grid
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
       ONLY: local_range
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
!EOP


! *** local variables ***
  INTEGER :: i          ! counter
  INTEGER :: ilow, iup  ! Index for domain range for observations
  INTEGER :: domain_g   ! Index of local domain in global grid


! ******************************
! *** Initialize observation ***
! ******************************

  ! We have already initialized the full observation in init_dim_obs_f_pdaf. 
  ! Here, we use this array through the module mod_model. 
  ! The motivation for this approach is that in a real case, one would already
  ! need to read the observation file to count the number of totally available
  ! observations. At this time, one can already initialize the vector of
  ! observatios. 


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

END SUBROUTINE init_obs_l_pdaf

