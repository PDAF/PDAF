!$Id: g2l_obs_pdaf.F90 1587 2015-06-11 09:52:11Z lnerger $
!BOP
!
! !ROUTINE: g2l_obs_pdaf --- Restrict an obs. vector to local analysis domain
!
! !INTERFACE:
SUBROUTINE g2l_obs_pdaf(domain, step, dim_obs_f, dim_obs_l, mstate_f, &
     mstate_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the analysis step
! on each of the local analysis domains.
! It has to restrict the full vector of all 
! observations required for the loop of localized 
! analyses on the PE-local domain to the current 
! local analysis domain.
!
! This routine is called by all filter processes.
!
! Implementation for the dummy model with domain
! decomposition. In this variant a local observation 
! domain is used that is defined by the cut-off 
! distance lseik\_range around the current grid
! point that is updated.
!
! !REVISION HISTORY:
! 2007-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &
       ONLY: dim_state, local_dims
  USE mod_assimilation, &
       ONLY: local_range
  USE mod_parallel, &
       ONLY: mype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain     ! Current local analysis domain
  INTEGER, INTENT(in) :: step       ! Current time step
  INTEGER, INTENT(in) :: dim_obs_f  ! Dimension of full PE-local obs. vector
  INTEGER, INTENT(in) :: dim_obs_l  ! Local dimension of observation vector
  REAL, INTENT(in)    :: mstate_f(dim_obs_f)   ! Full PE-local obs. vector
  REAL, INTENT(out)   :: mstate_l(dim_obs_l)   ! Obs. vector on local domain

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_analysis   (as U_g2l_obs)
!EOP


! *** local variables ***
  INTEGER :: i, ilow, iup  ! Counters
  INTEGER :: domain_g      ! Global domain index


! *******************************************************
! *** Perform localization of some observation vector *** 
! *** to the current local analysis domain.           ***
! *******************************************************

  ! Get domain index in global grid
  domain_g = domain
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
     mstate_l(i - ilow + 1) = mstate_f(i)
  END DO

END SUBROUTINE g2l_obs_pdaf
