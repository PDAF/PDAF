!$Id: g2l_obs_pdaf.F90 2136 2019-11-22 18:56:35Z lnerger $
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
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mod_assim_pdaf, &
       ONLY: local_obs_nod2d, offset, ivariance_obs, ivariance_obs_l

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
! Called by: PDAF_lestkf_analysis  (as U_g2l_obs)
! Called by: PDAF_letkf_analysis   (as U_g2l_obs)
!EOP

! *** Local variables ***
  INTEGER :: i  ! Counter


! *******************************************************
! *** Perform localization of some observation vector *** 
! *** to the current local analysis domain.           ***
! *******************************************************

  ! Restrict vector on local observation domain according to index vector
  DO i = 1, dim_obs_l
     mstate_l(i) = mstate_f(local_obs_nod2d(i))
  END DO

  ! Store local inverse observation error variances
  IF (ALLOCATED(ivariance_obs_l)) DEALLOCATE(ivariance_obs_l)
  ALLOCATE(ivariance_obs_l(dim_obs_l))
 
  DO i = 1, dim_obs_l
     ivariance_obs_l(i) = ivariance_obs(local_obs_nod2d(i))
  END DO

END SUBROUTINE g2l_obs_pdaf
