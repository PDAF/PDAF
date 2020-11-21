!$Id$
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
! Generic implementation using index vector 
! ID_LOBS_IN_FOBS
!
! This routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2005-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: id_lobs_in_fobs

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


! *** local variables ***
  INTEGER :: i             ! Counter


! *******************************************************
! *** Perform localization of some observation vector *** 
! *** to the current local analysis domain.           ***
! *******************************************************

  ! Generic implementation with ID_LOBS_IN_FOBS from INIT_DIM_OBS_L_PDAF
  DO i = 1, dim_obs_l
     mstate_l(i) = mstate_f(id_lobs_in_fobs(i))
  END DO


END SUBROUTINE g2l_obs_pdaf
