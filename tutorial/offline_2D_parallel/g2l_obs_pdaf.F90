!$Id: g2l_obs_pdaf.F90 1366 2013-04-24 16:25:05Z lnerger $
!BOP
!
! !ROUTINE: g2l_obs_pdaf --- Restrict an obs. vector to local analysis domain
!
! !INTERFACE:
SUBROUTINE g2l_obs_pdaf(domain_p, step, dim_obs_f, dim_obs_l, mstate_f, &
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
! Implementation for the 2D offline example
! with or without parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: obs_index_l

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain_p   ! Current local analysis domain
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

  DO i = 1, dim_obs_l
     mstate_l(i) = mstate_f(obs_index_l(i))
  END DO

END SUBROUTINE g2l_obs_pdaf
