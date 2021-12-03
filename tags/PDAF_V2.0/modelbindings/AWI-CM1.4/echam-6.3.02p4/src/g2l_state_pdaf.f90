!$Id$
!> Routine to restrict a model state to a local analysis domain
!!
!! User-supplied call-back routine for PDAF.
!!
!! The routine is called during the loop over all
!! local analysis domains in the domain local filters
!! before the analysis on a single local analysis 
!! domain. It has to project the full process-local 
!! model state onto the current local analysis 
!! domain.
!!
!! The routine is called by each filter process.
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
SUBROUTINE g2l_state_pdaf(step, domain_p, dim_p, state_p, dim_l, state_l)

  USE mod_assim_pdaf, &           ! Variables for assimilation
       ONLY: id_lstate_in_pstate
  USE mod_assim_atm_pdaf, ONLY: dp

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step             !< Current time step
  INTEGER, INTENT(in) :: domain_p         !< Current local analysis domain
  INTEGER, INTENT(in) :: dim_p            !< PE-local full state dimension
  INTEGER, INTENT(in) :: dim_l            !< Local state dimension
  REAL(dp), INTENT(in)  :: state_p(dim_p) !< PE-local full state vector 
  REAL(dp), INTENT(out) :: state_l(dim_l) !< State vector on local analysis domain

! *** Local variables *** 
  INTEGER :: i


! *************************************
! *** Initialize local state vector ***
! *************************************

  DO i = 1, dim_l
     state_l(i) = state_p(id_lstate_in_pstate(i))
  ENDDO

END SUBROUTINE g2l_state_pdaf
