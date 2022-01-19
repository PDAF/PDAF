!$Id: l2g_state_pdaf.F90 2233 2020-04-03 14:17:38Z lnerger $
!>  Routine to initialize full state from local analysis
!!
!! User-supplied call-back routine for PDAF.
!!
!! The routine is called during the loop over all
!! local analysis domains in the domain local filters
!! after the analysis and ensemble transformation 
!! on a single local analysis domain. It has to 
!! initialize elements of the PE-local full state 
!! vector from the provided analysis state vector 
!! on the local analysis domain.
!!
!! The routine is called by each filter process.
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
SUBROUTINE l2g_state_pdaf(step, domain, dim_l, state_l, dim_p, state_p)

  USE mod_assim_pdaf, &           ! Variables for assimilation
       ONLY: id_lstate_in_pstate

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step           !< Current time step
  INTEGER, INTENT(in) :: domain         !< Current local analysis domain
  INTEGER, INTENT(in) :: dim_l          !< Local state dimension
  INTEGER, INTENT(in) :: dim_p          !< Process-local full state dimension
  REAL, INTENT(in)    :: state_l(dim_l) !< State vector on local analysis domain
  REAL, INTENT(inout) :: state_p(dim_p) !< Process-local full state vector 
  
! *** Local variables *** 
  INTEGER :: i


! **************************************************
! *** Initialize elements of global state vector ***
! **************************************************

  DO i = 1, dim_l
     state_p(id_lstate_in_pstate(i)) = state_l(i)
  END DO

END SUBROUTINE l2g_state_pdaf
