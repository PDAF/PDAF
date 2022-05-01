!$Id: l2g_state_pdaf.F90 1369 2013-04-24 16:38:17Z lnerger $
!>  Initialize full state from local analysis
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! The routine is called during the loop over all
!! local analysis domains in PDAF_X_update 
!! after the analysis and ensemble transformation 
!! on a single local analysis domain. It has to 
!! initialize elements of the PE-local full state 
!! vector from the provided analysis state vector 
!! on the local analysis domain.
!!
!! Generic implementation using index vector 
!! ID_LSTATE_IN_PSTATE.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE l2g_state_pdaf(step, domain_p, dim_l, state_l, dim_p, state_p)

  USE mod_assimilation, &
       ONLY: id_lstate_in_pstate

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step           !< Current time step
  INTEGER, INTENT(in) :: domain_p       !< Current local analysis domain
  INTEGER, INTENT(in) :: dim_l          !< Local state dimension
  INTEGER, INTENT(in) :: dim_p          !< PE-local full state dimension
  REAL, INTENT(in)    :: state_l(dim_l) !< State vector on local analysis domain
  REAL, INTENT(inout) :: state_p(dim_p) !< PE-local full state vector 

! *** local variables ***
  INTEGER :: i                          ! Counter


! **************************************************
! *** Initialize elements of global state vector ***
! **************************************************

  ! Generic initialization using ID_LSTATE_IN_PSTATE set in INIT_DIM_L_PDAF
  DO i = 1, dim_l
     state_p(id_lstate_in_pstate(i)) = state_l(i)
  END DO

END SUBROUTINE l2g_state_pdaf
