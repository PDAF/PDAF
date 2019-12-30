!$Id$
!>  Restrict a model state to a local analysis domain
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! The routine is called during the loop over all
!! local analysis domains in PDAF_X_update
!! before the analysis on a single local analysis 
!! domain.  It has to initialize elements of the 
!! state vector for the local analysis domains from
!! the PE-local full state vector.
!!
!! Implementation for the 2D online example
!! with or without parallelization.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE g2l_state_pdaf(step, domain_p, dim_p, state_p, dim_l, state_l)

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step           !< Current time step
  INTEGER, INTENT(in) :: domain_p       !< Current local analysis domain
  INTEGER, INTENT(in) :: dim_p          !< PE-local full state dimension
  INTEGER, INTENT(in) :: dim_l          !< Local state dimension
  REAL, INTENT(in)    :: state_p(dim_p) !< PE-local full state vector 
  REAL, INTENT(out)   :: state_l(dim_l) !< State vector on local analysis domain


! *************************************
! *** Initialize local state vector ***
! *************************************
  
  ! Here, the local domain is a single grid point and variable
  state_l(1) = state_p(domain_p)

END SUBROUTINE g2l_state_pdaf
