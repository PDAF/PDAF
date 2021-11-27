!$Id$
!> Set number of local analysis domains
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! The routine is called in PDAF_X_update 
!! at the beginning of the analysis step before 
!! the loop through all local analysis domains. 
!! It has to set the number of local analysis 
!! domains for the PE-local domain.
!!
!! __Revision history:__
!! * 2005-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_n_domains_pdaf(step, n_domains_p)

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step        !< Current time step
  INTEGER, INTENT(out) :: n_domains_p !< PE-local number of analysis domains


! ************************************
! *** Initialize number of domains ***
! ************************************
 
  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE init_n_domains_pdaf.F90: Set number of local analysis domains here!'

!   n_domains_p = ?

  ! SET n_domains_p FOR TESTING
  n_domains_p = 1

END SUBROUTINE init_n_domains_pdaf
