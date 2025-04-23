!> Set number of local analysis domains
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF/LKNETF
!! and the 3DEnVar and hybrid 3DVAR
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

  USE mod_assimilation, &
       ONLY: dim_state_p

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step        !< Current time step
  INTEGER, INTENT(out) :: n_domains_p !< PE-local number of analysis domains


! ************************************
! *** Initialize number of domains ***
! ************************************
  
  ! Here simply the state dimension
  n_domains_p = dim_state_p

END SUBROUTINE init_n_domains_pdaf
