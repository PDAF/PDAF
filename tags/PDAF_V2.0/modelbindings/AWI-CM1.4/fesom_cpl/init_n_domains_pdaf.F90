!$Id: init_n_domains_pdaf.F90 2339 2020-06-09 17:11:58Z lnerger $
!>  Routine to set number of local analysis domains
!!
!! User-supplied call-back routine for PDAF.
!!
!! The routine is called in PDAF_X_update 
!! or PDAF_assimialte_X at the beginning of the 
!! analysis step before the loop through all local
!! analysis domains. It has to set the number of 
!! local analysis domains for the process-local domain.
!!
!! The routine is called by all filter processes.
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
SUBROUTINE init_n_domains_pdaf(step, n_domains_p)

  USE g_parfe, ONLY: mydim_nod2d

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step        !< Current time step
  INTEGER, INTENT(out) :: n_domains_p !< Process-local number of analysis domains


! ************************************
! *** Initialize number of domains ***
! ************************************

  n_domains_p = myDim_nod2d

END SUBROUTINE init_n_domains_pdaf
