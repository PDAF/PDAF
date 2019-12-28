!$Id$
!>  Implementation of observation operator
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! The routine is called in PDAF_X_update
!! before the loop over all local analysis domains
!! is entered.  The routine has to perform the 
!! operation of the observation operator acting on 
!! a state vector.  The full vector of all 
!! observations required for the localized analysis
!! on the PE-local domain has to be initialized.
!! This is usually data on the PE-local domain plus 
!! some region surrounding the PE-local domain. 
!! This data is gathered by MPI operations. The 
!! gathering has to be done here, since in the loop 
!! through all local analysis domains, no global
!! MPI operations can be performed, because the 
!! number of local analysis domains can vary from 
!! PE to PE. 
!! 
!! Implementation for the 2D online example
!! with or without parallelization.
!!
!! \date 2019-06 - Lars Nerger - Initial code for PDAF_OMI
!! \date Later revisions - see repository log
!!
SUBROUTINE obs_op_f_pdaf(step, dim_p, dim_obs_f, state_p, m_state_f)

  USE interface_pdafomi, &     ! PDAF-OMI interface routine
       ONLY: obs_op_f_pdafomi

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                 !< Current time step
  INTEGER, INTENT(in) :: dim_p                !< PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs_f            !< Dimension of observed state
  REAL, INTENT(in)    :: state_p(dim_p)       !< PE-local model state
  REAL, INTENT(inout) :: m_state_f(dim_obs_f) !< PE-local observed state


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

  ! For PDAF-OMI we just call the interface routine
  ! than contains the observation-specific calls

  CALL obs_op_f_pdafomi(step, dim_p, dim_obs_f, state_p, m_state_f)

END SUBROUTINE obs_op_f_pdaf
