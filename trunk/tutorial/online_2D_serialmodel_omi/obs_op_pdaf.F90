!$Id$
!>  Apply observation operator to state vector
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: SEEK/SEIK/EnKF/ETKF/ESTKF
!!
!! The routine is called during the analysis step.
!! It has to perform the operation of the
!! observation operator acting on a state vector.
!!
!! For domain decomposition, the action is on the
!! PE-local sub-domain of the state and has to 
!! provide the observed sub-state for the PE-local 
!! domain.
!!
!! Implementation for the 2D online example
!! with or without parallelization.
!!
!! \date 2013-02 - Lars Nerger - Initial code
!! \date Later revisions - see repository log
!!
SUBROUTINE obs_op_pdaf(step, dim_p, dim_obs_p, state_p, ostate_p)

  USE mod_assimilation, &   ! Assimilation variables
       ONLY: obs_index_p

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                !< Currrent time step
  INTEGER, INTENT(in) :: dim_p               !< PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs_p           !< Dimension of observed state
  REAL, INTENT(in)    :: state_p(dim_p)      !< PE-local model state
  REAL, INTENT(out)   :: ostate_p(dim_obs_p) !< PE-local observed state

! *** local variables ***
  INTEGER :: i       ! Counter


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

  DO i = 1, dim_obs_p
     ostate_p(i) = state_p(obs_index_p(i))
  END DO

END SUBROUTINE obs_op_pdaf
