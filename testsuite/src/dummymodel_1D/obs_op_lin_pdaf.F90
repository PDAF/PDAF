!$Id$
!BOP
!
! !ROUTINE: obs_op_lin_pdaf --- Implementation of linearized observation operator
!
! !INTERFACE:
SUBROUTINE obs_op_lin_pdaf(step, dim_p, dim_obs_p, state_p, m_state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in: 3D-Var, ensemble 3D-Var, hybrid 3D-Var
!
! The routine is called during the analysis step.
! It has to perform the operation of the
! observation operator acting on a state vector.
! For domain decomposition, the action is on the
! PE-local sub-domain of the state and has to 
! provide the observed sub-state for the PE-local 
! domain.
!
! The routine is called by all filter processes.
!
! For the dummy-model and PDAF with domain
! decomposition the state is fully observed.
!
! !REVISION HISTORY:
! 2021-03 - Lars Nerger - Initial code based on obs_op_pdaf
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step               ! Currrent time step
  INTEGER, INTENT(in) :: dim_p              ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs_p          ! Dimension of observed state
  REAL, INTENT(in)    :: state_p(dim_p)     ! PE-local model state
  REAL, INTENT(out) :: m_state_p(dim_obs_p) ! PE-local observed state

! !CALLING SEQUENCE:
! Called by: PDAF_3dvar_costf_cvt
! Called by: PDAF_3dvar_costf_cg_cvt
!EOP


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

  ! For the dummy model the observation operator is the identity
  m_state_p(:) = state_p(:)
  
END SUBROUTINE obs_op_lin_pdaf
