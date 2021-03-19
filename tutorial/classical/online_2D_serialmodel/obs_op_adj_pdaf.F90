!$Id: obs_op_pdaf.F90 80 2019-02-09 15:30:15Z lnerger $
!BOP
!
! !ROUTINE: obs_op_adj_pdaf --- Implementation of adjoint observation operator
!
! !INTERFACE:
SUBROUTINE obs_op_adj_pdaf(step, dim_p, dim_obs_p, m_state_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in: 3D-Var, ensemble 3D-Var, hybrid 3D-Var
!
! The routine is called during the analysis step.
! It has to perform the operation of the
! adjoint observation operator acting on an
! observation vector.
! For domain decomposition, the action is on the
! PE-local part of observation vector and has to
! provide the sub-state for the PE-local domain.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2021-03 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: obs_index_p

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step               ! Currrent time step
  INTEGER, INTENT(in) :: dim_p              ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs_p          ! PE-local dimension of observed state
  REAL, INTENT(in) :: m_state_p(dim_obs_p)  ! PE-local observed state
  REAL, INTENT(out)   :: state_p(dim_p)     ! PE-local model state

! !CALLING SEQUENCE:
! Called by: PDAF_3dvar_costf_cvt
! Called by: PDAF_3dvar_costf_cg_cvt
!EOP

! *** local variables ***
  INTEGER :: i       ! Counter


! ***************************************************
! *** Perform application of adjoint observation  ***
! *** operator H^T on vector or matrix column     ***
! ***************************************************

  state_p = 0.0
  DO i = 1, dim_obs_p
     state_p(obs_index_p(i)) = m_state_p(i)
  END DO
  
END SUBROUTINE obs_op_adj_pdaf
