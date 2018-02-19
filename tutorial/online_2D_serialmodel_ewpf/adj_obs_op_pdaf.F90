!$Id: obs_op_pdaf.F90 3 2013-09-05 10:28:51Z lnerger $
!BOP
!
! !ROUTINE: adj_obs_op_pdaf --- Implementation of adjoint observation operator
!
! !INTERFACE:
SUBROUTINE adj_obs_op_pdaf(step, dim_obs_p, dim_p, m_state_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filter: EWPF
!
! The routine is called during the analysis step.
! It has to perform the operation of the adjoint
! observation operator acting on an observation vector.
! For domain decomposition, the action is on the
! PE-local sub-domain of the state and has to 
! provide the observed sub-state for the PE-local 
! domain.
!
! The routine is called by all filter processes.
!
! !REVISION HISTORY:
! 2014-05 - Paul Kirchgessner
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: obs_index_p

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step               ! Currrent time step
  INTEGER, INTENT(in) :: dim_p              ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs_p          ! Dimension of observed state
  REAL, INTENT(out)   :: state_p(dim_p)     ! PE-local model state
  REAL, INTENT(in) :: m_state_p(dim_obs_p)  ! PE-local observed state
!EOP

! *** local variables ***
  INTEGER :: i       ! Counter


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

  state_p = 0
  DO i = 1, dim_obs_p
     state_p(obs_index_p(i)) = m_state_p(i) 
  END DO

END SUBROUTINE adj_obs_op_pdaf




