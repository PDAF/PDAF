!$Id: obs_op.F90 1631 2016-08-14 07:41:45Z lnerger $
!BOP
!
! !ROUTINE: obs_op --- Implementation of observation operator
!
! !INTERFACE:
SUBROUTINE obs_op(step, dim, dim_obs, state, m_state)

! !DESCRIPTION:
! User-supplied routine for PDAF (SEEK/SEIK/EnKF):
!
! The routine is called during the analysis step.
! It has to perform the operation of the
! observation operator acting on a state vector.
! For domain decomposition, the action is on the
! PE-local sub-domain of the state and has to 
! provide the observed sub-state for the PE-local 
! domain.
!
! This variant is for the Lorenz96 model without
! parallelization. The state is fully observed.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: use_obs_mask, obsindx

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step           ! Currrent time step
  INTEGER, INTENT(in) :: dim            ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs        ! Dimension of observed state
  REAL, INTENT(in)    :: state(dim)     ! PE-local model state
  REAL, INTENT(out) :: m_state(dim_obs) ! PE-local observed state

! !CALLING SEQUENCE:
! Called by: PDAF_seek_analysis   (as U_obs_op)
! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT
! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
!EOP

! *** Local variables ***
  INTEGER :: i               ! Counter


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

  IF (.NOT. use_obs_mask) THEN
     ! Full state is observed
     m_state(:) = state(:)
  ELSE
     ! Use gappy observations
     DO i = 1, dim_obs
        m_state(i) = state(obsindx(i))
     END DO
  END IF
 
END SUBROUTINE obs_op
