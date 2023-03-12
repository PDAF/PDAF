!$Id: obs_op_full.F90 1631 2016-08-14 07:41:45Z lnerger $
!BOP
!
! !ROUTINE: obs_op_full --- Observation operator for full domain 
!
! !INTERFACE:
SUBROUTINE obs_op_full(step, dim, dim_obs, state, m_state)

! !DESCRIPTION:
! User-supplied routine for PDAF (LSEIK):
!
! The routine is called in PDAF\_lseik\_update
! before the loop over all local analysis domains
! is entered.  The routine has to perform the 
! operation of the observation operator acting on 
! a state vector.  The full vector of all 
! observations required for the localized analysis
! on the PE-local domain has to be initialized.
! This is usually data on the PE-local domain plus 
! some region surrounding the PE-local domain. 
! This data is gathered by MPI operations. The 
! gathering has to be done here, since in the loop 
! through all local analysis domains, no global
! MPI operations can be performed, because the 
! number of local analysis domains can vary from 
! PE to PE.
!
! This variant is for the Lorenz96 model without
! parallelization. The state is fully observed. We
! initialize here the global observation vector.
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
  INTEGER, INTENT(in) :: step      ! Currrent time step
  INTEGER, INTENT(in) :: dim     ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs ! Dimension of observed state
  REAL, INTENT(in)  :: state(dim)         ! PE-local model state
  REAL, INTENT(inout) :: m_state(dim_obs) ! PE-local observed state

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_obs_op)
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

END SUBROUTINE obs_op_full
