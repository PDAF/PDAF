!$Id: obs_op_f_pdaf.F90 1864 2017-12-20 19:53:30Z lnerger $
!BOP
!
! !ROUTINE: obs_op_f_pdaf --- Implementation of observation operator
!
! !INTERFACE:
SUBROUTINE obs_op_f_pdaf(step, dim_p, dim_obs_f, state_p, m_state_f)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called in PDAF\_X\_update
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
! Implementation for the 2D online example
! without parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: obs_index_p

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                 ! Current time step
  INTEGER, INTENT(in) :: dim_p                ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs_f            ! Dimension of observed state
  REAL, INTENT(in)    :: state_p(dim_p)       ! PE-local model state
  REAL, INTENT(inout) :: m_state_f(dim_obs_f) ! PE-local observed state

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_obs_op)
! Called by: PDAF_lestkf_update  (as U_obs_op)
! Called by: PDAF_letkf_update   (as U_obs_op)
!EOP

! *** local variables ***
  INTEGER :: i       ! Counter


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

  DO i = 1, dim_obs_f
     m_state_f(i) = state_p(obs_index_p(i))
  END DO

END SUBROUTINE obs_op_f_pdaf
