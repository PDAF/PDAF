!$Id: obs_op_f_pdaf.F90 1249 2012-01-27 13:44:54Z lnerger $
!BOP
!
! !ROUTINE: obs_op_f_pdaf --- Implementation of observation operator 
!
! !INTERFACE:
SUBROUTINE obs_op_f_pdaf(step, dim_p, dim_obs_f, state_p, m_state_f)

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
! The routine is called by all filter processes, 
! and the operation has to be performed by each 
! these processes for its PE-local domain.
!
! For the dummy-model and PDAF with domain
! decomposition the state is fully observed. We
! initialize here the global observation vector.
!
! !REVISION HISTORY:
! 2005-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mpi
  USE mod_parallel, &
       ONLY: npes_filter, COMM_filter, MPIerr, mype_filter
  USE mod_model, &
       ONLY: local_dims, dim_state
  USE mod_assimilation, &
       ONLY: dim_obs, local_dims_obs

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step      ! Currrent time step
  INTEGER, INTENT(in) :: dim_p     ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs_f ! Dimension of observed state
  REAL, INTENT(in)  :: state_p(dim_p)         ! PE-local model state
  REAL, INTENT(inout) :: m_state_f(dim_obs_f) ! PE-local observed state

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_obs_op)
!EOP

! *** local variables ***
  INTEGER  :: i                          ! Counter
  INTEGER, ALLOCATABLE :: local_dis(:)   ! Array of displacements for MPI
  REAL, ALLOCATABLE :: mstate_gather(:)  ! Gathered observation vector

! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

  ! For the dummy model the observation operator is the identity
  ! For using all avaible observations also with the local
  ! analysis we have to gather the state vectors 

  ALLOCATE(local_dis(npes_filter))
  ALLOCATE(mstate_gather(dim_obs))

  ! Init array of displacements
  local_dis(1) = 0
  DO i = 2, npes_filter
     local_dis(i) = local_dis(i - 1) + local_dims_obs(i - 1)
  END DO

  CALL MPI_AllGatherV(state_p(1:local_dims_obs(mype_filter+1)), local_dims_obs(mype_filter+1), &
       MPI_DOUBLE_PRECISION, &
       mstate_gather, local_dims_obs, local_dis, MPI_DOUBLE_PRECISION, &
       COMM_filter, MPIerr)
  
  m_state_f = mstate_gather


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(local_dis, mstate_gather)
  
END SUBROUTINE obs_op_f_pdaf
