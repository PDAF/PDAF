!$Id: obs_op_f_pdaf.F90 1253 2012-01-30 18:56:08Z lnerger $
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
  USE mod_parallel, &
       ONLY: npes_filter, COMM_filter, MPIerr, MPI_DOUBLE_PRECISION
  USE mod_model, &
       ONLY: local_dims, dim_state

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
  ALLOCATE(mstate_gather(dim_state))

  ! Init array of displacements
  local_dis(1) = 0
  DO i = 2, npes_filter
     local_dis(i) = local_dis(i - 1) + local_dims(i - 1)
  END DO

  CALL MPI_AllGatherV(state_p, dim_p, MPI_DOUBLE_PRECISION, &
       mstate_gather, local_dims, local_dis, MPI_DOUBLE_PRECISION, &
       COMM_filter, MPIerr)
  
  m_state_f = mstate_gather


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(local_dis, mstate_gather)
  
END SUBROUTINE obs_op_f_pdaf
