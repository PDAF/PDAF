!$Id: distribute_state_pdaf.F90 2196 2020-03-26 13:26:59Z lnerger $
!>   Routine to initialize model fields from state vector
!!
!! User-supplied call-back routine for PDAF.
!!
!! During the forecast phase of the filter this
!! subroutine is called from PDAF_get_state
!! supplying a model state which has to be evolved. 
!! The routine has to initialize the fields of the 
!! model (typically available through a module) from 
!! the state vector of PDAF. With parallelization, 
!! MPI communication might be required to 
!! initialize all subdomains on the model PEs.
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! This routine is particular for restarting: It 
!! does nothing, since the ensemble is read from
!! modle restart files.
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
SUBROUTINE distribute_state_restart_pdaf(dim_p, state_p)

  USE mod_parallel_pdaf, &
       ONLY: mype_submodel, task_id

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< process-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  !< local state vector


! **********************
! *** Initialization ***
! **********************

  ! For restarting we use the fields from the model restart files

  if (mype_submodel==0) write (*,*) 'FESOM-PDAF distribute_state_restart_pdaf, task: ', task_id

END SUBROUTINE distribute_state_restart_pdaf
