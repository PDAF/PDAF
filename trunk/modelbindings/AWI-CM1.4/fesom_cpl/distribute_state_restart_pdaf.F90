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
       ONLY: mype_submodel, mype_world, task_id
  USE mod_assim_pdaf, &
       ONLY: offset
  USE g_parfe, &
       ONLY: mydim_nod2d, mydim_nod3d, ToDim_nod2D, eDim_nod3D
  USE o_array, &
       ONLY: uf, ssh, tracer, Tsurf, Ssurf, w, ucori_back, vcori_back
  USE i_array, &
       ONLY: a_ice, m_ice, m_snow, u_ice, v_ice
  USE o_mesh, &
       ONLY: nod3D_below_nod2D

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< process-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  !< local state vector

! *** Local variables ***
  INTEGER :: i         ! Counter
  INTEGER :: node      ! Node index
  INTEGER :: cnt_aice, cnt_mice, cnt_msnow  ! Cound number of invalid points for ice


! **********************
! *** Initialization ***
! **********************

  if (mype_submodel==0) write (*,*) 'FESOM-PDAF distribute_state_restart_pdaf, task: ', task_id

END SUBROUTINE distribute_state_restart_pdaf
