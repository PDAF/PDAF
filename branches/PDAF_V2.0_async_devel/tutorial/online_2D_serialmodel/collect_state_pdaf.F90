!$Id$
!>  Initialize state vector from model fields
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all filters.
!!
!! This subroutine is called during the forecast 
!! phase from PDAF_put_state_X or PDAF_assimilate_X
!! after the propagation of each ensemble member. 
!! The supplied state vector has to be initialized
!! from the model fields (typically via a module). 
!! With parallelization, MPI communication might be 
!! required to initialize state vectors for all 
!! subdomains on the model PEs. 
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! For the 2D tutorial model the state vector and
!! the model field are identical. Hence, state vector
!! directly initialized from the model field by
!! each model PE.
!!
!! __Revision history:__
!! * 2013-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE collect_state_pdaf(dim_p, state_p)

  USE mod_model, &             ! Model variables
       ONLY: nx, ny, field
  USE mod_assimilation, &
       ONLY: async
  USE mod_parallel_pdaf, &
       ONLY: MPI_REAL8, MPI_COMM_WORLD, mpierr, mype_world, filterpe
  USE obs_A_pdafomi, ONLY: thisobs, ostate_A, oens_A, obs_op_A

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  !< local state vector

! *** local variables ***
  INTEGER :: j          ! Counters
  INTEGER :: assim_stat ! Flag whether is called inside the analysis step
  REAL, ALLOCATABLE :: ostate_tmp(:)

! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

  DO j = 1, nx
     state_p(1 + (j-1)*ny : j*ny) = field(1:ny, j)
  END DO

  ! Check whether we are inside the analysis step
  CALL PDAF_get_assim_flag(assim_stat)


! ***********************************************************
! *** Asynchronous DA: Initialize observed state ensemble ***
! ***********************************************************

  ! When collect_state is called with asynchronous DA
  ! inside the analysis, we need to gather the observed ensemble

  IF (async .AND. assim_stat==1) THEN
     CALL MPI_Gather(ostate_A, thisobs%dim_obs_p, MPI_REAL8, &
          oens_A, thisobs%dim_obs_p, MPI_REAL8, 0, MPI_COMM_WORLD, MPIerr)

     ! This call is required to set the OMI-internal pointer to the observations
     ALLOCATE(ostate_tmp(thisobs%dim_obs_f))
     IF (.not. filterpe) CALL obs_op_A(dim_p, thisobs%dim_obs_f, ostate_A, ostate_tmp)
     DEALLOCATE(ostate_tmp)
  END IF

END SUBROUTINE collect_state_pdaf
