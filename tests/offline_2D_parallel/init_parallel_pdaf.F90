!>  Interface routine to call initialization of parallelization for PDAF
!!
!! This routine calls the parallelization routine for PDAF, which 
!! initializes the communicators for handling the analysis step
!! of the data assimilation.
!!
!! The parallelization variables returned from the PDAF initialization
!! routines are stored in the module parallel_pdaf_mod so that they can
!! be used in the different user-provided routines.
!!
!! The routine is generic, but has to be part of the user code
!! because it uses the module parallel_pdaf_mod.
!!
!! __Revision history:__
!! * 2004-11 - Lars Nerger - Initial code
!! * 2026-02 - Lars Nerger - Revision for using PDAF3_init_forecast 
!! * Later revisions - see repository log
!!
SUBROUTINE init_parallel_pdaf(screen)

  USE mpi
  USE PDAF, &                     ! Command line parser
       ONLY: PDAF3_init_parallel
  USE mod_parallel_pdaf, &        ! PDAF parallelization variables
       ONLY: n_modeltasks, task_id, mype_world, npes_world, &
       mype_model, npes_model, COMM_model, mype_filter, npes_filter, COMM_filter

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)    :: screen           !< Whether screen information is shown

  ! Set number of model tasks for offline mode
  n_modeltasks = 1

  ! Initialize ensemble parallelization
  CALL PDAF3_init_parallel(screen, 0, 0, 0, n_modeltasks, &
     COMM_model, mype_model, npes_model, &
     COMM_filter, mype_filter, npes_filter, &
     task_id)

  ! Initialize variables for all processes as they are used in some routines
  mype_world = mype_model
  npes_world = npes_model

END SUBROUTINE init_parallel_pdaf
