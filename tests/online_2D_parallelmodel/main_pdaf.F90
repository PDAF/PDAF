!>  Main driver for PDAF tutorial
!!
!! This is a simple model program to demonstrate the
!! fully-parallel implementation of the online mode of PDAF. 
!!
!! The simple model has a 2-dimensional mesh. The initial state
!! is read from a file. The time stepping consists in shifting
!! the field vertically (in the direction of the first array index)
!! by one grid point per time step. A period boundary condition is
!! applied by inserting the field from the upper boundary into the
!! lower one. 
!!
!! In this code variant the coupling to PDAF is completed.
!!
!! __Revision history:__
!! * 2013-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
PROGRAM MAIN

  USE mod_parallel_model, &      ! Model parallelization variables
       ONLY: mype_world, init_parallel, finalize_parallel

  IMPLICIT NONE


! ********************************
! ***      INITIALIZATION      ***
! ********************************

  ! Initialize parallelization
  CALL init_parallel()

#ifdef USE_PDAF
  ! Revise parallelization for ensemble assimilation
  CALL init_parallel_pdaf(0, 1)
#endif

! *** Initial Screen output ***
  IF (mype_world==0) THEN
     WRITE (*, '(/19x, a/)') '+++++ PDAF online mode +++++'
     WRITE (*, '(17x, a)') 'Online data assimilation with PDAF'
     WRITE (*, '(/)')
  END IF

  ! Initialize model
  CALL initialize()

#ifdef USE_PDAF
  ! Initialize PDAF
  CALL init_pdaf()
#endif


! *****************************
! ***      Integration      ***
! *****************************

  ! *** Perform integration
  CALL integrate_pdaf()


! **************************
! ***      Clean up      ***
! **************************

#ifdef USE_PDAF
  ! End parallelization
  CALL finalize_pdaf()
#endif

  CALL finalize_parallel()

END PROGRAM MAIN
