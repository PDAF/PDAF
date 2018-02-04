!$Id: main_pdaf.F90 1411 2013-09-25 14:04:41Z lnerger $
!BOP
!
! !ROUTINE: main --- Main driver for PDAF testsuite
!
! !INTERFACE:
PROGRAM MAIN

! !DESCRIPTION:
! This is a simple model program to demonstrate the
! fully-parallel implementation of the online mode of PDAF. 
!
! The simple model has a 2-dimensional mesh. The initial state
! is read from a file. The time stepping consists in shifting
! the field vertically (in the direction of the first array index)
! by one grid point per time step. A period boundary condition is
! applied by inserting the field from the upper boundary into the
! lower one. 
!
! !REVISION HISTORY:
! 2013-09 - Lars Nerger - Initial code based on dummy model example
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_model, &
       ONLY: mype_world, init_parallel, finalize_parallel, modelpe

  IMPLICIT NONE
!EOP


! ********************************
! ***      INITIALIZATION      ***
! ********************************


  ! Initialize parallelization
  CALL init_parallel()

  ! Revise parallelization for ensemble assimilation
  CALL init_parallel_pdaf(0, 1)

! *** Initial Screen output ***
  IF (mype_world==0) THEN
     WRITE (*, '(/17x, a/)') '+++++ PDAF tutorial - online mode +++++'
     WRITE (*, '(17x, a)') 'Tutorial: 2D model with parallelization'
     WRITE (*, '(/)')
  END IF

  ! *** Distinguish between model processes and filter processes
  modelpes: IF (modelpe) THEN

     ! Initialize model
     CALL initialize()

     ! Initialize PDAF
     CALL init_pdaf()

     ! Perform integration
     CALL integrate_pdaf()

  ELSE modelpes

     ! Filtering routine is executed only on filter processes
     CALL filtering()

  END IF modelpes


! **************************
! ***      Clean up      ***
! **************************

  CALL finalize_parallel()

#ifdef USE_PDAF
  ! End parallelization
  CALL finalize_pdaf()
#endif


END PROGRAM MAIN
