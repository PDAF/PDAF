!$Id: initialize.F90 1589 2015-06-12 11:57:58Z lnerger $
!BOP
!
! !ROUTINE: initialize --- Initialize model
!
! !INTERFACE:
SUBROUTINE initialize()

! !DESCRIPTION:
! Initialization routine for the simple 2D model without
! parallelization of the model.
!
! The routine defines the size of the model grid and
! read the initial state from a file. 
!
! !REVISION HISTORY:
! 2013-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &
       ONLY: nx, ny, field, total_steps

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: main
!EOP

! *** local variables ***
  INTEGER :: i              ! Counters


! **********************
! *** INITIALIZATION ***
! **********************

! *** Model specifications ***
  nx = 36          ! Extent of grid in x-direction
  ny = 18          ! Extent of grid in y-direction
  total_steps = 18 ! Number of time steps to perform

! *** Screen output ***
  WRITE (*, '(1x, a)') 'INITIALIZE 2D TUTORIAL MODEL'
  WRITE (*, '(10x,a,i4,1x,a1,1x,i4)') 'Grid size:', nx, 'x', ny
  WRITE (*, '(10x,a,i4)') 'Time steps', total_steps

  ! allocate memory for temporary fields
  ALLOCATE(field(ny, nx))


! ************************************
! *** Read initial field from file ***
! ************************************

  OPEN(11, file = '../inputs_online/true_initial.txt', status='old')
 
  DO i = 1, ny
     READ (11, *) field(i, :)
  END DO

  CLOSE(11)
  
END SUBROUTINE initialize
