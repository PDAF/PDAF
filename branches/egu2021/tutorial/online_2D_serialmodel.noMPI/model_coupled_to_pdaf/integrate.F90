!$Id$
!BOP
!
! !ROUTINE: integrate --- Time stepping loop with adaption for assimilation
!
! !INTERFACE:
SUBROUTINE integrate(nsteps)

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

! !ARGUMENTS:
  INTEGER, INTENT(in) :: nsteps ! Number of time steps to be performed

! *** local variables ***
  INTEGER :: step, i, j        ! Counters
  CHARACTER(len=2) :: stepstr  ! String for time step
  REAL :: store                ! Store single field element


! ****************
! *** STEPPING ***
! ****************

  stepping: DO step = 1 , nsteps

     WRITE (*,*) 'step', step

! *** Time step: Shift field vertically ***
     DO j = 1, nx
        store = field(ny, j)

        DO i = ny-1,1,-1
           field(i+1, j) = field(i, j)
        END DO

        field(1, j) = store
     END DO


#ifndef USE_PDAF     
! *** Write new field into file ***
     WRITE (stepstr, '(i2.2)') step
     OPEN(11, file = 'true_step'//TRIM(stepstr)//'.txt', status = 'replace')

     DO i = 1, ny
        WRITE (11, *) field(i, :)
     END DO

     CLOSE(11)
#endif     

  END DO stepping

END SUBROUTINE integrate
