!$Id: integrate_pdaf.F90 1864 2017-12-20 19:53:30Z lnerger $
!BOP
!
! !ROUTINE: integrate_pdaf --- Time stepping loop with adaption for assimilation
!
! !INTERFACE:
SUBROUTINE integrate_pdaf()

! !DESCRIPTION:
! Initialization routine for the simple 2D model without
! parallelization of the model.
!
! The routine defines the size of the model grid and
! read the initial state from a file. 
!
! !REVISION HISTORY:
! 2013-09 - Lars Nerger - Initial code based on offline_1D
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &
       ONLY: nx, ny, field, total_steps
  USE mod_parallel_pdaf, &
       ONLY: mype_world

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: main
!EOP

! *** local variables ***
  INTEGER :: step, i, j        ! Counters
  CHARACTER(len=2) :: stepstr  ! String for time step
  REAL :: store                ! Store single field element


! ****************
! *** STEPPING ***
! ****************

  WRITE (*, '(1x, a)') 'START INTEGRATION'

  stepping: DO step = 1 , total_steps

     IF (mype_world==0) WRITE (*,*) 'step', step

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

#ifdef USE_PDAF
     CALL assimilate_pdaf()
#endif
     

  END DO stepping

END SUBROUTINE integrate_pdaf
