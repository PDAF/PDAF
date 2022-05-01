!$Id$
!>  Time stepping loop with adaption for assimilation
!!
!! Time integration for simple 2D tutorial model
!! without parallelization of the model. In this
!! code variant the coupling to PDAF for ensemble
!! assimilation is completed.
!!
!! Each time step the field is shifted by one grid 
!! point in the vertical direction (first array index).
!!
!! __Revision history:__
!! * 2013-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE integrate_pdaf()

  USE mod_model, &          ! Include model variables
       ONLY: nx, ny, field, fieldB, total_steps
  USE mod_parallel_pdaf, &  ! Include parallelization variables
       ONLY: mype_world

  IMPLICIT NONE

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

        ! Second field (fieldB)
        store = fieldB(ny, j)

        DO i = ny-1,1,-1
           fieldB(i+1, j) = fieldB(i, j)
        END DO

        fieldB(1, j) = store
     END DO


#ifndef USE_PDAF     
! *** Write new fields into files ***

     WRITE (stepstr, '(i2.2)') step
     OPEN(11, file = 'true_step'//TRIM(stepstr)//'.txt', status = 'replace')

     DO i = 1, ny
        WRITE (11, *) field(i, :)
     END DO

     CLOSE(11)

     OPEN(12, file = 'trueB_step'//TRIM(stepstr)//'.txt', status = 'replace')

     DO i = 1, ny
        WRITE (12, *) fieldB(i, :)
     END DO

     CLOSE(12)
#endif

#ifdef USE_PDAF
     CALL assimilate_pdaf()
#endif

  END DO stepping

END SUBROUTINE integrate_pdaf
