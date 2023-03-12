!$Id: integrate.F90 1411 2013-09-25 14:04:41Z lnerger $
!BOP
!
! !ROUTINE: integrate --- Time stepping loop of tutorial model
!
! !INTERFACE:
SUBROUTINE integrate()

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
       ONLY: nx, ny, nx_p, field_p, total_steps
  USE mod_parallel_model, &
       ONLY: mype_world, MPI_DOUBLE_PRECISION, MPIErr, COMM_model

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: main
  !EOP

! *** local variables ***
  INTEGER :: step, i, j        ! Counters
  CHARACTER(len=2) :: stepstr  ! String for time step
  REAL :: store                ! Store single field element
  REAL, ALLOCATABLE :: field(:,:) ! GLobal model field


! ****************
! *** STEPPING ***
! ****************

  IF (mype_world==0) WRITE (*, '(1x, a)') 'START INTEGRATION'

  stepping: DO step = 1 , total_steps

     IF (mype_world==0) WRITE (*,*) 'step', step

! *** Time step: Shift field vertically ***
     DO j = 1, nx_p
        store = field_p(ny, j)

        DO i = ny-1,1,-1
           field_p(i+1, j) = field_p(i, j)
        END DO

        field_p(1, j) = store

     END DO

! *** Write new field into file ***

     ! Gather global field on process 0
     ALLOCATE(field(ny, nx))

     CALL MPI_Gather(field_p, nx_p*ny, MPI_DOUBLE_PRECISION, field, nx_p*ny, &
          MPI_DOUBLE_PRECISION, 0, COMM_model, MPIerr)

     ! Write file from process 0
     IF (mype_world==0) THEN
        WRITE (stepstr, '(i2.2)') step
        OPEN(11, file = 'true_step'//TRIM(stepstr)//'.txt', status = 'replace')

        DO i = 1, ny
           WRITE (11, *) field(i, :)
        END DO

        CLOSE(11)     
     END IF

     DEALLOCATE(field)

  END DO stepping

END SUBROUTINE integrate
