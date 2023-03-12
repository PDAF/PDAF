!$Id: integrate_pdaf.F90 1411 2013-09-25 14:04:41Z lnerger $
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

  USE mpi                     ! MPI
  USE mod_model, &            ! Model variables
       ONLY: nx, ny, nx_p, field_p, total_steps
  USE mod_parallel_model, &   ! Model parallelization variables
       ONLY: mype_model, MPIErr, COMM_model
  USE mod_parallel_pdaf, &    ! PDAF parallelization variables
       ONLY: task_id

  IMPLICIT NONE

! *** local variables ***
  INTEGER :: step, i, j        ! Counters
  CHARACTER(len=2) :: stepstr  ! String for time step
  REAL :: store                ! Store single field element
  REAL, ALLOCATABLE :: field(:,:) ! GLobal model field


! ****************
! *** STEPPING ***
! ****************

  IF (task_id==1 .AND. mype_model==0) WRITE (*, '(1x, a)') 'MODEL: START INTEGRATION'

  stepping: DO step = 1 , total_steps

     IF (task_id==1 .AND. mype_model==0) WRITE (*,*) 'step', step

! *** Time step: Shift field vertically ***
     DO j = 1, nx_p
        store = field_p(ny, j)

        DO i = ny-1,1,-1
           field_p(i+1, j) = field_p(i, j)
        END DO

        field_p(1, j) = store

     END DO

#ifndef USE_PDAF     
! *** Write new field into file ***

     ! Gather global field on process 0
     ALLOCATE(field(ny, nx))

     CALL MPI_Gather(field_p, nx_p*ny, MPI_DOUBLE_PRECISION, field, nx_p*ny, &
          MPI_DOUBLE_PRECISION, 0, COMM_model, MPIerr)

     ! Write file from process 0
     IF (task_id==1 .AND. mype_model==0) THEN
        WRITE (stepstr, '(i2.2)') step
        OPEN(11, file = 'true_step'//TRIM(stepstr)//'.txt', status = 'replace')

        DO i = 1, ny
           WRITE (11, *) field(i, :)
        END DO

        CLOSE(11)     
     END IF

     DEALLOCATE(field)
#endif

#ifdef USE_PDAF
     CALL assimilate_pdaf()
#endif

  END DO stepping

END SUBROUTINE integrate_pdaf
