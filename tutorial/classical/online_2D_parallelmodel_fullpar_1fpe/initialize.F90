!$Id: initialize.F90 1426 2013-09-25 15:21:35Z lnerger $
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
       ONLY: nx, ny, nx_p, field_p, total_steps
  USE mod_parallel_model, &
       ONLY: mype_model, npes_model, abort_parallel
#ifdef USE_PDAF
  USE mod_parallel_pdaf, &
       ONLY: task_id
#endif

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: main
!EOP

! *** local variables ***
  INTEGER :: i, j                 ! Counters
  REAL, ALLOCATABLE :: field(:,:) ! Global model field


! **********************
! *** INITIALIZATION ***
! **********************

! *** Model specifications ***
  nx = 36          ! Extent of grid in x-direction
  ny = 18          ! Extent of grid in y-direction
  total_steps = 18 ! Number of time steps to perform

! *** Initialize size of local nx for parallelization ***
  IF (npes_model==1 .OR. npes_model==2 .OR. npes_model==3 .OR. npes_model==4 .OR. &
       npes_model==6 .OR.npes_model==9) THEN
     ! Split x-diection in chunks of equal size
     nx_p = nx / npes_model
  ELSE
     WRITE (*,*) 'ERROR: Invalid number of processes'
     CALL abort_parallel()
  END IF

! *** Screen output ***
#ifndef USE_PDAF
  IF (mype_model == 0) THEN
#else
  IF (task_id==1 .AND. mype_model == 0) THEN
#endif
     WRITE (*, '(1x, a)') 'MODEL-side: INITIALIZE PARALLELIZED 2D TUTORIAL MODEL'
     WRITE (*, '(10x,a,i4,1x,a1,1x,i4)') 'Grid size:', nx, 'x', ny
     WRITE (*, '(10x,a,i4)') 'Time steps', total_steps

     WRITE (*, '(/2x, a, i3, a)') &
          '-- Domain decomposition over', npes_model, ' PEs'
     WRITE (*, '(2x,a,i3,a,i3)') &
          '-- local domain sizes (nx_p x ny): ', nx_p, ' x', ny
  END IF

  ! allocate memory for temporary fields
  ALLOCATE(field_p(ny, nx_p))


! ************************************
! *** Read initial field from file ***
! ************************************

  ALLOCATE(field(ny, nx))

  ! Read global model field
  OPEN(11, file = '../../inputs_online/true_initial.txt', status='old')
 
  DO i = 1, ny
     READ (11, *) field(i, :)
  END DO

  CLOSE(11)

  ! Initialize local part of model field
  DO j = 1, nx_p
     DO i = 1, ny
        field_p(i,j) = field(i, nx_p*mype_model + j)
     END DO
  END DO

  DEALLOCATE(field)

END SUBROUTINE initialize
