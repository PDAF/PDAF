!>  Initialize model
!!
!! Initialization routine for the simple 2D model without
!! parallelization of the model.
!!
!! The routine defines the size of the model grid and
!! reads the initial state from a file. 
!!
!! __Revision history:__
!! * 2013-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE initialize()

  USE mod_model, &              ! Model variables
       ONLY: nx, ny, nx_p, ndim, field_p, total_steps
  USE mod_parallel_model, &     ! Model parallelzation variables
       ONLY: mype_world, mype_model, npes_model, abort_parallel
  USE parser, &                 ! Parser function
       ONLY: parse

  IMPLICIT NONE

! *** local variables ***
  INTEGER :: i, j                 ! Counters
  REAL, ALLOCATABLE :: field(:,:) ! GLobal model field
  CHARACTER(len=32) :: handle     ! handle for command line parser
  INTEGER :: gridsize=1           ! Setting for grid dimensions: (1) 18x36, (2) 512x2048


! **********************
! *** INITIALIZATION ***
! **********************

  ! Parse grid size index
  handle='gridsize'
  CALL parse(handle, gridsize)

! *** Model specifications ***

  ! Number of time steps to perform
  total_steps = 18

  ! Number of coordinate directions
  ndim = 2

  IF (gridsize==1) THEN
     nx = 36    ! Extent of grid in x-direction
     ny = 18    ! Extent of grid in y-direction
  ELSEIF (gridsize==2) THEN
     nx = 256    ! Extent of grid in x-direction
     ny = 128    ! Extent of grid in y-direction
  ELSE
     nx = 2048   ! Extent of grid in x-direction
     ny = 512    ! Extent of grid in y-direction
  END IF


! *** Screen output ***
  IF (mype_world == 0) THEN
     WRITE (*, '(1x, a)') 'INITIALIZE PARALLELIZED 2D TUTORIAL MODEL'
     WRITE (*, '(10x,a,i4,1x,a1,1x,i4)') 'Grid size:', nx, 'x', ny
     WRITE (*, '(10x,a,i4)') 'Time steps', total_steps
  END IF

! *** Initialize size of local nx for parallelization ***
  IF (npes_model==1 .OR. npes_model==2 .OR. npes_model==3 .OR. npes_model==4 .OR. &
       npes_model==6 .OR. npes_model==9 .OR. npes_model==12 .OR. npes_model==18) THEN
     ! Split x-direction in chunks of equal size
     nx_p = nx / npes_model
  ELSE
     WRITE (*,*) 'ERROR: Invalid number of processes'
     CALL abort_parallel()
  END IF

  IF (mype_world == 0) THEN
     WRITE (*, '(/2x, a, i3, a)') &
          '-- Domain decomposition over', npes_model, ' PEs'
     WRITE (*, '(2x,a,i3,a,i3)') &
          '-- local domain sizes (nx_p x ny): ', nx_p, ' x', ny
  END IF

  ! allocate memory for process-local part of field
  ALLOCATE(field_p(ny, nx_p))


! ************************************
! *** Read initial field from file ***
! ************************************

  ALLOCATE(field(ny, nx))

  ! Read global model field
  IF (nx==36) THEN
     OPEN(11, file = '../inputs_online.18x36/true_initial.txt', status='old')
  ELSEIF (nx==256) THEN
     OPEN(11, file = '../inputs_online.128x256/true_initial.txt', status='old')
  ELSE
     OPEN(11, file = '../inputs_online.512x2048/true_initial.txt', status='old')
  END IF
 
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
