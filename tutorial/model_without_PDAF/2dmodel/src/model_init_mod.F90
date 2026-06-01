!>  Initialize model
!!
!! Initialization routine for the simple 2D model with
!! parallelization of the model.
!!
!! The routine defines the size of the model grid and
!! reads the initial state from a file. 
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
module model_init_mod

contains

  subroutine initialize()

    use model_mod, &                ! Model variables
         only: n_dim, nx, ny, fieldA, fieldB, &
         coords_x, coords_y, total_steps
    use model_io_mod, &             ! File operations
         only: io_read_sngl

    implicit none

! *** local variables ***
    integer :: i, j                  ! Counters
    character(len=100) :: filename   ! Name of output file
    character(len=100) :: str1, str2 ! String for cammand line parsing
    character(len=100) :: path       ! Path for input files


! **********************
! *** INITIALIZATION ***
! **********************

! *** Model specifications ***
    nx = 36                    ! Extent of grid in x-direction
    ny = 18                    ! Extent of grid in y-direction
    n_dim = 2                  ! Number of model dimensions
    total_steps = 20           ! Number of time steps to perform
    path = '../../inputs_2fields' ! Path to input files


! *** Parse command line arguments ***

    ! Number of time steps
    IF (command_argument_count() > 0) THEN 
       DO i = 1, command_argument_count() - 1 
          CALL get_command_argument(i, str1)
          CALL get_command_argument(i+1, str2)
          IF (str1 == '-nsteps') READ(str2, *) total_steps
       ENDDO
    ENDIF

    ! Path for input files
    IF (command_argument_count() > 0) THEN 
       DO i = 1, command_argument_count() - 1 
          CALL get_command_argument(i, str1)
          CALL get_command_argument(i+1, str2)
          IF (str1 == '-path') READ(str2, '(a)') path
       ENDDO
    ENDIF

! *** Screen output ***
    write (*, '(1x, a)') 'INITIALIZE PARALLELIZED 2D TUTORIAL MODEL'
    write (*, '(10x,a,i4,1x,a1,1x,i4)') 'Grid size:', nx, 'x', ny
    write (*, '(10x,a,i4)') 'Time steps', total_steps
    write (*, '(10x,a,a)') 'Read inputs from ', path


! *** allocate memory for fields and coordinates
    allocate(fieldA(ny, nx))
    allocate(fieldB(ny, nx))
    allocate(coords_x(nx))
    allocate(coords_y(ny))


! *************************************
! *** Read initial fields from file ***
! *************************************

    filename = trim(path)//'/trueA_initial.txt'
    call io_read_sngl(filename, fieldA)

    filename = trim(path)//'/trueB_initial.txt'
    call io_read_sngl(filename, fieldB)


! *************************************
! *** Initialize coordinates        ***
! *************************************

    ! The model coordinates are the grid point indices
    ! stored as real values

    do i = 1, nx
       coords_x(i) = real(i)
    end do

    do j = 1, ny
       coords_y(j) = real(j)
    end do

  end subroutine initialize

end module model_init_mod
