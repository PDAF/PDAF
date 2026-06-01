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
         only: n_dim, nx, ny, nx_p, offset_x_p, fieldA_p, fieldB_p, &
         coords_x_p, coords_y_p, total_steps
    use model_parallel_mod, &       ! Model parallelzation variables
         only: myproc_world, myproc_2Dmodel, nproc_2Dmodel, abort_parallel
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
    if (myproc_world == 0) then
       write (*, '(1x, a)') 'INITIALIZE PARALLELIZED 2D TUTORIAL MODEL'
       write (*, '(10x,a,i4,1x,a1,1x,i4)') 'Grid size:', nx, 'x', ny
       write (*, '(10x,a,i4)') 'Time steps', total_steps
       write (*, '(10x,a,a)') 'Read inputs from ', path
    end if    


! *** Initialize size of local nx for parallelization ***
    if (nproc_2Dmodel==1 .or. nproc_2Dmodel==2 .or. nproc_2Dmodel==3 .or. nproc_2Dmodel==4 .or. &
         nproc_2Dmodel==6 .or. nproc_2Dmodel==9 .or. nproc_2Dmodel==12 .or. nproc_2Dmodel==18) then
       ! Split x-direction in chunks of equal size
       nx_p = nx / nproc_2Dmodel
    else
       write (*,*) 'ERROR: Invalid number of processes'
       call abort_parallel()
    end if

    if (myproc_world == 0 .and. nproc_2Dmodel > 1) then
       write (*, '(/2x, a, i3, a)') &
            '-- Domain decomposition over', nproc_2Dmodel, ' Processs'
       write (*, '(2x,a,i3,a,i3)') &
            '-- local domain sizes (nx_p x ny): ', nx_p, ' x', ny
    end if

    ! Set offset of process-local grid in global grid
    offset_x_p = nx_p*myproc_2Dmodel

    ! allocate memory for process-local part of fields
    allocate(fieldA_p(ny, nx_p))
    allocate(fieldB_p(ny, nx_p))
    allocate(coords_x_p(nx_p))
    allocate(coords_y_p(ny))


! *************************************
! *** Read initial fields from file ***
! *************************************

    filename = trim(path)//'/trueA_initial.txt'
    call io_read_sngl(filename, fieldA_p)

    filename = trim(path)//'/trueB_initial.txt'
    call io_read_sngl(filename, fieldB_p)


! *************************************
! *** Initialize coordinates        ***
! *************************************

    ! The model coordinates are the grid point indices
    ! stored as real values

    ! Account for decomposition in x-direction
    do i = 1, nx_p
       coords_x_p(i) = real(i + offset_x_p)
    end do

    ! We don't use decomposition in y-direction
    do j = 1, ny
       coords_y_p(j) = real(j)
    end do

  end subroutine initialize

end module model_init_mod
