!>  Time stepping loop of tutorial model
!!
!! Time stepping for simple 2D tutorial model
!! with domain decomposition.
!!
!! Each time step the field is shifted by one grid 
!! point in the vertical direction (first array index).
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
module model_step_mod

contains

  subroutine stepping()

    use mpi                         ! MPI
    use model_mod, &                ! Model variables
         only: ny, nx, fieldA, fieldB, total_steps
    use model_io_mod, &             ! File operations
         only: io_write_sngl

    implicit none

! *** local variables ***
    integer :: step, i, j           ! Counters
    real :: store                   ! Store single field element
    character(len=100) :: filename  ! Name of output file
    character(len=2) :: stepstr     ! String for time step


! ****************
! *** STEPPING ***
! ****************

    write (*, '(1x, a)') 'MODEL INTEGRATION'

    steps: do step = 1 , total_steps

       write (*,*) 'step', step

       ! *** Time step: Shift fields vertically ***
       
       do j = 1, nx
          ! Field A
          store = fieldA(ny, j)

          do i = ny-1,1,-1
             fieldA(i+1, j) = fieldA(i, j)
          end do

          fieldA(1, j) = store

          ! Field B
          store = fieldB(ny, j)

          do i = ny-1,1,-1
             fieldB(i+1, j) = fieldB(i, j)
          end do

          fieldB(1, j) = store
       end do

       ! *** Write current fields into files ***

       write (stepstr, '(i2.2)') step

       filename = 'fieldA_step'//trim(stepstr)//'.txt'
       call io_write_sngl(step, filename, fieldA)

       filename = 'fieldB_step'//trim(stepstr)//'.txt'
       call io_write_sngl(step, filename, fieldB)

    end do steps

  end subroutine stepping

end module model_step_mod
