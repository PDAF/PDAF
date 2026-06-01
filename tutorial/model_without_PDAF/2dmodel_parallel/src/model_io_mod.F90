!> Module for netcdf operations in 2D tutorial model
!!
!! This module provides functionality to read and write
!! netcdf files for the 2-dimensional tutorial model
!! with parallelization.
!!
!! __Revision history:__
!! * 2026-02 - Lars Nerger - Initial code for advanced tutorial revising tutorial case
!! * Later revisions - see repository log
!!
module model_io_mod

  implicit none
  save
  private

  public io_write_sngl, io_read_sngl

contains

!-------------------------------------------------------------------------------
!> Write a field into a netCDF output file
!!
!! Routine to write a specific field at one time step.
!! Each call generates a single file.
!!
  subroutine io_write_sngl(step, filename, field_p)

    use mpi
    use model_parallel_mod, &
         only: myproc_world, COMM_2Dmodel
    use model_mod, &
         only: nx_p, nx, ny

    implicit none

    ! Arguments
    integer, intent(in) :: step                 !< Model time step
    character(len=100), intent(in) :: filename  !< Name of output file
    real, intent(in) :: field_p(:,:)            !< Decomposed model field

    ! Local variables
    integer :: i                          ! Counter
    character(len=2) :: stepstr           ! String for time step
    real, allocatable :: field(:,:)       ! Array for global model field
    integer :: MPIerr                     ! MPI error flag


    ! *** Gather global field on process 0
    allocate(field(ny, nx))

    call MPI_Gather(field_p, nx_p*ny, MPI_DOUBLE_PRECISION, field, nx_p*ny, &
         MPI_DOUBLE_PRECISION, 0, COMM_2Dmodel, MPIerr)

    ! *** Write file on process 0

    if (myproc_world==0) then

        write (stepstr, '(i2.2)') step
        open(11, file = trim(filename), status = 'replace')

        do i = 1, ny
           write (11, *) field(i, :)
        end do

        close(11)     

    end if

    deallocate(field)

  end subroutine io_write_sngl

!-------------------------------------------------------------------------------
!> Read field on subdomain from a netCDF output file
!!
!! Routine to read the subdomain-part of a specific field
!! at one time step.
!!
  subroutine io_read_sngl(filename, field_p)

    use model_mod, &
         only: nx, ny, nx_p, offset_x_p

    implicit none

    ! Arguments
    character(len=100), intent(in) :: filename  !< Name of output file
    real, intent(inout) :: field_p(:,:)         !< Decomposed model field

    ! Local variables
    integer :: i,j                  ! Counters
    real, allocatable :: field(:,:) ! Global model field

    ! *** Read field

    allocate(field(ny, nx))

    ! Read global model field
    open(11, file = trim(filename), status='old')
 
    do i = 1, ny
       read (11, *) field(i, :)
    end do

    close(11)

    ! Initialize local part of model field
    do j = 1, nx_p
       do i = 1, ny
          field_p(i,j) = field(i, offset_x_p + j)
       end do
    end do

    deallocate(field)

  end subroutine io_read_sngl

end module model_io_mod
