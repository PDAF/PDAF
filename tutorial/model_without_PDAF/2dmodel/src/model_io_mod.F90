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
  subroutine io_write_sngl(step, filename, field)

    use mpi
    use model_mod, &
         only: nx, ny

    implicit none

    ! Arguments
    integer, intent(in) :: step                 !< Model time step
    character(len=100), intent(in) :: filename  !< Name of output file
    real, intent(in) :: field(:,:)              !< Model field

    ! Local variables
    integer :: i                          ! Counter
    character(len=2) :: stepstr           ! String for time step


    ! *** Write file

     write (stepstr, '(i2.2)') step
     open(11, file = trim(filename), status = 'replace')

     do i = 1, ny
        write (11, *) field(i, :)
     end do

     close(11)     

  end subroutine io_write_sngl

!-------------------------------------------------------------------------------
!> Read field on subdomain from a text file
!!
!! Routine to read the subdomain-part of a specific field
!! at one time step.
!!
  subroutine io_read_sngl(filename, field)

    use model_mod, &
         only: nx, ny

    implicit none

    ! Arguments
    character(len=100), intent(in) :: filename  !< Name of output file
    real, intent(inout) :: field(:,:)           !< Model field

    ! Local variables
    integer :: i                     ! Counter

    ! *** Read field

    open(11, file = trim(filename), status='old')
 
    do i = 1, ny
       read (11, *) field(i, :)
    end do

    close(11)

  end subroutine io_read_sngl

end module model_io_mod
