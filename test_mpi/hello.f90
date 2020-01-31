program main  ! hello world

  implicit none
  include 'mpif.h'
  integer :: size, rank, ierr


  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD, size, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

  write(6, "(2(a,i3),2a)") "Hello from MPI: size = ", size, &
                                " rank = ", rank
  call mpi_finalize(ierr)
end

