! This is a strongly reduced header file for those
! MPI variables used in PDAF MPI calls. This is only
! intended for compilation without a real MPI library.

      INTEGER MPI_STATUS_SIZE
      PARAMETER (MPI_STATUS_SIZE=4)

      INTEGER MPI_INTEGER, MPI_REAL, MPI_DOUBLE_PRECISION 
      PARAMETER (MPI_REAL=26,MPI_DOUBLE_PRECISION=27,MPI_INTEGER=28)
 
      INTEGER MPI_COMM_WORLD
      PARAMETER (MPI_COMM_WORLD=91)
 
  INTEGER MPI_SUM, MPI_MAX
