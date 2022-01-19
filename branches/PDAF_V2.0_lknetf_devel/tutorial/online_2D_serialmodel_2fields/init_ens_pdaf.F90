!$Id$
!>  Initialize ensemble
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all ensemble filters.
!!
!! The routine is called when the filter is
!! initialized in PDAF_filter_init.  It has
!! to initialize an ensemble of dim_ens states.
!!
!! The routine is called by all filter processes and 
!! initializes the ensemble for the PE-local domain.
!!
!! Implementation for the 2D online example
!! without parallelization. Here, the ensmeble is
!! directly read from files.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

  USE mod_model, &         ! Model variables
       ONLY: nx, ny
  USE mod_assimilation, &  ! Assimilation variables
       ONLY: off_fields, id

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: filtertype                !< Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                     !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                   !< Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)            !< PE-local model state
  !< (It is not necessary to initialize the array 'state_p' for ensemble filters.
  !< It is available here only for convenience and can be used freely.)
  REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) !< Array not referenced for ensemble filters
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)     !< PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                   !< PDAF status flag

! *** local variables ***
  INTEGER :: i, j, member             ! Counters
  REAL, ALLOCATABLE :: readfield(:,:) ! global model field read from file
  CHARACTER(len=2) :: ensstr          ! String for ensemble member


! **********************
! *** INITIALIZATION ***
! **********************

  ! *** Generate full ensemble on filter-PE 0 ***
  WRITE (*, '(/9x, a)') 'Initialize state ensemble'
  WRITE (*, '(9x, a)') '--- read ensemble from files'
  WRITE (*, '(9x, a, i5)') '--- Ensemble size:  ', dim_ens
  
  ! allocate memory for temporary fields
  ALLOCATE(readfield(ny, nx))


! ********************************
! *** Read ensemble from files ***
! ********************************

  DO member = 1, dim_ens
     WRITE (ensstr, '(i1)') member

     ! Read field
     OPEN(11, file = '../inputs_online_2fields/ens_'//TRIM(ensstr)//'.txt', status='old')
 
     DO i = 1, ny
        READ (11, *) readfield(i, :)
     END DO
     DO j = 1, nx
        ens_p(off_fields(id%fieldA) + 1 + (j-1)*ny : off_fields(id%fieldA) + j*ny, member) = readfield(1:ny, j)
     END DO

     CLOSE(11)

     ! Read fieldB
     OPEN(12, file = '../inputs_online_2fields/ensB_'//TRIM(ensstr)//'.txt', status='old')
 
     DO i = 1, ny
        READ (12, *) readfield(i, :)
     END DO
     DO j = 1, nx
        ens_p(off_fields(id%fieldB) + 1 + (j-1)*ny : off_fields(id%fieldB) + j*ny, member) = readfield(1:ny, j)
     END DO

     CLOSE(12)
  END DO


! ****************
! *** clean up ***
! ****************

  DEALLOCATE(readfield)

END SUBROUTINE init_ens_pdaf
