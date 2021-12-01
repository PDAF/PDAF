!$Id: init_ens_offline.F90 1589 2015-06-12 11:57:58Z lnerger $
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
!! Implementation for the 2D offline example
!! without parallelization.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_ens_offline(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

  USE mod_assimilation, &    ! Model variables
       ONLY: nx, ny

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
  REAL, ALLOCATABLE :: field(:,:)     ! global model field
  CHARACTER(len=2) :: ensstr          ! String for ensemble member


! **********************
! *** INITIALIZATION ***
! **********************

  ! *** Generate full ensemble on filter-PE 0 ***
  WRITE (*, '(/9x, a)') 'Initialize state ensemble'
  WRITE (*, '(9x, a)') '--- read ensemble from files'
  WRITE (*, '(9x, a, i5)') '--- Ensemble size:  ', dim_ens
  
  ! allocate memory for temporary fields
  ALLOCATE(field(ny, nx))


! ********************************
! *** Read ensemble from files ***
! ********************************

  DO member = 1, dim_ens
     WRITE (ensstr, '(i1)') member
     OPEN(11, file = '../inputs_offline/ens_'//TRIM(ensstr)//'.txt', status='old')

     DO i = 1, ny
        READ (11, *) field(i, :)
     END DO
     DO j = 1, nx
        ens_p(1 + (j-1)*ny : j*ny, member) = field(1:ny, j)
     END DO

     CLOSE(11)
  END DO


! ****************
! *** clean up ***
! ****************

  DEALLOCATE(field)

END SUBROUTINE init_ens_offline
