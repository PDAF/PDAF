!$Id$
!>  Initialize 3D-Var
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in parameterized 3D-Var
!!
!! The routine is called when the filter is
!! initialized in PDAF_filter_init.  It has
!! to initialize an ensemble of dim_ens states.
!!
!! The routine is called by all filter processes and 
!! initializes the ensemble for the PE-local domain.
!!
!! Implementation for the 2D online example
!! without parallelization. Here, the ensemble 
!! information is directly read from files. Then
!! it is used to define the square root of the
!! B-matrix for 3D-Var.
!!
!! __Revision history:__
!! * 2021-04 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_3dvar_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

  USE mod_model, &         ! Model variables
       ONLY: nx, ny
  USE mod_assimilation, &  ! Assimilation variables
       ONLY: ensgroup, Vmat_p, dim_cvec

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: filtertype                !< Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                     !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                   !< Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)            !< PE-local model state
  !< (It is not necessary to initialize the array 'state_p' for ensemble filters.
  !< It is available here only for convenience and can be used freely.)
  REAL, INTENT(inout) :: Uinv(1,1)                 !< Array not referenced for 3D-Var
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)     !< PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                   !< PDAF status flag

! *** local variables ***
  INTEGER :: i, j, member             ! Counters
  REAL, ALLOCATABLE :: field(:,:)     ! global model field
  CHARACTER(len=2) :: ensstr          ! String for ensemble member
  REAL :: invdim_ens                  ! Inverse ensemble size
  REAL :: fact                        ! Scaling factor


! **********************
! *** INITIALIZATION ***
! **********************

  ! *** Read initial state and generate square root of B ***
  ! *** by reading the full ensemble on filter-PE 0      ***
  WRITE (*, '(/9x, a)') 'Initialize state and B^1/2 for 3D-Var'
  WRITE (*, '(9x, a)') '--- read ensemble from files'
  WRITE (*, '(9x, a, i5)') '--- members in B^1/2:  ', dim_cvec

  ! Initialize numbers 
  invdim_ens = 1.0 / REAL(dim_cvec)

  ! allocate memory for temporary fields
  ALLOCATE(field(ny, nx))

  ! Allocate matrix holding B^1/2 (from mod_assimilation)
  ALLOCATE(Vmat_p(dim_p, dim_cvec))


! ********************************
! *** Read ensemble from files ***
! ********************************

  DO member = 1, dim_cvec
     WRITE (ensstr, '(i1)') member
     IF (ensgroup==1) THEN
        OPEN(11, file = '../../inputs_online/ens_'//TRIM(ensstr)//'.txt', status='old')
     ELSE
        OPEN(11, file = '../../inputs_online/ensB_'//TRIM(ensstr)//'.txt', status='old')
     END IF

     DO i = 1, ny
        READ (11, *) field(i, :)
     END DO
     DO j = 1, nx
        Vmat_p(1 + (j-1)*ny : j*ny, member) = field(1:ny, j)
     END DO

     CLOSE(11)
  END DO


! **********************************************
! *** Initialize square-root of P for 3D-Var ***
! **********************************************

  WRITE (*, '(9x, a)') 'Initialize B^1/2'

  ! Here, we simply use the scaled ensemble perturbations

  ! Compute ensemble mean
  state_p = 0.0
  DO member = 1, dim_cvec
     DO i = 1, dim_p
        state_p(i) = state_p(i) + Vmat_p(i, member)
     END DO
  END DO
  state_p(:) = invdim_ens * state_p(:)

  ! Initialize ensemble perturbations
  DO member = 1, dim_cvec
     Vmat_p(:,member) = Vmat_p(:,member) - state_p(:)
  END DO

  fact = 1.0/SQRT(REAL(dim_cvec-1))

  Vmat_p = Vmat_p * fact


! ******************************************
! *** Initialize ensemble array for PDAF ***
! ******************************************

  ! This is only a single state for parameterized 3D-Var
  ens_p(:,1) = state_p(:)


! ****************
! *** clean up ***
! ****************

  DEALLOCATE(field)

END SUBROUTINE init_3dvar_pdaf
