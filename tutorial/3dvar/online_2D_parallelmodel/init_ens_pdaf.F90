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

  USE mod_model, &           ! Model variables
       ONLY: nx, ny, nx_p
  USE mod_parallel_model, &  ! Model parallelization variables
       ONLY: mype_model
  USE mod_parallel_pdaf, &   ! Assimilation parallelization variables
       ONLY: mype_filter
  USE mod_assimilation, &    ! Assimilation variables 
       ONLY: Vmat_p, dim_cvec, subtype

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
  REAL :: invdim_ens                  ! Inverse ensemble size
  REAL :: fact                        ! Scaling factor


! **********************
! *** INITIALIZATION ***
! **********************

  ! *** Generate full ensemble on filter-PE 0 ***
  IF (mype_filter==0) THEN
     WRITE (*, '(/9x, a)') 'Initialize state ensemble'
     WRITE (*, '(9x, a)') '--- read ensemble from files'
     WRITE (*, '(9x, a, i5)') '--- Ensemble size:  ', dim_ens
  END IF

  ! allocate memory for temporary fields
  ALLOCATE(field(ny, nx))


! ********************************
! *** Read ensemble from files ***
! ********************************

  DO member = 1, dim_ens
     WRITE (ensstr, '(i1)') member
     OPEN(11, file = '../../inputs_online/ens_'//TRIM(ensstr)//'.txt', status='old')

     ! Read global field
     DO i = 1, ny
        READ (11, *) field(i, :)
     END DO

     ! Initialize process-local part of ensemble
     DO j = 1, nx_p
        ens_p(1 + (j-1)*ny : j*ny, member) = field(1:ny, nx_p*mype_model + j)
     END DO

     CLOSE(11)
  END DO


! **********************************************
! *** Initialize square-root of P for 3D-Var ***
! **********************************************

  IF (filtertype==200 .AND. (subtype==0 .OR. subtype==6 .OR. subtype==7)) THEN
     
     WRITE (*, '(9x, a)') 'Initialize B^1/2 for 3D-Var'

     ! Here, we simply use the scaled ensemble perturbations

     ! Initialize numbers 
     invdim_ens = 1.0 / REAL(dim_cvec)

     ! Compute ensemble mean
     state_p = 0.0
     DO member = 1, dim_cvec
        DO i = 1, dim_p
           state_p(i) = state_p(i) + ens_p(i, member)
        END DO
     END DO
     state_p(:) = invdim_ens * state_p(:)

     ALLOCATE(Vmat_p(dim_p, dim_cvec))
  
     DO member = 1, dim_ens
        Vmat_p(:,member) = ens_p(:,member) - state_p(:)
     END DO

     fact = 1.0/SQRT(REAL(dim_cvec-1))

     Vmat_p = Vmat_p * fact
  END IF


! ****************
! *** clean up ***
! ****************

  DEALLOCATE(field)

END SUBROUTINE init_ens_pdaf
