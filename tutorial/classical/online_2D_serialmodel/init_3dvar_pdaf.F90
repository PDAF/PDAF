!$Id: init_3dvar.F90 1589 2015-06-12 11:57:58Z lnerger $
!BOP
!
! !ROUTINE: init_3dvar_pdaf --- Initialize ensemble
!
! !INTERFACE:
SUBROUTINE init_3dvar_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in: 3D-Var
!
! The routine is called when the filter is
! initialized in PDAF\_filter\_init. It has
! to initialize the initial state estimate
! for the assimilation. In addition it 
! initializes the square-root of P that is 
! used for 3D-Var (in this example we use
! explicitly a matrix holding the square-root
! of P, which is given by the scaled ensemble
! perturbations.
!
! The routine is called by all filter processes and
! initializes the ensemble for the PE-local domain.
!
! Implementation for the 2D online example
! without parallelization.
!
! !REVISION HISTORY:
! 2021-03 - Lars Nerger - Initial code based on init_ens_pdaf
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &
       ONLY: nx, ny
  USE mod_assimilation, &
       ONLY: ensgroup, dim_cvec, Vmat_p

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state
  ! It is not necessary to initialize the array 'state_p' for 3D-Var
  ! It is available here only for convenience and can be used freely.
  REAL, INTENT(inout) :: Uinv(1,1)               ! Array not referenced for 3D-Var
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_filter_init    (as U_ens_init)
!EOP

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

  ens_p(:,1) = state_p(:)


! ****************
! *** clean up ***
! ****************

  DEALLOCATE(field)

END SUBROUTINE init_3dvar_pdaf
