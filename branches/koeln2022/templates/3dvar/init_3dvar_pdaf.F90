!$Id: init_3dvar_pdaf.F90 901 2021-11-30 13:43:16Z lnerger $
!>  Initialize 3D-Var
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in parameterized 3D-Var
!!
!! The routine is called when the filter is
!! initialized in PDAF_filter_init.
!!
!! This routine ha to fill the initial state vector. 
!! In addition one can initialize the parameterized
!! square-root of the background state covariance matrix.
!!
!! The routine is called by all filter processes and 
!! initializes the ensemble for the local domain.
!!
!! __Revision history:__
!! * 2021-12 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_3dvar_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

  USE mod_assimilation, &  ! Assimilation variables
       ONLY: Vmat_p, dim_cvec

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
  REAL :: invdim_ens                  ! Inverse ensemble size


! **********************
! *** INITIALIZATION ***
! **********************

  WRITE (*, *) 'TEMPLATE init_3dvar_pdaf.F90: Initialize state and covariance matrix'

  ! *** Read initial state and generate square root of B ***
  ! *** by reading the full ensemble on filter-PE 0      ***
  WRITE (*, '(/9x, a)') 'Initialize state and B^1/2 for 3D-Var'


! *************************************************
! *** Initialize initial state and covar matrix ***
! *************************************************

!  state_p = ??


  ! Allocate matrix holding B^1/2 (from mod_assimilation)
!  ALLOCATE(Vmat_p(dim_p, dim_cvec))

!  Vmat_p = ??  


! ******************************************
! *** Initialize ensemble array for PDAF ***
! ******************************************

  ! This is only a single state for parameterized 3D-Var
  ens_p(:,1) = state_p(:)


! ****************
! *** clean up ***
! ****************


END SUBROUTINE init_3dvar_pdaf
