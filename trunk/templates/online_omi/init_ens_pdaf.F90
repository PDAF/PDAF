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
!! __Revision history:__
!! * 2010-07 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

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


! *******************************************************
! *** Call initialization routine for selected filter ***
! *******************************************************

  IF (filtertype == 0) THEN
     ! EOF initialization for SEEK
     CALL init_seek(filtertype, dim_p, dim_ens, state_p, Uinv, &
          ens_p, flag)
  ELSE IF (filtertype == 2) THEN
     ! Use random sampling initialization
     CALL init_enkf(filtertype, dim_p, dim_ens, state_p, Uinv, &
          ens_p, flag)
  ELSE
     ! Use 2nd-order exact sampling
     CALL init_seik(filtertype, dim_p, dim_ens, state_p, Uinv, &
          ens_p, flag)
  END IF

END SUBROUTINE init_ens_pdaf
