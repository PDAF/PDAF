!$Id: init_ens_pdaf.f90 2135 2019-11-22 18:56:29Z lnerger $
!BOP
!
! !ROUTINE: init_ens_pdaf --- Initialize ensemble for filter
!
! !INTERFACE:
SUBROUTINE init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! If only a single filter algorithm is used, the 
! ensemble initialization can be performed directly
! in this routine. If a single filter is implemented,
! one can perform the initialization directly here.
!
! This variant is used with the simplified interface of
! PDAF. In this case, the name of the routine is defined
! within PDAF. This routine just calls the particular
! ensemble initialization routine for the selected filter.
!
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, ONLY: mype_world
  USe mo_kind_pdaf, ONLY: dp

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
  REAL(dp), INTENT(inout) :: state_p(dim_p)          ! PE-local model state
  REAL(dp), INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
  REAL(dp), INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_init       (as U_init_ens)
!EOP


! *******************************************************
! *** Call initialization routine for selected filter ***
! *******************************************************

  IF (mype_world==768) THEN
     WRITE (*,*) 'ECHAM-PDAF TEMPLATE init_ens_pdaf'
  END IF

  ens_p = 1.0

END SUBROUTINE init_ens_pdaf
  
