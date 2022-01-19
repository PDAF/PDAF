!$Id: init_ens_pdaf.F90 1253 2012-01-30 18:56:08Z lnerger $
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
! This routine simply calls the initialization
! routine for the specified filter. 
!
! If only a single filter algorithm is used, the 
! ensemble initialization can be performed directly
! in this routine.
!
! !REVISION HISTORY:
! 2010-07 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state
  REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_init_flx    (as U_ens_init)
! Calls: init_seik
!EOP


! *******************************************************
! *** Call initialization routine for selected filter ***
! *******************************************************

  IF (filtertype == 0) THEN
     CALL init_seek_offline(filtertype, dim_p, dim_ens, state_p, Uinv, &
          ens_p, flag)
  ELSE IF (filtertype == 2) THEN
     CALL init_enkf_offline(filtertype, dim_p, dim_ens, state_p, Uinv, &
          ens_p, flag)
  ELSE
     CALL init_seik_offline(filtertype, dim_p, dim_ens, state_p, Uinv, &
          ens_p, flag)
  END IF

END SUBROUTINE init_ens_pdaf
