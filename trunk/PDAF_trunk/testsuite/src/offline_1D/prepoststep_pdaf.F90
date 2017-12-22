!$Id: prepoststep_pdaf.F90 1253 2012-01-30 18:56:08Z lnerger $
!BOP
!
! !ROUTINE: prepoststep_pdaf - Routine controlling ensemble integration for PDAF
!
! !INTERFACE:
SUBROUTINE prepoststep_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! This routine is called directly by PDAF and subsequently calls the
! prepoststep routine corresponding to the selected filter algorithm.
!
! 
!
! !REVISION HISTORY:
! 2010-07 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, & ! Variables for assimilation
       ONLY: filtertype

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
     ! (When the routine is called before the analysis -step is provided.)
  INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   ! PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   ! PE-local dimension of observation vector
  REAL, INTENT(inout) :: state_p(dim_p) ! PE-local forecast/analysis state
  ! The array 'state_p' is not generally not initialized in the case of SEIK.
  ! It can be used freely here.
  REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) ! Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      ! PE-local state ensemble
  INTEGER, INTENT(in) :: flag        ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_get_state      (as U_prepoststep)
! Called by: PDAF_seek_update    (as U_prepoststep)
! Called by: PDAF_seik_update    (as U_prepoststep)
! Called by: PDAF_enkf_update    (as U_prepoststep)
! Called by: PDAF_lseik_update    (as U_prepoststep)
! Called by: PDAF_etkf_update    (as U_prepoststep)
! Called by: PDAF_letkf_update    (as U_prepoststep)
!EOP


! **************************************************************
! *** Call prepoststep routine according to filter algorithm ***
! **************************************************************

  IF (filtertype == 0) THEN
     ! Special call for SEEK
     CALL prepoststep_seek_offline(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
          state_p, Uinv, ens_p, flag)
  ELSE
     ! General call for ensemble-based KFs
     CALL prepoststep_ens_offline(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
          state_p, Uinv, ens_p, flag)
  END IF

END SUBROUTINE prepoststep_pdaf
