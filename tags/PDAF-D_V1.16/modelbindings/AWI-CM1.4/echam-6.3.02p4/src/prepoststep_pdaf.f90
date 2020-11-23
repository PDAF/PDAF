!$Id: prepoststep_pdaf.f90 2135 2019-11-22 18:56:29Z lnerger $
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
! This variant is used with the simplified interface of
! PDAF. In this case, the name of the routine is defined
! within PDAF. This routine just calls the prepoststep
! routine corresponding to the selected filter algorithm.
!
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, ONLY: mype_filter

  USE mod_assim_pdaf, & ! Variables for assimilation
       ONLY: filtertype
  USE mo_kind_pdaf, ONLY: dp

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
     ! (When the routine is called before the analysis -step is provided.)
  INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   ! PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   ! PE-local dimension of observation vector
  REAL(dp), INTENT(inout) :: state_p(dim_p) ! PE-local forecast/analysis state
  ! The array 'state_p' is not generally not initialized in the case of SEIK.
  ! It can be used freely here.
  REAL(dp), INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) ! Inverse of matrix U
  REAL(dp), INTENT(inout) :: ens_p(dim_p, dim_ens)      ! PE-local state ensemble
  INTEGER, INTENT(in) :: flag        ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_get_state      (as U_prepoststep)
! Called by: PDAF_X_update       (as U_prepoststep)
!EOP


! **************************************************************
! *** Call prepoststep routine according to filter algorithm ***
! **************************************************************

  if (mype_filter==0) write (*,*) 'ECHAM-PDAF TEMPLATE prepoststep_pdaf'

END SUBROUTINE prepoststep_pdaf
