!$Id: obs_op_pdaf.F90 1864 2017-12-20 19:53:30Z lnerger $
!BOP
!
! !ROUTINE: cvtmat_ens_pdaf --- Generate matrix of localized ensemble perturbations
!
! !INTERFACE:
SUBROUTINE cvtmat_ens_pdaf(step, dim_p, dim_ens, dim_obs_p, dim_cvec_ens, pens_p, HV_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in: ensemble 3D-Var and hybrid 3D-Var
!
! The routine is called during the analysis step.
! It has to initialize the transformation matrix
! HV_p, which transforms from the model state
! to the control vector. For ensemble 3D-Var
! HV_p is in general the matrix of observed ensemble
! perturbations scaled by 1/(dim_cvec_ens-1). The 
! number of columns is defined by the ensemble
! size times the localization factors.
!
! For domain decomposition, the action is on the
! PE-local sub-domain of the state and has to 
! provide the observed sub-state for the PE-local 
! domain.
!
! Implementation for the 2D online example
! with or without parallelization.
!
! !REVISION HISTORY:
! 2021-03 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: mcols_cvec_ens

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step               ! Currrent time step
  INTEGER, INTENT(in) :: dim_p              ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_ens            ! Ensemble size
  INTEGER, INTENT(in) :: dim_obs_p          ! Dimension of observed state
  INTEGER, INTENT(in) :: dim_cvec_ens       ! Number of columns in HV_p
  REAL, INTENT(in)    :: pens_p(dim_p, dim_ens)  ! PE-local model state
  REAL, INTENT(inout) :: HV_p(dim_obs_p, dim_cvec_ens) ! PE-local tranform matrix
!EOP

! *** local variables ***
  INTEGER :: i, member       ! Counters


! *******************************
! *** Initialize matrix HV_p  ***
! *******************************

  ! Here, we only apply the observation operator for dim_ens columns
  ! More columns would result from localization
  DO member = 1, dim_ens
     ! [Hx_1 ... Hx_N]
     CALL obs_op_pdaf(step, dim_p, dim_obs_p, pens_p(:, member), HV_p(:, member))
  END DO

  ! Here er simulate additional columns by repeating the ensemble perturbations
  DO i = 2, mcols_cvec_ens
     DO member = (i-1)*dim_ens+1, i*dim_ens
        ! [Hx_1 ... Hx_N]
        CALL obs_op_pdaf(step, dim_p, dim_obs_p, pens_p(:, member-(i-1)*dim_ens), HV_p(:, member))
     END DO
  END DO

  HV_p = HV_p / SQRT(REAL(dim_cvec_ens-1))

END SUBROUTINE cvtmat_ens_pdaf
