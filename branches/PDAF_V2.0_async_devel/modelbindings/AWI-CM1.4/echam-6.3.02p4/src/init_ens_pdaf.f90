!$Id$
!>  Routine to initialize ensemble for PDAF
!!
!! User-supplied routine for PDAF.
!!
!! If only a single filter algorithm is used, the 
!! ensemble initialization can be performed directly
!! in this routine. If a single filter is implemented,
!! one can perform the initialization directly here.
!!
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
SUBROUTINE init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

  USE mod_parallel_pdaf, ONLY: mype_world
  USE mod_assim_atm_pdaf, ONLY: dp

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: filtertype                    !< Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                         !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                       !< Size of ensemble
  REAL(dp), INTENT(inout) :: state_p(dim_p)            !< PE-local model state
  REAL(dp), INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) !< Array not referenced for SEIK
  REAL(dp), INTENT(out)   :: ens_p(dim_p, dim_ens)     !< PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                       !< PDAF status flag

! *** local variables ***
  INTEGER :: col         ! Counters


! ***************************
! *** Initialize ensemble ***
! ***************************

  ! Here, we simply use the state from the model restart
  ! i.e. zero ensemble spread, and let the model spin up

  CALL collect_state_pdaf(dim_p,state_p)

  DO col = 1, dim_ens
     ens_p(1:dim_p,col) = state_p(1:dim_p)
  END DO

END SUBROUTINE init_ens_pdaf
  
