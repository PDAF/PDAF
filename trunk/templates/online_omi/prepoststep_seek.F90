!$Id$
!> Used-defined Pre/Poststep routine for PDAF
!!
!! User-supplied routine for PDAF.
!!
!! Used in the filters: SEEK
!! 
!! The routine is called before the analysis
!! and after the re-diagonalization.  Also it 
!! is called once at the initial time before 
!! any forecasts are computed. 
!! The routine provides full access to the state 
!! estimate and the state covariance matrix 
!! to the user.  Thus, user-controlled pre- and 
!! poststep operations can be performed here. 
!! For example the forecast and the analysis 
!! states and error covariance matrix can be 
!! analized, e.g. by computing the estimated 
!! variance.  In addition, the estimates can be 
!! written to disk.  If a user considers to 
!! perform adjustments to the estimates (e.g. 
!! for balances), this routine is the right 
!! place for it.
!!
!! The routine is called by all filter processes.
!!
!! __Revision history:__
!! * 2004-11 - Lars Nerger - Initial code
!! * Later revisions - see svn log
!!
SUBROUTINE prepoststep_seek(step, dim_p, dim_eof, dim_eof_p, dim_obs_p, &
     state_p, Uinv, eofV_p, flag)

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step        !< Current time step (negative for call after forecast)
  INTEGER, INTENT(in) :: dim_p       !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     !< Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   !< PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   !< PE-local dimension of observation vector
  REAL, INTENT(inout) :: state_p(dim_p) !< PE-local forecast/analysis state
  ! *** The covariance P is decomposed as P = V U V^T ***
  REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) !< Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      !< PE-local state ensemble
  INTEGER, INTENT(in) :: flag        !< PDAF status flag


! ****************************
! *** Perform pre/poststep ***
! ****************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE prepoststep_seek.F90: Implement prepoststep for SEEK here!'


END SUBROUTINE prepoststep_seek
