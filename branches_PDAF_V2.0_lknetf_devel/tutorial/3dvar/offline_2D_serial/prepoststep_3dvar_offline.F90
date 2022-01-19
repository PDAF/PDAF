!$Id$
!>  Used-defined Pre/Poststep routine for PDAF
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in: 3D-Var
!! 
!! The routine is called for the 3D-Var with parameterized
!! covariances. It is called before and after the analysis.
!! The routine provides full access to the state 
!! estimate to the user.
!!
!! The routine provides full access to the state 
!! estimate and the state ensemble to the user.
!! Thus, user-controlled pre- and poststep 
!! operations can be performed here. For example 
!! the forecast and the analysis states and ensemble
!! covariance matrix can be analyzed, e.g. by 
!! computing the estimated variances. 
!! For the offline mode, this routine is the place
!! in which the writing of the analysis ensemble
!! can be performed.
!!
!! If a user considers to perform adjustments to the 
!! estimates (e.g. for balances), this routine is 
!! the right place for it.
!!
!! Implementation for the 2D example for 3D-Var
!! without model parallelization.
!!
!! __Revision history:__
!! * 2021-05 - Lars Nerger - Initial code based on prepoststep_ens_offline
!! * Later revisions - see repository log
!!
SUBROUTINE prepoststep_3dvar_offline(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

  USE mod_assimilation, &
       ONLY: nx, ny, dim_cvec, Vmat_p

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step        !< Current time step (negative for call after forecast)
  INTEGER, INTENT(in) :: dim_p       !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     !< Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   !< PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   !< PE-local dimension of observation vector
  REAL, INTENT(inout) :: state_p(dim_p) !< PE-local forecast/analysis state
  !< (The array 'state_p' is not generally not initialized in the case of SEIK.
  !< It can be used freely here.)
  REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) !< Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      !< PE-local state ensemble
  INTEGER, INTENT(in) :: flag        !< PDAF status flag


! *** local variables ***
  INTEGER :: i, j, member             ! Counters
  LOGICAL, SAVE :: firsttime = .TRUE. ! Routine is called for first time?
  REAL :: invdim_ens                  ! Inverse ensemble size
  REAL :: invdim_ensm1                ! Inverse of ensemble size minus 1
  REAL :: rmserror_est                ! estimated RMS error
  REAL, ALLOCATABLE :: variance(:)    ! model state variances
  REAL, ALLOCATABLE :: field(:,:)     ! global model field
  CHARACTER(len=2) :: ensstr          ! String for ensemble member
  REAL :: fact                        ! Scaling factor


! **********************
! *** INITIALIZATION ***
! **********************

  IF (firsttime) THEN
     WRITE (*, '(8x, a)') 'Analyze forecasted state for 3D-Var'
  ELSE
     WRITE (*, '(8x, a)') 'Analyze and write assimilated state for 3D-Var'
  END IF

  ! Allocate fields
  ALLOCATE(variance(dim_p))

  ! Initialize numbers
  rmserror_est  = 0.0
  invdim_ens    = 1.0 / REAL(dim_ens)  
  invdim_ensm1  = 1.0 / REAL(dim_ens - 1)


! **************************************************************
! *** Perform prepoststep for 3D-Var in which dim_ens=1      ***
! *** The sampled error is here computed from B^(1/2)        ***
! **************************************************************

  ! *** Initialize state estimate (here we only have a single state)
  state_p(:) = ens_p(:,1)

  ! *** Compute sampled variances ***
  variance(:) = 0.0
  DO member = 1, dim_cvec
     DO j = 1, dim_p
        variance(j) = variance(j) &
             + Vmat_p(j,member) * Vmat_p(j,member)
     END DO
  END DO


! ************************************************************
! *** Compute RMS errors according to sampled covar matrix ***
! ************************************************************

  ! total estimated RMS error
  DO i = 1, dim_p
     rmserror_est = rmserror_est + variance(i)
  ENDDO
  rmserror_est = SQRT(rmserror_est / dim_p)


! *****************
! *** Screen IO ***
! *****************

  ! Output RMS errors given by sampled covar matrix
  WRITE (*, '(12x, a, es12.4)') &
       'RMS error according to sampled variance: ', rmserror_est


! *******************
! *** File output ***
! *******************

  IF (.not. firsttime) THEN

     WRITE (*, '(8x, a)') '--- write ensemble and state estimate'

     ALLOCATE(field(ny, nx))

     ! Write analysis ensemble
     DO member = 1, dim_ens
        DO j = 1, nx
           field(1:ny, j) = ens_p(1 + (j-1)*ny : j*ny, member)
        END DO

        WRITE (ensstr, '(i2.2)') member
        OPEN(11, file = 'ens_'//TRIM(ensstr)//'_ana.txt', status = 'replace')
 
        DO i = 1, ny
           WRITE (11, *) field(i, :)
        END DO

        CLOSE(11)
     END DO

     ! Write analysis state
     DO j = 1, nx
        field(1:ny, j) = state_p(1 + (j-1)*ny : j*ny)
     END DO

     OPEN(11, file = 'state_ana.txt', status = 'replace')
 
     DO i = 1, ny
        WRITE (11, *) field(i, :)
     END DO

     CLOSE(11)


     DEALLOCATE(field)
  END IF


! ********************
! *** finishing up ***
! ********************

  DEALLOCATE(variance)

  firsttime = .FALSE.

END SUBROUTINE prepoststep_3dvar_offline
