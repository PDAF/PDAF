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
!! * 2021-05 - Lars Nerger - Initial code based on offline_2D
!! * Later revisions - see repository log
!!
SUBROUTINE prepoststep_3dvar_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

  USE mod_model, &             ! Model variables
       ONLY: nx, ny
  USE mod_assimilation, &      ! Assimilation variables
       ONLY: dim_cvec, Vmat_p
  USE PDAF, &                  ! PDAF diagnostic routine
       ONLY: PDAF_diag_stddev_nompi

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
  INTEGER :: pdaf_status              ! status flag
  REAL :: ens_stddev                  ! estimated RMS error
  REAL, ALLOCATABLE :: field(:,:)     ! global model field
  CHARACTER(len=2) :: ensstr          ! String for ensemble member
  CHARACTER(len=2) :: stepstr         ! String for time step
  CHARACTER(len=3) :: anastr          ! String for call type (initial, forecast, analysis)


! **********************
! *** INITIALIZATION ***
! **********************

  IF (firsttime) THEN
     WRITE (*, '(8x, a)') 'Analyze initial state for 3D-Var'
     anastr = 'ini'
  ELSE
     IF (step<0) THEN
        WRITE (*, '(8x, a)') 'Analyze and write forecasted state for 3D-Var'
        anastr = 'for'
     ELSE
        WRITE (*, '(8x, a)') 'Analyze and write assimilated state for 3D-Var'
        anastr = 'ana'
     END IF
  END IF


! **************************************************************
! *** Compute sample standard deviation for Vmat_p           ***
! *** We use here the routine PDAF_diag_stddev_nompi which   ***
! *** is designed for ensembles. Since for Vmat_p the        ***
! *** sample normalization by SQRT(dim_cvec-1) is not valid, ***
! *** we multiply the result of the routine by this factor.  ***
! **************************************************************

  ! Set mean to zero to handle modes in Vmat_p
  state_p(:) = 0

  CALL PDAF_diag_stddev_nompi(dim_p, dim_cvec, state_p, Vmat_p, &
        ens_stddev, 0, pdaf_status)

  ! Scale for correct standard deviation
  ens_stddev = ens_stddev*SQRT(REAL(dim_cvec-1))


  ! *** Initialize state estimate (here we only have a single state)
  state_p(:) = ens_p(:,1)


! *****************
! *** Screen IO ***
! *****************

  ! Output RMS errors given by sampled covar matrix
  WRITE (*, '(12x, a, es12.4)') &
       'RMS error according to modes Vmat_p: ', ens_stddev


! *******************
! *** File output ***
! *******************

  IF (.not. firsttime) THEN

     WRITE (*, '(8x, a)') '--- write ensemble and state estimate'

     ALLOCATE(field(ny, nx))

     ! Set string for time step
     IF (step>=0) THEN
        WRITE (stepstr, '(i2.2)') step
     ELSE
        WRITE (stepstr, '(i2.2)') -step
     END IF

     ! Write analysis ensemble
     DO member = 1, dim_ens
        DO j = 1, nx
           field(1:ny, j) = ens_p(1 + (j-1)*ny : j*ny, member)
        END DO

        WRITE (ensstr, '(i2.2)') member

        OPEN(11, file = 'ens_'//TRIM(ensstr)//'_step'//TRIM(stepstr)//'_'//TRIM(anastr)//'.txt', status = 'replace')
 
        DO i = 1, ny
           WRITE (11, *) field(i, :)
        END DO

        CLOSE(11)
     END DO

     ! Write analysis state
     DO j = 1, nx
        field(1:ny, j) = state_p(1 + (j-1)*ny : j*ny)
     END DO

     OPEN(11, file = 'state_step'//TRIM(stepstr)//'_'//TRIM(anastr)//'.txt', status = 'replace')
 
     DO i = 1, ny
        WRITE (11, *) field(i, :)
     END DO

     CLOSE(11)


     DEALLOCATE(field)
  END IF


! ********************
! *** finishing up ***
! ********************

  firsttime = .FALSE.

END SUBROUTINE prepoststep_3dvar_pdaf
