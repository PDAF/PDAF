!$Id$
!>  Used-defined Pre/Poststep routine for PDAF
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all ensemble filters.
!! 
!! The routine is called for global filters (e.g. ESTKF)
!! before the analysis and after the ensemble transformation.
!! For local filters (e.g. LESTKF) the routine is called
!! before and after the loop over all local analysis
!! domains.
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
!! Implementation for the 2D example
!! without model parallelization.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code based on offline_1D
!! * Later revisions - see repository log
!!
SUBROUTINE prepoststep_ens_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

  USE mod_model, &
       ONLY: nx, ny
  USE mod_assimilation, &
       ONLY: n_fields, dim_fields, off_fields, id

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
  REAL, ALLOCATABLE :: rmserror_est(:) ! estimated RMS errors
  REAL, ALLOCATABLE :: variance(:)    ! model state variances
  REAL, ALLOCATABLE :: field_tmp(:,:) ! global model field for file output
  CHARACTER(len=2) :: ensstr          ! String for ensemble member
  CHARACTER(len=2) :: stepstr         ! String for time step
  CHARACTER(len=3) :: anastr          ! String for call type (initial, forecast, analysis)


! **********************
! *** INITIALIZATION ***
! **********************

  IF (firsttime) THEN
     WRITE (*, '(8x, a)') 'Analyze initial state ensemble'
     anastr = 'ini'
  ELSE
     IF (step<0) THEN
        WRITE (*, '(8x, a)') 'Analyze and write forecasted state ensemble'
        anastr = 'for'
     ELSE
        WRITE (*, '(8x, a)') 'Analyze and write assimilated state ensemble'
        anastr = 'ana'
     END IF
  END IF

  ! Allocate fields
  ALLOCATE(variance(dim_p))
  ALLOCATE(rmserror_est(n_fields))

  ! Initialize numbers
  rmserror_est  = 0.0
  invdim_ens    = 1.0 / REAL(dim_ens)  
  invdim_ensm1  = 1.0 / REAL(dim_ens - 1)


! **************************************************************
! *** Perform prepoststep for SEIK with re-inititialization. ***
! *** The state and error information is completely in the   ***
! *** ensemble.                                              ***
! *** Also performed for SEIK without re-init at the initial ***
! *** time.                                                  ***
! **************************************************************

  ! *** Compute mean state
  WRITE (*, '(8x, a)') '--- compute ensemble mean'

  state_p = 0.0
  DO member = 1, dim_ens
     DO i = 1, dim_p
        state_p(i) = state_p(i) + ens_p(i, member)
     END DO
  END DO
  state_p(:) = invdim_ens * state_p(:)

  ! *** Compute sampled variances ***
  variance(:) = 0.0
  DO member = 1, dim_ens
     DO j = 1, dim_p
        variance(j) = variance(j) &
             + (ens_p(j, member) - state_p(j)) &
             * (ens_p(j, member) - state_p(j))
     END DO
  END DO
  variance(:) = invdim_ensm1 * variance(:)


! ************************************************************
! *** Compute RMS errors according to sampled covar matrix ***
! ************************************************************

  DO j = 1, n_fields
     ! total estimated RMS error per field
     DO i = 1+off_fields(j), dim_fields(j)+off_fields(j)
        rmserror_est(j) = rmserror_est(j) + variance(i)
     ENDDO
     rmserror_est(j) = SQRT(rmserror_est(j) / dim_fields(j))
  ENDDO

  DEALLOCATE(variance)


! *****************
! *** Screen IO ***
! *****************

  ! Output RMS errors given by sampled covar matrix
  WRITE (*, '(12x, a, 2es12.4)') &
       'RMS errors according to sampled variance: ', rmserror_est

  DEALLOCATE(rmserror_est)

  
! *******************
! *** File output ***
! *******************

  IF (.not. firsttime) THEN

     WRITE (*, '(8x, a)') '--- write ensemble and state estimate'

     ALLOCATE(field_tmp(ny, nx))

     ! Set string for time step
     IF (step>=0) THEN
        WRITE (stepstr, '(i2.2)') step
     ELSE
        WRITE (stepstr, '(i2.2)') -step
     END IF

     ! Write analysis ensemble fields
     DO member = 1, dim_ens

        ! Field
        DO j = 1, nx
           field_tmp(1:ny, j) = ens_p(off_fields(id%fieldA) + 1 + (j-1)*ny : off_fields(id%fieldA) + j*ny, member)
        END DO

        WRITE (ensstr, '(i2.2)') member

        OPEN(11, file = 'ens_'//TRIM(ensstr)//'_step'//TRIM(stepstr)//'_'//TRIM(anastr)//'.txt', status = 'replace')
 
        DO i = 1, ny
           WRITE (11, *) field_tmp(i, :)
        END DO

        CLOSE(11)

        ! FieldB
        DO j = 1, nx
           field_tmp(1:ny, j) = ens_p(off_fields(id%fieldB) + 1 + (j-1)*ny : off_fields(id%fieldB) + j*ny, member)
        END DO

        WRITE (ensstr, '(i2.2)') member

        OPEN(12, file = 'ensB_'//TRIM(ensstr)//'_step'//TRIM(stepstr)//'_'//TRIM(anastr)//'.txt', status = 'replace')
 
        DO i = 1, ny
           WRITE (12, *) field_tmp(i, :)
        END DO

        CLOSE(12)
     END DO

     ! Write analysis state fields

     ! Field
     DO j = 1, nx
        field_tmp(1:ny, j) = state_p(off_fields(id%fieldA) + 1 + (j-1)*ny : off_fields(id%fieldA) + j*ny)
     END DO

     OPEN(11, file = 'state_step'//TRIM(stepstr)//'_'//TRIM(anastr)//'.txt', status = 'replace')
 
     DO i = 1, ny
        WRITE (11, *) field_tmp(i, :)
     END DO

     CLOSE(11)

     ! FieldB
     DO j = 1, nx
        field_tmp(1:ny, j) = state_p(off_fields(id%fieldB) + 1 + (j-1)*ny : off_fields(id%fieldB) + j*ny)
     END DO

     OPEN(12, file = 'stateB_step'//TRIM(stepstr)//'_'//TRIM(anastr)//'.txt', status = 'replace')
 
     DO i = 1, ny
        WRITE (12, *) field_tmp(i, :)
     END DO

     CLOSE(12)

     DEALLOCATE(field_tmp)
  END IF


! ********************
! *** finishing up ***
! ********************

  firsttime = .FALSE.

END SUBROUTINE prepoststep_ens_pdaf
