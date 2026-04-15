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
SUBROUTINE prepoststep_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

  USE mod_model, &             ! Model variables
       ONLY: nx, ny
  USE mod_statevector_pdaf, &  ! Statevector variables
       ONLY: n_fields, sfields, id
  USE mod_parallel_pdaf, &     ! Parallelization variables
       ONLY: COMM_filter
  USE PDAF, &                  ! PDAF diagnostic routine
       ONLY: PDAF_diag_stddev

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
  INTEGER :: i, j, member, fid        ! Counters
  INTEGER :: istart, iend             ! stard and end index of a field in state vector
  INTEGER :: pdaf_status              ! status flag
  LOGICAL, SAVE :: firsttime = .TRUE. ! Routine is called for first time?
  REAL, ALLOCATABLE :: ens_stddev(:)  ! estimated RMS errors
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


! ************************************************************
! *** Compute ensemble mean and standard deviation         ***
! *** (=RMS errors according to sampled covar matrix)      ***
! ************************************************************

  ! Allocate fields
  ALLOCATE(ens_stddev(n_fields))

  ! Compute ensemble deviation and mean separately
  ! for each field in the state vector
  DO j = 1, n_fields
     ! Start and end index
     istart = 1 + sfields(j)%off
     iend = sfields(j)%dim + sfields(j)%off

     CALL PDAF_diag_stddev(sfields(j)%dim, dim_ens, &
          state_p(istart:iend), ens_p(istart:iend,:), &
          ens_stddev(j), 1, COMM_filter, pdaf_status)
  END DO


! *****************
! *** Screen IO ***
! *****************

  ! Output RMS errors given by sampled covar matrix
  WRITE (*, '(8x,a)') 'Ensemble standard deviation (estimated RMS error)'
  DO i = 1, n_fields
     WRITE (*,'(8x,a13,4x,a10,2x,es12.4)') &
          'stddev-'//anastr, TRIM(sfields(i)%name), ens_stddev(i)
  END DO

  DEALLOCATE(ens_stddev)

  
! *******************
! *** File output ***
! *******************

  IF (.NOT. firsttime) THEN

     WRITE (*, '(8x, a)') '--- write ensemble and state estimate'

     ALLOCATE(field_tmp(ny, nx))

     ! Set string for time step
     IF (step>=0) THEN
        WRITE (stepstr, '(i2.2)') step
     ELSE
        WRITE (stepstr, '(i2.2)') -step
     END IF

     ! Write analysis ensemble fields

     DO fid = 1, n_fields

        DO member = 1, dim_ens

           DO j = 1, nx
              field_tmp(1:ny, j) = ens_p(sfields(fid)%off + 1 + (j-1)*ny : sfields(fid)%off + j*ny, member)
           END DO

           WRITE (ensstr, '(i2.2)') member
           OPEN(11, file = 'ens'//TRIM(sfields(fid)%name)//'_'//TRIM(ensstr)//'_step'//TRIM(stepstr)//'_'//TRIM(anastr)//'.txt', &
                status = 'replace')

           DO i = 1, ny
              WRITE (11, *) field_tmp(i, :)
           END DO

           CLOSE(11)
        END DO
     END DO

     ! Write analysis state fields

     DO fid = 1, n_fields
        DO j = 1, nx
           field_tmp(1:ny, j) = state_p(sfields(fid)%off + 1 + (j-1)*ny : sfields(fid)%off + j*ny)
        END DO

        OPEN(11, file = 'state'//TRIM(sfields(fid)%name)//'_step'//TRIM(stepstr)//'_'//TRIM(anastr)//'.txt', status = 'replace')

        DO i = 1, ny
           WRITE (11, *) field_tmp(i, :)
        END DO

        CLOSE(11)
     END DO

     DEALLOCATE(field_tmp)
  END IF


! ********************
! *** finishing up ***
! ********************

  firsttime = .FALSE.

END SUBROUTINE prepoststep_pdaf
