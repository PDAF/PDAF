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
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code based on offline_1D
!! * Later revisions - see repository log
!!
SUBROUTINE prepoststep_ens_offline(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

  USE mpi                      ! MPI
  USE mod_parallel, &          ! Parallelization
       ONLY: mype_filter, npes_filter, COMM_filter, MPIerr, MPIstatus
  USE mod_assimilation, &      ! Assimilation variables
       ONLY: dim_state

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
  INTEGER :: i, j, member, domain     ! Counters
  LOGICAL, SAVE :: firstio = .TRUE.   ! File output is peformed for first time?
  LOGICAL, SAVE :: firsttime = .TRUE. ! Routine is called for first time?
  REAL :: invdim_ens                  ! Inverse ensemble size
  REAL :: invdim_ensm1                ! Inverse of ensemble size minus 1
  REAL :: rmserror_est                ! estimated RMS error
  REAL, ALLOCATABLE :: variance_p(:)  ! model state variances
  REAL, ALLOCATABLE :: field(:,:)     ! global model field
  CHARACTER(len=2) :: ensstr          ! String for ensemble member
  ! Variables for parallelization - global fields
  INTEGER :: offset   ! Row-offset according to domain decomposition
  REAL, ALLOCATABLE :: variance(:)    ! local variance
  REAL, ALLOCATABLE :: ens(:,:)       ! global ensemble
  REAL, ALLOCATABLE :: state(:)       ! global state vector
  REAL,ALLOCATABLE :: ens_p_tmp(:,:)  ! Temporary ensemble for some PE-domain
  REAL,ALLOCATABLE :: state_p_tmp(:)  ! Temporary state for some PE-domain


! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype_filter == 0) THEN
     IF (firsttime) THEN
        WRITE (*, '(8x, a)') 'Analyze forecasted state ensemble'
     ELSE
        WRITE (*, '(8x, a)') 'Analyze and write assimilated state ensemble'
     END IF
  END IF
  ! Allocate fields
  ALLOCATE(variance_p(dim_p))
  ALLOCATE(variance(dim_state))

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
  IF (mype_filter == 0) WRITE (*, '(8x, a)') '--- compute ensemble mean'

  state_p = 0.0
  DO member = 1, dim_ens
     DO i = 1, dim_p
        state_p(i) = state_p(i) + ens_p(i, member)
     END DO
  END DO
  state_p(:) = invdim_ens * state_p(:)

  ! *** Compute sampled variances ***
  variance_p(:) = 0.0
  DO member = 1, dim_ens
     DO j = 1, dim_p
        variance_p(j) = variance_p(j) &
             + (ens_p(j, member) - state_p(j)) &
             * (ens_p(j, member) - state_p(j))
     END DO
  END DO
  variance_p(:) = invdim_ensm1 * variance_p(:)


! ******************************************************
! *** Assemble global variance vector on filter PE 0 ***
! ******************************************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE prepoststep_ens_offline.F90: Initialize variance, either directly or with MPI'

!   PE0_a: IF (mype_filter /= 0) THEN
! 
!      ! send sub-fields from PEs /=0
!      CALL MPI_send(variance_p(1 : dim_p), dim_p, &
!           MPI_DOUBLE_PRECISION,0, mype_filter, COMM_filter, MPIerr)
! 
!   ELSE PE0_a
!      ! receive and assemble variance field
! 
!      ! On PE 0 init variance directly
!      variance(1 : dim_p) = variance_p(1 : dim_p)
! 
!      ! Receive part of variance field from PEs > 0 into 
!      ! correct part of global variance
! 
!      offset = 0
! 
!      DO i = 2, npes_filter
!         ! Increment offset
!         offset = offset + local_dims(i - 1)
! 
!         ! Receive variance part
!         CALL MPI_recv(variance(1 + offset), local_dims(i), &
!              MPI_DOUBLE_PRECISION, i - 1, i - 1, COMM_filter, MPIstatus, MPIerr)
!      END DO
!       
!   END IF PE0_a

  DEALLOCATE(variance_p)


! ************************************************************
! *** Compute RMS errors according to sampled covar matrix ***
! ************************************************************

  ! Total estimated RMS error
  ! This example is univariate - one should distinguish different fields

  DO i = 1, dim_state
     rmserror_est = rmserror_est + variance(i)
  ENDDO
  rmserror_est = SQRT(rmserror_est / dim_state)

  DEALLOCATE(variance)


! *****************
! *** Screen IO ***
! *****************

  ! Output RMS errors given by sampled covar matrix
  IF (mype_filter == 0) THEN
     WRITE (*, '(12x, a, es12.4)') &
       'sampled ensemble standard deviation: ', rmserror_est
  END IF

 
! *******************
! *** File output ***
! *******************

  notfirst: IF (.not. firsttime) THEN

     ! Template reminder - delete when implementing functionality
     WRITE (*,*) 'TEMPLATE prepoststep_ens_offline.F90: Implement writing of output files here!'

  END IF notfirst



! ********************
! *** finishing up ***
! ********************

  firsttime = .FALSE.

END SUBROUTINE prepoststep_ens_offline
