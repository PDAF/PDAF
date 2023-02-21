!$Id: prepoststep_3dvar_pdaf.F90 1864 2017-12-20 19:53:30Z lnerger $
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
!! with model parallelization.
!!
!! __Revision history:__
!! * 2013-04 - Lars Nerger - Initial code based on offline_1D
!! * Later revisions - see repository log
!!
SUBROUTINE prepoststep_3dvar_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

  USE mpi                  ! MPI
  USE mod_model, &         ! Model variables
       ONLY: nx, ny, nx_p
  USE mod_assimilation, &  ! Assimilation variables
       ONLY: dim_state, Vmat_p, dim_cvec
  USE mod_parallel_pdaf, & ! Model parallelization variables
       ONLY: mype_filter, npes_filter, COMM_filter, &
       MPIerr, MPIstatus

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
  LOGICAL, SAVE :: firsttime = .TRUE. ! Routine is called for first time?
  REAL :: invdim_ens                  ! Inverse ensemble size
  REAL :: invdim_ensm1                ! Inverse of ensemble size minus 1
  REAL :: rmserror_est                ! estimated RMS error
  REAL, ALLOCATABLE :: variance_p(:)  ! model state variances
  REAL, ALLOCATABLE :: field(:,:)     ! global model field
  CHARACTER(len=2) :: ensstr          ! String for ensemble member
  CHARACTER(len=2) :: stepstr         ! String for time step
  CHARACTER(len=3) :: anastr          ! String for call type (initial, forecast, analysis)
  ! Variables for parallelization - global fields
  INTEGER :: off_p   ! Row-offset according to domain decomposition
  REAL, ALLOCATABLE :: variance(:)    ! local variance
  REAL, ALLOCATABLE :: ens(:,:)       ! global ensemble
  REAL, ALLOCATABLE :: state(:)       ! global state vector
  REAL,ALLOCATABLE :: ens_p_tmp(:,:)  ! Temporary ensemble for some PE-domain


! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype_filter == 0) THEN
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
  END IF

  ! Allocate fields
  ALLOCATE(variance_p(dim_p))
  ALLOCATE(variance(dim_state))

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
  variance_p(:) = 0.0
  DO member = 1, dim_cvec
     DO j = 1, dim_p
        variance_p(j) = variance_p(j) &
             + Vmat_p(j,member) * Vmat_p(j,member)
     END DO
  END DO


! ******************************************************
! *** Assemble global variance vector on filter PE 0 ***
! ******************************************************

  PE0_a: IF (mype_filter /= 0) THEN

     ! send sub-fields from PEs /=0
     CALL MPI_send(variance_p(1 : dim_p), dim_p, &
          MPI_DOUBLE_PRECISION,0, mype_filter, COMM_filter, MPIerr)

  ELSE PE0_a
     ! receive and assemble variance field

     ! On PE 0 init variance directly
     variance(1 : dim_p) = variance_p(1 : dim_p)

     ! Receive part of variance field from PEs > 0 into 
     ! correct part of global variance

     off_p = 0

     DO i = 2, npes_filter
        ! Increment offset
        off_p = off_p + nx_p*ny

        ! Receive variance part
        CALL MPI_recv(variance(1 + off_p), nx_p*ny, &
             MPI_DOUBLE_PRECISION, i - 1, i - 1, COMM_filter, MPIstatus, MPIerr)
     END DO
      
  END IF PE0_a

  DEALLOCATE(variance_p)


! ************************************************************
! *** Compute RMS errors according to sampled covar matrix ***
! ************************************************************

  ! total estimated RMS error
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
       'RMS error according to sampled variance: ', rmserror_est
  END IF

 
! *******************
! *** File output ***
! *******************

  notfirst: IF (.not. firsttime) THEN

     ALLOCATE(ens(dim_state, dim_ens))
     ALLOCATE(state(dim_state))

     ! Gather full ensemble on process with rank 0 and write file
     mype0b: IF (mype_filter /= 0) THEN

        ! *** Send ensemble substates on filter-PEs with rank > 0 ***

        CALL MPI_send(ens_p, dim_ens * dim_p, &
             MPI_DOUBLE_PRECISION, 0, 1, COMM_filter, MPIerr)

     ELSE mype0b

        ! *** Initialize and receive sub-states on PE 0 ***

        ! Initialize sub-ensemble for PE 0
        DO member = 1, dim_ens
           DO i=1, dim_p
              ens(i, member) = ens_p(i, member)
           END DO
        END DO

        ! Define offset in state vectors
        off_p = dim_p

        DO domain = 2, npes_filter
           ! Initialize sub-ensemble for other PEs and send sub-arrays

           ! Allocate temporary buffer array
           ALLOCATE(ens_p_tmp(nx_p*ny, dim_ens))

           ! Receive sub-arrays
           CALL MPI_recv(ens_p_tmp, nx_p*ny * dim_ens, MPI_DOUBLE_PRECISION, &
                domain - 1, 1, COMM_filter, MPIstatus, MPIerr)

           ! Initialize MPI buffer for local ensemble
           DO member = 1, dim_ens
              DO i = 1, nx_p*ny
                 ens(i + off_p, member) = ens_p_tmp(i, member)
              END DO
           END DO

           DEALLOCATE(ens_p_tmp)

           ! Increment offset
           off_p = off_p + nx_p*ny

        END DO


        ! *** Now write analysis ensemble ***

        WRITE (*, '(8x, a)') '--- write ensemble and state estimate'

        ! Set string for time step
        IF (step>=0) THEN
           WRITE (stepstr, '(i2.2)') step
        ELSE
           WRITE (stepstr, '(i2.2)') -step
        END IF

        ALLOCATE(field(ny, nx))

        DO member = 1, dim_ens
           DO j = 1, nx
              field(1:ny, j) = ens(1 + (j-1)*ny : j*ny, member)
           END DO

           WRITE (ensstr, '(i2.2)') member

           OPEN(11, file = 'ens_'//TRIM(ensstr)//'_step'//TRIM(stepstr)//'_'//TRIM(anastr)//'.txt', status = 'replace')
 
           DO i = 1, ny
              WRITE (11, *) field(i, :)
           END DO

           CLOSE(11)
        END DO

     END IF mype0b

     ! Gather full state vector on process with rank 0 and write to file
     mype0c: IF (mype_filter /= 0) THEN

        ! *** Send ensemble substates on filter-PEs with rank > 0 ***

        CALL MPI_send(state_p, dim_p, &
             MPI_DOUBLE_PRECISION, 0, 1, COMM_filter, MPIerr)

     ELSE mype0c

        ! *** Initialize and receive sub-states on PE 0 ***

        ! Initialize sub-state for PE 0
        DO i = 1, dim_p
           state(i) = state_p(i)
        END DO

        ! Define offset in state vectors
        off_p = dim_p

        DO domain = 2, npes_filter
           ! Initialize sub-ensemble for other PEs and send sub-arrays

           ! Receive sub-arrays
           CALL MPI_recv(state(1+off_p), nx_p*ny, MPI_DOUBLE_PRECISION, &
                domain - 1, 1, COMM_filter, MPIstatus, MPIerr)

           ! Increment offset
           off_p = off_p + nx_p*ny

        END DO
     
        ! *** Now write analysis state estimate ***

        DO j = 1, nx
           field(1:ny, j) = state(1 + (j-1)*ny : j*ny)
        END DO

        OPEN(11, file = 'state_step'//TRIM(stepstr)//'_'//TRIM(anastr)//'.txt', status = 'replace')
 
        DO i = 1, ny
           WRITE (11, *) field(i, :)
        END DO

        CLOSE(11)

        DEALLOCATE(field)
     END IF mype0c

     DEALLOCATE(ens, state)

  END IF notfirst


! ********************
! *** finishing up ***
! ********************

  firsttime = .FALSE.

END SUBROUTINE prepoststep_3dvar_pdaf
