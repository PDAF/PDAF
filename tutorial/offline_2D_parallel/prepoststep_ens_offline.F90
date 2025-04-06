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
!! For the offline mode, this routine is the place
!! in which the writing of the analysis ensemble
!! can be performed.
!!
!! If a user considers to perform adjustments to the 
!! estimates (e.g. for balances), this routine is 
!! the right place for it.
!!
!! Implementation for the 2D offline example
!! with model parallelization.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code based on offline_1D
!! * Later revisions - see repository log
!!
SUBROUTINE prepoststep_ens_offline(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

  USE mpi                      ! MPI
  USE mod_parallel_pdaf, &     ! Parallelization
       ONLY: mype_filter, npes_filter, COMM_filter, MPIerr, MPIstatus
  USE mod_assimilation, &      ! Assimilation variables
       ONLY: nx, ny, dim_state, local_dims
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
  INTEGER :: i, j, member, domain     ! Counters
  INTEGER :: pdaf_status              ! status flag
  LOGICAL, SAVE :: firsttime = .TRUE. ! Routine is called for first time?
  REAL :: ens_stddev                  ! ensemble STDDEV = estimated RMS error
  INTEGER :: nobs                     ! Number of observations in diagnostics
  REAL, POINTER :: obsRMSD(:)         ! Array of observation RMS deviations
  REAL, ALLOCATABLE :: field(:,:)     ! global model field
  CHARACTER(len=2) :: ensstr          ! String for ensemble member
  ! Variables for parallelization - global fields
  INTEGER :: offset   ! Row-offset according to domain decomposition
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


! ************************************************************
! *** Compute ensemble mean and standard deviation         ***
! *** (=RMS errors according to sampled covar matrix)      ***
! ************************************************************

  CALL PDAF_diag_stddev(dim_p, dim_ens, state_p, ens_p, &
        ens_stddev, 1, COMM_filter, pdaf_status)


! *****************
! *** Screen IO ***
! *****************

  ! Output RMS errors given by sampled covar matrix
  IF (mype_filter == 0) THEN
     WRITE (*, '(12x, a, es12.4)') &
          'RMS error according to sampled standard deviation: ', ens_stddev
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
        offset = dim_p

        DO domain = 2, npes_filter
           ! Initialize sub-ensemble for other PEs and send sub-arrays

           ! Allocate temporary buffer array
           ALLOCATE(ens_p_tmp(local_dims(domain), dim_ens))

           ! Receive sub-arrays
           CALL MPI_recv(ens_p_tmp, local_dims(domain) * dim_ens, MPI_DOUBLE_PRECISION, &
                domain - 1, 1, COMM_filter, MPIstatus, MPIerr)

           ! Initialize MPI buffer for local ensemble
           DO member = 1, dim_ens
              DO i = 1, local_dims(domain)
                 ens(i + offset, member) = ens_p_tmp(i, member)
              END DO
           END DO

           DEALLOCATE(ens_p_tmp)

           ! Increment offset
           offset = offset + local_dims(domain)

        END DO


        ! *** Now write analysis ensemble ***

        WRITE (*, '(8x, a)') '--- write ensemble and state estimate'

        ALLOCATE(field(ny, nx))

        DO member = 1, dim_ens
           DO j = 1, nx
              field(1:ny, j) = ens(1 + (j-1)*ny : j*ny, member)
           END DO

           WRITE (ensstr, '(i2.2)') member
           OPEN(11, file = 'ens_'//TRIM(ensstr)//'_ana.txt', status = 'replace')
 
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
        offset = dim_p

        DO domain = 2, npes_filter
           ! Initialize sub-ensemble for other PEs and send sub-arrays

           ! Receive sub-arrays
           CALL MPI_recv(state(1+offset), local_dims(domain), MPI_DOUBLE_PRECISION, &
                domain - 1, 1, COMM_filter, MPIstatus, MPIerr)

           ! Increment offset
           offset = offset + local_dims(domain)

        END DO
     
        ! *** Now write analysis state estimate ***

        DO j = 1, nx
           field(1:ny, j) = state(1 + (j-1)*ny : j*ny)
        END DO

        OPEN(11, file = 'state_ana.txt', status = 'replace')
 
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

END SUBROUTINE prepoststep_ens_offline
