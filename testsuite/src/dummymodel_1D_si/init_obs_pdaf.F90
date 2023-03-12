!$Id: init_obs_pdaf.F90 1253 2012-01-30 18:56:08Z lnerger $
!BOP
!
! !ROUTINE: init_obs_pdaf --- Initialize observation vector
!
! !INTERFACE:
SUBROUTINE init_obs_pdaf(step, dim_obs_p, observation_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/ETKF/ESTKF
!
! The routine is called during the analysis step. 
! It has to provide the PE-local observation vector 
! for the current time step.
!
! The routine is called by all filter processes.
!
! Version for the dummy model with domain decomposition.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel, &
       ONLY: mype_filter, npes_filter, COMM_filter, MPI_DOUBLE_PRECISION, &
       MPIerr, MPIstatus
  USE mod_assimilation, &
       ONLY: rms_obs
  USE mod_model, &
       ONLY: local_dims, dim_state, dt

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step             ! Current time step
  INTEGER, INTENT(in) :: dim_obs_p        ! PE-local dimension of obs. vector
  REAL, INTENT(out)   :: observation_p(dim_obs_p) ! PE-local observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_seek_analysis    (as U_init_obs)
! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT
! Called by: PDAF_enkf_obs_ensemble
! Calls: dlarnv (LAPACK)
! Calls: MPI_send (MPI)
! Calls: MPI_recv (MPI)
!EOP

! *** local variables ***
  INTEGER :: i                         ! Counter
  REAL, ALLOCATABLE :: observation(:)  ! Global observation vector
  REAL, ALLOCATABLE :: obs_errors(:)   ! Global array holding obs. errors
  INTEGER, SAVE :: first = 1 ! Flag for initialization of random number seed
  INTEGER, SAVE :: iseed(4)  ! Seed array for random number generator
  ! variables and arrays for domain decomposition
  INTEGER :: offset          ! Row-offset according to domain decomposition
  INTEGER :: domain          ! Domain counter
  REAL, ALLOCATABLE :: observation_temp_p(:)  ! temporary sub-array


! ******************************
! *** Initialize observation ***
! ******************************

  ! For the dummy model the observation is just the true state
  ! plus some Gaussian error.
  ! We generate first the global observation vector which is then
  ! distributed. This is motivated by the consistency with the 
  ! mode-decomposition variant with regard to the random-number
  ! generation.
  mype0: IF (mype_filter == 0) THEN
     WRITE (*, '(8x, a, i7)') &
          '--- initialize observation at time step', step

     ! Allocate global observation vector
     ALLOCATE(observation(dim_state))

     ! *** Compute global true state
     observation(:) = 1.0
     DO i = 1, step
        observation(:) = observation(:) + 1.0 * dt
     END DO

     ! *** generate array of observation errors
     ALLOCATE(obs_errors(dim_state))

     ! Initialized seed for random number routine
     IF (first == 1) THEN
        iseed(1) = 1000
        iseed(2) = 2034
        iseed(3) = 0
        iseed(4) = 3
        first = 2
     END IF

     ! generate random number with Gaussian distribution variance 1
     CALL dlarnv(3, iseed, dim_state, obs_errors(1 : dim_state))

     ! disturb true state
     DO i = 1, dim_state
        observation(i) = observation(i) + rms_obs * obs_errors(i)
     END DO

  END IF mype0


! ****************************
! *** Distribute substates ***
! ****************************

  mype0b: IF (mype_filter == 0) THEN
     ! *** Initialize and send sub-state on PE 0 ***

     offset = 0

     DO domain = 1, npes_filter

        whichdomain: IF (domain == 1) THEN
           ! Initialize sub-state and sub_ensemble for PE 0

           ! perform reordering of state for PE 0
           DO i = 1, local_dims(1)
              observation_p(i) = observation(i)
           END DO
    
        ELSE whichdomain
           ! Initialize sub-observation vectors for other PEs
           ! and send sub-vectors

           ! allocate temporary sub-arrays
           ALLOCATE(observation_temp_p(local_dims(domain)))

           ! perform reordering of observation vector
           DO i = 1, local_dims(domain)
              observation_temp_p(i) = observation(i + offset)
           END DO

           ! Send sub-arrays
           CALL MPI_send(observation_temp_p, local_dims(domain), &
             MPI_DOUBLE_PRECISION, domain - 1, 2, COMM_filter, MPIerr)

           DEALLOCATE(observation_temp_p)

        END IF whichdomain

        ! Increment offset
        offset = offset + local_dims(domain)

     END DO

  ELSE mype0b
     ! *** Receive substate on filter-PEs with rank > 0 ***

     CALL MPI_recv(observation_p, local_dims(mype_filter + 1), &
          MPI_DOUBLE_PRECISION, 0, 2, COMM_filter, MPIstatus, MPIerr)

  END IF mype0b


! ********************
! *** Finishing up ***
! ********************

  IF (mype_filter == 0) DEALLOCATE(obs_errors, observation)

END SUBROUTINE init_obs_pdaf

