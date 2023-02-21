!$Id: assimilation_pdaf_si.F90 1678 2016-12-11 12:32:53Z lnerger $
!BOP
!
! !ROUTINE: assimilation_pdaf - Routine controlling ensemble integration for PDAF
!
! !INTERFACE:
SUBROUTINE assimilation_pdaf(time)

! !DESCRIPTION:
! This routine performs the ensemble forecasts and includes
! the calls to PDAF to perform the analysis steps of the 
! data assimilation process.
!
! The model gets state vectors to be evolved as well as
! the number of time steps and the current model time 
! from PDAF by calling PDAF\_get\_state.
! Each forecasted state is written back into the ensemble 
! matrix of PDAF by calling a filter-specific routine
! PDAF\_put\_state\_X. When all ensemble members are 
! evolved and hence the forecast phase is completed, 
! PDAF\_put\_state\_X executes the analysis step of the
! chosen filter algorithm.
!
! This variant uses the simplified interface of PDAF. That is,
! no names of user-supplied subroutines have to be specified.
! However, this variant required that the user-supplied
! have specific names defined inside PDAF.
!
! This example is for the 1D dummy model. All 6 filters
! are implemented. If one only wants a single filter, one
! can delete all calls to PDAF_put_state_X except for
! those for the filter of choice. 
! In this implementation we consider the case that the time
! stepper of the model is implemented as a subroutine. If this
! is not the case, the operations of this routine have to be
! added to the model code around the time stepping loop of
! the model.
!
! !REVISION HISTORY:
! 2011-03 - Lars Nerger - Initial code based on assimilation_pdaf
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel, &     ! Parallelization
       ONLY: Comm_model, MPIerr, mype_world, abort_parallel
  USE mod_assimilation, & ! Variables for assimilation
       ONLY: filtertype

  IMPLICIT NONE

! !ARGUMENTS:
  REAL, INTENT(INOUT) :: time  ! Model time

! !CALLING SEQUENCE:
! Called by: main
! Calls: PDAF_get_state_si
! Calls: integration
! Calls: PDAF_put_state_seek_si
! Calls: PDAF_put_state_seik_si
! Calls: PDAF_put_state_enkf_si
! Calls: PDAF_put_state_lseik_si
! Calls: MPI_barrier (MPI)
!EOP

! local variables
  INTEGER :: nsteps    ! Number of time steps to be performed in current forecast
  INTEGER :: doexit    ! Whether to exit forecasting (1=true)
  INTEGER :: status    ! Status flag for filter routines
  REAL    :: timenow   ! Current model time


! *************************
! *** Perform forecasts ***
! *************************

  ! PDAF: External loop around model time stepper loop
  pdaf_ensembleloop: DO  

     ! *** PDAF: Get state and forecast information (nsteps,time) ***
     CALL PDAF_get_state_si(nsteps, timenow, doexit, status)

     ! Check whether forecast has to be performed
     checkforecast: IF (doexit /= 1 .AND. status == 0) THEN

        ! *** Forecast ensemble state ***
      
        IF (nsteps > 0) THEN

           ! Initialize current model time
           time = timenow

           ! *** call time stepper ***  
           CALL integration(time, nsteps)

        END IF

        ! *** PDAF: Send state forecast to filter;                         ***
        ! *** PDAF: Perform assimilation if ensemble forecast is completed ***
        ! *** PDAF: Distinct calls for each filter                         ***
        IF (filtertype == 0) THEN
           CALL PDAF_put_state_seek_si(status)
        ELSE IF (filtertype == 1) THEN
           CALL PDAF_put_state_seik_si(status)
        ELSE IF (filtertype == 2) THEN
           CALL PDAF_put_state_enkf_si(status)
        ELSE IF (filtertype == 3) THEN
           CALL PDAF_put_state_lseik_si(status)
        ELSE IF (filtertype == 4) THEN
           CALL PDAF_put_state_etkf_si(status)
        ELSE IF (filtertype == 5) THEN
           CALL PDAF_put_state_letkf_si(status)
        ELSE IF (filtertype == 6) THEN
           CALL PDAF_put_state_estkf_si(status)
        ELSE IF (filtertype == 7) THEN
           CALL PDAF_put_state_lestkf_si(status)
        ELSE IF (filtertype == 8) THEN
           CALL PDAF_put_state_lenkf_si(status)
        ELSE IF (filtertype == 9) THEN
           CALL PDAF_put_state_netf_si(status)
        ELSE IF (filtertype == 10) THEN
           CALL PDAF_put_state_lnetf_si(status)
        END IF

        CALL MPI_barrier(COMM_model, MPIERR)

     ELSE checkforecast

        ! *** No more work, exit modeling loop
        EXIT pdaf_ensembleloop

     END IF checkforecast

  END DO pdaf_ensembleloop


! ************************
! *** Check error flag ***
! ************************

  IF (status /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a47,i4,a1/)') &
          'ERROR ', status, &
          ' during assimilation with PDAF - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF

END SUBROUTINE assimilation_pdaf
