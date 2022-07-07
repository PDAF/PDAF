!$Id: integration.F90 1135 2011-08-24 08:35:17Z lnerger $
!BOP
!
! !ROUTINE: integration  --- integration routine for the 1D dummy model
!
! !INTERFACE:
SUBROUTINE integration(time, nsteps)

! !DESCRIPTION:
! Routine to perform model integration with the 1D dummy model.
!
! For simplicity of the implementation with PDAF,
! the time stepping is separated into a single routine.
! This allows to simply put the assimilation routine
! assimilation\_pdaf() in between the main program and
! the integration routine. If the time stepping is not
! available as a separate routine, a different 
! implementation style is required adding the 
! functionality of assimilation\_pdaf() to the model
! code before enclosing the time stepping.
!
! !REVISION HISTORY:
! 2004-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE timer, &            ! Timing
       ONLY: timeit      
  USE mod_parallel, &     ! Parallelization
       ONLY: mype_model  
  USE mod_model, &        ! Model variables
       ONLY: field, dt
#ifdef USE_PDAF
  USE mod_assimilation, & ! Variables for assimilation
       ONLY: filtertype, incremental
#endif

  IMPLICIT NONE

! !ARGUMENTS:
  REAL, INTENT(inout) :: time   ! Model time
  INTEGER, INTENT(in) :: nsteps ! Number of time steps to be performed
!EOP

! local variables
  INTEGER :: step               ! Time step counter
  INTEGER :: doexit             ! Whether to exit forecasting (1=true)


! *********************************
! *** Perform model integration ***
! *********************************

  CALL timeit(5, 'new')

#ifndef USE_PDAF
!  Show model state element 1 if no assimilation is done
  IF (mype_model == 0) &
       WRITE (*,'(1x, i3, 1x, a, i6, 1x, f10.3)') &
       mype_model, 'step, field(1)', 0, field(1)
#endif


! *** time stepping loop ***
  integrate: DO step = 1, nsteps
     
#ifdef USE_PDAF
     ! For incremental updating (SEEK, SEIK, and LSEIK)
     IF (incremental == 1 &
          .AND. (filtertype==0 .OR. filtertype == 1 .OR. filtertype == 3)) THEN
        CALL PDAF_incremental_si(nsteps)
     END IF

  ! *** PDAF: Add model error ***      
  ! This is the place where model error can be simulated.
  ! This routine is not implemented here!
!   call add_model_noise(step)
#endif

! *** model time step ***
     field(:) = field(:) + 1.0 * dt


#ifndef USE_PDAF
!  Show model state element 1 if no assimilation is done
     IF (mype_model == 0) &
          WRITE (*,'(1x, i3, 1x, a, i6, 1x, f10.3)') &
          mype_model, 'step, field(1)', step, field(1)
#endif

     ! Increment time
     time = time + dt

  END DO integrate

  CALL timeit(5, 'old')

END SUBROUTINE integration
