!$Id: main.F90 731 2009-06-17 08:39:07Z lnerger $
!BOP
!
! !ROUTINE: main --- Main program for example of PDAF offline implementation
!
! !INTERFACE:
PROGRAM MAIN_OFFLINE

! !DESCRIPTION:
! This is the main program for an example implementation of
! PDAF with domain-decomposition and offline configuration.
! That is, we assume that the ensemble integrations are 
! performed in a different program and the forecasted 
! ensemble can be read-in from files. After initializing 
! the ensemble information by reading model outputs,
! a single analysis step is performed. Subsequently, the 
! analysis ensemble can be written to files that can be 
! used to initialize another ensemble forecast.
!
! Using PDAF for domain-decomposition, the offline
! mode can be used to perform assimilation with domain-
! decomposed models. If the models write results for each 
! sub-domain, these can be read here using the same 
! parallelization. Then, the filter analysis can be 
! performed utilizing this parallelization. If the files
! contain the full model state, PDAF in offline mode
! can be used either on a single processor, or the 
! fields can be distributed in this program to utilize
! the parallelization of the filters.
!
! Parameters can be set in the code, or - preferably -
! by command line arguments that are parsed by the 
! module PARSER. The format for this is
! EXECUTABLE -HANDLE1 VALUE1 -HANDLE2 VALUE2 ...
! The handles are defined in the code before the calls
! to the routine PARSE.
!
! !REVISION HISTORY:
! 2008-07 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel, &     ! Parallelization variables
       ONLY: MPI_COMM_WORLD, MPIerr, npes_world, mype_world, &
       init_parallel, finalize_parallel
  USE timer, &            ! Timing
       ONLY: timeit, time_tot
  USE mod_memcount, &     ! Counting allocated memory
       ONLY: memcount, memcount_ini, memcount_get

  IMPLICIT NONE
!EOP


! Local variables
  INTEGER :: i                 ! Counter


! **********************
! *** Initialize MPI ***
! **********************

  CALL init_parallel() ! initializes MPI


! **********************************************************
! ***               PROGRAM CONTROL                      ***
! **********************************************************


! ********************************
! ***      INITIALIZATION      ***
! ********************************

! *** Initial Screen output ***
  initscreen: IF (mype_world == 0) THEN

     WRITE (*, '(/8x, a/)') '+++++ PDAF test suite - offline mode +++++'
     WRITE (*, '(16x, a)') 'Data assimilation with PDAF'

     IF (npes_world > 1) THEN
        WRITE (*, '(/21x, a, i3, a/)') 'Running on ', npes_world, ' PEs'
     ELSE
        WRITE (*, '(/21x, a/)') 'Running on 1 PE'
     END IF
     WRITE (*, '(/)')
     
  END IF initscreen

! *** set number of timers ***
  CALL timeit(6, 'ini')

! *** set first timer ***
  CALL timeit(1, 'new')

! *** set number of memory counters ***
  CALL memcount_ini(3)

  
! *** Initialize MPI communicators for PDAF (model and filter) ***
! *** NOTE: It is always n_modeltasks=1 for offline mode       ***
  CALL init_parallel_pdaf(0, 1)

! *** Initialize model information ***
! *** This should only be information on the model dimension
! *** Generally, this could be joined with init_pdaf.
  CALL timeit(2, 'new')
  CALL initialize()
  CALL timeit(2, 'old')


! *******************************
! ***      ASSIMILATION       ***
! *******************************

  ! *** Initialize PDAF ***
  CALL timeit(4, 'new')
  CALL init_pdaf()
  CALL timeit(4, 'old')

  ! *** Perform analysis ***
  CALL timeit(3, 'new')
  IF (mype_world == 0) &
       WRITE (*, '(/2x, a)') 'PDAF test suite - offline mode: START ASSIMILATION'
  CALL assimilation_pdaf_offline()

  ! Syncronize at barrier for exit
  CALL MPI_Barrier(MPI_COMM_WORLD, MPIerr) 
  WRITE (*,*) 'model PE exited: mype ', mype_world

  CALL timeit(3, 'old')


! ********************
! *** Finishing up ***
! ********************

  CALL timeit(1, 'old')

! *** Final screen output ***
  screen3: IF (mype_world == 0) THEN
     WRITE (*, '(/1x, a)') 'PDAF test suite: EXITED ASSIMILATION'

     ! *** Show allocated memory for the model ***
     WRITE (*, '(/18x, a)') 'Model - Memory overview'
     WRITE (*, '(10x, 45a)') ('-', i=1, 45)
     WRITE (*, '(21x, a, f10.3, a)') 'Allocated memory  (MB)'
     WRITE (*, '(14x, a, f10.5, a)') &
          'Model field:', memcount_get(1, 'M'), ' MB (persistent)'
     WRITE (*, '(12x, a, f10.5, a)') &
          'ensemble init:', memcount_get(2, 'M'), ' MB (temporary)'
     WRITE (*, '(13x, a, f10.5, a)') &
          'Pre-Poststep:', memcount_get(3, 'M'), ' MB (temporary)'

     ! Show allocated memory for PDAF
     CALL PDAF_print_info(2)

     ! *** Print timings onto screen ***

     ! Show timings for PDAF
     CALL PDAF_print_info(1)

     WRITE (*, '(/17x, a)') 'Offline - Timing information'
     WRITE (*, '(10x, 45a)') ('-', i=1, 45)
     ! Timing summary for assimilation
     WRITE (*, '(19x, a, F11.3, 1x, a)') 'initialize model:', time_tot(2), 's'
     WRITE (*, '(18x, a, F11.3, 1x, a)') 'initialize filter:', time_tot(4), 's'
     WRITE (*, '(23x, a, F11.3, 1x, a)') 'assimilation:', time_tot(3), 's'
     WRITE (*, '(19x, a, F11.3, 1x, a)') 'total run time:', time_tot(1), 's'

     WRITE (*, '(/1x, a)') 'PDAF test suite: END'
  END IF screen3

! *** deallocate timers ***
  CALL timeit(6, 'fin')

! *** Terminate MPI
  CALL finalize_parallel()

END PROGRAM MAIN_OFFLINE
