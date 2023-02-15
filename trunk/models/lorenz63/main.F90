!$Id$
!BOP
!
! !ROUTINE: main --- Main driver for PDAF testsuite
!
! !INTERFACE:
PROGRAM MAIN

! !DESCRIPTION:
! This is the main program for the PDAF testsuite for
! PDAF with domain-decomposition. 
! The additions to the model code when PDAF is attached
! are visible from the preprocessor directives for USE_PDAF.
!
! Parameters can be set in the code, or - preferably -
! by command line arguments that are parsed by the 
! module PARSER. The format for this is
! EXECUTABLE -HANDLE1 VALUE1 -HANDLE2 VALUE2 ...
! The handles are defined in the code before the calls to
! the routine PARSE.
!
! !REVISION HISTORY:
! 2009-05 - Lars Nerger - Initial code based on dummy model example
! Later revisions - see svn log
!
! !USES:
  USE mpi
  USE mod_parallel, &     ! Parallelization variables
       ONLY: MPIerr, npes_world, mype_world, n_modeltasks, &
       init_parallel, finalize_parallel
  USE mod_modeltime, &    ! Model time information
       ONLY: time, total_steps
  USE timer, &            ! Timing
       ONLY: timeit, time_tot
  USE mod_memcount, &     ! Counting allocated memory
       ONLY: memcount, memcount_ini, memcount_get
#if defined USE_PDAF
  USE pdaf_interfaces_module, &   ! Interface definitions to PDAF core routines
       ONLY: PDAF_print_info, PDAF_deallocate
  USE parser, &           ! Parse command lines
       ONLY: handle, parse
#endif
#if defined (_OPENMP)
  USE omp_lib, &
       ONLY: omp_get_num_threads, omp_get_thread_num
#endif

  IMPLICIT NONE
!EOP


! Local variables
  INTEGER :: i                         ! Counter
  INTEGER, SAVE :: mythread, nthreads  ! Thread variables for OpenMP
!$OMP THREADPRIVATE(mythread, nthreads)


! ********************************
! ***      INITIALIZATION      ***
! ********************************

! *** Initialize MPI parallelization
  CALL init_parallel()

#if defined (_OPENMP)
! Initialize threading information in case of OpenMP
!$OMP PARALLEL DEFAULT(none)
  nthreads = omp_get_num_threads()
  mythread = omp_get_thread_num()
!$OMP END PARALLEL
#else
  nthreads = 1
  mythread = 0
#endif

! *** Initial Screen output ***
  initscreen: IF (mype_world == 0) THEN

     WRITE (*, '(/16x, a/)') '+++++ PDAF test suite +++++'
#ifndef USE_PDAF
     WRITE (*, '(18x, a)') 'Run only model forecast'
#else
     WRITE (*, '(16x, a)') 'Data assimilation with PDAF'
#endif

     IF (npes_world > 1) THEN
        WRITE (*, '(/21x, a, i3, a/)') 'Running on ', npes_world, ' PEs'
     ELSE
        WRITE (*, '(/21x, a/)') 'Running on 1 PE'
     END IF
     IF (nthreads > 1) THEN
        WRITE (*, '(15x, a, i3, a/)') 'Use OpenMP with ', nthreads, ' threads'
     END IF
     WRITE (*, '(/)')
     
  END IF initscreen

! *** set number of timers ***
  CALL timeit(6, 'ini')

! *** set first timer ***
  CALL timeit(1, 'new')

! *** set number of memory counters ***
  CALL memcount_ini(3)

#ifdef USE_PDAF
  ! *** Parse settings for parallelization with PDAF ***  
  handle = 'tasks'                   ! Set number of concurrent model tasks
  CALL parse(handle, n_modeltasks)

! *** Initialize MPI communicators for PDAF (model and filter) ***
  CALL init_parallel_pdaf(0, 1)
#endif

! *** Initialize model ***
  CALL timeit(2, 'new')
  CALL initialize()
  CALL timeit(2, 'old')

#ifdef USE_PDAF
! *** Initialize PDAF ***
  CALL timeit(4, 'new')
  CALL init_pdaf()
  CALL timeit(4, 'old')
#endif


! *****************************
! ***      Integration      ***
! *****************************

  CALL timeit(3, 'new')

  ! *** Perform integration
#ifndef USE_PDAF
  ! Normal integration without assimilation
  IF (mype_world == 0) WRITE (*, '(/1x, a)') 'PDAF test suite: START INTEGRATION'
  CALL integration(time, total_steps)
#else
  ! With PDAF perform assimilation with ensemble integration
  ! Note: As we have the model time stepper as a subroutine, we simply push the
  ! routine assimilation_pdaf between the main program and the integration routine.
  IF (mype_world == 0) WRITE (*, '(/1x, a)') 'PDAF test suite: START ASSIMILATION'
  CALL assimilation_pdaf(time)
#endif

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
#ifndef USE_PDAF
     WRITE (*, '(/1x, a)') 'PDAF test suite: EXITED INTEGRATIONS'
#else
     WRITE (*, '(/1x, a)') 'PDAF test suite: EXITED ASSIMILATION'
#endif

     ! *** Show allocated memory for the model ***
     WRITE (*, '(/18x, a)') 'Model - Memory overview'
     WRITE (*, '(10x, 45a)') ('-', i=1, 45)
     WRITE (*, '(21x, a, f10.3, a)') 'Allocated memory  (MB)'
     WRITE (*, '(14x, a, f10.5, a)') &
          'Model field:', memcount_get(1, 'M'), ' MB (persistent)'
#ifdef USE_PDAF
     WRITE (*, '(12x, a, f10.5, a)') &
          'ensemble init:', memcount_get(2, 'M'), ' MB (temporary)'
     WRITE (*, '(13x, a, f10.5, a)') &
          'Pre-Poststep:', memcount_get(3, 'M'), ' MB (temporary)'

     ! Show allocated memory for PDAF
     CALL PDAF_print_info(10)
#endif

     ! *** Print timings onto screen ***
#ifdef USE_PDAF
     ! Show timings for PDAF
     CALL PDAF_print_info(1)
     CALL PDAF_print_info(3)
#endif
     WRITE (*, '(/17x, a)') 'Model - Timing information'
     WRITE (*, '(10x, 45a)') ('-', i=1, 45)
#ifndef USE_PDAF
     ! Timing summary for pure model integration
     WRITE (*, '(19x, a, F11.3, 1x, a)') 'initialize model:', time_tot(2), 's'
     WRITE (*, '(24x, a, F11.3, 1x, a)') 'integration:', time_tot(3), 's'
     WRITE (*, '(19x, a, F11.3, 1x, a)') 'total run time:', time_tot(1), 's'
#else
     ! Timing summary for assimilation
     WRITE (*, '(19x, a, F10.3, 1x, a)') 'Time of forecasts:', time_tot(5), 's'
     WRITE (*, '(17x, a)') '------------ Summary ------------'
     WRITE (*, '(19x, a, F11.3, 1x, a)') 'initialize model:', time_tot(2), 's'
     WRITE (*, '(18x, a, F11.3, 1x, a)') 'initialize filter:', time_tot(4), 's'
     WRITE (*, '(23x, a, F11.3, 1x, a)') 'assimilation:', time_tot(3), 's'
     WRITE (*, '(19x, a, F11.3, 1x, a)') 'total run time:', time_tot(1), 's'
#endif

     WRITE (*, '(/1x, a)') 'PDAF test suite: END'
  END IF screen3

! *** deallocate timers ***
  CALL timeit(6, 'fin')

#ifdef USE_PDAF
  CALL PDAF_deallocate()
#endif

! *** Terminate MPI
  CALL finalize_parallel()

END PROGRAM MAIN
