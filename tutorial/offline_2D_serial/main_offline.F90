!$Id: main_offline.F90 1864 2017-12-20 19:53:30Z lnerger $
!>  Main driver for PDAF offline tutorial
!!
!! This is the main program for an example implementation of
!! PDAF with domain-decomposition and offline configuration.
!!
!! In the offline mode, we assume that the ensemble
!! integrations are performed in a separate program (model)
!! and the forecasted ensemble can be read from files. After
!! initializing the ensemble information by reading model
!! outputs, a single analysis step is performed. Subsequently,
!! the analysis ensemble can be written to files that can be 
!! used to initialize another ensemble forecast.
!!
!! Using PDAF for domain-decomposition, the offline
!! mode can be used to perform assimilation with domain-
!! decomposed models. If the models write results for each 
!! sub-domain, these can be read here using the same 
!! parallelization. Then, the filter analysis can be 
!! performed utilizing this parallelization. If the files
!! contain the full model state, PDAF in offline mode
!! can be used either on a single processor, or the 
!! fields can be distributed in this program to utilize
!! the parallelization of the filters.
!!
!! Parameters can be set in the code, or - preferably -
!! by command line arguments that are parsed by the 
!! module PARSER. The format for this is
!! EXECUTABLE -HANDLE1 VALUE1 -HANDLE2 VALUE2 ...
!! The handles are defined in the code before the calls
!! to the routine PARSE.
!!
!! __Revision history:__
!! * 2008-07 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
PROGRAM MAIN_OFFLINE

  USE mpi                 ! MPI
  USE mod_parallel, &     ! Parallelization
       ONLY: MPIerr, npes_world, mype_world, &
       init_parallel, finalize_parallel

  IMPLICIT NONE


! **********************
! *** Initialize MPI ***
! **********************

  CALL init_parallel() ! initializes MPI


! ********************************
! ***      INITIALIZATION      ***
! ********************************

! *** Initial Screen output ***
  initscreen: IF (mype_world == 0) THEN

     WRITE (*, '(/8x, a/)') '+++++ PDAF offline mode +++++'
     WRITE (*, '(9x, a)') 'Data assimilation with PDAF'

     IF (npes_world > 1) THEN
        WRITE (*, '(/21x, a, i3, a/)') 'Running on ', npes_world, ' PEs'
     ELSE
        WRITE (*, '(/21x, a/)') 'Running on 1 PE'
     END IF
     WRITE (*, '(/)')
     
  END IF initscreen

  
! *** Initialize MPI communicators for PDAF (model and filter) ***
! *** NOTE: It is always n_modeltasks=1 for offline mode       ***

  CALL init_parallel_pdaf(0, 1)

! *** Initialize model information ***
! *** This should only be information on the model dimension
! *** Generally, this could be joined with init_pdaf.

  CALL initialize()


! *******************************
! ***      ASSIMILATION       ***
! *******************************

  ! *** Initialize PDAF ***

  CALL init_pdaf()


  ! *** Perform analysis ***

  IF (mype_world == 0) &
       WRITE (*, '(/2x, a)') 'PDAF offline mode: START ASSIMILATION'

  CALL assimilation_pdaf_offline()


  ! Synchronize at barrier for exit
  CALL MPI_Barrier(MPI_COMM_WORLD, MPIerr) 
  WRITE (*,*) 'model PE exited: mype ', mype_world


! ********************
! *** Finishing up ***
! ********************

! *** Final screen output ***
  IF (mype_world == 0) &
       WRITE (*, '(/1x, a)') 'PDAF offline mode: EXITED ASSIMILATION'

  ! *** Finalize PDAF - print memory and timing information
  CALL finalize_pdaf(0)

  IF (mype_world == 0) &
       WRITE (*, '(/1x, a)') 'PDAF offline mode: END'


! *** Terminate MPI
  CALL finalize_parallel()

END PROGRAM MAIN_OFFLINE
