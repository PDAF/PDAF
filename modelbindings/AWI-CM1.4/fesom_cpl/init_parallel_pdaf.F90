!$Id: init_parallel_pdaf.F90 2271 2020-04-08 13:04:09Z lnerger $
!>  Routine to initialize communicators for PDAF
!!
!! Parallelization routine for a model with attached PDAF.
!! The subroutine is called in the main program subsequently
!! to the initialization of MPI. It initializes MPI
!! communicators for the model tasks, filter tasks and the
!! coupling between model and filter tasks. In addition some
!! other variables for the parallelization are initialized.
!! The communicators and variables are handed over to PDAF
!! in the call to PDAF_init.
!!
!! 3 Communicators are generated:
!! - COMM_filter: Communicator in which the filter itself operates
!! - COMM_model: Communicators for parallel model forecasts
!! - COMM_couple: Communicator for coupling between models and filter
!! Other variables that have to be initialized are:
!! - filterpe - Logical: Does the process execute the filter?
!! - my_ensemble - Integer: The index of the model task of this process
!! - local_npes_model - Integer array holding numbers of processes 
!!     per model task
!!
!! For COMM_filter and COMM_model also the size of the
!! communicators (npes_filter and npes_model) and the rank
!! of each process (mype\_filter, mype\_model) are initialized. 
!! These variables can be used in the model part of the program, 
!! but are not handed over to PDAF.
!!
!! For AWI-CM we initialize the model-specific variables
!! comm_cplmod, cplmodid, and strcplmodid 
!! after the generic initializations.
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
SUBROUTINE init_parallel_pdaf(dim_ens_in, screen)

  USE mod_parallel_pdaf, &        ! Parallelization variables
       ONLY: mype_filter, npes_filter, COMM_filter, filterpe, n_modeltasks, &
       local_npes_model, task_id, COMM_couple, MPI_UNDEFINED,  &
       mype_world, npes_world, MPI_COMM_WORLD, mype_model, npes_model, &
       COMM_model, MPIerr, abort_parallel, pairs
  USE mod_assim_pdaf, &           ! Variables for assimilation
       ONLY: dim_ens, da_couple_type
  use mod_oasis_data, only: comm_cplmod, cplmodid, strcplmodid

  IMPLICIT NONE    
  
! *** Arguments ***
  INTEGER, INTENT(inout) :: dim_ens_in !< Ensemble size or number of EOFs (only SEEK)
  ! Often dim_ens=0 when calling this routine, because the real ensemble size
  ! is initialized later in the program. For dim_ens=0 no consistency check
  ! for ensemble size with number of model tasks is performed.
  INTEGER, INTENT(in)    :: screen     !< Whether screen information is shown

! *** Local variables ***
  INTEGER :: i, j               ! Counters
  INTEGER :: COMM_ensemble      ! Communicator of all PEs doing model tasks
  INTEGER :: mype_ens, npes_ens ! rank and size in COMM_ensemble
  INTEGER :: mype_couple, npes_couple ! Rank and size in COMM_couple
  INTEGER :: pe_index           ! Index of PE
  INTEGER :: my_color, color_couple ! Variables for communicator-splitting 
  LOGICAL :: iniflag            ! Flag whether MPI is initialized
  CHARACTER(len=32) :: handle   ! handle for command line parser
  INTEGER :: ncpus_atm          ! number of processes for each atmosphere task
  INTEGER :: ncpus_ocn          ! number of processes for each ocean task
  
  namelist /pdaf_parallel/ dim_ens, ncpus_atm, ncpus_ocn, pairs, da_couple_type


  ! *** Initialize MPI if not yet initialized ***
  CALL MPI_Initialized(iniflag, MPIerr)
  IF (.not.iniflag) THEN
     CALL MPI_Init(MPIerr)
  END IF

  ! *** Initialize PE information on COMM_world ***
  CALL MPI_Comm_size(MPI_COMM_WORLD, npes_world, MPIerr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, mype_world, MPIerr)

  ! *** Get paralelization information for ensemble runs ***
  OPEN (10, file='namelist.pdaf')
  READ (10, NML=pdaf_parallel)
  CLOSE (10)
  n_modeltasks = dim_ens

  ! *** Initialize communicators for ensemble evaluations ***
  IF (mype_world == 0) &
       WRITE (*, '(/1x, a)') 'FESOM-PDAF: Initialize communicators for assimilation with PDAF'

  IF (mype_world == 0) THEN
     WRITE (*,'(a,3x,a,i6)') 'FESOM-PDAF', 'model tasks: ', n_modeltasks
     WRITE (*,'(a,3x,a,i6)') 'FESOM-PDAF', 'CPUs per ocean task:      ', ncpus_ocn
     WRITE (*,'(a,3x,a,i6)') 'FESOM-PDAF', 'CPUs per atmosphere task: ', ncpus_atm
  END IF


  ! *** Check consistency of number of parallel ensemble tasks ***
  consist1: IF (n_modeltasks > npes_world) THEN
     ! *** # parallel tasks is set larger than available PEs ***
     n_modeltasks = npes_world
     IF (mype_world == 0) WRITE (*, '(3x, a)') &
          '!!! Resetting number of parallel ensemble tasks to total number of PEs!'
  END IF consist1
  IF (dim_ens > 0) THEN
     ! Check consistency with ensemble size
     consist2: IF (n_modeltasks > dim_ens) THEN
        ! # parallel ensemble tasks is set larger than ensemble size
        n_modeltasks = dim_ens
        IF (mype_world == 0) WRITE (*, '(5x, a)') &
             '!!! Resetting number of parallel ensemble tasks to number of ensemble states!'
     END IF consist2
  END IF


  ! ***              COMM_ENSEMBLE                ***
  ! *** Generate communicator for ensemble runs   ***
  ! *** only used to generate model communicators ***
  COMM_ensemble = MPI_COMM_WORLD

  npes_ens = npes_world
  mype_ens = mype_world


  ! *** Store # PEs per ensemble                 ***
  ! *** used for info on PE 0 and for generation ***
  ! *** of model communicators on other Pes      ***

  ALLOCATE(local_npes_model(n_modeltasks))

  IF (pairs) THEN
    local_npes_model = ncpus_ocn + ncpus_atm 
  ELSE
    local_npes_model = ncpus_ocn
  END IF

  IF (n_modeltasks*(ncpus_ocn+ncpus_atm) /= npes_world) THEN
     WRITE (*,*) 'ERROR: Inconsistent setting of processes for ensemble integrations!'
     call abort_parallel()
  END IF


  ! ***              COMM_MODEL               ***
  ! *** Generate communicators for model runs ***
  ! *** (Split COMM_ENSEMBLE)                 ***
  pe_index = 0
  doens1: DO i = 1, n_modeltasks
     DO j = 1, local_npes_model(i)
        IF (mype_ens == pe_index) THEN
           task_id = i
           EXIT doens1
        END IF
        pe_index = pe_index + 1
     END DO
  END DO doens1


  CALL MPI_Comm_split(COMM_ensemble, task_id, mype_ens, &
       COMM_model, MPIerr)

  ! *** Re-initialize PE informations   ***
  ! *** according to model communicator ***
  CALL MPI_Comm_Size(COMM_model, npes_model, MPIerr)
  CALL MPI_Comm_Rank(COMM_model, mype_model, MPIerr)

  if (screen > 1) then
    write (*,'(a,i5,a,i5,a,i5,a,i5)') 'FESOM-PDAF: mype(w)= ', mype_world, &
         '; model task: ', task_id, '; mype(m)= ', mype_model, '; npes(m)= ', npes_model
  end if


  ! Init flag FILTERPE (all PEs of model task 1)
  IF (task_id == 1) THEN
     filterpe = .TRUE.
  ELSE
     filterpe = .FALSE.
  END IF

  ! ***         COMM_FILTER                 ***
  ! *** Generate communicator for filter    ***
  ! *** For simplicity equal to COMM_couple ***
  IF (filterpe) THEN
     my_color = task_id
  ELSE
     my_color = MPI_UNDEFINED
  ENDIF

  CALL MPI_Comm_split(MPI_COMM_WORLD, my_color, mype_world, &
       COMM_filter, MPIerr)

  ! *** Initialize PE informations         ***
  ! *** according to coupling communicator ***
  IF (filterpe) THEN
     CALL MPI_Comm_Size(COMM_filter, npes_filter, MPIerr)
     CALL MPI_Comm_Rank(COMM_filter, mype_filter, MPIerr)
  ENDIF


  ! ***              COMM_COUPLE                 ***
  ! *** Generate communicators for communication ***
  ! *** between model and filter PEs             ***
  ! *** (Split COMM_ENSEMBLE)                    ***

  color_couple = mype_model + 1

  CALL MPI_Comm_split(MPI_COMM_WORLD, color_couple, mype_world, &
       COMM_couple, MPIerr)

  ! *** Initialize PE informations         ***
  ! *** according to coupling communicator ***
  CALL MPI_Comm_Size(COMM_couple, npes_couple, MPIerr)
  CALL MPI_Comm_Rank(COMM_couple, mype_couple, MPIerr)

  IF (screen > 0) THEN
     IF (mype_world == 0) THEN
        WRITE (*, '(/a, 18x,a)') 'FESOM-PDAF Pconf','PE configuration:'
        WRITE (*, '(a, 2x, a6, a9, a10, a14, a13, /a, 2x, a5, a9, a7, a7, a7, a7, a7, /a, 2x, a)') &
             'FESOM-PDAF Pconf', 'world', 'filter', 'model', 'couple', 'filterPE', &
             'FESOM-PDAF Pconf', 'rank', 'rank', 'task', 'rank', 'task', 'rank', 'T/F', &
             'FESOM-PDAF Pconf', '----------------------------------------------------------'
     END IF

     CALL MPI_Barrier(MPI_COMM_WORLD, MPIerr)
     IF (task_id == 1) THEN
        WRITE (*, '(a, 2x, i6, 3x, i4, 4x, i3, 4x, i3, 4x, i3, 4x, i3, 4x, l3)') &
             'FESOM-PDAF Pconf', mype_world, mype_filter, task_id, mype_model, color_couple, &
             mype_couple, filterpe
     ENDIF
     IF (task_id > 1) THEN
        WRITE (*,'(a, 2x, i6, 11x, i3, 4x, i3, 4x, i3, 4x, i3, 4x, l3)') &
         'FESOM-PDAF Pconf', mype_world, task_id, mype_model, color_couple, mype_couple, filterpe
     END IF
     CALL MPI_Barrier(MPI_COMM_WORLD, MPIerr)

     IF (mype_world == 0) WRITE (*, '(/a)') ''

  END IF


! ******************************************************************************
! *** Initialize model equivalents to COMM_model, npes_model, and mype_model ***
! ******************************************************************************

  ! Set communicator for each coupled model task
  comm_cplmod = comm_model

  ! Store task ID
  cplmodid = task_id

  ! Store string of coupled task ID
  write (strcplmodid,'(a1,i5.5)') '_',cplmodid

END SUBROUTINE init_parallel_pdaf
