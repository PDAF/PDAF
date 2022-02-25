!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz and Guowei He (Forschungszentrum Juelich GmbH)
!
!This file is part of TerrSysMP-PDAF
!
!TerrSysMP-PDAF is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TerrSysMP-PDAF is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU LesserGeneral Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public License
!along with TerrSysMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------
!init_parallel_pdaf.F90: TerrSysMP-PDAF implementation of routine
!                        'init_parallel_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: init_parallel_pdaf.F90 1442 2013-10-04 10:35:19Z lnerger $
!BOP
!
! !ROUTINE: init_parallel_pdaf --- Initialize communicators for PDAF
!
! !INTERFACE:
SUBROUTINE init_parallel_pdaf(dim_ens, screen)

! !DESCRIPTION:
! Parallelization routine for a model with 
! attached PDAF. The subroutine is called in 
! the main program subsequently to the 
! initialization of MPI. It initializes
! MPI communicators for the model tasks, filter 
! tasks and the coupling between model and
! filter tasks. In addition some other variables 
! for the parallelization are initialized.
! The communicators and variables are handed
! over to PDAF in the call to 
! PDAF\_filter\_init.
!
! 3 Communicators are generated:\\
! - COMM\_filter: Communicator in which the
!   filter itself operates\\
! - COMM\_model: Communicators for parallel
!   model forecasts\\
! - COMM\_couple: Communicator for coupling
!   between models and filter\\
! Other variables that have to be initialized are:\\
! - filterpe - Logical: Does the PE execute the 
! filter?\\
! - my\_ensemble - Integer: The index of the PE's 
! model task\\
! - local\_npes\_model - Integer array holding 
! numbers of PEs per model task
!
! For COMM\_filter and COMM\_model also
! the size of the communicators (npes\_filter and 
! npes\_model) and the rank of each PE 
! (mype\_filter, mype\_model) are initialized. 
! These variables can be used in the model part 
! of the program, but are not handed over to PDAF.
!
! This variant is for a domain decomposed 
! model.
!
! This is a template that is expected to work 
! with many domain-decomposed models. However, 
! it might be necessary to adapt the routine 
! for a particular model. Inportant is that the
! communicator COMM_model equals the communicator 
! used in the model. If one plans to run the
! ensemble forecast in parallel COMM_model cannot 
! be MPI_COMM_WORLD! Thus, if the model uses 
! MPI_COMM_WORLD it has to be replaced by an 
! alternative communicator named, e.g., COMM_model.
!
! !REVISION HISTORY:
! 2004-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! TEMPLATE: The include from mod_parallel_model has to be adapted to a particular model
 ! gh, universal pe world
   USE mod_parallel_model, &
        ONLY: MPI_COMM_WORLD, mype_model, npes_model, COMM_model, &
            mype_world, npes_world
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, npes_filter, COMM_filter, filterpe, n_modeltasks, &
       local_npes_model, task_id, COMM_couple, MPIerr
  ! Un-comment the following 2 lines in case a serial model is used 
!   USE mod_parallel_pdaf, &
!        ONLY: MPI_COMM_WORLD, mype_model, npes_model, COMM_model, npes_world, mype_world
!   then remove local declaration of npes_world and mype_world
  USE parser, &
       ONLY: parse
#if defined use_comm_da
  USE mod_oasis_data, &
       only: da_comm
#endif

#if (defined CLMSA)
  use enkf_clm_mod, only: da_comm_clm
#endif

#if (defined COUP_OAS_COS)
  use mod_parallel_pdaf, only : task_id
  use mod_tsmp, only: nprocpf, nprocclm
  use data_parallel, only: cosmo_input_suffix
#endif

  IMPLICIT NONE    
  
! !ARGUMENTS:
  INTEGER, INTENT(inout) :: dim_ens ! Ensemble size or number of EOFs (only SEEK)
  ! Often dim_ens=0 when calling this routine, because the real ensemble size
  ! is initialized later in the program. For dim_ens=0 no consistency check
  ! for ensemble size with number of model tasks is performed.
  INTEGER, INTENT(in)    :: screen ! Whether screen information is shown

! !CALLING SEQUENCE:
! Called by: main program
! Calls: MPI_Comm_size
! Calls: MPI_Comm_rank
! Calls: MPI_Comm_split
! Calls: MPI_Barrier
!EOP

  ! local variables
  INTEGER :: i, j               ! Counters
  INTEGER :: COMM_ensemble      ! Communicator of all PEs doing model tasks
  INTEGER :: mype_ens, npes_ens ! rank and size in COMM_ensemble
  INTEGER :: mype_couple, npes_couple ! Rank and size in COMM_couple
  INTEGER :: pe_index           ! Index of PE
  INTEGER :: my_color, color_couple ! Variables for communicator-splitting 
  LOGICAL :: iniflag            ! Flag whether MPI is initialized
  CHARACTER(len=32) :: handle   ! handle for command line parser
  ! gh, universal pe world
  !INTEGER :: npes_world, mype_world   ! Rank and size on MPI_COMM_WORLD


  ! *** Initialize MPI if not yet initialized ***
  CALL MPI_Initialized(iniflag, MPIerr)
  IF (.not.iniflag) THEN
     CALL MPI_Init(MPIerr)
  END IF

  ! *** Initialize PE information on COMM_world ***
  CALL MPI_Comm_size(MPI_COMM_WORLD, npes_world, MPIerr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, mype_world, MPIerr)

  ! *** Parse number of model tasks ***
  handle = 'n_modeltasks'
  CALL parse(handle, n_modeltasks)


  ! *** Initialize communicators for ensemble evaluations ***
  IF (mype_world == 0) &
       WRITE (*, '(/1x, a)') 'Initialize communicators for assimilation with PDAF'


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
  consist3: IF (modulo(npes_world, n_modeltasks) \= 0) THEN
    ! *** # parallel tasks is a divisor of # available PEs ***
    IF (mype_world == 0) WRITE (*, '(3x, a)') &
        'Error: n_modeltasks should be a divisor of npes_world.'
    stop
  END IF consist3


  ! ***              COMM_ENSEMBLE                ***
  ! *** Generate communicator for ensemble runs   ***
  ! *** only used to generate model communicators ***
  COMM_ensemble = MPI_COMM_WORLD

  npes_ens = npes_world
  mype_ens = mype_world


  ! *** Store # PEs per ensemble                 ***
  ! *** used for info on PE 0 and for generation ***
  ! *** of model communicators on other Pes      ***
  IF (ALLOCATED(local_npes_model)) DEALLOCATE(local_npes_model)
  ALLOCATE(local_npes_model(n_modeltasks))

  local_npes_model = FLOOR(REAL(npes_world) / REAL(n_modeltasks))
  DO i = 1, (npes_world - n_modeltasks * local_npes_model(1))
     local_npes_model(i) = local_npes_model(i) + 1
  END DO
  

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
    write (*,*) 'MODEL: mype(w)= ', mype_world, '; model task: ', task_id, &
         '; mype(m)= ', mype_model, '; npes(m)= ', npes_model
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
  my_color = task_id

  CALL MPI_Comm_split(MPI_COMM_WORLD, my_color, mype_world, &
       COMM_filter, MPIerr)

  ! *** Initialize PE informations         ***
  ! *** according to coupling communicator ***
  CALL MPI_Comm_Size(COMM_filter, npes_filter, MPIerr)
  CALL MPI_Comm_Rank(COMM_filter, mype_filter, MPIerr)

  ! ***              COMM_COUPLE                 ***
  ! *** Generate communicators for communication ***
  ! *** between model and filter PEs             ***
  ! *** (Split COMM_ENSEMBLE)                    ***

  color_couple = mype_filter + 1

  CALL MPI_Comm_split(MPI_COMM_WORLD, color_couple, mype_world, &
       COMM_couple, MPIerr)

  ! *** Initialize PE informations         ***
  ! *** according to coupling communicator ***
  CALL MPI_Comm_Size(COMM_couple, npes_couple, MPIerr)
  CALL MPI_Comm_Rank(COMM_couple, mype_couple, MPIerr)

  IF (screen > 0) THEN
     IF (mype_world == 0) THEN
        WRITE (*, '(/18x, a)') 'PE configuration:'
        WRITE (*, '(2x, a6, a9, a10, a14, a13, /2x, a5, a9, a7, a7, a7, a7, a7, /2x, a)') &
             'world', 'filter', 'model', 'couple', 'filterPE', &
             'rank', 'rank', 'task', 'rank', 'task', 'rank', 'T/F', &
             '----------------------------------------------------------'
     END IF
     CALL MPI_Barrier(MPI_COMM_WORLD, MPIerr)
     IF (task_id == 1) THEN
        WRITE (*, '(2x, i4, 4x, i4, 4x, i3, 4x, i3, 4x, i3, 4x, i3, 5x, l3)') &
             mype_world, mype_filter, task_id, mype_model, color_couple, &
             mype_couple, filterpe
     ENDIF
     IF (task_id > 1) THEN
        WRITE (*,'(2x, i4, 12x, i3, 4x, i3, 4x, i3, 4x, i3, 5x, l3)') &
         mype_world, task_id, mype_model, color_couple, mype_couple, filterpe
     END IF
     CALL MPI_Barrier(MPI_COMM_WORLD, MPIerr)

     IF (mype_world == 0) WRITE (*, '(/a)') ''

  END IF



! ******************************************************************************
! *** Initialize model equivalents to COMM_model, npes_model, and mype_model ***
! ******************************************************************************

  ! If the names of the variables for COMM_model, npes_model, and 
  ! mype_model are different in the numerical model, the 
  ! model-internal variables should be initialized at this point.
!
#if defined use_comm_da
  da_comm = comm_model
#endif

#if (defined CLMSA)
    da_comm_clm = comm_model
#endif

    ! set certain variables in component models
#if (defined COUP_OAS_COS)
    if(mype_model.ge.(nprocpf+nprocclm)) then
        cosmo_input_suffix = task_id-1
    end if
#endif

END SUBROUTINE init_parallel_pdaf
