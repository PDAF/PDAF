! Copyright (c) 2004-2021 Lars Nerger
!
! This file is part of PDAF.
!
! PDAF is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version
! 3 of the License, or (at your option) any later version.
!
! PDAF is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with PDAF.  If not, see <http://www.gnu.org/licenses/>.
!
!$Id$
!BOP
!
! !MODULE:
MODULE PDAF_mod_filterMPI

! !DESCRIPTION:
! This module provides variables for the parallelization of PDAF. 
! In addition, the routine PDAF\_init\_par to intialize 
! internal parallelization variables is provided.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2003-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  use mpi

  IMPLICIT NONE
  SAVE

!  INCLUDE 'mpif.h'

! !PUBLIC DATA MEMBERS:
  INTEGER :: mype_world, npes_world     ! PE information for MPI_COMM_world
  INTEGER :: mype_filter, npes_filter   ! PE information for COMM_filter
  INTEGER :: mype, npes                 ! Aliases to mype_filter, npes_filter
  INTEGER :: mype_couple, npes_couple   ! PE information for COMM_couple
  INTEGER :: mype_model, npes_model     ! PE rank in COMM_model
  INTEGER :: error, MPIerr              ! MPI error flags
  INTEGER :: dim_ens_l                  ! Ensemble size of my task
  INTEGER :: dim_eof_l                  ! Number of EOFs in my task (SEEK only)
  INTEGER :: COMM_filter, COMM_couple   ! MPI communicators
  INTEGER :: COMM_pdaf                  ! MPI communicator for all PEs involved in PDAF
  LOGICAL :: isset_comm_pdaf = .false.  ! Whether COMM_pdaf was set externally
  INTEGER :: task_id                    ! Which ensemble task I am belonging to
  LOGICAL :: filterpe                   ! Whether PE belongs to the filter PEs
  LOGICAL :: modelpe                    ! Whether PE belongs to the model PEs
  INTEGER :: n_modeltasks               ! Number of parallel model tasks
  INTEGER :: MPIstatus(MPI_STATUS_SIZE) ! Status array for MPI
  INTEGER, ALLOCATABLE :: all_dim_eof_l(:)    ! Number of EOFs per task
  INTEGER, ALLOCATABLE :: all_dis_eof_l(:)    ! Displacements
  INTEGER, ALLOCATABLE :: all_dim_ens_l(:)    ! Size of ensembles per task
  INTEGER, ALLOCATABLE :: all_dis_ens_l(:)    ! Displacements
  INTEGER, ALLOCATABLE :: all_npes_model_l(:) ! # PEs per model task
  INTEGER :: statetask = -1 ! SEEK: Index of model task holding the state forecast
  LOGICAL :: filter_no_model = .FALSE.
  INTEGER, ALLOCATABLE :: all_dim_obs_p(:)    ! PE-Local observation dimensions
  INTEGER, ALLOCATABLE :: all_dis_obs_p(:)    ! PE-Local observation displacements
  INTEGER :: dimobs_p, dimobs_f             ! PE-local and global observation dimension
!EOP

CONTAINS
!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: PDAF_init_parallel - Initialize parallelization variables for PDAF
!
! !INTERFACE: PDAF_init_parallel
  SUBROUTINE PDAF_init_parallel(dim_ens, ensemblefilter, fixedbasis, &
       COMM_model, in_COMM_filter, in_COMM_couple, &
       in_n_modeltasks, in_task_id, screen, flag)

! !DESCRIPTION:
! This subroutine initializes internal parallelization 
! information for PDAF.
!
! !USES:
    IMPLICIT NONE    

! !ARGUMENTS:
    INTEGER, INTENT(inout) :: dim_ens        ! Rank of covar matrix/ensemble size
    LOGICAL, INTENT(in) :: ensemblefilter    ! Is the filter ensemble-based?
    LOGICAL, INTENT(in) :: fixedbasis        ! Run with fixed error-space basis?
    INTEGER, INTENT(in) :: COMM_model        ! Model communicator (not shared)
    INTEGER, INTENT(in) :: in_COMM_filter    ! Filter communicator
    INTEGER, INTENT(in) :: in_COMM_couple    ! Coupling communicator
    INTEGER, INTENT(in) :: in_task_id        ! Task ID of current PE
    INTEGER, INTENT(in) :: in_n_modeltasks   ! Number of model tasks
    INTEGER, INTENT(in) :: screen            ! Whether screen information is shown
    INTEGER, INTENT(inout):: flag            ! Status flag

! !CALLING SEQUENCE:
! Called by: PDAF_filter_init
!EOP


! local variables
    INTEGER :: i, task                       ! Counters
    CHARACTER(len=200) :: tskstr1, tskstr2   ! Strings for ensemble overview
    CHARACTER(len=200) :: tskstr3, tskstrtmp ! Strings for ensemble overview


    ! Check parallelization setting
    IF (fixedbasis .AND. in_n_modeltasks > 1) THEN
       IF (mype_filter == 0) THEN
          WRITE (*, '(/5x, a)') 'PDAF-ERROR: Fixed basis filters can only be run'
          WRITE (*, '(5x, a)') 'PDAF-ERROR: with a single model task!'
          WRITE (*, '(5x, a)') 'PDAF-ERROR: STOPPING PROGRAM !!!'
       END IF
       CALL  MPI_Finalize(MPIerr)
       STOP
    END IF

    ! Initialize values
    task_id = in_task_id
    n_modeltasks = in_n_modeltasks

    ! Initialize model variables for communicators
    COMM_filter = in_COMM_filter
    COMM_couple = in_COMM_couple

    ! *** Initialize PE information on COMM_world ***
    CALL MPI_Comm_size(COMM_pdaf, npes_world, MPIerr)
    CALL MPI_Comm_rank(COMM_pdaf, mype_world, MPIerr)


    ! *** Initialize PE information on COMM_filter ***
    IF (COMM_filter /= MPI_COMM_NULL) THEN
       CALL MPI_Comm_Size(COMM_filter, npes_filter, MPIerr)
       CALL MPI_Comm_Rank(COMM_filter, mype_filter, MPIerr)

       ! For short in the filter routines
       mype = mype_filter
       npes = npes_filter
    END IF

    ! *** Initialize PE information on COMM_couple ***
    IF (COMM_couple /= MPI_COMM_NULL) THEN
       CALL MPI_Comm_Size(COMM_couple, npes_couple, MPIerr)
       CALL MPI_Comm_Rank(COMM_couple, mype_couple, MPIerr)
    ELSE
       npes_couple = 0
       mype_couple = -1
    END IF

    ! *** Initialize PE information on COMM_model ***
    IF (COMM_model /= MPI_COMM_NULL) THEN
       CALL MPI_Comm_Size(COMM_model, npes_model, MPIerr)
       CALL MPI_Comm_Rank(COMM_model, mype_model, MPIerr)
       modelpe = .TRUE.
    ELSE
       filter_no_model = .TRUE.
       modelpe = .FALSE.
    END IF

    ! Distribute information whether filter pes are also model pes
    CALL MPI_BCAST(filter_no_model, 1, MPI_LOGICAL, 0, COMM_pdaf, MPIerr)

    ! *** Get # PEs per ensemble ***
    ! *** used only for info     ***
    ALLOCATE(all_npes_model_l(n_modeltasks))

    IF (filter_no_model) THEN
       all_npes_model_l = FLOOR(REAL(npes_world-npes_filter) / REAL(n_modeltasks))
       DO i = 1, (npes_world - npes_filter - n_modeltasks * all_npes_model_l(1))
          all_npes_model_l(i) = all_npes_model_l(i) + 1
       END DO
    ELSE
       all_npes_model_l = FLOOR(REAL(npes_world) / REAL(n_modeltasks))
       DO i = 1, (npes_world - n_modeltasks * all_npes_model_l(1))
          all_npes_model_l(i) = all_npes_model_l(i) + 1
       END DO
    END IF

    IF (screen > 2 .AND. filterpe) THEN
       IF (filter_no_model) THEN
          WRITE (*,*) 'PDAF: FILTER-PE - no model task: ', &
               'mype(f)= ', mype_filter, '; npes(f)= ', npes_filter
       ELSE
          WRITE (*,*) 'PDAF: FILTER-PE - model task: ', task_id, &
               '; mype(w)= ', mype_world, '; mype(m)= ', mype_model, &
               '; mype(f)= ', mype_filter, '; npes(f)= ', npes_filter
       END IF
    ELSE IF (screen > 2) THEN
       WRITE (*,*) 'PDAF: MODEL-PE - model task: ', task_id,'; mype(w)= ', &
            mype_world, '; mype(m)= ', mype_model
    END IF

    ! *** store local ensemble sizes and displacements on filter PE ***
    fclass: IF (ensemblefilter) THEN
       ! *** Ensemble filters SEIK/EnKF/LSEIK

       ALLOCATE(all_dim_ens_l(n_modeltasks))

       all_dim_ens_l = FLOOR(REAL(dim_ens) / REAL(n_modeltasks))
       DO i = 1, (dim_ens - n_modeltasks * all_dim_ens_l(1))
          all_dim_ens_l(i) = all_dim_ens_l(i) + 1
       END DO

       ! Initialize PE-local ensemble sizes
       IF (task_id>0) THEN
          dim_ens_l = all_dim_ens_l(task_id)
       ELSE
          dim_ens_l = dim_ens
       END IF

       IF (screen > 2 .AND. modelpe) &
            WRITE (*,*) 'PDAF: model task ', task_id, &
            ' mype(m)', mype_model, '; local ensemble size=', dim_ens_l

       ! *** Initialize array of displacements (for GATHER) ***
       ALLOCATE(all_dis_ens_l(n_modeltasks))

       all_dis_ens_l = 0
       DO i = 1, n_modeltasks - 1
          all_dis_ens_l(i + 1) = all_dis_ens_l(i) + all_dim_ens_l(i)
       END DO

       IF (screen > 2 .AND. filterpe) &
            WRITE (*,*) 'PDAF: PE(filter) ', mype_filter, &
            '; local displacements for ensemble=', all_dis_ens_l

    ELSE fclass
       ! *** Mode-based filter (SEEK) ***
       
       ALLOCATE(all_dim_ens_l(n_modeltasks))
       ALLOCATE(all_dim_eof_l(n_modeltasks))

       ! +1 required for SEEK to joint evolution of ens and state
       all_dim_eof_l = FLOOR( REAL(dim_ens) / REAL(n_modeltasks))
       DO i = 1, (dim_ens - n_modeltasks * all_dim_eof_l(1))
          all_dim_eof_l(i) = all_dim_eof_l(i) + 1
          statetask = i + 1
       END DO
       ! determine which task evolves the central state
       IF (statetask > n_modeltasks .OR. statetask == -1) statetask = 1

       ! Initialize PE-local numbers of EOFs
       IF (task_id>0) THEN 
          dim_eof_l = all_dim_eof_l(task_id)
       ELSE
          dim_eof_l = dim_ens
       END IF

       IF (screen > 2) &
            WRITE (*,*) 'PDAF: model task ', task_id, &
            ' mype(m)', mype_model, '; local number of modes=', dim_eof_l

       ! *** Initialize array of displacements (for GATHER) ***
       ALLOCATE(all_dis_ens_l(n_modeltasks))
       ALLOCATE(all_dis_eof_l(n_modeltasks))

       all_dis_eof_l = 0
       DO i = 1, n_modeltasks - 1
          all_dis_eof_l(i + 1) = all_dis_eof_l(i) + all_dim_eof_l(i)
       END DO

       IF (screen > 2 .AND. filterpe) &
            WRITE (*,*) 'PDAF: PE(filter) ', mype_filter, &
            '; local displacements for EOFs=', all_dis_eof_l

       ! Initialize _ens_ variables for unified use
       all_dim_ens_l = all_dim_eof_l
       dim_ens_l = all_dim_eof_l(task_id)
       all_dis_ens_l = all_dis_eof_l

    END IF fclass



! *********************
! *** Screen output ***
! *********************

    filter_pe: IF (filterpe .AND. mype == 0 .AND. screen > 0 .AND. flag == 0) THEN

       WRITE (*, '(/a)') 'PDAF: Initialize Parallelization'

       ! *** Parallelization information ***
       WRITE (*, '(a)') 'PDAF     Parallelization - Filter on model PEs:'
       WRITE (*, '(a, i6)') 'PDAF                 Total number of PEs: ', npes_world
       WRITE (*, '(a, i6)') 'PDAF      Number of parallel model tasks: ', n_modeltasks
       WRITE (*, '(a, i6)') 'PDAF                      PEs for Filter: ', npes_filter
       WRITE (*, '(a)') 'PDAF     # PEs per ensemble task and local ensemble sizes: '

       tskstr1 = ''
       DO task = 1, n_modeltasks
          WRITE (tskstrtmp, '(i6)') task
          tskstr1 = TRIM(tskstr1) // tskstrtmp
       ENDDO

       tskstr2 = ''
       DO task = 1, n_modeltasks
          IF (all_npes_model_l(task) < 10) THEN
             WRITE (tskstrtmp, '(5x, i1)') all_npes_model_l(task)
          ELSEIF (all_npes_model_l(task) < 100) THEN
             WRITE (tskstrtmp, '(4x, i2)') all_npes_model_l(task)
          ELSEIF (all_npes_model_l(task) < 1000) THEN
             WRITE (tskstrtmp, '(3x, i3)') all_npes_model_l(task)
          ELSEIF (all_npes_model_l(task) < 10000) THEN
             WRITE (tskstrtmp, '(2x, i4)') all_npes_model_l(task)
          ELSE
             WRITE (tskstrtmp, '(1x, i5)') all_npes_model_l(task)
          ENDIF
          tskstr2 = TRIM(tskstr2) // tskstrtmp
       ENDDO
       tskstr3 = ''
       DO task = 1, n_modeltasks
          IF (all_dim_ens_l(task) < 10) THEN
             WRITE (tskstrtmp, '(5x, i1)') all_dim_ens_l(task)
          ELSEIF (all_dim_ens_l(task) < 100) THEN
             WRITE (tskstrtmp, '(4x, i2)') all_dim_ens_l(task)
          ELSEIF (all_dim_ens_l(task) < 1000) THEN
             WRITE (tskstrtmp, '(3x, i3)') all_dim_ens_l(task)
          ELSEIF (all_dim_ens_l(task) < 10000) THEN
             WRITE (tskstrtmp, '(2x, i4)') all_dim_ens_l(task)
          ELSE
             WRITE (tskstrtmp, '(1x, i5)') all_dim_ens_l(task)
          ENDIF
          tskstr3 = TRIM(tskstr3) // tskstrtmp
       ENDDO
       WRITE (*, '(a13, a)') 'PDAF     Task', TRIM(tskstr1)
       WRITE (*, '(a13, a)') 'PDAF     #PEs', TRIM(tskstr2)
       WRITE (*, '(a13, a)') 'PDAF        N', TRIM(tskstr3)
       IF (.NOT. ensemblefilter) &
            WRITE (*, '(a, i3)') 'PDAF   Evolve central state on task: ', statetask

    END IF filter_pe

  END SUBROUTINE PDAF_init_parallel

END MODULE PDAF_mod_filterMPI
