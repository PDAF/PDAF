!$Id: initialize.F90 1590 2015-06-12 11:59:35Z lnerger $
!BOP
!
! !ROUTINE: initialize  --- initialize the 1D dummy model for PDAF offline
!
! !INTERFACE:
SUBROUTINE initialize()

! !DESCRIPTION:
! Routine to perform initialization of the 1D dummy model for the
! offline configuration of PDAF. Here, only the global size of the
! model state vector as well as sizes for the domain decomposition
! of the model domain need to be initialized.
! Generally, this could also be joined with the routine init_pdaf().
!
! For the 1D dummy model, the state vector size is simply defined by
! dim_state. The local dimensions are obtained by evenly splitting
! the global dimension.
!
! !REVISION HISTORY:
! 2009-05 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &        ! Model variables
       ONLY: dim_state, dim_state_p, local_dims
  USE mod_parallel, &     ! Parallelization variables
       ONLY: MPI_COMM_WORLD, mype_world, mype_model, npes_model, task_id, &
       init_parallel, finalize_parallel
  USE mod_memcount, &     ! Counting allocated memory
       ONLY: memcount, memcount_ini, memcount_get
  USE parser, &           ! Parse command lines
       ONLY: handle, parse

  IMPLICIT NONE

! !ARGUMENTS:
!EOP

! local variables
  INTEGER :: i   ! Counter

! *** Model specifications ***
  dim_state   = 300 ! State dimension (shared via MOD_OFFLINE)

! *** Parse command line options - optional ***
  handle = 'dim_state'               ! state dimension of dummy model
  CALL parse(handle, dim_state)

! *** Screen output ***
  screen2: IF (mype_world == 0) THEN
     WRITE (*, '(1x, a)') 'INITIALIZE MODEL INFORMATION FOR PDAF OFFLINE MODE'
     WRITE (*, '(22x,a)') 'MODEL: 1D Dummy Model'
     WRITE (*, '(5x, a, i7)') &
          'Global model state dimension:', dim_state
  END IF screen2

! *** Initialize dimensions and fields with domain decompsition

  ! Determine dimensions of local domains
  ALLOCATE (local_dims(npes_model))
  ! count allocated memory
  CALL memcount(1, 'i', npes_model)

  local_dims = FLOOR(REAL(dim_state) / REAL(npes_model))
  DO i = 1, (dim_state - npes_model * local_dims(1))
     local_dims(i) = local_dims(i) + 1
  END DO
  IF (mype_world == 0) THEN
     WRITE (*, '(/2x, a, i3, a)') &
          '-- Domain decomposition over each', npes_model, ' PEs'
     DO i = 1, npes_model
        WRITE (*, '(5x, a, i3, a, i3, a, i5)') &
             'task ', task_id, ' PE(model) ', i-1, &
             ' dim_local(state): ', local_dims(i)
     END DO
  END IF
  
  ! State dimension for my PE-local domain
  dim_state_p = local_dims(mype_model + 1)

END SUBROUTINE initialize
