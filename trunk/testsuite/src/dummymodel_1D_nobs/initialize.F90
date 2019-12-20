!$Id: initialize.F90 1548 2015-01-17 09:22:02Z lnerger $
!BOP
!
! !ROUTINE: initialize  --- initialize the 1D dummy model
!
! !INTERFACE:
SUBROUTINE initialize()

! !DESCRIPTION:
! Routine to perform initialization of the 1D dummy model.
!
! !REVISION HISTORY:
! 2009-05 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &        ! Model variables
       ONLY: dim_state_p, local_dims, dt, step_null, &
       step_final, field
  USE mod_modeltime, &    ! Time information for model integration
       ONLY: time, total_steps
  USE mod_parallel, &     ! Parallelization variables
       ONLY: mype_world, npes_model
  USE mod_memcount, &     ! Counting allocated memory
       ONLY: memcount, memcount_ini, memcount_get
  USE parser, &           ! Parse command lines
       ONLY: handle, parse

  IMPLICIT NONE

! !ARGUMENTS:
!EOP

! local variables
  REAL :: rdim_state
  
! *** Model specifications ***
  dim_state_p = 300 ! State dimension per process
  dt          = 1.0 ! Size of time step
  step_null   = 0   ! initial time step
  total_steps = 20  ! number of time steps

! *** Parse command line options - optional ***
  handle = 'dim_state_p'               ! state dimension of dummy model
  CALL parse(handle, dim_state_p)
  handle = 'dt'                      ! Time step size
  CALL parse(handle, dt)
  handle = 'step_null'               ! Initial time step
  CALL parse(handle, step_null)
  handle = 'total_steps'             ! total number of time steps
  CALL parse(handle, total_steps)

! *** Define final time step ***  
  step_final = step_null + total_steps

! *** Compute global state dimension (double precision integer)
  rdim_state = REAL(dim_state_p) * REAL(npes_model)

! *** Screen output ***
  screen2: IF (mype_world == 0) THEN
     WRITE (*, '(1x, a)') 'INITIALIZE MODEL'
     WRITE (*, '(22x,a)') 'MODEL: 1D Dummy Model'
     WRITE (*, '(5x, a, es14.5)') &
          'Global model dimension:', rdim_state
     WRITE (*, '(17x, a, i7, a, i7, a1)') &
          'Time steps:', total_steps, '  (final step:', step_final, ')'
     IF (step_null /= 0) WRITE (*, '(10x, a, i7)') 'Initial time step:', step_null
     WRITE (*, '(13x, a, f10.3)') 'Time step size:', dt
  END IF screen2


! *** Initialize parallelization of the model by domain decomposition ***

  ! Determine dimensions of local domains
  ! - split domain into equally sized sub-domains
  ALLOCATE (local_dims(npes_model))
  ! count allocated memory
  CALL memcount(1, 'i', npes_model)

!   local_dims = FLOOR(REAL(dim_state) / REAL(npes_model))
!   DO i = 1, (dim_state - npes_model * local_dims(1))
!      local_dims(i) = local_dims(i) + 1
!   END DO
  local_dims(:) = dim_state_p
  IF (mype_world == 0) THEN
     WRITE (*, '(/2x, a, i6, a)') &
          '-- Domain decomposition over ', npes_model, ' PEs'
     WRITE (*, '(5x, a, i9)') &
          'state dimension per process: ', dim_state_p
  END IF
  
!   ! Set state dimension for my PE-local domain
!   dim_state_p = local_dims(mype_model + 1)

  ! Allocate the model field for my PE-local domain
  ALLOCATE(field(dim_state_p))
  ! count allocated memory
  CALL memcount(1, 'r', dim_state_p)
  
  ! Initialize model field
  field(:) = 1.0

  ! initialize model time 
  time = REAL(step_null)

END SUBROUTINE initialize
