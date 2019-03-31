!$Id: initialize.F90 784 2009-12-07 10:30:03Z lnerger $
!BOP
!
! !ROUTINE: initialize  --- initialize the Lorenz 63 model
!
! !INTERFACE:
SUBROUTINE initialize()

! !DESCRIPTION:
! Routine to perform initialization of the Lorenz 63 model.
! The model is introduced in E. N. Lorenz, Deterministic non-periodic
! flows. J. Atmos. Sci. 20 (1963) 130-141. Here, the model is 
! implemented using a 4th order Runge-Kutta scheme.
! Because the model has only 3 variables, it is implemented without 
! parallelization.
!
! !REVISION HISTORY:
! 2009-06 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &        ! Model variables
       ONLY: dt, step_null, step_final, x, gamma, &
       rho, beta
  USE mod_modeltime, &    ! Time information for model integration
       ONLY: time, total_steps
  USE mod_parallel, &     ! Parallelization variables
       ONLY: MPI_COMM_WORLD, MPIerr, npes_world, mype_world, &
       n_modeltasks, mype_model, npes_model, task_id, &
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
  REAL :: x0, y0, z0 ! Initial state
  
! *** Model specifications ***
  x0          = 1.508870  ! Initial state
  y0          = -1.531271 ! Initial state
  z0          = 25.46091  ! Initial state
  gamma       = 10.0      ! Model parameter gamma
  rho         = 28.0      ! Model parameter r
  beta        = 8.0/3.0   ! Model parameter beta
  dt          = 0.005     ! Size of time step
  step_null   = 0         ! initial time step
  total_steps = 20        ! number of time steps


! *** Parse command line options ***
  handle = 'x0'                      ! Initial state - first variable
  CALL parse(handle, x0)
  handle = 'y0'                      ! Initial state - 2nd variable
  CALL parse(handle, y0)
  handle = 'z0'                      ! Initial state - 3rd variable
  CALL parse(handle, z0)
  handle = 'gamma'                   ! Model parameter gamma
  CALL parse(handle, gamma)
  handle = 'rho'                     ! Model parameter r
  CALL parse(handle, rho)
  handle = 'beta'                    ! Model parameter beta
  CALL parse(handle, beta)
  handle = 'dt'                      ! Time step size
  CALL parse(handle, dt)
  handle = 'step_null'               ! Initial time step
  CALL parse(handle, step_null)
  handle = 'total_steps'             ! total number of time steps
  CALL parse(handle, total_steps)

! *** Define final time step ***  
  step_final = step_null + total_steps

! *** Screen output ***
  screen2: IF (mype_world == 0) THEN
     WRITE (*, '(1x, a)') 'INITIALIZE MODEL'
     WRITE (*, '(22x,a)') 'MODEL: Lorenz 63'
     WRITE (*, '(3x,a)') 'Initial state'
     WRITE (*, '(25x,a3,f12.8,/25x,a3,f12.8,/25x,a3,f12.8)') &
          'x0:', x0, 'y0:', y0, 'z0:', z0
     WRITE (*, '(3x,a)') 'Model parameters:'
     WRITE (*, '(22x,a6,f12.8,/26x,a2,f12.8,/23x,a5,f12.8)') &
          'Gamma:',gamma,'rho:',rho,'beta:',beta
     WRITE (*, '(17x, a, i7, a, i7, a1)') &
          'Time steps:', total_steps, '  (final step:', step_final, ')'
     IF (step_null /= 0) WRITE (*, '(10x, a, i7)') 'Initial time step:', step_null
     WRITE (*, '(13x, a, f10.3)') 'Time step size:', dt
  END IF screen2


! *** Initialize dimensions and fields with domain decompsition

  ! Allocate a model field
  ALLOCATE(x(3))
  ! count allocated memory
  CALL memcount(1, 'r', 3)
  
  ! Initialize model field
  x(1) = x0
  x(2) = y0
  x(3) = z0

  ! initialize model time 
  time = REAL(step_null)

END SUBROUTINE initialize
