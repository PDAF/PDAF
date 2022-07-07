!$Id: initialize.F90 783 2009-12-07 10:28:43Z lnerger $
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
       ONLY: dims, dim_l, dims_l_all, dt, step_null, &
       step_final, field
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
  INTEGER :: i                 ! Counter
  INTEGER :: dim_decomp(2)     ! Number of processes in 2D domain decomposition
  INTEGER :: dim_decomp_tmp(2) ! Temporary array to compute dim_decomp

! *** Model specifications ***
  dims(1)  = 10      ! State dimension - first direction (x)
  dims(2)  = dims(1) ! State dimension - second direction (y)
  dims(3)  = 10      ! State dimension - third direction (z)
  dt          = 1.0  ! Size of time step
  step_null   = 0    ! initial time step
  total_steps = 2   ! number of time steps

! *** Parse command line options ***
  handle = 'dim1'               ! first dimension of dummy model
  CALL parse(handle, dims(1))
  dims(2) = dims(1)  ! Use quadratic domain, if not overwritten by command line

  handle = 'dim2'               ! second dimension of dummy model
  CALL parse(handle, dims(2))
  handle = 'dim3'               ! third dimension of dummy model
  CALL parse(handle, dims(3))

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
     WRITE (*, '(22x,a)') 'MODEL: 3D Dummy Model'
     WRITE (*, '(//10x,a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     WRITE (*, '(10x,a)') '!!! This model is not yet fully tested    !!!'
     WRITE (*, '(10x,a)') '!!! Implementation for PDAF is incomplete !!!'
     WRITE (*, '(10x,a//)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     WRITE (*, '(5x, a, 3i7)') &
          'Global model dimensions (x, y, z):', dims
     WRITE (*, '(17x, a, i7, a, i7, a1)') &
          'Time steps:', total_steps, '  (final step:', step_final, ')'
     IF (step_null /= 0) WRITE (*, '(10x, a, i7)') 'Initial time step:', step_null
     WRITE (*, '(13x, a, f10.3)') 'Time step size:', dt
  END IF screen2


! *** Initialize dimensions and fields with domain decompsition

  ! Determine dimensions of local domains
  ALLOCATE (dims_l_all(3, npes_model))
  ! count allocated memory
  CALL memcount(1, 'i', 3*npes_model)

  ! Determine domain decomposition
  dim_decomp(1) = FLOOR( SQRT( REAL(npes_model)))
  dim_decomp(2) = FLOOR( REAL(npes_model) / REAL(dim_decomp(1))) 
  IF (dim_decomp(1) * dim_decomp(2) /= npes_model) THEN
     DO i = 1, dim_decomp(1) - 1
        dim_decomp_tmp(1) = dim_decomp(1) -i
        dim_decomp_tmp(2) = FLOOR( REAL(npes_model) / REAL(dim_decomp_tmp(1)))
        IF (dim_decomp_tmp(1) * dim_decomp_tmp(2) == npes_model) THEN
           dim_decomp = dim_decomp_tmp
           EXIT
        END IF
     END DO
  END IF

  screendd: IF (mype_world == 0 .AND. npes_model > 1) THEN
     WRITE (*,'(/5x,a,i7,a)') &
          'Initialize domain decomposition for ', npes_model,' model processes.'
     IF (dim_decomp(1) <= dims(1) .AND. dim_decomp(2) <= dims(2)) THEN 
        WRITE (*,'(10x,a,2i5)') 'PEs per direction:', dim_decomp
     ELSE
        WRITE (*,'(/5x,a/)') 'ERROR: Number of processors incompatible with mesh dimensions!'
     END IF
  END IF screendd

  ! Initalize dimensions of local domains
  dims_l_all(1, :) = FLOOR( REAL(dims(1)) / REAL(dim_decomp(1)))
  DO i = 1, (dims(1) - dim_decomp(1) * dims_l_all(1, 1))
     dims_l_all(1, i) = dims_l_all(1, i) + 1
  END DO
  dims_l_all(2, :) = FLOOR( REAL(dims(2)) / REAL(dim_decomp(2)))
  DO i = 1, (dims(2) - dim_decomp(2) * dims_l_all(2, 1))
     dims_l_all(2, i) = dims_l_all(2, i) + 1
  END DO
  dims_l_all(3, :) = dims(3)
  IF (mype_world == 0) THEN
     WRITE (*, '(2x, a, i3, a)') &
          '-- Domain decomposition over each', npes_model, ' PEs'
     DO i = 1, npes_model
        WRITE (*, '(5x, a, i3, a, i3, a, 3i5)') &
             'task ', task_id, ' PE(model) ', i - 1, &
             ' dims_local(state): ', dims_l_all(1, i), dims_l_all(2, i), dims_l_all(3, i)
     END DO
  END IF

  ! Initialize array of PE-local dimensions
  dim_l(1) = dims_l_all(1, mype_model + 1)
  dim_l(2) = dims_l_all(2, mype_model + 1)
  dim_l(3) = dims_l_all(3, mype_model + 1)

  ! Allocate a model field for my PE-local domain
  ALLOCATE(field(dim_l(1), dim_l(2), dim_l(3)))
  ! count allocated memory
  CALL memcount(1, 'r', dim_l(1) * dim_l(2) * dim_l(3))
  
  ! Initialize model field
  field(:,:,:) = 1.0

  ! initialize model time 
  time = REAL(step_null)

END SUBROUTINE initialize
