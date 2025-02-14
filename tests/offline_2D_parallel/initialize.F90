!>  Initialize model
!!
!! Routine to perform initialization of the 2D offline example for
!! PDAF. Implementation with parallelization.
!! Here, only the global size of the model domain, the global size
!! of the model state vector and the sizes for decomposition of the 
!! state vector need to be initialized.
!! Generally, this could also be joined with the routine init_pdaf().
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE initialize()

  USE mod_assimilation, &   ! Assimilation variables
       ONLY: nx, ny, dim_state, dim_state_p, local_dims
  USE mod_parallel_pdaf, &  ! Parallelization variables
       ONLY: mype_world, mype_model, npes_model, task_id
  USE parser, &             ! Parser function
       ONLY: parse

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: i                 ! Counter
  CHARACTER(len=32) :: handle  ! handle for command line parser
  INTEGER :: gridsize=1        ! Setting for grid dimensions: (1) 18x36, (2) 512x2048


! **********************
! *** INITIALIZATION ***
! **********************

  ! Parse grid size index
  handle='gridsize'
  CALL parse(handle, gridsize)

! *** Model specifications ***
  IF (gridsize==1) THEN
     nx = 36    ! Extent of grid in x-direction
     ny = 18    ! Extent of grid in y-direction
  ELSEIF (gridsize==2) THEN
     nx = 256    ! Extent of grid in x-direction
     ny = 128    ! Extent of grid in y-direction
  ELSE
     nx = 2048   ! Extent of grid in x-direction
     ny = 512    ! Extent of grid in y-direction
  END IF

  dim_state   = nx * ny ! State dimension (shared via MOD_OFFLINE)

! *** Screen output ***
  screen2: IF (mype_world == 0) THEN
     WRITE (*, '(1x, a)') 'INITIALIZE MODEL INFORMATION FOR PDAF OFFLINE MODE'
     WRITE (*, '(22x,a)') '2D Offline Test case for verification'
     WRITE (*, '(24x,a,i6,1x,a1,1x,i6)') 'Grid size:',nx,'x',ny
     WRITE (*, '(5x, a, i9)') &
          'Global model state dimension:', dim_state
  END IF screen2

! *** Initialize dimensions and fields with domain decompsition

  ! Determine dimensions of local domains
  ALLOCATE (local_dims(npes_model))

  local_dims = FLOOR(REAL(dim_state) / REAL(npes_model))
  DO i = 1, (dim_state - npes_model * local_dims(1))
     local_dims(i) = local_dims(i) + 1
  END DO
  IF (mype_world == 0) THEN
     WRITE (*, '(/2x, a, i3, a)') &
          '-- Domain decomposition over each', npes_model, ' PEs'
     DO i = 1, npes_model
        WRITE (*, '(5x, a, i3, a, i3, a, i9)') &
             'task ', task_id, ' PE(model) ', i-1, &
             ' dim_local(state): ', local_dims(i)
     END DO
  END IF
  
  ! State dimension for my PE-local domain
  dim_state_p = local_dims(mype_model + 1)

END SUBROUTINE initialize
