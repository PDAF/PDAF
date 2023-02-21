!$Id: initialize.F90 1366 2013-04-24 16:25:05Z lnerger $
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
  USE mod_parallel, &       ! Parallelization variables
       ONLY: mype_world, mype_model, npes_model, task_id

  IMPLICIT NONE

! *** Local variables ***
  INTEGER :: i   ! Counter


! **********************
! *** INITIALIZATION ***
! **********************

! *** Model specifications ***
  nx = 36    ! Extent of grid in x-direction
  ny = 18    ! Extent of grid in y-direction

  dim_state   = nx * ny ! State dimension (shared via MOD_OFFLINE)

! *** Screen output ***
  screen2: IF (mype_world == 0) THEN
     WRITE (*, '(1x, a)') 'INITIALIZE MODEL INFORMATION FOR PDAF OFFLINE MODE'
     WRITE (*, '(22x,a)') 'MODEL: 2D Offline Example'
     WRITE (*, '(24x,a,i4,1x,a1,1x,i4)') 'Grid size:',nx,'x',ny
     WRITE (*, '(5x, a, i7)') &
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
        WRITE (*, '(5x, a, i3, a, i3, a, i5)') &
             'task ', task_id, ' PE(model) ', i-1, &
             ' dim_local(state): ', local_dims(i)
     END DO
  END IF
  
  ! State dimension for my PE-local domain
  dim_state_p = local_dims(mype_model + 1)

END SUBROUTINE initialize
