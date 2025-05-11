!>  Initialize model
!!
!! Routine to perform initialization of model dimensions used to
!! set the dimension of the state vector for PDAF and for use
!! in call-back routines.
!! Here, the global size of the model domain and the 
!! dimension of the global state vector are initialized.
!! With MPI-parallelization one will usually initialize both the
!! global and process-local sized of the state vector.
!! Generally, this could also be joined with the routine init_pdaf().
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE initialize()

  USE mod_assimilation, &   ! Variables for state vector dimension and model grid
        ONLY: dim_state_p !, dim_state, local_dims
  USE mod_parallel_pdaf, &  ! Parallelization variables
       ONLY: mype_world, mype_model, npes_model, task_id

  IMPLICIT NONE

! *** Local variables ***


! **********************
! *** INITIALIZATION ***
! **********************

! *** Model specifications ***

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE initialize.F90: Initialize mesh dimensions here!'

! nx, ny = ?

  dim_state_p = 10    ! TEMPLATE: This is just a dummy value so that the program runs to the end


! *** Screen output ***
  screen2: IF (mype_world == 0) THEN
     WRITE (*, '(1x, a)') 'INITIALIZE MODEL INFORMATION FOR PDAF OFFLINE MODE'
     WRITE (*, '(22x,a)') 'MODEL: ...'
!      WRITE (*, '(24x,a,i4,1x,a1,1x,i4)') 'Grid size:',nx,'x',ny
     WRITE (*, '(5x, a, i7)') &
          'Global model state dimension:', dim_state_p
  END IF screen2

! *** Initialize dimensions and fields with domain decomposition

  ! Determine dimensions of local domains
!  ALLOCATE (local_dims(npes_model))

!   local_dims = ...

END SUBROUTINE initialize
