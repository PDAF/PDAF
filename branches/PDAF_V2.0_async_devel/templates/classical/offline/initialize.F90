!$Id$
!BOP
!
! !ROUTINE: initialize  --- initialize the model for PDAF offline
!
! !INTERFACE:
SUBROUTINE initialize()

! !DESCRIPTION:
! Routine to perform initialization of the model infromation for
! PDAF. Here, the global size of the model domain, the global size
! of the model state vector and the sizes for decomposition of the 
! state vector need to be initialized.
! Generally, this could also be joined with the routine init_pdaf().
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &        ! Model variables
        ONLY: dim_state_p
  USE mod_parallel, &     ! Parallelization variables
       ONLY: mype_world, mype_model, npes_model, task_id

  IMPLICIT NONE

! !ARGUMENTS:
!EOP

! local variables


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



END SUBROUTINE initialize
