!$Id$
!>  Initialize model
!!
!! Routine to perform initialization of the 2D offline example for
!! PDAF. Here, only the global size of the model domain and the 
!! dimension of the global state vector are initialized.
!! Generally, this could also be joined with the routine init_pdaf().
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE initialize()

  USE mod_assimilation, &   ! Model variables
       ONLY: dim_state_p, nx, ny

  IMPLICIT NONE


! **********************
! *** INITIALIZATION ***
! **********************

! *** Model specifications ***
  nx = 36    ! Extent of grid in x-direction
  ny = 18    ! Extent of grid in y-direction

  ! *** Define state dimension ***
  dim_state_p = nx * ny

! *** Screen output ***
  WRITE (*, '(1x, a)') 'INITIALIZE MODEL INFORMATION FOR PDAF OFFLINE MODE'
  WRITE (*, '(22x,a)') 'MODEL: 2D Offline Example for Tutorial'
  WRITE (*, '(24x,a,i4,1x,a1,1x,i4)') 'Grid size:',nx,'x',ny
  WRITE (*, '(5x, a, i7)') &
       'Global model state dimension:', dim_state_p

END SUBROUTINE initialize
