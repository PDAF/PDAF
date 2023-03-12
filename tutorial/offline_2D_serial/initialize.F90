!$Id: initialize.F90 1589 2015-06-12 11:57:58Z lnerger $
!BOP
!
! !ROUTINE: initialize  --- initialize the 2D offline example for PDAF
!
! !INTERFACE:
SUBROUTINE initialize()

! !DESCRIPTION:
! Routine to perform initialization of the 2D offline example for
! PDAF. Here, only the global size of the model domain and the
! global size of the model state vector need to be initialized.
! Generally, this could also be joined with the routine init_pdaf().
!
! For the 2D offline tutorial example, the domain is defined by the
! dimensions nx and ny. The state vector size is nx*ny.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, & ! Model variables
       ONLY: dim_state_p, nx, ny
  USE mod_parallel, &     ! Parallelization variables
       ONLY: MPI_COMM_WORLD, init_parallel, finalize_parallel

  IMPLICIT NONE

!EOP

! *** Model specifications ***
  nx = 36    ! Extent of grid in x-direction
  ny = 18    ! Extent of grid in y-direction

  dim_state_p   = nx * ny ! State dimension (shared via MOD_OFFLINE)



! *** Screen output ***
  WRITE (*, '(1x, a)') 'INITIALIZE MODEL INFORMATION FOR PDAF OFFLINE MODE'
  WRITE (*, '(22x,a)') 'MODEL: 2D Offline Example for Tutorial'
  WRITE (*, '(24x,a,i4,1x,a1,1x,i4)') 'Grid size:',nx,'x',ny
  WRITE (*, '(5x, a, i7)') &
       'Global model state dimension:', dim_state_p

END SUBROUTINE initialize
