!$Id: add_model_noise.F90 1666 2016-11-28 17:52:28Z lnerger $
!BOP
!
! !ROUTINE: add_model_noise  --- integration routine for the Lorenz96 model
!
! !INTERFACE:
SUBROUTINE add_model_noise(dt, dim_state, x)

! !DESCRIPTION:
! Routine to add model error to the model state of the
! Lorenz96 model. The variance of the model error is
! given by model_err_amp times the time step size.
!
! !REVISION HISTORY:
! 2016-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, & ! Variables for assimilation
       ONLY: model_err_amp

  IMPLICIT NONE

! !ARGUMENTS:
  REAL, INTENT(in) :: dt              ! Time step size
  INTEGER, INTENT(in) :: dim_state    ! State dimension
  REAL, INTENT(inout) :: x(dim_state) ! Model state
!EOP

! local variables
  REAL, ALLOCATABLE :: noise(:)      ! Random noise
  INTEGER, SAVE :: iseed(4)          ! Seed for random number generator
  INTEGER, SAVE :: firststep=1       ! Flag for first call


! ***********************
! *** Add model error ***
! ***********************

  IF (firststep==1) THEN
     WRITE (*,'(9x, a)') '--- Initialize seed for model error'
     iseed(1)=2*220+1
     iseed(2)=2*100+5
     iseed(3)=2*10+7
     iseed(4)=2*30+9
     firststep=0
  ENDIF

  ! Generate random Gaussian noise
  ALLOCATE(noise(dim_state))

  CALL dlarnv(3, iseed, dim_state, noise)

  ! Add noise to generate observations
  x = x + SQRT(model_err_amp) * SQRT(dt) * noise

  DEALLOCATE(noise)

END SUBROUTINE add_model_noise
