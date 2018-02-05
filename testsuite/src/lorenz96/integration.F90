!$Id: integration.F90 1666 2016-11-28 17:52:28Z lnerger $
!BOP
!
! !ROUTINE: integration  --- integration routine for the Lorenz96 model
!
! !INTERFACE:
SUBROUTINE integration(time, nsteps)

! !DESCRIPTION:
! Routine to perform model integration with the Lorenz96 model.
!
! For simplicity of the implementation with PDAF,
! the time stepping is separated into a single routine.
! This allows to simply put the assimilation routine
! assimilation\_pdaf() in between the main program and
! the integration routine. If the time stepping is not
! available as a separate routine, a different 
! implementation style is required.
!
! !REVISION HISTORY:
! 2009-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE timer, &            ! Timing
       ONLY: timeit      
  USE mod_model, &        ! Model variables
       ONLY: x, dt, dim_state
#ifdef USE_PDAF
  USE mod_assimilation, & ! Variables for assimilation
       ONLY: filtertype, incremental, model_error
#endif
  USE output_netcdf, &    ! NetCDF output
       ONLY: write_netcdf, close_netcdf

  IMPLICIT NONE

! !ARGUMENTS:
  REAL, INTENT(inout) :: time   ! Model time
  INTEGER, INTENT(in) :: nsteps ! Number of time steps to be performed
!EOP

! local variables
  INTEGER :: step               ! Time step counter
  REAL, ALLOCATABLE :: x1(:), x2(:), x3(:), x4(:) ! Temporary arrays for RK4

#ifdef USE_PDAF
  EXTERNAL :: distribute_stateinc ! Routine to add state increment for IAU
#endif


! **********************
! *** Initialization ***
! **********************

  ! Allocate temporary arrays for RK4
  ALLOCATE(x1(dim_state))
  ALLOCATE(x2(dim_state))
  ALLOCATE(x3(dim_state))
  ALLOCATE(x4(dim_state))


! *********************************
! *** Perform model integration ***
! *********************************

  CALL timeit(5, 'new')

! *** time stepping loop ***
  integrate: DO step = 1, nsteps
     
#ifdef USE_PDAF
     ! For incremental updating (SEEK, SEIK, and LSEIK)
     IF (incremental == 1 &
          .AND. (filtertype==0 .OR. filtertype == 1 .OR. filtertype == 3)) THEN
        CALL PDAF_incremental(nsteps, distribute_stateinc)
     END IF
#endif

! *** model time step - RK4 ***

     ! Intermediate steps
     CALL lorenz96_dxdt(dim_state, x, x1)
     x1 = dt * x1
     CALL lorenz96_dxdt(dim_state, x + x1/2.0, x2)
     x2 = dt * x2
     CALL lorenz96_dxdt(dim_state, x + x2/2.0, x3)
     x3 = dt * x3
     CALL lorenz96_dxdt(dim_state, x + x3, x4)
     x4 = dt * x4

     ! New value of x
     x = x + x1/6.0 + x2/3.0 + x3/3.0 + x4/6.0

     ! Increment time
     time = time + dt

     
#ifdef USE_PDAF
     ! *** PDAF: Add model error ***      
     IF (model_error) CALL add_model_noise(dt, dim_state, x)
#endif

#ifndef USE_PDAF
     ! Write NetCDF output
     CALL write_netcdf(step, time, dim_state, x)
#endif

  END DO integrate

#ifndef USE_PDAF
  ! Close NetCDF file
  CALL close_netcdf()
#endif

  DEALLOCATE(x1, x2, x3, x4)

  CALL timeit(5, 'old')

END SUBROUTINE integration

! !ROUTINE: lorenz96_dxdt  --- compute dx/dt for Lorenz96 model
!
! !INTERFACE:
SUBROUTINE lorenz96_dxdt(dim_state, x, dxdt)

! !DESCRIPTION:
! This function computes the time derivate for the Lorenz96 model
!
! !REVISION HISTORY:
! 2009-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &        ! Model variables
       ONLY: forcing

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_state     ! Dimension of model state
  REAL, INTENT(in)  :: x(dim_state)    ! Model state
  REAL, INTENT(out) :: dxdt(dim_state) ! Time derivate
!EOP

! local variables
  INTEGER :: i, im1, im2, ip1

! Compute derivate
  DO i = 1, dim_state

     ! Cyclic boundary conditions:
     ip1 = i + 1
     IF(ip1 > dim_state) ip1 = 1
     im1 = i - 1
     IF(im1 < 1) im1 = dim_state
     im2 = i - 2
     IF(im2 < 1) im2 = dim_state + im2
   
     ! Compute derivative
     dxdt(i) = (x(ip1) - x(im2)) * x(im1) - x(i) + forcing

  END DO

END SUBROUTINE lorenz96_dxdt
