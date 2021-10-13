!$Id$
!BOP
!
! !ROUTINE: integration  --- integration routine for the Lorenz 63 model
!
! !INTERFACE:
SUBROUTINE integration(time, nsteps)

! !DESCRIPTION:
! Routine to perform model integration with the Lorenz63 model.
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
! 2009-06 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE timer, &            ! Timing
       ONLY: timeit      
  USE mod_parallel, &     ! Parallelization
       ONLY: mype_model  
  USE mod_model, &        ! Model variables
       ONLY: x, dt
#ifdef USE_PDAF
  USE mod_assimilation, & ! Variables for assimilation
       ONLY: filtertype, incremental
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
  REAL :: x1(3), x2(3), x3(3), x4(3) ! Temporary arrays for RK4

#ifdef USE_PDAF
  EXTERNAL :: distribute_stateinc_pdaf ! Routine to add state increment for IAU
#endif


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
        CALL PDAF_incremental(nsteps, distribute_stateinc_pdaf)
     END IF

  ! *** PDAF: Add model error ***      
  ! This is the place where model error can be simulated.
  ! This routine is not implemented here!
!   call add_model_noise(step)
#endif

! *** model time step - RK4 ***

     ! Intermediate steps
     call lorenz_dxdt(x, x1)
     x1 = dt * x1
     call lorenz_dxdt(x + x1/2.0, x2)
     x2 = dt * x2
     call lorenz_dxdt(x + x2/2.0, x3)
     x3 = dt * x3
     call lorenz_dxdt(x + x3, x4)
     x4 = dt * x4

     ! New value of x
     x = x + x1/6.0 + x2/3.0 + x3/3.0 + x4/6.0

     ! Increment time
     time = time + dt

#ifndef USE_PDAF
     ! Write NetCDF output
     CALL write_netcdf(step, time, 3, x)
#endif

  END DO integrate

#ifndef USE_PDAF
  ! Close NetCDF file
  CALL close_netcdf()
#endif

  CALL timeit(5, 'old')

END SUBROUTINE integration

! !ROUTINE: lorenz_dxdt  --- compute dx/dt for Lorenz 63 model
!
! !INTERFACE:
SUBROUTINE lorenz_dxdt(x, dxdt)

! !DESCRIPTION:
! This function computes the time derivate for the Lorenz63 model
!
! !REVISION HISTORY:
! 2009-06 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &        ! Model variables
       ONLY: gamma, rho, beta

  IMPLICIT NONE

! !ARGUMENTS:
  REAL, INTENT(in)  :: x(3)    ! Model state
  REAL, INTENT(out) :: dxdt(3) ! Time derivate
!EOP

! local variables

! Compute derivate
  dxdt(1) = gamma * (x(2)-x(1))
  dxdt(2) = rho*x(1) - x(2) - x(1)*x(3)
  dxdt(3) = x(1)*x(2) - beta*x(3)

END SUBROUTINE lorenz_dxdt
