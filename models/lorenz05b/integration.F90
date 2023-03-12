!$Id: integration.F90 261 2019-11-28 11:36:49Z lnerger $
!BOP
!
! !ROUTINE: integration  --- integration routine for the Lorenz2005 model
!
! !INTERFACE:
SUBROUTINE integration(time, nsteps)

! !DESCRIPTION:
! Routine to perform model integration with the Lorenz2005 model.
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
  EXTERNAL :: distribute_stateinc_pdaf ! Routine to add state increment for IAU
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
        CALL PDAF_incremental(nsteps, distribute_stateinc_pdaf)
     END IF
#endif

! *** model time step - RK4 ***

     ! Intermediate steps
     CALL lorenz05b_dxdt(dim_state, x, x1)
     x1 = dt * x1
     CALL lorenz05b_dxdt(dim_state, x + x1/2.0, x2)
     x2 = dt * x2
     CALL lorenz05b_dxdt(dim_state, x + x2/2.0, x3)
     x3 = dt * x3
     CALL lorenz05b_dxdt(dim_state, x + x3, x4)
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

! !ROUTINE: lorenz05b_dxdt  --- compute dx/dt for Lorenz2005 model II
!
! !INTERFACE:
SUBROUTINE lorenz05b_dxdt(dim_state, x, dxdt)

! !DESCRIPTION:
! This function computes the time derivate for the Lorenz2005 model II.
!
! !REVISION HISTORY:
! 2020-11 - Nabir Mamnun - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &        ! Model variables
       ONLY: forcing, k_avg

  IMPLICIT NONE

  ! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_state        ! Dimension of model state
  REAL, INTENT(in)    :: x(dim_state)     ! Model state
  REAL, INTENT(out)   :: dxdt(dim_state)  ! Time derivate


  !EOP
  ! local variables
  INTEGER :: i, j, imj, hk                ! Indices
  INTEGER :: imkpj, ipkpj, imk, im2k      ! Indices
  INTEGER :: j_low, j_up                  ! Counters
  REAL    :: w(dim_state)                 ! Temporary array w
  REAL    :: xx(dim_state)                ! [x,x] of the model equation
  INTEGER :: odd_check


  ! *** Initialization ***

  odd_check = mod(k_avg, 2)

  ! Lower and upper loop limits for averaging

  IF (odd_check==0) THEN ! k_avg is even
     hk = k_avg/2
     j_low = -hk + 1
     j_up  = hk - 1
  ELSE ! k_avg is odd
    j_low = - (k_avg - 1)/2
    j_up  = (k_avg - 1)/2
  ENDIF

    ! *** Compute the vector W ***

  w = 0.0
  DO i = 1, dim_state
     DO j= j_low, j_up
        ! Cyclic boundary conditions:
        imj = i - j
        IF(imj > dim_state) imj = imj - dim_state
        IF(imj < 1) imj = imj + dim_state
        w(i) = w(i) + x(imj)
      END DO
   END DO

   IF (odd_check==0) THEN
      DO i = 1, dim_state
         ! First and last term of the sum for even k_avg
         imj = i + hk
         IF(imj > dim_state) imj = imj - dim_state
         w(i) = w(i) + x(imj)/2.0
         imj = i-hk
         IF(imj < 1) imj = imj + dim_state
         w(i) = w(i) + x(imj)/2.0
      END DO
   ENDIF

   w = w/REAL(k_avg)

    ! *** Compute the xx ***

   xx = 0.0
   DO i = 1, dim_state
      DO j= j_low, j_up
         imkpj = i - k_avg + j
         IF(imkpj < 1) imkpj = imkpj + dim_state
         ipkpj = i + k_avg + j
         IF(ipkpj > dim_state) ipkpj = ipkpj - dim_state
         xx(i) = xx(i) + w(imkpj)*x(ipkpj)
      END DO
   END DO

   IF (odd_check==0) THEN
      DO i = 1, dim_state
         ! First and last term of the sum for even k_avg
         imkpj = i - hk*3  ! k_avg + k_avg/2 = (k_avg/2)*3 = hk*3
         ipkpj = i + hk    !-k_avg + k_avg/2 = -(k_avg/2) = -hk
         IF(imkpj < 1) imkpj = imkpj + dim_state
         IF(ipkpj > dim_state) ipkpj = ipkpj - dim_state
         xx(i) = xx(i) + w(imkpj)*x(ipkpj)/2.0

         imkpj = i-hk      !-k_avg + k_avg/2 = -(k_avg/2) = -hk
         ipkpj = i + hk*3  ! k_avg + k_avg/2 = (k_avg/2)*3 = hk*3
         IF(imkpj < 1) imkpj = imkpj + dim_state
         IF(ipkpj > dim_state) ipkpj = ipkpj - dim_state
         xx(i) = xx(i) + w(imkpj)*x(ipkpj)/2.0
      END DO
   ENDIF

   xx = xx/REAL(k_avg)

   ! now add the both terms
   DO i = 1, dim_state
      im2k = i - 2*k_avg
      IF(im2k < 1) im2k = im2k + dim_state
      imk = i-k_avg
      IF(imk < 1) imk = imk + dim_state
      xx(i) = - w(im2k)*w(imk) + xx(i)
    END DO

    ! *** Compute derivate ***
    DO i = 1, dim_state
      dxdt(i) = xx(i) - x(i) + forcing
    END DO

END SUBROUTINE lorenz05b_dxdt
