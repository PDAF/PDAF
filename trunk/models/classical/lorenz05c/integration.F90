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
     CALL lorenz05c_dzdt(dim_state, x, x1)
     x1 = dt * x1
     CALL lorenz05c_dzdt(dim_state, x + x1/2.0, x2)
     x2 = dt * x2
     CALL lorenz05c_dzdt(dim_state, x + x2/2.0, x3)
     x3 = dt * x3
     CALL lorenz05c_dzdt(dim_state, x + x3, x4)
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

! !ROUTINE: lorenz05c_dzdt  --- compute dx/dt for Lorenz2005 model III
!
! !INTERFACE:
SUBROUTINE lorenz05c_dzdt(dim_state, z, dzdt)

! !DESCRIPTION:
! This function computes the time derivate for the Lorenz2005 model III.
!
! !REVISION HISTORY:
! 2020-11 - Nabir Mamnun - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &        ! Model variables
       ONLY: forcing, k_avg, fluctuate_coef, coupling_coef, smoothing_scale

  IMPLICIT NONE

  ! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_state        ! Dimension of model state
  REAL, INTENT(in)    :: z(dim_state)     ! Model State
  REAL, INTENT(out)   :: dzdt(dim_state)  ! Time derivate


  !EOP
  ! local variables
  INTEGER   :: i, j, imj, ipj, hk            ! Indices
  INTEGER  :: im2, im1, ip1                  ! Indices
  INTEGER   :: imkpj, ipkpj, imk, im2k       ! Indices
  INTEGER   :: j_low, j_up                   ! Counters
  REAL      :: alpha, beta
  REAL      :: amb                           ! alpha - beta
  REAL      :: x(dim_state)                  ! slow variable
  REAL      :: y(dim_state)                  ! fast variable
  REAL      :: w(dim_state)                  ! Temporary array w
  REAL      :: xx(dim_state)                 ! [x,x] of the model equation
  REAL      :: yy(dim_state)                 ! [Y,Y] of the model equation
  REAL      :: yx(dim_state)                 ! [Y,X] of the model equation
  INTEGER   :: k_odd_check


  ! *** Initialization ***

  alpha = (3.0*(smoothing_scale**2.0) + 3.0) &
      / (2.0*(smoothing_scale**3.0) + 4.0*smoothing_scale)
  beta  = (2.0*(smoothing_scale**2.0) + 1.0) &
      / (1.0*(smoothing_scale**4.0) + 2.0*(smoothing_scale**2.0))


   ! *** construction of x and y through z ***

   x = 0.0
   DO i = 1, dim_state
      amb = alpha - beta*smoothing_scale    ! abs(-smoothing scale)
      ipj = i - smoothing_scale             ! i + (- smoothing_scale)
      IF(ipj < 1) ipj = dim_state + ipj
      x(i) = x(i) + amb*z(ipj)/2.0
      DO j = - (smoothing_scale - 1), (smoothing_scale - 1)
         amb = alpha - beta*abs(j)
         ! Cyclic boundary conditions
         ipj = i + j
         IF(ipj < 1) ipj = dim_state + ipj
         IF(ipj > dim_state) ipj = ipj - dim_state
         x(i) = x(i) + amb*z(ipj)
      END DO
      amb = alpha - beta*smoothing_scale    ! abs(smoothing_scale)
      ipj = i + smoothing_scale
      IF(ipj > dim_state) ipj = ipj - dim_state
      x(i) = x(i) + amb*z(ipj)/2.0
   END DO

   DO i = 1, dim_state
      y(i) = z(i) - x(i)
   END DO


  ! Lower and upper loop limits for averaging

  k_odd_check = mod(k_avg, 2)

  IF (k_odd_check==0) THEN ! k_avg is even
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

   IF (k_odd_check==0) THEN
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

   IF (k_odd_check==0) THEN
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


   ! *** Compute the yy ***
   ! [Y,Y]_1,n = -Y_n-2 Y_n-1 + Y_n-1 Y_n+1
   DO i = 1, dim_state
      ! Cyclic boundary conditions:
      ip1 = i + 1
      IF(ip1 > dim_state) ip1 = 1
      im1 = i - 1
      IF(im1 < 1) im1 = dim_state
      im2 = i - 2
      IF(im2 < 1) im2 = dim_state + im2
      yy(i) = - y(im2)*y(im1) + y(im1)*y(ip1)
   END DO


   ! *** Compute the xy ***
   ! [Y,X]_1,n = -Y_n-2 X_n-1 + Y_n-1 X_n+1
   DO i = 1, dim_state
      ! Cyclic boundary conditions:
      ip1 = i + 1
      IF(ip1 > dim_state) ip1 = 1
      im1 = i - 1
      IF(im1 < 1) im1 = dim_state
      im2 = i - 2
      IF(im2 < 1) im2 = dim_state + im2
      yx(i) = - y(im2)*x(im1) + y(im1)*x(ip1)
   END DO


  ! *** Compute derivate ***
  DO i = 1, dim_state
     dzdt(i) = xx(i) + (fluctuate_coef**2)*yy(i) + coupling_coef*yx(i) - x(i) &
     - fluctuate_coef*y(i) + forcing
  END DO

END SUBROUTINE lorenz05c_dzdt
