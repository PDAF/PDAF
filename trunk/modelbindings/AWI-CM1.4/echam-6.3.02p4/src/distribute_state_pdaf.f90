!$Id$
!>   Routine to initialize model fields from state vector
!!
!! User-supplied call-back routine for PDAF.
!!
!! During the forecast phase of the filter this
!! subroutine is called from PDAF_get_state
!! supplying a model state which has to be evolved. 
!! The routine has to initialize the fields of the 
!! model (typically available through a module) from 
!! the state vector of PDAF. With parallelization, 
!! MPI communication might be required to 
!! initialize all subdomains on the model PEs.
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
SUBROUTINE distribute_state_pdaf(dim_p, state_p)

  USE mod_parallel_pdaf, ONLY: mype_submodel, task_id, mype_filter, mype_world
  USE mod_assim_atm_pdaf, ONLY: dp, wp
  USE mod_assim_pdaf,    ONLY: offset
  USE mo_memory_g3b,     ONLY: aps
  USE mo_decomposition,  ONLY: dc=>local_decomposition
  USE mo_scan_buffer,    ONLY: t, alps, vo, d, u, v, alpha, alnpr
  USE mo_memory_g1a,     ONLY: tm1, alpsm1, vom1, dm1, qm1
  USE mo_memory_g2a,     ONLY: um1, vm1
  USE mo_memory_g3a,     ONLY: geospm
  USE mo_memory_gl,      ONLY: q
  USE mo_control,        ONLY: nvclev, nlev, vct


  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p               !< PE-local state dimension
  REAL(dp), INTENT(inout) :: state_p(dim_p)  !< PE-local state vector

! *** Local variables ***
  INTEGER :: i, k, jk, jrow, jl          ! Counters
  INTEGER :: nproma, ngpblks, nbdim
  REAL(dp) :: increment  
  REAL(wp) :: hyam(nlev), hybm(nlev)
  REAL(dp), DIMENSION(nvclev) :: hyai, hybi
  REAL(wp), ALLOCATABLE :: zp(:,:,:)
  REAL(dp), ALLOCATABLE :: geom1(:,:)
  LOGICAL :: loc_vert
  REAL, ALLOCATABLE :: fac_incre(:,:,:)  
  REAL(wp) :: pressure_low, pressure_high, pressure_diff


! *******************************************
! *** Initialize model fields from state  ***
! *** Each model PE knows his sub-state   ***
!********************************************

  if (mype_submodel==0) write (*,*) 'ECHAM-PDAF distribute_state_pdaf, task: ', task_id

  ngpblks=dc%ngpblks
  !nlev=dc%nlev

  ! Define two pressure levels used for vertical localization 
  ! Unit: Pa

  ! pressure_low: between surface and pressure_low full increment is added to
  ! the state vector
  ! pressure_high: between pressure_low and pressure_high increment by DA is
  ! linearly reduced added to the state vector
  ! Above pressure_high: no update
  ! If pressure_low is set equal to pressure_high: stepwise update until
  ! pressure_low
  
  pressure_low = 80000.0
  pressure_high = 60000.0
 
  pressure_diff = pressure_low - pressure_high

  ! Calculate the pressure
  hyai = vct(1:nvclev)
  hybi = vct(nvclev+1:2*nvclev)
  DO i=1,nlev
    hyam(i)=0.5_wp*(hyai(i)+hyai(i+1))
    hybm(i)=0.5_wp*(hybi(i)+hybi(i+1))
  END DO

  ALLOCATE(zp(dc%nproma,nlev,ngpblks))
  DO jk=1,nlev
    DO jrow = 1, ngpblks

      IF ( jrow == ngpblks ) THEN
        nproma = dc%npromz
      ELSE
        nproma = dc%nproma
      END IF

      DO jl = 1, nproma
        zp(jl,jk,jrow)=hyam(jk)+hybm(jk)*aps(jl,jrow)
      END DO
    END DO
  END DO

  
  ! Define an increment factor
  ALLOCATE(fac_incre(dc%nproma,nlev,ngpblks))  
  fac_incre=1.0

! Vertical localisation
  loc_vert=.true.

  IF (loc_vert) THEN
    if (mype_submodel==0) write (*,*) 'ECHAM-PDAF vertical localisation'
  !*****************************
  ! vertical localisation 
  !*****************************


  DO jk=1,nlev
    DO jrow = 1, ngpblks

      IF ( jrow == ngpblks ) THEN
        nproma = dc%npromz
      ELSE
        nproma = dc%nproma
      END IF

      DO jl = 1, nproma
        IF (zp(jl,jk,jrow) >= pressure_low) THEN
          fac_incre(jl,jk,jrow)=1.0
        ELSE IF (zp(jl,jk,jrow) < pressure_high) THEN
          fac_incre(jl,jk,jrow)=0.0
        ELSE IF (pressure_high /= pressure_low) THEN
        ! calculate fac_incre depending on the pressure
          fac_incre(jl,jk,jrow)=zp(jl,jk,jrow)/pressure_diff - pressure_high/pressure_diff
        END IF
      END DO
    END DO
  END DO


  END IF 

 
  nbdim = dc% nproma
  ALLOCATE(geom1(dc%nproma,nlev))
  DO jrow=1,ngpblks
    CALL geopot(geom1,tm1(1:nproma,:,jrow),alnpr(:,:,jrow),alpha(:,:,jrow),geospm,nbdim,nproma)
  END DO
  DEALLOCATE(geom1)
 

! 3D temperature
  k = 1
  DO jk = nlev,1,-1
    DO jrow = 1, ngpblks

      IF ( jrow == ngpblks ) THEN
        nproma = dc%npromz
      ELSE
        nproma = dc%nproma
      END IF

      DO jl = 1, nproma
        increment = state_p(k+offset(1)) - tm1(jl,jk,jrow)
        tm1(jl,jk,jrow) = tm1(jl,jk,jrow) + fac_incre(jl,jk,jrow) * increment
        t(jl,jk,jrow) = t(jl,jk,jrow) + fac_incre(jl,jk,jrow) * increment
        k = k + 1
      END DO

    END DO
  END DO


! 2D log surface pressure
  k = 1
  DO jrow = 1, ngpblks
    IF ( jrow == ngpblks ) THEN
       nproma = dc%npromz
    ELSE
       nproma = dc%nproma
    END IF

    DO jl = 1, nproma
      increment = state_p(k+offset(2)) - alpsm1(jl,jrow)
      alpsm1(jl,jrow) = alpsm1(jl,jrow) + fac_incre(jl,nlev,jrow) * increment
      alps(jl,jrow) = alps(jl,jrow) + fac_incre(jl,nlev,jrow) * increment
      k = k + 1
    END DO

  END DO

! 3D vorticity
  k = 1
  DO jk = nlev,1,-1
    DO jrow = 1, ngpblks

      IF ( jrow == ngpblks ) THEN
        nproma = dc%npromz
      ELSE
        nproma = dc%nproma
      END IF

      DO jl = 1, nproma
        increment = state_p(k+offset(3)) - vom1(jl,jk,jrow)
        vom1(jl,jk,jrow) = vom1(jl,jk,jrow) + fac_incre(jl,jk,jrow) * increment
        vo(jl,jk,jrow) = vo(jl,jk,jrow) + fac_incre(jl,jk,jrow) * increment
        k = k + 1
      END DO

    END DO
  END DO

! 3D divergence
  k = 1
  DO jk = nlev,1,-1
    DO jrow = 1, ngpblks

      IF ( jrow == ngpblks ) THEN
        nproma = dc%npromz
      ELSE
        nproma = dc%nproma
      END IF

      DO jl = 1, nproma
        increment = state_p(k+offset(4)) - dm1(jl,jk,jrow)
        dm1(jl,jk,jrow) = dm1(jl,jk,jrow) + fac_incre(jl,jk,jrow) * increment
        d(jl,jk,jrow) = d(jl,jk,jrow) + fac_incre(jl,jk,jrow) * increment
        k = k + 1
      END DO

    END DO
  END DO

! 3D specific humidity
  k = 1
  DO jk = nlev,1,-1
    DO jrow = 1, ngpblks

      IF ( jrow == ngpblks ) THEN
        nproma = dc%npromz
      ELSE
        nproma = dc%nproma
      END IF

      DO jl = 1, nproma
        increment = state_p(k+offset(5)) - qm1(jl,jk,jrow)
        qm1(jl,jk,jrow) = qm1(jl,jk,jrow) + fac_incre(jl,jk,jrow) * increment
        q(jl,jk,jrow) = q(jl,jk,jrow) + fac_incre(jl,jk,jrow) * increment
        k = k + 1
      END DO

    END DO
  END DO

! 3D u
  k = 1
  DO jk = nlev,1,-1
    DO jrow = 1, ngpblks

      IF ( jrow == ngpblks ) THEN
        nproma = dc%npromz
      ELSE
        nproma = dc%nproma
      END IF

      DO jl = 1, nproma
        increment = state_p(k+offset(6)) - um1(jl,jk,jrow)
        um1(jl,jk,jrow) = um1(jl,jk,jrow) + fac_incre(jl,jk,jrow) * increment
        u(jl,jk,jrow) = u(jl,jk,jrow) + fac_incre(jl,jk,jrow) * increment
        k = k + 1
      END DO

    END DO
  END DO

! 3D v
  k = 1
  DO jk = nlev,1,-1
    DO jrow = 1, ngpblks

      IF ( jrow == ngpblks ) THEN
        nproma = dc%npromz
      ELSE
        nproma = dc%nproma
      END IF

      DO jl = 1, nproma
        increment = state_p(k+offset(7)) - vm1(jl,jk,jrow)
        vm1(jl,jk,jrow) = vm1(jl,jk,jrow) + fac_incre(jl,jk,jrow) * increment
        v(jl,jk,jrow) = v(jl,jk,jrow) + fac_incre(jl,jk,jrow) * increment
        k = k + 1
      END DO

    END DO
  END DO

  DEALLOCATE(zp)
  DEALLOCATE(fac_incre)

END SUBROUTINE distribute_state_pdaf
