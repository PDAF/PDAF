!$Id$
!>  Routine to initialize state vector from model fields
!!
!! User-supplied call-back routine for PDAF.
!!
!! This subroutine is called during the forecast 
!! phase from PDAF_put_state_X or PDAF_assimilate_X
!! after the  propagation of each ensemble member. 
!! The supplied state vector has to be initialized
!! from the model fields (typically via a module). 
!! With parallelization, MPI communication might be 
!! required to initialize state vectors for all 
!! subdomains on the model PEs. 
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
SUBROUTINE collect_state_pdaf(dim_p, state_p)

  USE mod_assim_pdaf,     ONLY: off_fields_p
  USE mod_assim_atm_pdaf, ONLY: dp
  USE mo_decomposition,   ONLY: dc=>local_decomposition 
  USE mo_scan_buffer,     ONLY: t
  USE mo_memory_g1a,      ONLY: tm1, alpsm1, vom1, dm1, qm1
  USE mo_memory_g2a,      ONLY: um1, vm1

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p               !< PE-local state dimension
  REAL(dp), INTENT(inout) :: state_p(dim_p)  !< local state vector
  
! *** Local variables ***
  INTEGER :: k                        ! Counter
  INTEGER :: jrow, jk, jl             ! Counters
  INTEGER :: nproma, ngpblks, nlev    ! Indices


! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

  ngpblks = dc%ngpblks
  nlev    = dc%nlev

! 3D temperature
  k = 1
  DO jk = nlev, 1, -1
     DO jrow = 1, ngpblks

        IF ( jrow == ngpblks ) THEN
           nproma = dc%npromz
        ELSE
           nproma = dc%nproma
        END IF

        DO jl = 1, nproma
           state_p(k+off_fields_p(1)) = tm1(jl,jk,jrow)
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
        state_p(k+off_fields_p(2)) = alpsm1(jl,jrow)
        k = k + 1
     END DO

  END DO

! 3D vorticity
  k = 1
  DO jk = nlev, 1, -1
     DO jrow = 1, ngpblks

        IF ( jrow == ngpblks ) THEN
           nproma = dc%npromz
        ELSE
           nproma = dc%nproma
        END IF

        DO jl = 1, nproma
           state_p(k+off_fields_p(3)) = vom1(jl,jk,jrow)
           k = k + 1
        END DO

     END DO
  END DO

! 3D divergence
  k = 1
  DO jk = nlev, 1, -1
     DO jrow = 1, ngpblks

        IF ( jrow == ngpblks ) THEN
           nproma = dc%npromz
        ELSE
           nproma = dc%nproma
        END IF

        DO jl = 1, nproma
           state_p(k+off_fields_p(4)) = dm1(jl,jk,jrow)
           k = k + 1
        END DO

     END DO
  END DO

! 3D specific humidity
  k = 1
  DO jk = nlev, 1, -1
     DO jrow = 1, ngpblks

        IF ( jrow == ngpblks ) THEN
           nproma = dc%npromz
        ELSE
           nproma = dc%nproma
        END IF

        DO jl = 1, nproma
           state_p(k+off_fields_p(5)) = qm1(jl,jk,jrow)
           k = k + 1
        END DO

     END DO
  END DO

! 3D u
  k = 1
  DO jk = nlev, 1, -1
     DO jrow = 1, ngpblks

        IF ( jrow == ngpblks ) THEN
           nproma = dc%npromz
        ELSE
           nproma = dc%nproma
        END IF

        DO jl = 1, nproma
           state_p(k+off_fields_p(6)) = um1(jl,jk,jrow)
           k = k + 1
        END DO

     END DO
  END DO

! 3D v
  k = 1
  DO jk = nlev, 1, -1
     DO jrow = 1, ngpblks

        IF ( jrow == ngpblks ) THEN
           nproma = dc%npromz
        ELSE
           nproma = dc%nproma
        END IF

        DO jl = 1, nproma
           state_p(k+off_fields_p(7)) = vm1(jl,jk,jrow)
           k = k + 1
        END DO

     END DO
  END DO

END SUBROUTINE collect_state_pdaf
