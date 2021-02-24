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
!! Since ECHAM uses a leap-frog time stepping we
!! here update both time levels using the same
!! increment.
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
SUBROUTINE distribute_state_pdaf(dim_p, state_p)

  USE mod_parallel_pdaf, ONLY: mype_submodel, task_id
  USE mod_assim_pdaf,    ONLY: off_fields_p
  USE mod_assim_atm_pdaf, ONLY: dp
  USE mo_memory_g3b,     ONLY: aps
  USE mo_decomposition,  ONLY: dc=>local_decomposition
  USE mo_scan_buffer,    ONLY: t, alps, vo, d, u, v, alpha, alnpr
  USE mo_memory_g1a,     ONLY: tm1, alpsm1, vom1, dm1, qm1
  USE mo_memory_g2a,     ONLY: um1, vm1
  USE mo_memory_gl,      ONLY: q
  USE mo_control,        ONLY: nlev


  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p               !< PE-local state dimension
  REAL(dp), INTENT(inout) :: state_p(dim_p)  !< PE-local state vector

! *** Local variables ***
  INTEGER :: i, k, jk, jrow, jl          ! Counters
  INTEGER :: nproma, ngpblks, nbdim
  REAL(dp) :: increment  


! *******************************************
! *** Initialize model fields from state  ***
! *** Each model PE knows his sub-state   ***
!********************************************

  if (mype_submodel==0) write (*,*) 'ECHAM-PDAF distribute_state_pdaf, task: ', task_id

  ngpblks=dc%ngpblks
 

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
           increment = state_p(k+off_fields_p(1)) - tm1(jl,jk,jrow)
           tm1(jl,jk,jrow) = tm1(jl,jk,jrow) + increment
           t(jl,jk,jrow) = t(jl,jk,jrow) + increment
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
        increment = state_p(k+off_fields_p(2)) - alpsm1(jl,jrow)
        alpsm1(jl,jrow) = alpsm1(jl,jrow) + increment
        alps(jl,jrow) = alps(jl,jrow) + increment
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
           increment = state_p(k+off_fields_p(3)) - vom1(jl,jk,jrow)
           vom1(jl,jk,jrow) = vom1(jl,jk,jrow) + increment
           vo(jl,jk,jrow) = vo(jl,jk,jrow) + increment
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
           increment = state_p(k+off_fields_p(4)) - dm1(jl,jk,jrow)
           dm1(jl,jk,jrow) = dm1(jl,jk,jrow) + increment
           d(jl,jk,jrow) = d(jl,jk,jrow) + increment
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
           increment = state_p(k+off_fields_p(5)) - qm1(jl,jk,jrow)
           qm1(jl,jk,jrow) = qm1(jl,jk,jrow) + increment
           q(jl,jk,jrow) = q(jl,jk,jrow) + increment
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
           increment = state_p(k+off_fields_p(6)) - um1(jl,jk,jrow)
           um1(jl,jk,jrow) = um1(jl,jk,jrow) + increment
           u(jl,jk,jrow) = u(jl,jk,jrow) + increment
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
           increment = state_p(k+off_fields_p(7)) - vm1(jl,jk,jrow)
           vm1(jl,jk,jrow) = vm1(jl,jk,jrow) + increment
           v(jl,jk,jrow) = v(jl,jk,jrow) + increment
           k = k + 1
        END DO

     END DO
  END DO

END SUBROUTINE distribute_state_pdaf
