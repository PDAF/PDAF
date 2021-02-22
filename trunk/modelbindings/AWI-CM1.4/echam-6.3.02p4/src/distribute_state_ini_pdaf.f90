!$Id$
!BOP
!
! !ROUTINE: distribute_state_pdaf --- Initialize model fields from state vector
!
! !INTERFACE:
SUBROUTINE distribute_state_ini_pdaf(dim_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/EnKF/SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! During the forecast phase of the filter this
! subroutine is called from PDAF\_get\_state
! supplying a model state which has to be evolved. 
! The routine has to initialize the fields of the 
! model (typically available through a module) from 
! the state vector of PDAF. With parallelization, 
! MPI communication might be required to 
! initialize all subdomains on the model PEs.
!
! The routine is executed by each process that is
! participating in the model integrations.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, ONLY: mype_submodel, task_id,mype_filter,mype_world
  USE mod_assim_pdaf,   ONLY: offset
  USE mod_assim_atm_pdaf, ONLY: dp
  USE mo_memory_g3b,    ONLY: aps
  USE mo_decomposition, ONLY: dc=>local_decomposition
  USE mo_scan_buffer,   ONLY: t
  USE mo_memory_g1a,    ONLY: tm1, alpsm1, vom1, dm1, qm1
  USE mo_memory_g2a,    ONLY: um1, vm1


  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL(dp), INTENT(inout) :: state_p(dim_p)  ! PE-local state vector

! !CALLING SEQUENCE:
! Called by: PDAF_get_state   (as U_dist_state)
!EOP

! Local variables
  INTEGER :: i, j, k, jk, jrow, jl          ! Counter
  INTEGER :: nproma, ngpblks, nlev
  REAL :: increment  

! *******************************************
! *** Initialize model fields from state  ***
! *** Each model PE knows his sub-state   ***
!********************************************

  if (mype_submodel==0) write (*,*) 'ECHAM-PDAF distribute_state_pdaf, task: ', task_id

  ngpblks=dc%ngpblks
  nlev=dc%nlev

  k = 1
  DO jk = nlev,1,-1
    DO jrow = 1, ngpblks

      IF ( jrow == ngpblks ) THEN
        nproma = dc%npromz
      ELSE
        nproma = dc%nproma
      END IF

      DO jl = 1, nproma
        tm1(jl,jk,jrow) = state_p(k+offset(1))
        k = k + 1
      END DO

    END DO
  END DO

  k = 1
  DO jrow = 1, ngpblks
    IF ( jrow == ngpblks ) THEN
       nproma = dc%npromz
    ELSE
       nproma = dc%nproma
    END IF

    DO jl = 1, nproma
      alpsm1(jl,jrow) = state_p(k+offset(2))
        k = k + 1
    END DO

  END DO

  k = 1
  DO jk = nlev,1,-1
    DO jrow = 1, ngpblks

      IF ( jrow == ngpblks ) THEN
        nproma = dc%npromz
      ELSE
        nproma = dc%nproma
      END IF

      DO jl = 1, nproma
        vom1(jl,jk,jrow) = state_p(k+offset(3))
        k = k + 1
      END DO

    END DO
  END DO

  k = 1
  DO jk = nlev,1,-1
    DO jrow = 1, ngpblks

      IF ( jrow == ngpblks ) THEN
        nproma = dc%npromz
      ELSE
        nproma = dc%nproma
      END IF

      DO jl = 1, nproma
        dm1(jl,jk,jrow) = state_p(k+offset(4))
        k = k + 1
      END DO

    END DO
  END DO

  k = 1
  DO jk = nlev,1,-1
    DO jrow = 1, ngpblks

      IF ( jrow == ngpblks ) THEN
        nproma = dc%npromz
      ELSE
        nproma = dc%nproma
      END IF

      DO jl = 1, nproma
        qm1(jl,jk,jrow) = state_p(k+offset(5))
        k = k + 1
      END DO

    END DO
  END DO

  k = 1
  DO jk = nlev,1,-1
    DO jrow = 1, ngpblks

      IF ( jrow == ngpblks ) THEN
        nproma = dc%npromz
      ELSE
        nproma = dc%nproma
      END IF

      DO jl = 1, nproma
        um1(jl,jk,jrow) = state_p(k+offset(6))
        k = k + 1
      END DO

    END DO
  END DO

  k = 1
  DO jk = nlev,1,-1
    DO jrow = 1, ngpblks

      IF ( jrow == ngpblks ) THEN
        nproma = dc%npromz
      ELSE
        nproma = dc%nproma
      END IF

      DO jl = 1, nproma
        vm1(jl,jk,jrow) = state_p(k+offset(7))
        k = k + 1
      END DO

    END DO
  END DO

END SUBROUTINE distribute_state_ini_pdaf
