!$Id: collect_state_pdaf.F90 2136 2019-11-22 18:56:35Z lnerger $
!BOP
!
! !ROUTINE: collect_state_pdaf --- Initialize state vector from model fields
!
! !INTERFACE:
SUBROUTINE collect_state_pdaf(dim_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/EnKF/SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! This subroutine is called during the forecast 
! phase from PDAF\_put\_state\_X after the 
! propagation of each ensemble member. 
! The supplied state vector has to be initialized
! from the model fields (typically via a module). 
! With parallelization, MPI communication might be 
! required to initialize state vectors for all 
! subdomains on the model PEs. 
!
! The routine is executed by each process that is
! participating in the model integrations.
!
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, ONLY: mype_world
  USE mod_assim_pdaf, ONLY: offset
  USE g_parfe, &
       ONLY: mydim_nod2d, mydim_nod3d, eDim_nod3D
  USE o_array, &
       ONLY: uf, ssh, tracer, Tsurf, Ssurf, w
  USE i_array, &
       ONLY: a_ice, m_ice, m_snow, u_ice, v_ice

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! local state vector

! !CALLING SEQUENCE:
! Called by: PDAF_put_state_X   (as U_coll_state)
!EOP

! Local variables
  INTEGER :: i         ! Counter
  

! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

  DO i = 1, myDim_nod2D
     state_p(i + offset(1)) = ssh(i)
  END DO

  DO i = 1, myDim_nod3D
     state_p(i + offset(2)) = uf(i)
  END DO

  DO i = 1, myDim_nod3D
     state_p(i + offset(3)) = uf(i + myDim_nod3d + eDim_nod3D)
  END DO

  DO i = 1, myDim_nod3D
     state_p(i + offset(4)) = w(i)
  END DO

  DO i = 1, myDim_nod3D
     state_p(i + offset(5)) = tracer(i, 1)
  END DO

  DO i = 1, myDim_nod3D
     state_p(i + offset(6)) = tracer(i, 2)
  END DO

  DO i = 1, myDim_nod2D
     state_p(i + offset(7)) = a_ice(i)
  END DO

  DO i = 1, myDim_nod2D
     state_p(i + offset(8)) = m_ice(i)
  END DO

  DO i = 1, myDim_nod2D
     state_p(i + offset(9)) = m_snow(i)
  END DO

  DO i = 1, myDim_nod2D
     state_p(i + offset(10)) = u_ice(i)
  END DO

  DO i = 1, myDim_nod2D
     state_p(i + offset(11)) = v_ice(i)
  END DO
  
END SUBROUTINE collect_state_pdaf
