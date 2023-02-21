!$Id: collect_state_pdaf.F90 1950 2018-11-12 15:02:24Z qtang $
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

  USE mod_parallel_pdaf, ONLY: mype_world
  USE mod_assim_pdaf, ONLY: off_fields_p
  USE g_parfe, &
       ONLY: mydim_nod2d, mydim_nod3d, eDim_nod3D
  USE o_array, &
       ONLY: uf, ssh, tracer, Tsurf, Ssurf, w
  USE i_array, &
       ONLY: a_ice, m_ice, m_snow, u_ice, v_ice

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< process-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  !< local state vector

! *** Local variables ***
  INTEGER :: i         ! Counter
  

! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

   DO i = 1, myDim_nod2D
      state_p(i + off_fields_p(1)) = ssh(i)
   END DO

   DO i = 1, myDim_nod3D
      state_p(i + off_fields_p(2)) = uf(i)
   END DO

   DO i = 1, myDim_nod3D
      state_p(i + off_fields_p(3)) = uf(i + myDim_nod3d + eDim_nod3D)
   END DO

   DO i = 1, myDim_nod3D
      state_p(i + off_fields_p(4)) = w(i)
   END DO

   DO i = 1, myDim_nod3D
      state_p(i + off_fields_p(5)) = tracer(i, 1)
   END DO
   
   DO i = 1, myDim_nod3D
      state_p(i + off_fields_p(6)) = tracer(i, 2)
   END DO

   DO i = 1, myDim_nod2D
      state_p(i + off_fields_p(7)) = a_ice(i)
   END DO

   DO i = 1, myDim_nod2D
      state_p(i + off_fields_p(8)) = m_ice(i)
   END DO

   DO i = 1, myDim_nod2D
      state_p(i + off_fields_p(9)) = m_snow(i)
   END DO

   DO i = 1, myDim_nod2D
      state_p(i + off_fields_p(10)) = u_ice(i)
   END DO

   DO i = 1, myDim_nod2D
      state_p(i + off_fields_p(11)) = v_ice(i)
   END DO
  
END SUBROUTINE collect_state_pdaf
