!$Id: distribute_state_pdaf.F90 2196 2020-03-26 13:26:59Z lnerger $
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

  USE mod_parallel_pdaf, &
       ONLY: mype_submodel, mype_world, task_id
  USE mod_assim_pdaf, &
       ONLY: off_fields_p
  USE g_parfe, &
       ONLY: mydim_nod2d, mydim_nod3d, ToDim_nod2D, eDim_nod3D
  USE o_array, &
       ONLY: uf, ssh, tracer, Tsurf, Ssurf, w, ucori_back, vcori_back
  USE i_array, &
       ONLY: a_ice, m_ice, m_snow, u_ice, v_ice
  USE o_mesh, &
       ONLY: nod3D_below_nod2D

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< process-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  !< local state vector

! *** Local variables ***
  INTEGER :: i         ! Counter
  INTEGER :: node      ! Node index
  INTEGER :: cnt_aice, cnt_mice, cnt_msnow  ! Count number of invalid points for ice


! **********************
! *** Initialization ***
! **********************

  if (mype_submodel==0) write (*,*) 'FESOM-PDAF distribute_state_pdaf, task: ', task_id

  ! Initialize counters

  cnt_aice = 0
  cnt_mice = 0
  cnt_msnow = 0


! *******************************************
! *** Initialize model fields from state  ***
! *** Each model PE knows his sub-state   ***
!********************************************

  DO i = 1, myDim_nod2D
     ssh(i) = state_p(i + off_fields_p(1))
  END DO

  DO i = 1, myDim_nod3D
     uf(i) = state_p(i + off_fields_p(2))
  END DO

  DO i = 1, myDim_nod3D
     uf(i + myDim_nod3D + eDim_nod3D) = state_p(i + off_fields_p(3))
  END DO

  DO i = 1, myDim_nod3D
     w(i) = state_p(i + off_fields_p(4))
  END DO

  DO i = 1, myDim_nod3D
     tracer(i, 1) = state_p(i + off_fields_p(5))
  END DO

  DO i = 1, myDim_nod3D
     tracer(i, 2) = state_p(i + off_fields_p(6))
  END DO

  DO i = 1, myDim_nod2D
     a_ice(i) = state_p(i + off_fields_p(7))
     IF (a_ice(i) < 0.0) THEN
        a_ice(i) = 0.0
        cnt_aice = cnt_aice + 1
     END IF
  END DO

  DO i = 1, myDim_nod2D
     m_ice(i) = state_p(i + off_fields_p(8))
     IF (m_ice(i) < 0.0) THEN
        m_ice(i) = 0.0
        cnt_mice = cnt_mice + 1
     END IF
  END DO

  DO i = 1, myDim_nod2D
     m_snow(i) = state_p(i + off_fields_p(9))
     IF (m_snow(i) < 0.0) THEN
        m_snow(i) = 0.0
        cnt_msnow = cnt_msnow + 1
     END IF
  END DO

  DO i = 1, myDim_nod2D
     u_ice(i) = state_p(i + off_fields_p(10))
  END DO

  DO i = 1, myDim_nod2D
     v_ice(i) = state_p(i + off_fields_p(11))
  END DO


! *********************************
! *** Initialize external nodes ***
!**********************************

  CALL com_2d(ssh)
  CALL com_3d(uf(1 : myDim_nod3d+eDim_nod3D))
  CALL com_3d(uf(myDim_nod3d+eDim_nod3D+1 : 2*(eDim_nod3D+myDim_nod3d)))
  CALL com_3d(tracer(:, 1))
  CALL com_3d(tracer(:, 2))
  CALL com_3d(w)
  CALL com_2d(a_ice)
  CALL com_2d(m_ice)
  CALL com_2d(m_snow)
  CALL com_2d(u_ice)
  CALL com_2d(v_ice)
 
 DO i=1, ToDim_nod2D     
    node=nod3D_below_nod2D(1,i)       
    Tsurf(i)=tracer(node,1)          
    Ssurf(i)=tracer(node,2)          
 END DO
 
 ucori_back = 0.0
 vcori_back = 0.0

  IF (cnt_aice>0) WRITE (*,*) 'PE ',mype_world, &
       ' number of points with ice conc. below 0: ',cnt_aice
  IF (cnt_mice>0) WRITE (*,*) 'PE ',mype_world,&
       ' number of points with ice thickness below 0: ',cnt_mice
  IF (cnt_msnow>0) WRITE (*,*) 'PE ',mype_world,&
       ' number of points with show height below 0: ',cnt_msnow

END SUBROUTINE distribute_state_pdaf
