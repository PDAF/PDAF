! Copyright (c) 2004-2021 Lars Nerger
!
! This file is part of PDAF.
!
! PDAF is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version
! 3 of the License, or (at your option) any later version.
!
! PDAF is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with PDAF.  If not, see <http://www.gnu.org/licenses/>.
!
!$Id$
!BOP
!
! !ROUTINE: PDAF_en3dvar_optim_lbfgs --- Optimization with LBFGS for En3dVar
!
! !INTERFACE:
SUBROUTINE PDAF_en3dvar_optim_lbfgs(step, dim_p, dim_ens, dim_cvec_p, dim_obs_p, &
     ens_p, obs_p, dy_p, v_p, U_prodRinvA, &
     U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
     opt_parallel, screen)

! !DESCRIPTION:
! Optimization routine for ensemble 3D-Var using the LBFGS solver
!
! Variant for domain decomposed states.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2021-03 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: mype

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                  ! Current time step
  INTEGER, INTENT(in) :: dim_p                 ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens               ! ensemble size
  INTEGER, INTENT(in) :: dim_cvec_p            ! Size of control vector
  INTEGER, INTENT(in) :: dim_obs_p             ! PE-local dimension of observation vector
  REAL, INTENT(in) :: ens_p(dim_p, dim_ens)    ! PE-local state ensemble
  REAL, INTENT(in) :: obs_p(dim_obs_p)         ! Vector of observations
  REAL, INTENT(in) :: dy_p(dim_obs_p)          ! Background innovation
  REAL, INTENT(inout) :: v_p(dim_cvec_p)       ! Control vector
  INTEGER, INTENT(in) :: opt_parallel          ! Whether to use a decomposed control vector
  INTEGER, INTENT(in) :: screen                ! Verbosity flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA, &              ! Provide product R^-1 A
       U_cvt_ens, &                       ! Apply control vector transform matrix to control vector
       U_cvt_adj_ens, &                   ! Apply adjoint control vector transform matrix
       U_obs_op_lin, &                    ! Linearized observation operator
       U_obs_op_adj                       ! Adjoint observation operator

! !CALLING SEQUENCE:
! Called by: PDAF_3dvar_analysis_cvt
! Calls: PDAF_timeit
! Calls: PDAF_memcount
!EOP

! *** local variables ***
  INTEGER :: iter                      ! Counter
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  REAL :: J_tot                        ! Cost function
  REAL, ALLOCATABLE :: gradJ_p(:)      ! PE-local part of gradient of J

  ! Variables for LFBGS
  INTEGER, PARAMETER :: m = 5
  INTEGER            :: iprint
  CHARACTER(len=60)  :: task, csave
  LOGICAL            :: lsave(4)
  INTEGER            :: isave(44)
  REAL, PARAMETER    :: factr  = 1.0e+7, pgtol  = 1.0e-5
  REAL               :: dsave(29)
  INTEGER, ALLOCATABLE :: nbd(:), iwa(:)
  REAL, ALLOCATABLE  :: lvec(:), uvec(:), wa(:)



! **********************
! *** INITIALIZATION ***
! **********************

  ! Set verbosity of solver
  IF (screen>0 .AND. screen<2) THEN
     iprint = -1
  ELSEIF (screen<3) THEN
     iprint = 0
     IF (mype>0) iprint = -1
  ELSE
     iprint = 99
  END IF

  ! Allocate arrays
  ALLOCATE(nbd(dim_cvec_p), lvec(dim_cvec_p), uvec(dim_cvec_p))
  ALLOCATE (iwa(3*dim_cvec_p))
  ALLOCATE (wa(2*m*dim_cvec_p + 5*dim_cvec_p + 11*m*m + 8*m))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 11*dim_cvec_p + 2*m*dim_cvec_p + 11*m*m + 8*m)

  ! Settings for LBGFS
  nbd = 0  ! Values are unbounded
  task = 'START'
  iter = 1

  ! Initialize control vector
  v_p = 0.0
  

! ***************************
! ***   Iterative solving ***
! ***************************

  ! Prepare arrays for iterations
  ALLOCATE(gradJ_p(dim_cvec_p))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_cvec_p)

  IF (mype==0 .AND. screen > 0) &
       WRITE (*, '(a, 5x, a)') 'PDAF', '--- OPTIMIZE' 

  minloop: DO

     IF (.NOT.(task(1:2)=='FG' .OR. TRIM(task)=='NEW_X'.OR. &
          TRIM(task)=='START') ) THEN
        IF (mype==0 .AND. screen > 0) &
             WRITE (*,'(a, 5x, a, a)') 'PDAF', '--- Exit optimization, status ', task
        EXIT minloop
     END IF

     ! LBFGS
     CALL PDAF_timeit(54, 'new')
     CALL setulb(dim_cvec_p, m, v_p, lvec, uvec, nbd, &
          J_tot, gradJ_p, factr, pgtol, &
          wa, iwa, task, iprint,&
          csave, lsave, isave, dsave )
     CALL PDAF_timeit(54, 'old')


! ********************************
! ***   Evaluate cost function ***
! ********************************

     CALL PDAF_timeit(53, 'new')
     CALL PDAF_en3dvar_costf_cvt(step, iter, dim_p, dim_ens, dim_cvec_p, dim_obs_p, &
          ens_p, obs_p, dy_p, v_p, J_tot, gradJ_p, &
          U_prodRinvA, U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
          opt_parallel)
     CALL PDAF_timeit(53, 'old')

     IF (mype==0 .AND. screen >2) &
          WRITE (*,'(a, 8x, a, i5, es12.4)') 'PDAF', '--- iter, J: ', iter, J_tot

     iter = iter + 1

  END DO minloop


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(gradJ_p)
  DEALLOCATE(nbd, lvec, uvec, iwa, wa)

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_en3dvar_optim_lbfgs
