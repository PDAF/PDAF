! Copyright (c) 2004-2020 Lars Nerger
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
! !ROUTINE: PDAF_3dvar_analysis --- 3DVAR
!
! !INTERFACE:
SUBROUTINE PDAF_3dvar_analysis(step, dim_p, dim_obs_p, dim_ens, &
     state_p, Ainv, ens_p, state_inc_p, forget, forget_ana, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_init_obsvar, U_prodRinvA, &
     screen, incremental, type_forget, flag)

! !DESCRIPTION:
! Analysis step of 3DVAR.
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
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: mype, MPIerr, COMM_filter, MPI_SUM, MPI_REALTYPE
  USE PDAF_mod_filter, &
       ONLY: type_trans, filterstr, obs_member, observe_ens

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step         ! Current time step
  INTEGER, INTENT(in) :: dim_p        ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p   ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens      ! Size of ensemble
  REAL, INTENT(out)   :: state_p(dim_p)          ! on exit: PE-local forecast state
  REAL, INTENT(out)   :: Ainv(dim_ens, dim_ens)  ! on entry: uninitialized
                                      ! on exit: weight matrix for ensemble transformation
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  REAL, INTENT(inout) :: state_inc_p(dim_p)      ! PE-local state analysis increment
  REAL, INTENT(in)    :: forget       ! Forgetting factor
  REAL, INTENT(out)   :: forget_ana   ! Forgetting factor actually used in analysis
  INTEGER, INTENT(in) :: screen       ! Verbosity flag
  INTEGER, INTENT(in) :: incremental  ! Control incremental updating
  INTEGER, INTENT(in) :: type_forget  ! Type of forgetting factor
  INTEGER, INTENT(inout) :: flag      ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, & ! Initialize dimension of observation vector
       U_obs_op, &              ! Observation operator
       U_init_obsvar, &         ! Initialize mean observation error variance
       U_init_obs, &            ! Initialize observation vector
       U_prodRinvA              ! Provide product R^-1 A

! !CALLING SEQUENCE:
! Called by: PDAF_3dvar_update
! Calls: U_init_dim_obs
! Calls: U_obs_op
! Calls: U_init_obs
! Calls: U_prodRinvA
! Calls: PDAF_timeit
! Calls: PDAF_memcount
! Calls: PDAF_set_forget
! Calls: PDAF_3dvar_Tright
! Calls: PDAF_3dvar_Tleft
! Calls: PDAF_generate_rndmat
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: gemvTYPE (BLAS; dgemv or sgemv dependent on precision)
! Calls: syevTYPE (LAPACK; dsyev or ssyev dependent on precision)
! Calls: MPI_allreduce (MPI)
!EOP

! *** local variables ***
  INTEGER :: i, member, col, row      ! Counters
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  INTEGER :: syev_info                ! Status flag for SYEV
  INTEGER :: ldwork                   ! Size of work array for SYEV
  INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
  REAL :: invdimens                   ! Inverse global ensemble size
  REAL :: sqrtNm1                     ! Temporary variable: sqrt(dim_ens-1)
!  REAL, ALLOCATABLE :: HZ_p(:,:)      ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: Riresid_p(:,:)    ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: resid_p(:)     ! PE-local observation residual
  REAL, ALLOCATABLE :: obs_p(:)       ! PE-local observation vector
  REAL, ALLOCATABLE :: HXbar_p(:)     ! PE-local observed state
  REAL, ALLOCATABLE :: RiHZd(:)       ! Temporary vector for analysis 
  REAL, ALLOCATABLE :: RiHZd_p(:)     ! PE-local RiHZd
  REAL, ALLOCATABLE :: VRiHZd(:)      ! Temporary vector for analysis
  REAL, ALLOCATABLE :: tmp_Ainv(:,:)  ! Temporary storage of Ainv
  REAL, ALLOCATABLE :: Asqrt(:, :)    ! Square-root of matrix Ainv
  REAL, ALLOCATABLE :: rndmat(:,:)    ! Temporary random matrix
  REAL, ALLOCATABLE :: ens_blk(:,:)   ! Temporary block of state ensemble
  REAL, ALLOCATABLE :: svals(:)       ! Singular values of Ainv
  REAL, ALLOCATABLE :: work(:)        ! Work array for SYEV
  INTEGER :: incremental_dummy        ! Dummy variable to avoid compiler warning
  REAL :: state_inc_p_dummy(1)        ! Dummy variable to avoid compiler warning
  REAL, ALLOCATABLE :: z_p(:,:)       ! ensemble perturbation matrix
  REAL, ALLOCATABLE :: svdV(:,:)      ! Transpose of right singular vectors
  REAL :: svec(1)
  REAL :: J_B, J_obs


! **********************
! *** INITIALIZATION ***
! **********************

  CALL PDAF_timeit(51, 'new')

  ! Initialize variable to prevent compiler warning
  incremental_dummy = incremental
  state_inc_p_dummy(1) = state_inc_p(1)


  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 1x, i7, 3x, a)') &
          'PDAF', step, 'Assimilating observations - 3DVAR'
  END IF


! ***********************************
! *** Compute mean forecast state ***
! ***********************************

  CALL PDAF_timeit(11, 'new')

  state_p = 0.0
  invdimens = 1.0 / REAL(dim_ens)
  DO member = 1, dim_ens
     DO row = 1, dim_p
        state_p(row) = state_p(row) + invdimens * ens_p(row, member)
     END DO
  END DO
  
  CALL PDAF_timeit(11, 'old')
  CALL PDAF_timeit(51, 'new')


! *********************************
! *** Get observation dimension ***
! *********************************

  CALL PDAF_timeit(15, 'new')
  CALL U_init_dim_obs(step, dim_obs_p)
  CALL PDAF_timeit(15, 'old')
  
  IF (screen > 2) THEN
     WRITE (*, '(a, 5x, a13, 1x, i6, 1x, a, i10)') &
          'PDAF', '--- PE-domain', mype, 'dimension of observation vector', dim_obs_p
  END IF


! ************************
! *** Compute residual ***
! ***   d = H x - y    ***
! ************************

  CALL PDAF_timeit(12, 'new')
  
  haveobsB: IF (dim_obs_p > 0) THEN
     ! *** The residual only exists for domains with observations ***

     ALLOCATE(resid_p(dim_obs_p))
     ALLOCATE(obs_p(dim_obs_p))
     ALLOCATE(HXbar_p(dim_obs_p))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3 * dim_obs_p)

     ! Project state onto observation space
!      IF (.NOT.observe_ens) THEN
        obs_member = 0 ! Store member index (0 for central state)
        CALL PDAF_timeit(44, 'new')
        CALL U_obs_op(step, dim_p, dim_obs_p, state_p, HXbar_p)
        CALL PDAF_timeit(44, 'old')
!      ELSE
!         ! For nonlinear H: apply H to each ensemble state; then average
!         ALLOCATE(HZ_p(dim_obs_p, dim_ens))
!         IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * dim_ens)
! 
!         CALL PDAF_timeit(44, 'new')
!         ENS1: DO member = 1, dim_ens
!            ! Store member index to make it accessible with PDAF_get_obsmemberid
!            obs_member = member
! 
!            ! [Hx_1 ... Hx_N]
!            CALL U_obs_op(step, dim_p, dim_obs_p, ens_p(:, member), HZ_p(:, member))
!         END DO ENS1
!         CALL PDAF_timeit(44, 'old')
! 
!         CALL PDAF_timeit(51, 'new')
!         HXbar_p = 0.0
!         DO member = 1, dim_ens
!            DO row = 1, dim_obs_p
!               HXbar_p(row) = HXbar_p(row) + invdimens * HZ_p(row, member)
!            END DO
!         END DO
!         CALL PDAF_timeit(51, 'old')
!      END IF

     ! get observation vector
     CALL PDAF_timeit(50, 'new')
     CALL U_init_obs(step, dim_obs_p, obs_p)
     CALL PDAF_timeit(50, 'old')

     ! Get residual as difference of observation and observed state
     CALL PDAF_timeit(51, 'new')
     resid_p = HXbar_p - obs_p
     CALL PDAF_timeit(51, 'old')

     DEALLOCATE(HXbar_p)

  END IF haveobsB

  CALL PDAF_timeit(12, 'old')


! **********************************************
! ***   Compute analyzed matrix Ainv         ***
! ***                                        ***
! ***     -1                 T  -1           ***
! ***    A  = forget I + (HZ)  R   HZ        ***
! **********************************************

  CALL PDAF_timeit(10, 'new')

!  ALLOCATE(Asqrt(dim_ens, dim_ens))
!  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

  haveobsA: IF (dim_obs_p > 0) THEN
     ! *** The contribution of observation matrix ist only ***
     ! *** computed for domains with observations          ***


     CALL PDAF_timeit(31, 'new')


     ! ***                RiHZ = Rinv HZ                
     ! *** This is implemented as a subroutine thus that
     ! *** Rinv does not need to be allocated explicitly.
     ALLOCATE(Riresid_p(dim_obs_p, 1))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * dim_ens)

     CALL PDAF_timeit(48, 'new')
     CALL U_prodRinvA(step, dim_obs_p, 1, obs_p, resid_p, Riresid_p)
     CALL PDAF_timeit(48, 'old')

     ! ***  Compute  J_obs ***

     CALL PDAF_timeit(51, 'new')

     J_obs = 0.0
     DO i = 1, dim_obs_p
        J_obs = J_obs + resid_p(i)*Riresid_p(i,1)
     END DO

     J_obs = 0.5*J_obs

     CALL PDAF_timeit(51, 'old')

  END IF haveobsA

!   ! get total sum on all filter PEs
!   ALLOCATE(tmp_Ainv(dim_ens, dim_ens))
!   IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)
! 
!   CALL MPI_allreduce(Asqrt, tmp_Ainv, dim_ens**2, &
!        MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

  CALL PDAF_timeit(31, 'old')
  CALL PDAF_timeit(10, 'old')


! ***********************************************
! *** Compute weight for model state update   ***
! ***                                         ***
! ***              T                    f     ***
! ***    w = A RiHZ d  with d = (y - H x )    ***
! ***                                         ***
! ***********************************************

  CALL PDAF_timeit(51, 'new')
  CALL PDAF_timeit(13, 'new')

  ALLOCATE(RiHZd(dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * dim_ens)



  ! *** Compute weight vector for state analysis:        ***
  ! ***          w = A RiHZd                             ***
  ! *** Use singular value decomposition of Ainv         ***
  ! ***        Ainv = USV^T                              ***
  ! *** Then: A = U S^(-1) V                             ***
  ! *** The decomposition is also used for the symmetric ***
  ! *** square-root for the ensemble transformation.     ***

  ALLOCATE(z_p(dim_p, dim_ens))

  ! Get perturbation matrix X'=X-xmean
  DO col = 1, dim_ens
     DO row = 1, dim_p
        z_p(row, col) = ens_p(row, col) - state_p(row)
     END DO
  END DO


  ! *** Invert B using SVD
  ALLOCATE(svals(dim_p))
  ALLOCATE(svdV(dim_ens, dim_p))
  ALLOCATE(work(MAX(3 * MIN(dim_p, dim_ens) + &
       MAX(dim_p, dim_ens), 5 * MIN(dim_p, dim_ens))))
  ldwork = MAX(3 * MIN(dim_p, dim_ens) + &
          MAX(dim_p, dim_ens), 5 * MIN(dim_p, dim_ens))
!  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3 * dim_ens)
  
  
  ! Compute SVD of ens_p
  CALL gesvdTYPE('n', 'S', dim_p, dim_ens, z_p, &
       dim_p, svals, svec, dim_p, svdV, &
       dim_p, work, ldwork, flag)

  DEALLOCATE(work)

  ! Check if SVD was successful
  IF (syev_info == 0) THEN
     flag = 0
  ELSE
     WRITE (*, '(/5x, a/)') 'PDAF-ERROR(1): Problem in SVD of inverse of A !!!'
     flag = 1
  END IF

  ! *** Compute w = A RiHZd stored in RiHZd
  check0: IF (flag == 0) THEN

     ALLOCATE(VRiHZd(dim_p))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)

     CALL gemvTYPE('n', dim_ens, dim_p, 1.0, svdV, &
          dim_ens, state_p, 1, 0.0, VRiHZd, 1)

     DO row = 1, dim_ens
        VRiHZd(row) = VRiHZd(row) / svals(row)
     END DO

     J_B = 0.0
     DO i=1, dim_p
        J_B = J_B + VRiHZd(i)*VRiHZd(i)
     END DO
     J_B = 0.5*J_B

     DEALLOCATE(VRiHZd)

  END IF check0

write (*,*) 'J: ', J_B, J_obs

  CALL PDAF_timeit(13, 'old')

  CALL PDAF_timeit(20, 'old')


  CALL PDAF_timeit(51, 'old')


! ********************
! *** Finishing up ***
! ********************
!endif
  IF (allocated(obs_p)) DEALLOCATE(obs_p)
  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_3dvar_analysis
