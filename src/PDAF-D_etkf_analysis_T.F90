! Copyright (c) 2004-2016 Lars Nerger
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
!$Id: PDAF-D_etkf_analysis_T.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_etkf_analysis_T --- ETKF using T-matrix
!
! !INTERFACE:
SUBROUTINE PDAF_etkf_analysis_T(step, dim_p, dim_obs_p, dim_ens, &
     state_p, Ainv, ens_p, state_inc_p, forget, forget_ana, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_init_obsvar, U_prodRinvA, &
     screen, incremental, type_forget, flag)

! !DESCRIPTION:
! Analysis step of the ETKF using a matrix T analogous
! to the SEIK filter. The consistent use of the operation
! of removing the mean from an ensemble matrix in form of
! a linear transformation (T-matrix) reduces the computational
! complexity.
!
! The implementation also supports an adaptive forgetting factor.
!
! Variant for domain decomposed states.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2009-07 - Lars Nerger - Initial code
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
       ONLY: type_trans, filterstr, obs_member

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
! Called by: PDAF_etkf_update
! Calls: U_init_dim_obs
! Calls: U_obs_op
! Calls: U_init_obs
! Calls: U_prodRinvA
! Calls: PDAF_timeit
! Calls: PDAF_memcount
! Calls: PDAF_set_forget
! Calls: PDAF_etkf_Tright
! Calls: PDAF_etkf_Tleft
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
  REAL, ALLOCATABLE :: HZ_p(:,:)      ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: RiHZ_p(:,:)    ! Temporary matrices for analysis
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

  
! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 1x, i7, 3x, a)') &
          'PDAF', step, 'Assimilating observations - ETKF using T-matrix'
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


! *********************************
! *** Get observation dimension ***
! *********************************

  CALL PDAF_timeit(15, 'new')
  CALL U_init_dim_obs(step, dim_obs_p)
  CALL PDAF_timeit(15, 'old')
  
  IF (screen > 0) THEN
     WRITE (*, '(a, 5x, a13, 1x, i6, 1x, a, i10)') &
          'PDAF', '--- PE-domain', mype, 'dimension of observation vector', dim_obs_p
  END IF


! ************************
! *** Compute residual ***
! ***   d = y - H x    ***
! ************************

  CALL PDAF_timeit(12, 'new')
  
  haveobsB: IF (dim_obs_p > 0) THEN
     ! *** The residual only exists for domains with observations ***

     ALLOCATE(resid_p(dim_obs_p))
     ALLOCATE(obs_p(dim_obs_p))
     ALLOCATE(HXbar_p(dim_obs_p))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3 * dim_obs_p)

     ! Project state onto observation space
     obs_member = 0 ! Store member index (0 for central state)
     CALL U_obs_op(step, dim_p, dim_obs_p, state_p, HXbar_p)

     ! get observation vector
     CALL U_init_obs(step, dim_obs_p, obs_p)

     ! Get residual as difference of observation and observed state
     resid_p = obs_p - HXbar_p

  END IF haveobsB

  CALL PDAF_timeit(12, 'old')


! **********************************************
! ***   Compute analyzed matrix Ainv         ***
! ***                                        ***
! ***     -1                 T  -1           ***
! ***    A  = forget I + (HZ)  R   HZ        ***
! **********************************************

  CALL PDAF_timeit(10, 'new')

  ALLOCATE(Asqrt(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

  haveobsA: IF (dim_obs_p > 0) THEN
     ! *** The contribution of observation matrix ist only ***
     ! *** computed for domains with observations          ***

     CALL PDAF_timeit(30, 'new')

     ! *** Compute HZ = [Hx_1 ... Hx_N] T
     ALLOCATE(HZ_p(dim_obs_p, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * dim_ens)

     ENS: DO member = 1, dim_ens
        ! Store member index to make it accessible with PDAF_get_obsmemberid
        obs_member = member

        ! [Hx_1 ... Hx_N]
        CALL U_obs_op(step, dim_p, dim_obs_p, ens_p(:, member), HZ_p(:, member))
     END DO ENS

     ! Set forgetting factor
     forget_ana = forget
     IF (type_forget == 1) THEN
        CALL PDAF_set_forget(step, filterstr, dim_obs_p, dim_ens, HZ_p, &
             HXbar_p, obs_p, U_init_obsvar, forget, forget_ana)
     ENDIF
     DEALLOCATE(HXbar_p)

     ! Subtract ensemble mean: HZ = [Hx_1 ... Hx_N] T
     CALL PDAF_etkf_Tright(dim_obs_p, dim_ens, HZ_p)

     CALL PDAF_timeit(30, 'old')
     CALL PDAF_timeit(31, 'new')


     ! ***                RiHZ = Rinv HZ                
     ! *** This is implemented as a subroutine thus that
     ! *** Rinv does not need to be allocated explicitly.
     ALLOCATE(RiHZ_p(dim_obs_p, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * dim_ens)

     CALL U_prodRinvA(step, dim_obs_p, dim_ens, obs_p, HZ_p, RiHZ_p)
     DEALLOCATE(obs_p)

     ! *** Initialize Ainv = (N-1) I ***
     Ainv = 0.0
     DO i = 1, dim_ens
        Ainv(i, i) = REAL(dim_ens - 1)
     END DO

     ! ***             T        ***
     ! ***  Compute  HZ  RiHZ   ***

     CALL gemmTYPE('t', 'n', dim_ens, dim_ens, dim_obs_p, &
          1.0, HZ_p, dim_obs_p, RiHZ_p, dim_obs_p, &
          0.0, Asqrt, dim_ens)

     DEALLOCATE(HZ_p)

  ELSE haveobsA
     ! *** For domains with dim_obs_p=0 there is no ***
     ! *** direct observation-contribution to Ainv  ***

     CALL PDAF_timeit(31, 'new')
    
     ! Set forgetting factor
     forget_ana = forget
    
     ! *** Initialize Ainv = (N-1) I ***
     Ainv = 0.0
     DO i = 1, dim_ens
        Ainv(i, i) = REAL(dim_ens - 1)
     END DO

     ! No observation-contribution to Ainv from this domain
     Asqrt = 0.0

  END IF haveobsA

  ! get total sum on all filter PEs
  ALLOCATE(tmp_Ainv(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

  CALL MPI_allreduce(Asqrt, tmp_Ainv, dim_ens**2, &
       MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

  ! *** Complete computation of Ainv ***
  ! ***   -1          -1    T        ***
  ! ***  A  = forget A  + HZ RiHZ    ***
  Ainv = forget_ana * Ainv + tmp_Ainv

  CALL PDAF_timeit(31, 'old')
  CALL PDAF_timeit(10, 'old')


! ***********************************************
! *** Compute weight for model state update   ***
! ***                                         ***
! ***              T                    f     ***
! ***    w = A RiHZ d  with d = (y - H x )    ***
! ***                                         ***
! ***********************************************

  CALL PDAF_timeit(13, 'new')

  ! *** Compute RiHZd = RiHZ^T d ***
  ALLOCATE(RiHZd_p(dim_ens))
  ALLOCATE(RiHZd(dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * dim_ens)

  haveobsC: IF (dim_obs_p > 0) THEN
     ! *** RiHZd_p/=0 only with observations
    
     ! local products (partial sum)
     CALL gemvTYPE('t', dim_obs_p, dim_ens, 1.0, RiHZ_p, &
          dim_obs_p, resid_p, 1, 0.0, RiHZd_p, 1)

     DEALLOCATE(RiHZ_p, resid_p)

  ELSE haveobsC

     RiHZd_p = 0.0

  END IF haveobsC

  ! get total sum on all filter PEs
  CALL MPI_allreduce(RiHZd_p, RiHZd, dim_ens, &
       MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

  DEALLOCATE(RiHZd_p)


  ! *** Compute weight vector for state analysis:        ***
  ! ***          w = A RiHZd                             ***
  ! *** Use singular value decomposition of Ainv         ***
  ! ***        Ainv = USV^T                              ***
  ! *** Then: A = U S^(-1) V                             ***
  ! *** The decomposition is also used for the symmetric ***
  ! *** square-root for the ensemble transformation.     ***

  ! *** Invert Ainv using SVD
  ALLOCATE(svals(dim_ens))
  ALLOCATE(work(3 * dim_ens))
  ldwork = 3 * dim_ens
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3 * dim_ens)
    
  ! Compute SVD of Ainv
  CALL syevTYPE('v', 'l', dim_ens, Ainv, dim_ens, svals, work, ldwork, syev_info)

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

     ALLOCATE(VRiHZd(dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)

     CALL gemvTYPE('t', dim_ens, dim_ens, 1.0, Ainv, &
          dim_ens, RiHZd, 1, 0.0, VRiHZd, 1)

     DO row = 1, dim_ens
        VRiHZd(row) = VRiHZd(row) / svals(row)
     END DO
  
     CALL gemvTYPE('n', dim_ens, dim_ens, 1.0, Ainv, &
          dim_ens, VRiHZd, 1, 0.0, RiHZd, 1)

     DEALLOCATE(VRiHZd)

  END IF check0

  CALL PDAF_timeit(13, 'old')


! ************************************************
! ***     Transform state ensemble             ***
! ***              a   _f   f                  ***
! ***             X  = X + X  W                ***
! *** The weight matrix W is stored in Ainv.   ***
! ************************************************

! *** Prepare weight matrix for ensemble transformation ***

  check1: IF (flag == 0) THEN

     IF (mype == 0 .AND. screen > 0) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', 'Perform ensemble transformation'
     END IF

     CALL PDAF_timeit(20, 'new')

     ! Part 1: square-root of A
     DO col = 1, dim_ens
        DO row = 1, dim_ens
           tmp_Ainv(row, col) = Ainv(row, col) / SQRT(svals(col))
        END DO
     END DO
  
     sqrtNm1 = SQRT(REAL(dim_ens-1))
     CALL gemmTYPE('n', 't', dim_ens, dim_ens, dim_ens, &
          sqrtNm1, tmp_Ainv, dim_ens, Ainv, dim_ens, &
          0.0, Asqrt, dim_ens)


     ! Part 2 - Optional 
     ! Multiply by orthogonal random matrix with eigenvector (1,...,1)^T
     multrnd: IF (type_trans == 2) THEN
        WRITE (*,'(a, 5x, a)') 'PDAF', '--- Apply random rotation to ensemble'

        ALLOCATE(rndmat(dim_ens, dim_ens))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)
        
        ! Initialize random matrix
        CALL PDAF_generate_rndmat(dim_ens, rndmat, 2)

        CALL gemmTYPE('n', 'n', dim_ens, dim_ens, dim_ens, &
             1.0, Asqrt, dim_ens, rndmat, dim_ens, &
             0.0, tmp_Ainv, dim_ens)

        DEALLOCATE(rndmat)
     ELSE
        ! Non-random case
        tmp_Ainv = Asqrt
     END IF multrnd


     ! Part 3: W = sqrt(A) + w
     DO col = 1, dim_ens
        DO row = 1, dim_ens
           Ainv(row, col) = tmp_Ainv(row, col) + RiHZd(row)
        END DO
     END DO

     DEALLOCATE(tmp_Ainv, svals)
     DEALLOCATE(RiHZd, Asqrt)

     ! Part 4: T W
     CALL PDAF_etkf_Tleft(dim_ens, dim_ens, Ainv)

     CALL PDAF_timeit(20, 'old')

! *** Perform ensemble transformation ***

     ! Use block formulation for transformation
     maxblksize = 200
     IF (mype == 0 .AND. screen > 0) &
          WRITE (*, '(a, 5x, a, i5)') &
          'PDAF', '--- use blocking with size ', maxblksize
        
     ALLOCATE(ens_blk(maxblksize, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', maxblksize * dim_ens)

     blocking: DO blklower = 1, dim_p, maxblksize
           
        blkupper = MIN(blklower + maxblksize - 1, dim_p)

        ! Store forecast ensemble
        CALL PDAF_timeit(21, 'new')
        DO col = 1, dim_ens
           ens_blk(1 : blkupper - blklower + 1, col) &
                = ens_p(blklower : blkupper, col)
        END DO

        ! Store mean forecast in ensemble matrix
        DO col = 1,dim_ens
           ens_p(blklower : blkupper, col) = state_p(blklower : blkupper)
        END DO
        CALL PDAF_timeit(21, 'old')

        !                        a  _f   f
        ! Transform ensemble:   X = X + X  W
        CALL PDAF_timeit(22, 'new')

        CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
             1.0, ens_blk(1, 1), maxblksize, Ainv(1, 1), dim_ens, &
             1.0, ens_p(blklower, 1), dim_p)

        CALL PDAF_timeit(22, 'old')

     END DO blocking

     DEALLOCATE(ens_blk)

  END IF check1


! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_etkf_analysis_T
