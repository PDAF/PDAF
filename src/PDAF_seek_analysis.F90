! Copyright (c) 2004-2025 Lars Nerger
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
!
!> Perform SEEK analysis step
!!
!! Analysis step of the SEEK filter
!! with forgetting factor.
!!
!! Variant for domain decomposed states.
!!
!! !  This is a core routine of PDAF and 
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2003-10 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAF_seek_analysis

CONTAINS
SUBROUTINE PDAF_seek_ana(step, dim_p, dim_obs_p, dim_eof, state_p, &
     Ainv, V_p, forget, U_init_dim_obs, U_obs_op, &
     U_init_obs, U_prodRinvA, screen, incremental, flag)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE mpi
  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filter, &
       ONLY: state_inc, obs_member
  USE PDAF_mod_filtermpi, &
       ONLY: mype, MPIerr, COMM_filter

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step         !< Current time step
  INTEGER, INTENT(in) :: dim_p        !< PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p   !< PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_eof      !< Number of EOFs
  REAL, INTENT(inout) :: state_p(dim_p)        !< PE-local model state
  !< ! *** The covariance P is decomposed as P = V U V^T ***
  REAL, INTENT(inout) :: Ainv(dim_eof,dim_eof) !< Inverse of matrix U
  REAL, INTENT(inout) :: V_p(dim_p,dim_eof)    !< PE-local matrix V
  REAL, INTENT(in)    :: forget       !< Forgetting factor
  INTEGER, INTENT(in) :: screen       !< Verbosity flag
  INTEGER, INTENT(in) :: incremental  !< Control incremental updating
  INTEGER, INTENT(inout) :: flag      !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, &       !< Initialize dimension of observation vector
       U_obs_op, &                    !< Observation operator
       U_init_obs, &                  !< Initialize observation vector
       U_prodRinvA                    !< Provide product Rinv A for SEEK analysis

! *** local variables ***
  INTEGER :: j                        ! Counters
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  REAL, ALLOCATABLE :: HV_p(:,:)      ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: RiHV_p(:,:)    ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: Ainv_p(:,:)    ! local Ainv
  REAL, ALLOCATABLE :: Ainv_inc(:,:)  ! increment for Ainv
  REAL, ALLOCATABLE :: resid_p(:)     ! observation residual
  REAL, ALLOCATABLE :: obs_p(:)       ! observation vector
  REAL, ALLOCATABLE :: m_state_p(:)   ! state projected onto obs. space
  REAL, ALLOCATABLE :: RiHVd(:)       ! Temporary vector for analysis 
  REAL, ALLOCATABLE :: RiHVd_p(:)     ! local RiHVd
  REAL, ALLOCATABLE :: temp_Ainv(:,:) ! Temporary storage of Ainv
  INTEGER, ALLOCATABLE :: ipiv(:)     ! vector of pivot indices for GESV
  INTEGER :: gesv_info                ! Control flag for GESV

  
! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, i7, 3x, a)') 'PDAF ', step, 'Assimilating observations - SEEK'
  END IF


! *********************************
! *** Get observation dimension ***
! *********************************

  CALL U_init_dim_obs(step, dim_obs_p)

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a13, 1x, i6, 1x, a, i10)') &
          'PDAF', '--- PE-domain', mype, 'dimension of observation vector', dim_obs_p
  END IF


! ************************
! *** Compute residual ***
! ***   d = y - H x    ***
! ************************

  CALL PDAF_timeit(12, 'new')

  haveobsB: IF (dim_obs_p > 0) THEN
     ! The residual only exists for domains with observations

     ALLOCATE(resid_p(dim_obs_p))
     ALLOCATE(obs_p(dim_obs_p))
     ALLOCATE(m_state_p(dim_obs_p))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3 * dim_obs_p)

     ! *** Project state onto observation space and    ***
     ! *** compute observation residual (innovation) d ***

     ! Project state onto observation space
     obs_member = 0 ! Store member index (0 for central state)
     CALL U_obs_op(step, dim_p, dim_obs_p, state_p, m_state_p)

     ! get observation vector
     CALL U_init_obs(step, dim_obs_p, obs_p)

     ! get residual as difference of observation an projected state
     resid_p = obs_p - m_state_p

     DEALLOCATE(m_state_p)
  END IF haveobsB

  CALL PDAF_timeit(12, 'old')


! *****************************************
! ***   Compute analyzed matrix Ainv    ***
! ***                                   ***
! ***  -1          -1    T  T  -1       ***
! *** U  = forget*U   + V  H  R   H  V  ***
! ***  i           i-1   i  i  i   i  i ***
! *****************************************

  CALL PDAF_timeit(10, 'new')

  haveobsA: IF (dim_obs_p > 0) THEN
     ! *** The contribution of observation matrix ist only ***
     ! *** computed for domains with observations          ***

     CALL PDAF_timeit(30, 'new')

     ! HV = H V
     ALLOCATE(HV_p(dim_obs_p, dim_eof))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * dim_eof)

     DO j = 1, dim_eof
        ! Store member index to make it accessible with PDAF_get_obsmemberid
        obs_member = j

        ! Call observation operator
        CALL U_obs_op(step, dim_p, dim_obs_p, V_p(:, j), HV_p(:, j))
     ENDDO

     CALL PDAF_timeit(30,'old')

     ! ***                RiHV = Rinv HV                 ***
     ! *** this is implemented as a subroutine thus that ***
     ! *** Rinv does not need to be allocated explicitly ***
     ALLOCATE(RiHV_p(dim_obs_p, dim_eof))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * dim_eof)

     CALL U_prodRinvA(step, dim_obs_p, dim_eof, obs_p, HV_p, RiHV_p)

     ! *** Finish computation of Ainv  ***
     ! ***   -1          -1    T       ***
     ! ***  U  = forget U  + HV  RiHV  ***
     ALLOCATE(Ainv_p(dim_eof, dim_eof))
     ALLOCATE(Ainv_inc(dim_eof, dim_eof))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * dim_eof * dim_eof)

     CALL PDAF_timeit(31, 'new')
     ! partial increment
     CALL gemmTYPE('t', 'n', dim_eof, dim_eof, dim_obs_p, &
          1.0, HV_p, dim_obs_p, RiHV_p, dim_obs_p, &
          0.0, Ainv_p, dim_eof)

     DEALLOCATE(HV_p)

  ELSE haveobsA
     ! *** For domains with dim_obs_p=0 there is no ***
     ! *** direct observation-contribution to Ainv  ***
 
     CALL PDAF_timeit(31, 'new')

     ALLOCATE(Ainv_p(dim_eof, dim_eof))
     ALLOCATE(Ainv_inc(dim_eof, dim_eof))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * dim_eof * dim_eof)

     ! No observation-contribution to Ainv from this domain
     Ainv_p = 0.0
  END IF haveobsA

  ! get total increment on all filter PEs
  CALL MPI_allreduce(Ainv_p, Ainv_inc, dim_eof * dim_eof, &
       MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

  Ainv = forget * Ainv + Ainv_inc

  DEALLOCATE(Ainv_p, Ainv_inc)

  CALL PDAF_timeit(31, 'old')
  CALL PDAF_timeit(10, 'old')


! *******************************************************
! ***             update model state                  ***
! ***                                                 ***
! ***  a   f            f     f          T        f   ***
! *** x = x + K (y - H x ) = x + V U RiHV (y - H x )  ***
! ***                                                 ***
! *******************************************************


  ! ************************
  ! *** RiHVd = RiHV^T d ***
  ! ************************

  CALL PDAF_timeit(13, 'new')

  ALLOCATE(RiHVd(dim_eof))
  ALLOCATE(RiHVd_p(dim_eof))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * dim_eof)

  haveobsC: IF (dim_obs_p > 0) THEN
     ! *** RiHLd_p/=0 only with observations ***

     ! local products (partial sum)
     CALL gemvTYPE('t', dim_obs_p, dim_eof, 1.0, RiHV_p, &
          dim_obs_p, resid_p, 1, 0.0, RiHVd_p, 1)

     DEALLOCATE(RiHV_p, resid_p)
  ELSE haveobsC
     RiHVd_p = 0.0
  END IF haveobsC

  ! get total sum on all filter PEs
  CALL MPI_allreduce(RiHVd_p, RiHVd, dim_eof, &
       MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

  DEALLOCATE(RiHVd_p)


  ! ****************************************
  ! *** Compute  w = U RiHLd  by solving ***
  ! ***           -1                     ***
  ! ***          U  w = RiHLd            ***
  ! *** for w. We use the LAPACK         ***
  ! *** routine GESV.                    ***
  ! ****************************************

  ALLOCATE(temp_Ainv(dim_eof, dim_eof))
  ALLOCATE(ipiv(dim_eof))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_eof * dim_eof)
  IF (allocflag == 0) CALL PDAF_memcount(3, 'i', dim_eof)

  ! save matrix Ainv
  temp_Ainv = Ainv

  ! call solver (GESV - LU solver)
  CALL gesvTYPE(dim_eof, 1, temp_Ainv, dim_eof, ipiv, RiHVd, dim_eof, gesv_info)
  DEALLOCATE(temp_Ainv, ipiv)

  CALL PDAF_timeit(13, 'old')

  ! *** check whether solve was successful
  update: IF (gesv_info /= 0) THEN
     WRITE (*, '(/5x, a/)') 'PDAF-ERROR(1): Problem in solve for state analysis !!!'
     flag = 1
  ELSE

     CALL PDAF_timeit(14, 'new')

     ! **************************
     ! *** Update model state ***
     ! ***     a   f          ***
     ! ***   x = x + V RiHVd  ***
     ! **************************

     IF (incremental == 0) THEN
        ! Allocate only if no incremental updating is used. 
        ! With incremental STATE_INC is allocated in filter_init.
        ALLOCATE(state_inc(dim_p))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_p)
     END IF

     CALL gemvTYPE('n', dim_p, dim_eof, 1.0, V_p, &
          dim_p, RiHVd, 1, 0.0, state_inc, 1)
     DEALLOCATE(RiHVd)

     IF (incremental == 0) THEN
        ! update state only if incremental updating is not used
        state_p = state_p + state_inc
        DEALLOCATE(state_inc)
     END IF

     CALL PDAF_timeit(14, 'old')
        
  END IF update

    
! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_seek_ana


!-------------------------------------------------------------------------------
!> Perform rediagonalization of P in SEEK
!!
!! Re-orthogonalization of the modes V of the
!! low-rank approximated covariance matrix in
!! its decomposed form P = V U V$^T$.
!!
!! Compute eigenmodes of the matrix B = L$^T$ L = C D C$^T$
!! where L = V U$^{1/2}$ (from Cholesky decomposition)
!! and get new modes V as V = L C D$^{-1/2}$,
!! $D = diag(\lambda_1,...,\lambda_r),\ \lambda_i > \lambda_i+1$
!! and C matrix of corresponding eigenvectors.
!! The new U is given by the matrix D.
!!
!! Variant for domain decomposed states.
!!
!! New version to compute matrix B. More efficient for
!! dim $>>$ dim\_eof
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2003-10 - Lars Nerger - Initial code
!! *Other revisions - see repository log
!!
SUBROUTINE PDAF_seek_rediag(dim_p, dim_eof, Ainv, ens_p, subtype, &
     screen, flag)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE mpi
  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: mype, MPIerr, COMM_filter

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p    !< PE-Local state dimension
  INTEGER, INTENT(in) :: dim_eof  !< Number of EOFs
  REAL, INTENT(inout) :: Ainv(dim_eof,dim_eof) !< Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p,dim_eof)  !< PE-local matrix V
  INTEGER, INTENT(in) :: subtype  !< Filter subtype
  INTEGER, INTENT(in) :: screen   !< Verbosity flag
  INTEGER, INTENT(inout) :: flag  !< Status Flag
       
! *** local variables ***
  INTEGER :: row, col              ! counters
  INTEGER, SAVE :: allocflag = 0   ! Flag whether first time allocation is done
  INTEGER, ALLOCATABLE :: ipiv(:)  ! vector of pivot indices for SGESV
  REAL, ALLOCATABLE :: ev(:)       ! vector of eigenvalues
  REAL, ALLOCATABLE :: rwork(:)    ! workarray for eigenproblem
  REAL, ALLOCATABLE :: U(:,:)      ! temporary for covar matrix
  REAL, ALLOCATABLE :: L(:,:)      ! covariance matrix
  REAL, ALLOCATABLE :: Temp1(:,:)  ! temporary for covar matrix
  REAL, ALLOCATABLE :: B(:,:)      ! temporary for covar matrix
  INTEGER :: syev_info, gesv_info, potrf_info  ! Info flags for LAPACK calls


! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a)') 'PDAF', 'Re-orthogonalize covariance matrix modes'
  END IF


! **************************************
! *** Compute matrix B = A^T V^T V A ***
! **************************************

  CALL PDAF_timeit(20, 'new')
  ! *** Get U by inversion of Ainv ***
  ALLOCATE(U(dim_eof, dim_eof))
  ALLOCATE(ipiv(dim_eof))
  IF (allocflag == 0) CALL PDAF_memcount(4, 'r', dim_eof * dim_eof)
  IF (allocflag == 0) CALL PDAF_memcount(4, 'i', dim_eof)

  ! Initialize U
  U = 0.0
  DO row = 1, dim_eof
     U(row, row) = 1.0
  END DO

  ! call solver
  CALL PDAF_timeit(32, 'new')
  CALL gesvTYPE(dim_eof, dim_eof, Ainv, dim_eof, ipiv, &
       U, dim_eof, gesv_info)
  CALL PDAF_timeit(32, 'old')

  DEALLOCATE(ipiv)

  ! Check if inversion was successful
  INVok: IF (gesv_info == 0) THEN

     ! *** Cholesky decomposition of U: U = A A^T ***
     CALL PDAF_timeit(33, 'new')
     CALL potrfTYPE('l', dim_eof, U, dim_eof, potrf_info)
     CALL PDAF_timeit(33, 'old')

     ! *** set upper elements to zero ***
     DO col = 2, dim_eof
        DO row = 1, col - 1
           U(row, col) = 0.0
        END DO
     END DO

     !*** Compute B = A^T V^T V A ***
     ALLOCATE(Temp1(dim_eof, dim_eof))
     ALLOCATE(B(dim_eof, dim_eof))
     IF (allocflag == 0) CALL PDAF_memcount(4, 'r', 2 * dim_eof * dim_eof)

     ! local V^T V
     CALL gemmTYPE('t', 'n', dim_eof, dim_eof, dim_p, &
          1.0, ens_p, dim_p, ens_p, dim_p, &
          0.0, Temp1, dim_eof)
     CALL PDAF_timeit(20, 'old')

     CALL PDAF_timeit(21, 'new')
     CALL MPI_allreduce(Temp1, B, dim_eof * dim_eof, MPI_REALTYPE, &
          MPI_SUM, COMM_filter, MPIerr)
     CALL PDAF_timeit(21, 'old')

     CALL PDAF_timeit(20, 'new')
     ! (V^T V) A (A stored in U)
     CALL gemmTYPE('n', 'n', dim_eof, dim_eof, dim_eof, &
          1.0, B, dim_eof, U, dim_eof, &
          0.0, Temp1, dim_eof)

     ! B = A^T (V^T V A) (A stored in U)
     CALL gemmTYPE('t', 'n', dim_eof, dim_eof, dim_eof, &
          1.0, U, dim_eof, Temp1, dim_eof, &
          0.0, B, dim_eof)
     CALL PDAF_timeit(31, 'old')


! *******************************
! *** Eigendecomposition of B ***
! ***       B = C D C^T       ***
! *******************************

     ALLOCATE(ev(dim_eof))
     ALLOCATE(rwork(3 * dim_eof))
     IF (allocflag == 0) CALL PDAF_memcount(4, 'r', 4 * dim_eof)
  
     CALL syevTYPE('v', 'u', dim_eof, B, dim_eof, &
          ev, rwork, 3 * dim_eof, syev_info)
     CALL PDAF_timeit(20, 'old')
     DEALLOCATE(rwork)

     ! check if eigendecomposition was successful
     EVPok: IF (syev_info == 0) THEN
        ! Eigendecomposition OK, continue
      
      ! *** Reorder matrix of eigenvectors ***
!       eof_mid = floor(real(dim_eof)/2.0)

!       do col=1,eof_mid
!         rwork(1:dim_eof) = U(:,col)
!         U(:,col) = U(:,dim_eof-col+1)
!         U(:,dim_eof-col+1) = rwork(1:dim_eof)
!       end do


! *****************************************************
! *** Initialize covar with re-orthogonalized modes ***
! *****************************************************

        CALL PDAF_timeit(22, 'new')
        ! AC = A C (AC stored in Temp1, A stored in U, C stored in B) 
        CALL gemmTYPE('n', 'n', dim_eof, dim_eof, dim_eof, &
             1.0, U, dim_eof, B, dim_eof, &
             0.0, Temp1, dim_eof)

        ! initialize L from V
        ALLOCATE(L(dim_p, dim_eof))
        IF (allocflag == 0) CALL PDAF_memcount(4, 'r', dim_eof * dim_p)

        L = ens_p

        ! V = L AC (AC stored in Temp1)
        CALL gemmTYPE('n', 'n', dim_p, dim_eof, dim_eof, &
             1.0, L, dim_p, Temp1, dim_eof, &
             0.0, ens_p, dim_p)
        DEALLOCATE(L, U, Temp1, B)

        unitmodes: IF (subtype /= 1) THEN
           ! use eigenvalues in U with unit modes in V
           !    U = diag(lambda_1,...,lambda_r)

           IF (mype == 0 .AND. screen > 0) THEN
              WRITE (*, '(a, 5x, a)') 'PDAF', '--- Use normalized modes'
           END IF

           Ainv = 0.0
           DO row = 1, dim_eof
              Ainv(row, row) = 1.0 / ev(row)
              ev(row) = SQRT(1.0 / ev(row))
           END DO

           ! Rescale modes V
           DO col = 1, dim_eof
              DO row = 1, dim_p
                 ens_p(row, col) = ens_p(row, col) * ev(col)
              END DO
           END DO
        ELSE unitmodes
           ! use unit U with modes in V scaled by eigenvalues

           IF (mype==0 .AND. screen > 0) THEN
              WRITE (*, '(a, 5x, a)') 'PDAF', '--- Use unit U'
           END IF

           Ainv = 0.0
           DO row = 1, dim_eof
              Ainv(row, row) = 1.0
           END DO
        END IF unitmodes
        CALL PDAF_timeit(22, 'old')
      
     ELSE
        ! eigendecomposition failed
        IF (mype == 0) WRITE (*, '(/5x, a/)') &
             'PDAF-ERROR(1): Problem with EOF decomposition of matrix Lt L !!!'
        flag = 1
        
     ENDIF EVPok

     DEALLOCATE(ev)

  ELSE INVok
     ! Inversion failed
     IF (mype == 0) WRITE (*, '(/5x, a/)') &
          'PDAF-ERROR(2): Problem with inversion of Ainv !!!'
     flag = 2

  ENDIF INVok


! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1
  
END SUBROUTINE PDAF_seek_rediag


END MODULE PDAF_seek_analysis
