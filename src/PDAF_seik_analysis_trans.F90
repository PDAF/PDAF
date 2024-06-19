! Copyright (c) 2004-2024 Lars Nerger
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
! !ROUTINE: PDAF_seik_analysis_trans --- SEIK analysis/ensemble transformation
!
! !INTERFACE:
SUBROUTINE PDAF_seik_analysis_trans(step, dim_p, dim_obs_p, dim_ens, rank, &
     state_p, Uinv, ens_p, state_inc_p, forget, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_init_obsvar, U_prodRinvA, &
     screen, incremental, type_forget, type_sqrt, flag)

! !DESCRIPTION:
! Analysis step of the SEIK filter with direct
! transformation of the forecast into the 
! analysis ensemble. This variant does not
! compute the analysis state, but only the
! analysis ensemble, whose mean is the analysis
! state.
! Supported is also the adaptive forgetting factor.
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

  USE mpi
  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filter, &
       ONLY: Nm1vsN, type_trans, filterstr, obs_member, observe_ens, debug
  USE PDAF_mod_filtermpi, &
       ONLY: mype, MPIerr, COMM_filter
  USE PDAFomi, &
       ONLY: omi_n_obstypes => n_obstypes, omi_omit_obs => omit_obs

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step         ! Current time step
  INTEGER, INTENT(in) :: dim_p        ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p   ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens      ! Size of ensemble
  INTEGER, INTENT(in) :: rank         ! Rank of initial covariance matrix
  REAL, INTENT(inout) :: state_p(dim_p) ! on exit: PE-local forecast mean state
  REAL, INTENT(inout) :: Uinv(rank, rank)      ! Inverse of matrix U - temporary use only
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) ! PE-local state ensemble
  REAL, INTENT(inout) :: state_inc_p(dim_p)    ! PE-local state analysis increment
  REAL, INTENT(in)    :: forget       ! Forgetting factor
  INTEGER, INTENT(in) :: screen       ! Verbosity flag
  INTEGER, INTENT(in) :: incremental  ! Control incremental updating
  INTEGER, INTENT(in) :: type_forget  ! Type of forgetting factor
  INTEGER, INTENT(in) :: type_sqrt    ! Type of square-root of A
                                      ! (0): symmetric sqrt; (1): Cholesky decomposition
  INTEGER, INTENT(inout) :: flag      ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, & ! Initialize dimension of observation vector
       U_obs_op, &              ! Observation operator
       U_init_obsvar, &         ! Initialize mean observation error variance
       U_init_obs, &            ! Initialize observation vector
       U_prodRinvA              ! Provide product R^-1 A

! !CALLING SEQUENCE:
! Called by: PDAF_seek_update
! Calls: U_init_dim_obs
! Calls: U_obs_op
! Calls: U_init_obs
! Calls: U_prodRinvA
! Calls: PDAF_timeit
! Calls: PDAF_memcount
! Calls: PDAF_set_forget
! Calls: PDAF_seik_matrixT
! Calls: PDAF_seik_Uinv
! Calls: PDAF_seik_Omega
! Calls: PDAF_seik_TtimesA
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: gemvTYPE (BLAS; dgemv or sgemv dependent on precision)
! Calls: gesvTYPE (LAPACK; dgesv or sgesv dependent on precision)
! Calls: syevTYPE (LAPACK; dsyev or ssyev dependent on precision)
! Calls: potrfTYPE (LAPACK; dpotrf or spotrf dependent on precision)
! Calls: trtrsTYPE (LAPACK; dtrtrs or strtrs dependent on precision)
! Calls: MPI_allreduce (MPI)
!EOP

! *** local variables ***
  INTEGER :: i, j, member, col, row  ! counters
  INTEGER, SAVE :: allocflag = 0     ! Flag whether first time allocation is done
  INTEGER :: lib_info                ! Status flag for LAPACK calls
  INTEGER :: ldwork                  ! Size of work array for syev
  INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
  REAL    :: invdimens               ! Inverse global ensemble size
  REAL    :: fac                     ! Temporary variable sqrt(dim_ens) or sqrt(rank)
  REAL    :: forget_ana              ! Forgetting factor used for analysis
  LOGICAL :: storeOmega = .FALSE.    ! Store matrix Omega instead of recomputing it
  LOGICAL, SAVE :: firsttime = .TRUE.! Indicates first call to resampling
  REAL, ALLOCATABLE :: HL_p(:,:)     ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: RiHL_p(:,:)   ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: resid_p(:)    ! PE-local observation residual
  REAL, ALLOCATABLE :: obs_p(:)      ! PE-local observation vector
  REAL, ALLOCATABLE :: HXbar_p(:)    ! PE-local observed state
  REAL, ALLOCATABLE :: RiHLd(:)      ! Temporary vector for analysis 
  REAL, ALLOCATABLE :: RiHLd_p(:)    ! PE-local RiHLd
  REAL, ALLOCATABLE :: VRiHLd(:)     ! Temporary vector for analysis
  REAL, ALLOCATABLE :: Usqrt(:, :)   ! Square-root of matrix U
  REAL, ALLOCATABLE :: Uinv_p(:,:)   ! Uinv for PE-local domain
  REAL, ALLOCATABLE :: tmp_Uinv(:,:) ! Temporary storage of Uinv
  REAL, ALLOCATABLE :: Omega(:,:)    ! Orthogonal matrix Omega
  REAL, ALLOCATABLE :: OmegaT(:,:)   ! Transpose of Omega
  REAL, SAVE, ALLOCATABLE :: OmegaTsave(:,:) ! Saved transpose of Omega
  REAL, ALLOCATABLE :: TA(:,:)       ! Temporary matrix
  REAL, ALLOCATABLE :: ens_blk(:,:)  ! Temporary blocked state ensemble
  REAL, ALLOCATABLE :: svals(:)      ! Singular values of Uinv
  REAL, ALLOCATABLE :: work(:)       ! Work array for SYEV
  INTEGER, ALLOCATABLE :: ipiv(:)    ! vector of pivot indices for GESV
  INTEGER :: incremental_dummy       ! Dummy variable to avoid compiler warning
  REAL :: state_inc_p_dummy(1)       ! Dummy variable to avoid compiler warning

  
! **********************
! *** INITIALIZATION ***
! **********************

  CALL PDAF_timeit(51, 'new')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_seik_analysis -- START'

  ! Initialize variable to prevent compiler warning
  incremental_dummy = incremental
  state_inc_p_dummy(1) = state_inc_p(1)

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 1x, i7, 3x, a)') &
          'PDAF', step, 'Assimilating observations - SEIK with ensemble transformation'
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
  CALL PDAF_timeit(51, 'old')


! *********************************
! *** Get observation dimension ***
! *********************************

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_seik_analysis:', debug, '  dim_p', dim_p
     WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_seik_analysis -- call init_dim_obs'
  END IF

  CALL PDAF_timeit(15, 'new')
  CALL U_init_dim_obs(step, dim_obs_p)
  CALL PDAF_timeit(15, 'old')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_seik_analysis:', debug, '  dim_obs_p', dim_obs_p

  IF (screen > 2) THEN
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
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_etkf_analysis -- observe_ens', observe_ens
     IF (.NOT.observe_ens) THEN
        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, &
             'PDAF_etkf_analysis -- call obs_op for ensemble mean'

        obs_member = 0 ! Store member index (0 for central state)
        CALL PDAF_timeit(44, 'new')
        CALL U_obs_op(step, dim_p, dim_obs_p, state_p, HXbar_p)
        CALL PDAF_timeit(44, 'old')
     ELSE
        ! For nonlinear H: apply H to each ensemble state; then average
        ALLOCATE(HL_p(dim_obs_p, dim_ens))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * dim_ens)

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_seik_analysis -- call obs_op', dim_ens, 'times'

        CALL PDAF_timeit(44, 'new')
        ENS1: DO member = 1, dim_ens
           ! Store member index to make it accessible with PDAF_get_obsmemberid
           obs_member = member

           ! [Hx_1 ... Hx_(r+1)]
           CALL U_obs_op(step, dim_p, dim_obs_p, ens_p(:, member), HL_p(:, member))
        END DO ENS1
        CALL PDAF_timeit(44, 'old')

        CALL PDAF_timeit(51, 'new')
        HXbar_p = 0.0
        DO member = 1, dim_ens
           DO row = 1, dim_obs_p
              HXbar_p(row) = HXbar_p(row) + invdimens * HL_p(row, member)
           END DO
        END DO
        CALL PDAF_timeit(51, 'old')
     END IF

     ! get observation vector
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_seik_analysis -- call init_obs'

     CALL PDAF_timeit(50, 'new')
     CALL U_init_obs(step, dim_obs_p, obs_p)
     CALL PDAF_timeit(50, 'old')

     ! get residual as difference of observation and
     ! projected state
     CALL PDAF_timeit(51, 'new')
     resid_p = obs_p - HXbar_p
     CALL PDAF_timeit(51, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_seik_analysis:', debug, &
             'innovation d(1:min(dim_obs_p,10))', resid_p(1:min(dim_obs_p,10))
        WRITE (*,*) '++ PDAF-debug PDAF_seik_analysis:', debug, &
             'MIN/MAX of innovation', MINVAL(resid_p), MAXVAL(resid_p)
     END IF

     ! Omit observations with too high innovation
     IF (omi_omit_obs)  THEN
        CALL PDAF_timeit(51, 'new')
        CALL PDAFomi_omit_by_inno_cb(dim_obs_p, resid_p, obs_p)
        CALL PDAF_timeit(51, 'old')
     END IF

  ELSE IF (dim_obs_p == 0) THEN

     ! For OMI we need to call observation operator also for dim_obs_p=0
     ! in order to initialize the pointer to the observation types
     ! Further the observation operator has to be executed in cases
     ! in which the operation include a global communication
     IF (.NOT.observe_ens) THEN
        IF (omi_n_obstypes>0) THEN
           ALLOCATE(HXbar_p(1))
           obs_member = 0

           ! [Hx_1 ... Hx_N]
           CALL U_obs_op(step, dim_p, dim_obs_p, state_p, HXbar_p)

           DEALLOCATE(HXbar_p)
        ELSE
           ALLOCATE(HL_p(1,1))
           DO member = 1, dim_ens
              ! Store member index to make it accessible with PDAF_get_obsmemberid
              obs_member = member

              ! [Hx_1 ... Hx_N]
              CALL U_obs_op(step, dim_p, dim_obs_p, ens_p(:, member), HL_p(:, member))
           END DO
           DEALLOCATE(HL_p)
        END IF
     END IF

  END IF haveobsB

  CALL PDAF_timeit(12, 'old')


! *************************************************
! ***   Compute analyzed matrix Uinv            ***
! ***                                           ***
! ***  -1              T        T  -1           ***
! *** U  = forget*fac T T + (HL)  R  (HL)       ***
! ***                                           ***
! *** Here FAC is a scaling factor according    ***
! *** to the definition of the ensemble         ***
! *** covariance scaled by N^-1 (original SEIK) ***
! *** (N-1)^-1 (SEIK as ensemble KF)            ***
! *************************************************

  CALL PDAF_timeit(10, 'new')

  haveobsA: IF (dim_obs_p > 0) THEN
     ! *** The contribution of observation matrix ist only ***
     ! *** computed for domains with observations          ***

     CALL PDAF_timeit(30, 'new')

     IF (.NOT.observe_ens) THEN
        ! This part is only required if H is applied to the ensemble mean before

        ! *** Compute HL = [Hx_1 ... Hx_N] T
        ALLOCATE(HL_p(dim_obs_p, dim_ens))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * dim_ens)

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_seik_analysis -- call obs_op', dim_ens, 'times'

        CALL PDAF_timeit(44, 'new')
        ENS: DO member = 1, dim_ens
           ! Store member index to make it accessible with PDAF_get_obsmemberid
           obs_member = member

           ! [Hx_1 ... Hx_N]
           CALL U_obs_op(step, dim_p, dim_obs_p, ens_p(:, member), HL_p(:, member))
        END DO ENS
        CALL PDAF_timeit(44, 'old')
     END IF

     ! Set forgetting factor
     forget_ana = forget
     IF (type_forget == 1) THEN
        CALL PDAF_set_forget(step, filterstr, dim_obs_p, dim_ens, HL_p, &
             HXbar_p, obs_p, U_init_obsvar, forget, forget_ana)
     ENDIF
     DEALLOCATE(HXbar_p)

     ! Complete HL = [Hx_1 ... Hx_N] T
     CALL PDAF_timeit(51, 'new')
     CALL PDAF_seik_matrixT(dim_obs_p, dim_ens, HL_p)
     CALL PDAF_timeit(51, 'old')

     CALL PDAF_timeit(30, 'old')
     CALL PDAF_timeit(31, 'new')


     ! ***                RiHL = Rinv HL                 ***
     ! *** this is implemented as a subroutine thus that ***
     ! *** Rinv does not need to be allocated explicitly ***
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_seik_analysis -- call prodRinvA_l'

     ALLOCATE(RiHL_p(dim_obs_p, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * rank)

     CALL PDAF_timeit(48, 'new')
     CALL U_prodRinvA(step, dim_obs_p, rank, obs_p, HL_p, RiHL_p)
     CALL PDAF_timeit(48, 'old')
     DEALLOCATE(obs_p)
 
     CALL PDAF_timeit(51, 'new')

     ! *** Initialize Uinv = fac T^T T ***
     CALL PDAF_seik_Uinv(rank, Uinv)

     ! ***             T        ***
     ! ***  Compute  HL  RiHL   ***
     ALLOCATE(Uinv_p(rank, rank))
     ALLOCATE(tmp_Uinv(rank, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * rank**2)

     CALL gemmTYPE('t', 'n', rank, rank, dim_obs_p, &
          1.0, HL_p, dim_obs_p, RiHL_p, dim_obs_p, &
          0.0, Uinv_p, rank)

     DEALLOCATE(HL_p)

  ELSE haveobsA
     ! *** For domains with dim_obs_p=0 there is no ***
     ! *** direct observation-contribution to Uinv  ***

     CALL PDAF_timeit(31, 'new')
     CALL PDAF_timeit(51, 'new')
    
     ! Set forgetting factor
     forget_ana = forget

     ! Initialize Uinv = fac T^T T 
     CALL PDAF_seik_Uinv(rank, Uinv)

     ALLOCATE(Uinv_p(rank, rank))
     ALLOCATE(tmp_Uinv(rank, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * rank**2)

     ! No observation-contribution to Uinv from this domain
     Uinv_p = 0.0

     ! For OMI we need to call observation operator also for dim_obs_p=0
     ! in order to initialize the pointer to the observation types
     ! Further the observation operator has to be executed in cases
     ! in which the operation includes a global communication
     IF (omi_n_obstypes>0) THEN
        IF (.NOT.observe_ens) THEN
           ALLOCATE(HL_p(1,1))
           DO member = 1, dim_ens
              ! Store member index to make it accessible with PDAF_get_obsmemberid
              obs_member = member

              ! [Hx_1 ... Hx_N]
              CALL U_obs_op(step, dim_p, dim_obs_p, ens_p(:, member), HL_p(:, member))
           END DO
           DEALLOCATE(HL_p)
        END IF
     END IF

  END IF haveobsA

  ! get total sum on all filter PEs
  CALL MPI_allreduce(Uinv_p, tmp_Uinv, rank * rank, &
       MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

  ! *** Complete computation of Uinv  ***
  ! ***   -1          -1    T         ***
  ! ***  U  = forget U  + HL RiHL     ***

  Uinv = forget_ana * Uinv + tmp_Uinv

  DEALLOCATE(Uinv_p)

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_seik_analysis:', debug, '  U^-1', Uinv

  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(31, 'old')
  CALL PDAF_timeit(10, 'old')


! ***********************************************
! *** Compute weight for model state update   ***
! ***                                         ***
! ***              T                    f     ***
! ***    w = U RiHL d  with d = (y - H x )    ***
! ***********************************************

  CALL PDAF_timeit(51, 'new')
  CALL PDAF_timeit(13, 'new')

  ! *** RiHLd = RiHL^T d ***
  ALLOCATE(RiHLd_p(rank))
  ALLOCATE(RiHLd(rank))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * rank)

  haveobsC: IF (dim_obs_p > 0) THEN
     ! *** RiHLd_p/=0 only with observations ***
    
     ! local products (partial sum)
     CALL gemvTYPE('t', dim_obs_p, rank, 1.0, RiHL_p, &
          dim_obs_p, resid_p, 1, 0.0, RiHLd_p, 1)

     DEALLOCATE(RiHL_p, resid_p)

  ELSE haveobsC

     RiHLd_p = 0.0

  END IF haveobsC

  ! get total sum on all filter PEs
  CALL MPI_allreduce(RiHLd_p, RiHLd, rank, &
       MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

  DEALLOCATE(RiHLd_p)


  ! *** Compute weight vector for state analysis:     ***
  ! ***          w = U RiHLd                          ***
  ! *** For this, two variants are implemented:       ***
  ! *** 1. solve for w in:                            ***
  ! ***           -1                                  ***
  ! ***          U  w = RiHLd                         ***
  ! ***   We use the LAPACK routine GESV             ***
  ! *** 2. Compute singular value decomposition       ***
  ! ***   of Uinv: Uinv = ASB^T                       ***
  ! ***   Then: U = A S^(-1) B                        ***
  ! ***   This is combined with a symmetric           ***
  ! ***   square-root for the ensemble transformation ***

  typeuinv1: IF (type_sqrt==1) THEN
     ! *** Variant 1: Solve Uinv w= RiHLd for w

     ALLOCATE(ipiv(rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'i', rank)

     ! save matrix Uinv
     tmp_Uinv = Uinv

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_seik_analysis:', debug, &
          '  Invert U^-1 using solver GESV'

     ! call solver (GESV - LU solver)
     CALL gesvTYPE(rank, 1, tmp_Uinv, rank, ipiv, &
          RiHLd, rank, lib_info)

     DEALLOCATE(ipiv)

  ELSE typeuinv1
     ! *** Variant 2: Invert Uinv using SVD

     ALLOCATE(svals(rank))
     ALLOCATE(work(3 * rank))
     ldwork = 3 * rank
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3 * rank)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_seik_analysis:', debug, &
          '  Compute eigenvalue decomposition of U^-1'

     ! Compute SVD of Uinv
     CALL syevTYPE('v', 'l', rank, Uinv, rank, svals, work, ldwork, lib_info)

     DEALLOCATE(work)

     ! Compute product RiHLd U
     IF (lib_info==0) THEN
        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_seik_resample:', debug, '  eigenvalues', svals

        ALLOCATE(VRiHLd(rank))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank)

        CALL gemvTYPE('t', rank, rank, 1.0, Uinv, &
             rank, RiHLd, 1, 0.0, VRiHLd, 1)
     
        DO row = 1,rank
           VRiHLd(row) = VRiHLd(row) / svals(row)
        END DO
  
        CALL gemvTYPE('n', rank, rank, 1.0, Uinv, &
             rank, VRiHLd, 1, 0.0, RiHLd, 1)

        DEALLOCATE(VRiHLd)
     END IF
  END IF typeuinv1

  CALL PDAF_timeit(13, 'old')

  ! *** check if solve was successful
  IF (lib_info == 0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_seik_analysis:', debug, '  U(HL R^-1)^T d', RiHLd

     flag = 0
  ELSE
     WRITE (*, '(/5x, a/)') 'PDAF-ERROR(1): Problem in computation of analysis weights!!!'
     flag = 1
  END IF


! ****************************************************************
! ***     Transform state ensemble                             ***
! ***              a   _f   f                                  ***
! ***             X  = X + X  W                                ***
! *** with                               -T      T             ***
! ***          W = T (RiHLd + sqrt(FAC) C   Omega )            ***
! *** Here FAC depends on the use definition of the covariance ***
! *** matrix P using a factor (r+1)^-1 or r^-1.                ***
! ****************************************************************

! *** Prepare weight matrix for ensemble transformation ***

  check1: IF (flag == 0) THEN

     IF (mype == 0 .AND. screen > 0) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', 'Transform state ensemble'
        IF (type_sqrt == 1) THEN
           WRITE (*, '(a, 5x, a)') 'PDAF', '--- use Cholesky square-root of U'
        ELSE
           WRITE (*, '(a, 5x, a)') 'PDAF', '--- use symmetric square-root of U'
        END IF
     END IF

     CALL PDAF_timeit(20, 'new')
     CALL PDAF_timeit(32, 'new')

     ! Usqrt is allocated with dim_ens cols, because this is 
     ! required further below. Now only rank columns are used
     ALLOCATE(Usqrt(rank, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens * rank)

     ! Part 1: square-root of U
     typeuinv2: IF (type_sqrt == 1) THEN
        ! Variant, if Uinv has been inverted above by solving

        ! *** Store Uinv for temporary use
        tmp_Uinv(:,:) = Uinv(:,:)

        ! Cholesky decomposition of tmp_Uinv = C C^T
        CALL potrfTYPE('l', rank, tmp_Uinv, rank, lib_info)

        ! check if Cholesky decomposition was successful
        CholeskyOK: IF (lib_info == 0) THEN
           ! Decomposition OK, continue
           flag = 0
        ELSE
           ! Decomposition failed
           WRITE (*, '(/5x, a/)') &
                'PDAF-ERROR(3): Problem with Cholesky decomposition of Uinv !!!'
           flag = 3
        ENDIF CholeskyOK

     ELSE typeuinv2
        ! Variant, if SVD inversion of Uinv has been performed

        DO col = 1, rank
           DO row = 1, rank
              Usqrt(row, col) = Uinv(row, col) / SQRT(svals(col))
           END DO
        END DO

        CALL gemmTYPE('n', 't', rank, rank, rank, &
             1.0, Usqrt, rank, Uinv, rank, &
             0.0, tmp_Uinv, rank)
        DEALLOCATE(svals)
        
        ! Set flag
        flag = 0

     END IF typeuinv2

     CALL PDAF_timeit(32, 'old')

  END IF check1

  
  ! *** Part 2: Initialize Omega ***

  check2: IF (flag == 0) THEN
     ! allocate fields
     ALLOCATE(OmegaT(rank, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank * dim_ens)
     
     IF (storeOmega .AND. allocflag == 0) THEN
        ALLOCATE(OmegaTsave(rank, dim_ens))
        CALL PDAF_memcount(3, 'r', dim_ens * rank)
     END IF

     CALL PDAF_timeit(33, 'new')
     Omega_store: IF (storeOmega) THEN
        first: IF (firsttime) THEN
           ! *** At first call to SEIK_RESAMPLE initialize   ***
           ! *** the matrix Omega in SEIK_Omega and store it ***

           ALLOCATE(Omega(dim_ens, rank))
           IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens * rank)

           ! *** Generate uniform orthogonal matrix OMEGA ***
           CALL PDAF_seik_Omega(rank, Omega, type_trans, screen)
        
           ! transpose Omega
           IF (type_sqrt == 1) THEN
              OmegaT = TRANSPOSE(Omega)
              ! store transposed Omega
              OmegaTsave = OmegaT
           ELSE
              Usqrt = TRANSPOSE(Omega)
              ! store transposed Omega
              OmegaTsave = Usqrt
           END IF

           firsttime  = .FALSE.
      
           DEALLOCATE(Omega)

        ELSE first
           IF (mype == 0 .AND. screen > 0) &
                WRITE (*, '(a, 5x, a)') 'PDAF', '--- use stored Omega'
           IF (type_sqrt == 1) THEN
              OmegaT = OmegaTsave
           ELSE
              Usqrt = OmegaTsave
           END IF
        END IF first

     ELSE Omega_store

        ! *** Initialize the matrix Omega in SEIK_Omega ***
        ! *** each time SEIK_RESAMPLE is called         ***

        ALLOCATE(Omega(dim_ens, rank))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens * rank)
      
        ! *** Generate uniform orthogonal matrix OMEGA ***
        CALL PDAF_seik_Omega(rank, Omega, type_trans, screen)
        
        ! transpose Omega
        IF (type_sqrt == 1) THEN
           OmegaT = TRANSPOSE(Omega)
        ELSE
           Usqrt = TRANSPOSE(Omega)
        END IF

        DEALLOCATE(Omega)

     END IF Omega_store
     IF (debug>0) THEN
        IF (type_sqrt == 1) THEN
           WRITE (*,*) '++ PDAF-debug PDAF_seik_analysis:', debug, '  Omega^T', OmegaT
        ELSE
           WRITE (*,*) '++ PDAF-debug PDAF_seik_analysis:', debug, '  Omega^T', Usqrt
        END IF
     END IF

     CALL PDAF_timeit(33, 'old')


     ! *** Part 3: Product U^(1/2) Omega ***
    
     CALL PDAF_timeit(34, 'new')
     IF (type_sqrt == 1) THEN
        ! A = (Omega C^(-1)) by solving Ct A = OmegaT for A
        CALL trtrsTYPE('L', 'T', 'N', rank, dim_ens, &
             tmp_Uinv, rank, OmegaT, rank, lib_info)
     ELSE
        ! TMP_UINV already contains matrix C (no more inversion)

        CALL gemmTYPE('n', 'n', rank, dim_ens, rank, &
             1.0, tmp_Uinv, rank, Usqrt, rank, &
             0.0, OmegaT, rank)
     END IF
     CALL PDAF_timeit(34, 'old')

     ! check if solve was successful
     solveOK: IF (lib_info == 0) THEN
        ! Solve for A OK, continue
        flag = 0
     ELSE
        ! Solve for A failed
        WRITE (*, '(/5x, a/)') &
             'PDAF-ERROR(2): Problem in computation of transformation matrix !!!'
        flag = 2

        CALL PDAF_timeit(20, 'old')

     END IF solveOK

     DEALLOCATE(Usqrt)
        
  END IF check2

  check3: IF (flag == 0) THEN

     ! *** Part 4: Add RiHLd and multiply by scaling factor

     CALL PDAF_timeit(35, 'new')

     IF (Nm1vsN == 1) THEN
        ! Use factor (N-1)^-1
        fac = SQRT(REAL(dim_ens - 1))
     ELSE
        ! Use factor N^-1
        fac = SQRT(REAL(dim_ens))
     END IF

     ! *** Add RiHLd to At
     DO j = 1, dim_ens
        DO i = 1, rank
           OmegaT(i,j) = fac * OmegaT(i,j) + RiHLd(i)
        END DO
     END DO
     DEALLOCATE(RiHLd)

     ! *** T A^T (A^T stored in OmegaT) ***
     ALLOCATE(TA(dim_ens, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

     CALL PDAF_seik_TtimesA(rank, dim_ens, OmegaT, TA)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_seik_resample:', debug, '  transform', TA

     CALL PDAF_timeit(35, 'old')
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

        ! Store old state ensemble
        CALL PDAF_timeit(21, 'new')
        DO col = 1, dim_ens
           ens_blk(1 : blkupper - blklower + 1, col) &
                = ens_p(blklower : blkupper, col)
        END DO

        DO col = 1,dim_ens
           ens_p(blklower : blkupper, col) = state_p(blklower : blkupper)
        END DO

        CALL PDAF_timeit(21, 'old')

        !                        a  _f   f    T
        ! Transform ensemble:   X = X + X  T(A )
        CALL PDAF_timeit(22, 'new')

        CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
             1.0, ens_blk(1, 1), maxblksize, TA(1, 1), dim_ens, &
             1.0, ens_p(blklower, 1), dim_p)

        CALL PDAF_timeit(22, 'old')

     END DO blocking

     DEALLOCATE(ens_blk, TA)
     DEALLOCATE(OmegaT)

  END IF check3

  CALL PDAF_timeit(51, 'old')


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(tmp_Uinv)

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_seik_analysis -- END'

END SUBROUTINE PDAF_seik_analysis_trans
