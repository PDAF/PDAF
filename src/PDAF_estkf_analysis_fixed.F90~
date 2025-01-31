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
! !ROUTINE: PDAF_estkf_analysis_fixed --- ESTKF analysis without ensemble transformation
!
! !INTERFACE:
SUBROUTINE PDAF_estkf_analysis_fixed(step, dim_p, dim_obs_p, dim_ens, rank, &
     state_p, Ainv, ens_p, state_inc_p, forget, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_init_obsvar, U_prodRinvA, &
     screen, incremental, type_forget, type_sqrt, flag)

! !DESCRIPTION:
! Analysis step of the ESTKF with direct update
! of the state estimate, but no ensemble 
! transformation. The ensemble is only shifted
! to represent the analysis state. This variant
! is used for the filter variant with a fixed
! covariance matrix.
! Supported is also the adaptive forgetting factor.
!
! Variant for domain decomposed states.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2012-03 - Lars Nerger - Initial code based on dynamic ESTKF
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
       ONLY: filterstr, obs_member, observe_ens, debug
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
  REAL, INTENT(inout) :: Ainv(rank, rank)      ! Inverse of matrix A - temporary use only
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) ! PE-local state ensemble
  REAL, INTENT(inout) :: state_inc_p(dim_p)    ! PE-local state analysis increment
  REAL, INTENT(in)    :: forget       ! Forgetting factor
  INTEGER, INTENT(in) :: screen       ! Verbosity flag
  INTEGER, INTENT(in) :: incremental  ! Control incremental updating
  INTEGER, INTENT(in) :: type_forget  ! Type of forgetting factor
  INTEGER, INTENT(in) :: type_sqrt   ! Type of square-root of A
                                     ! (0): symmetric sqrt; (1): Cholesky decomposition
  INTEGER, INTENT(inout) :: flag      ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, & ! Initialize dimension of observation vector
       U_obs_op, &              ! Observation operator
       U_init_obsvar, &         ! Initialize mean observation error variance
       U_init_obs, &            ! Initialize observation vector
       U_prodRinvA              ! Provide product R^-1 with some matrix

! !CALLING SEQUENCE:
! Called by: PDAF_estkf_update
! Calls: U_init_dim_obs
! Calls: U_obs_op
! Calls: U_init_obs
! Calls: U_prodRinvA
! Calls: PDAF_timeit
! Calls: PDAF_memcount
! Calls: PDAF_set_forget
! Calls: PDAF_estkf_AOmega
! Calls: PDAF_estkf_OmegaA
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: gemvTYPE (BLAS; dgemv or sgemv dependent on precision)
! Calls: gesvTYPE (LAPACK; dgesv or sgesv dependent on precision)
! Calls: syevTYPE (LAPACK; dsyev or ssyev dependent on precision)
! Calls: MPI_allreduce (MPI)
!EOP

! *** local variables ***
  INTEGER :: i, member, col, row     ! counters
  INTEGER, SAVE :: allocflag = 0     ! Flag whether first time allocation is done
  INTEGER :: lib_info                ! Status flag for LAPACK calls
  INTEGER :: ldwork                  ! Size of work array for syevTYPE
  REAL    :: invdimens               ! Inverse global ensemble size
  REAL    :: forget_ana              ! Forgetting factor used for analysis
  REAL, ALLOCATABLE :: HL_p(:,:)     ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: RiHL_p(:,:)   ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: resid_p(:)    ! PE-local observation residual
  REAL, ALLOCATABLE :: obs_p(:)      ! PE-local observation vector
  REAL, ALLOCATABLE :: HXbar_p(:)    ! PE-local observed state
  REAL, ALLOCATABLE :: RiHLd(:)      ! Temporary vector for analysis 
  REAL, ALLOCATABLE :: RiHLd_p(:)    ! PE-local RiHLd
  REAL, ALLOCATABLE :: VRiHLd(:)     ! Temporary vector for analysis
  REAL, ALLOCATABLE :: Ainv_p(:,:)   ! Ainv for PE-local domain
  REAL, ALLOCATABLE :: tmp_Ainv(:,:) ! Temporary storage of Ainv
  REAL, ALLOCATABLE :: TRiHLd(:,:)   ! Temporary vector for analysis 
  REAL, ALLOCATABLE :: svals(:)      ! Singular values of Ainv
  REAL, ALLOCATABLE :: work(:)       ! Work array for SYEVTYPE
  INTEGER, ALLOCATABLE :: ipiv(:)    ! vector of pivot indices for GESVTYPE

  
! **********************
! *** INITIALIZATION ***
! **********************

  CALL PDAF_timeit(51, 'new')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_estkf_analysis -- START'

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, i7, 3x, a)') &
          'PDAF ', step, 'Assimilating observations - ESTKF with fixed ensemble'
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
     WRITE (*,*) '++ PDAF-debug PDAF_estkf_analysis:', debug, '  dim_p', dim_p
     WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_estkf_analysis -- call init_dim_obs'
  END IF

  CALL PDAF_timeit(15, 'new')
  CALL U_init_dim_obs(step, dim_obs_p)
  CALL PDAF_timeit(15, 'old')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_estkf_analysis:', debug, '  dim_obs_p', dim_obs_p

  IF (screen > 2) THEN
     WRITE (*, '(a4, 5x, a13, 1x, i6, 1x, a, i10)') &
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
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_estkf_analysis -- observe_ens', observe_ens
     IF (.NOT.observe_ens) THEN
        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, &
             'PDAF_estkf_analysis -- call obs_op for ensemble mean'

        obs_member = 0 ! Store member index (0 for central state)
        CALL PDAF_timeit(44, 'new')
        CALL U_obs_op(step, dim_p, dim_obs_p, state_p, HXbar_p)
        CALL PDAF_timeit(44, 'old')
     ELSE
        ! For nonlinear H: apply H to each ensemble state; then average
        ALLOCATE(HL_p(dim_obs_p, dim_ens))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * dim_ens)

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_estkf_analysis -- call obs_op', dim_ens, 'times'

        CALL PDAF_timeit(44, 'new')
        ENS1: DO member = 1, dim_ens
           ! Store member index to make it accessible with PDAF_get_obsmemberid
           obs_member = member

           ! [Hx_1 ... Hx_N]
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
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_estkf_analysis -- call init_obs'

     CALL PDAF_timeit(50, 'new')
     CALL U_init_obs(step, dim_obs_p, obs_p)
     CALL PDAF_timeit(50, 'old')

     ! get residual as difference of observation and
     ! projected state
     CALL PDAF_timeit(51, 'new')
     resid_p = obs_p - HXbar_p
     CALL PDAF_timeit(51, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_estkf_analysis:', debug, &
             'innovation d(1:min(dim_obs_p,10))', resid_p(1:min(dim_obs_p,10))
        WRITE (*,*) '++ PDAF-debug PDAF_estkf_analysis:', debug, &
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
! ***   Compute analyzed matrix Ainv            ***
! ***                                           ***
! ***  -1                       T  -1           ***
! *** A  = forget*(N-1) I + (HL)  R  (HL)       ***
! ***                                           ***
! *************************************************

  CALL PDAF_timeit(10, 'new')

  haveobsA: IF (dim_obs_p > 0) THEN
     ! *** The contribution of observation matrix ist only ***
     ! *** computed for domains with observations          ***

     CALL PDAF_timeit(30, 'new')

     IF (.NOT.observe_ens) THEN
        ! This part is only required if H is applied to the ensemble mean before

        ! *** Compute HL = [Hx_1 ... Hx_N] Omega
        ALLOCATE(HL_p(dim_obs_p, dim_ens))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * dim_ens)

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_estkf_analysis -- call obs_op', dim_ens, 'times'

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
     CALL PDAF_estkf_AOmega(dim_obs_p, dim_ens, HL_p)
     CALL PDAF_timeit(51, 'old')

     CALL PDAF_timeit(30, 'old')
     CALL PDAF_timeit(31, 'new')


     ! ***                RiHL = Rinv HL                 ***
     ! *** this is implemented as a subroutine thus that ***
     ! *** Rinv does not need to be allocated explicitly ***
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_estkf_analysis -- call prodRinvA_l'

     ALLOCATE(RiHL_p(dim_obs_p, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * rank)

     CALL PDAF_timeit(48, 'new')
     CALL U_prodRinvA(step, dim_obs_p, rank, obs_p, HL_p, RiHL_p)
     CALL PDAF_timeit(48, 'old')
     DEALLOCATE(obs_p)
 
     CALL PDAF_timeit(51, 'new')

     ! *** Initialize Ainv = (N-1) I ***
     Ainv = 0.0
     DO i = 1, rank
        Ainv(i,i) = REAL(dim_ens - 1)
     END DO

     ! ***             T        ***
     ! ***  Compute  HL  RiHL   ***
     ALLOCATE(Ainv_p(rank, rank))
     ALLOCATE(tmp_Ainv(rank, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * rank**2)

     CALL gemmTYPE('t', 'n', rank, rank, dim_obs_p, &
          1.0, HL_p, dim_obs_p, RiHL_p, dim_obs_p, &
          0.0, Ainv_p, rank)

     DEALLOCATE(HL_p)

  ELSE haveobsA
     ! *** For domains with dim_obs_p=0 there is no ***
     ! *** direct observation-contribution to Ainv  ***

     CALL PDAF_timeit(31, 'new')
     CALL PDAF_timeit(51, 'new')
    
     ! Set forgetting factor
     forget_ana = forget

     ! *** Initialize Ainv = (N-1) I ***
     Ainv = 0.0
     DO i = 1, rank
        Ainv(i,i) = REAL(dim_ens - 1)
     END DO

     ALLOCATE(Ainv_p(rank, rank))
     ALLOCATE(tmp_Ainv(rank, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * rank**2)

     ! No observation-contribution to Ainv from this domain
     Ainv_p = 0.0

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
  CALL MPI_allreduce(Ainv_p, tmp_Ainv, rank * rank, &
       MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

  ! *** Complete computation of Ainv  ***
  ! ***   -1                T         ***
  ! ***  A  = forget I  + HL RiHL     ***

  Ainv = forget_ana * Ainv + tmp_Ainv

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_estkf_analysis:', debug, '  A^-1', Ainv

  DEALLOCATE(Ainv_p)

  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(31, 'old')
  CALL PDAF_timeit(10, 'old')


! ***********************************************
! *** Compute weight for model state update   ***
! ***                                         ***
! ***              T                    f     ***
! ***    w = A RiHL d  with d = (y - H x )    ***
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

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_estkf_analysis:', debug, '  (HXT R^-1)^T d', RiHLd

  DEALLOCATE(RiHLd_p)


  ! *** Compute weight vector for state analysis:     ***
  ! ***          w = A RiHLd                          ***
  ! *** For this, two variants are implemented:       ***
  ! *** 1. solve for w in:                            ***
  ! ***           -1                                  ***
  ! ***          A  w = RiHLd                         ***
  ! ***   We use the LAPACK routine gesvTYPE          ***
  ! *** 2. Compute singular value decomposition       ***
  ! ***   of Ainv: Ainv = USV^T                       ***
  ! ***   Then: A = U S^(-1) V^T                      ***
  ! ***   This is combined with a symmetric           ***
  ! ***   square-root for the ensemble transformation ***

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_estkf_analysis -- type_sqrt', type_sqrt
  typeainv1: IF (type_sqrt==1) THEN
     ! *** Variant 1: Solve Ainv w= RiHLd for w

     ALLOCATE(ipiv(rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'i', rank)

     ! save matrix Ainv
     tmp_Ainv = Ainv

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_estkf_analysis:', debug, &
          '  Invert A^-1 using solver GESV'

     ! call solver (GESVTYPE - LU solver)
     CALL gesvTYPE(rank, 1, tmp_Ainv, rank, ipiv, &
          RiHLd, rank, lib_info)

     DEALLOCATE(ipiv)

  ELSE typeainv1
     ! *** Variant 2: Invert Ainv using SVD

     ALLOCATE(svals(rank))
     ALLOCATE(work(3 * rank))
     ldwork = 3 * rank
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3 * rank)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_estkf_analysis:', debug, &
          '  Compute eigenvalue decomposition of A^-1'

     ! save matrix Ainv
     tmp_Ainv = Ainv

     ! Compute SVD of Ainv
     CALL syevTYPE('v', 'l', rank, tmp_Ainv, rank, svals, work, ldwork, lib_info)

     DEALLOCATE(work)

     ! Compute product RiHLd A
     IF (lib_info==0) THEN
        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_estkf_analysis:', debug, '  eigenvalues', svals

        ALLOCATE(VRiHLd(rank))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank)

        CALL gemvTYPE('t', rank, rank, 1.0, tmp_Ainv, &
             rank, RiHLd, 1, 0.0, VRiHLd, 1)
     
        DO row = 1,rank
           VRiHLd(row) = VRiHLd(row) / svals(row)
        END DO
  
        CALL gemvTYPE('n', rank, rank, 1.0, tmp_Ainv, &
             rank, VRiHLd, 1, 0.0, RiHLd, 1)

        DEALLOCATE(svals, VRiHLd)
     END IF
  END IF typeainv1

  ! *** check whether solve was successful
  IF (lib_info == 0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_estkf_analysis:', debug, '  A(HXT R^-1)^T d', RiHLd

     flag = 0
  ELSE
     WRITE (*, '(/5x, a/)') 'PDAF-ERROR(1): Problem in computation of analysis weights!!!'
     flag = 1
  END IF


! ************************************
! ***      update model state      ***
! ***                              ***
! ***     a   f   f                ***
! ***    x = x + X  Omega RiHLd    ***
! ***                              ***
! ************************************

  check1: IF (flag == 0) THEN

     ! ******************************
     ! *** Compute vector Omega w ***
     ! ******************************

     ALLOCATE(TRiHLd(dim_ens, 1))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)

     CALL PDAF_estkf_OmegaA(rank, 1, RiHLd, TRiHLd)
     DEALLOCATE(RiHLd)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_estkf_analysis:', debug, '  transform vector', TRiHLd

     CALL PDAF_timeit(13, 'old')

     ! *****************************
     ! *** Update state estimate ***
     ! *****************************

     CALL PDAF_timeit(14, 'new')

     CALL gemvTYPE('n', dim_p, dim_ens, 1.0, ens_p, &
          dim_p, TRiHLd, 1, 0.0, state_inc_p, 1)
     DEALLOCATE(TRiHLd)
     
     ! Shift ensemble
     DO col = 1, dim_ens
        DO row = 1, dim_p
           ens_p(row, col) = ens_p(row, col) + state_inc_p(row)
        END DO
     END DO

     IF (incremental == 0) THEN
        ! update state here if incremental updating is not used
        state_p = state_p + state_inc_p
     END IF

     CALL PDAF_timeit(14, 'old')

  ELSE check1

     CALL PDAF_timeit(13, 'old')

  END IF check1
  
  CALL PDAF_timeit(51, 'old')


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(tmp_Ainv)

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_estkf_analysis -- END'
  
END SUBROUTINE PDAF_estkf_analysis_fixed
