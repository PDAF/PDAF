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
!$Id: PDAF-D_lestkf_analysis_fixed.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_lestkf_analysis_fixed --- LESTKF analysis without ens transformation
!
! !INTERFACE:
SUBROUTINE PDAF_lestkf_analysis_fixed(domain_p, step, dim_l, dim_obs_f, dim_obs_l, &
     dim_ens, rank, state_l, Ainv_l, ens_l, HX_f, &
     HXbar_f, state_inc_l, OmegaT_in, forget, U_g2l_obs, &
     U_init_obs_l, U_prodRinvA_l, U_init_obsvar_l, U_init_n_domains_p, &
     screen, incremental, type_forget, type_sqrt, flag)

! !DESCRIPTION:
! Analysis step of the LESTKF with direct update
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
! 2012-03 - Lars Nerger - Initial code
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
#if defined (_OPENMP)
  USE omp_lib, &
       ONLY: omp_get_num_threads, omp_get_thread_num
#endif

  IMPLICIT NONE

! !ARGUMENTS:
! ! Variable naming scheme:
! !   suffix _p: Denotes a full variable on the PE-local domain
! !   suffix _l: Denotes a local variable on the current analysis domain
  INTEGER, INTENT(in) :: domain_p    ! Current local analysis domain
  INTEGER, INTENT(in) :: step        ! Current time step
  INTEGER, INTENT(in) :: dim_l       ! State dimension on local analysis domain
  INTEGER, INTENT(in) :: dim_obs_f   ! PE-local dimension of full observation vector
  INTEGER, INTENT(in) :: dim_obs_l   ! Size of obs. vector on local ana. domain
  INTEGER, INTENT(in) :: dim_ens     ! Size of ensemble 
  INTEGER, INTENT(in) :: rank        ! Rank of initial covariance matrix
  REAL, INTENT(inout) :: state_l(dim_l)        ! on exit: state on local analysis domain
  REAL, INTENT(inout) :: Ainv_l(rank, rank)    ! Inverse of matrix U - temporary use only
  REAL, INTENT(inout) :: ens_l(dim_l, dim_ens) ! Local state ensemble
  REAL, INTENT(in) :: HX_f(dim_obs_f, dim_ens) ! PE-local full observed state ens.
  REAL, INTENT(in) :: HXbar_f(dim_obs_f)       ! PE-local full observed ens. mean
  REAL, INTENT(in) :: state_inc_l(dim_l)       ! Local state increment
  REAL, INTENT(inout) :: OmegaT_in(rank, dim_ens) ! Matrix omega
  REAL, INTENT(inout) :: forget      ! Forgetting factor
  INTEGER, INTENT(in) :: screen      ! Verbosity flag
  INTEGER, INTENT(in) :: incremental ! Control incremental updating
  INTEGER, INTENT(in) :: type_forget ! Type of forgetting factor
  INTEGER, INTENT(in) :: type_sqrt   ! Type of square-root of A
                                     ! (0): symmetric sqrt; (1): Cholesky decomposition
  INTEGER, INTENT(inout) :: flag     ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_g2l_obs, &   ! Restrict full obs. vector to local analysis domain
       U_init_obs_l, &       ! Init. observation vector on local analysis domain
       U_init_obsvar_l, &    ! Initialize local mean observation error variance
       U_init_n_domains_p, & ! Provide PE-local number of local analysis domains
       U_prodRinvA_l         ! Provide product R^-1 A for local analysis domain

! !CALLING SEQUENCE:
! Called by: PDAF_lestkf_update
! Calls: U_g2l_obs
! Calls: U_init_obs_l
! Calls: U_prodRinvA_l
! Calls: PDAF_set_forget_local
! Calls: PDAF_estkf_AOmega
! Calls: PDAF_estkf_OmegaA
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: gemvTYPE (BLAS; dgemv or sgemv dependent on precision)
! Calls: gesvTYPE (LAPACK; dgesv or sgesv dependent on precision)
! Calls: syevTYPE (LAPACK; dsyev or ssyev dependent on precision)
!EOP
       
! *** local variables ***
  INTEGER :: i, member, col, row       ! Counters
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  INTEGER :: lib_info                  ! Status flag for LAPACK calls
  INTEGER :: ldwork                    ! Size of work array for SYEVTYPE
  INTEGER, SAVE :: lastdomain = -1     ! store domain index
  LOGICAL :: screenout = .true.        ! Whether to print information to stdout
  REAL, ALLOCATABLE :: HL_l(:,:)       ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: RiHL_l(:,:)     ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: resid_l(:)      ! observation residual
  REAL, ALLOCATABLE :: obs_l(:)        ! local observation vector
  REAL, ALLOCATABLE :: HXbar_l(:)      ! state projected onto obs. space
  REAL, ALLOCATABLE :: RiHLd_l(:)      ! local RiHLd
  REAL, ALLOCATABLE :: VRiHLd_l(:)     ! Temporary vector for analysis
  REAL, ALLOCATABLE :: tmp_Ainv_l(:,:) ! Temporary storage of Ainv
  REAL, ALLOCATABLE :: TRiHLd_l(:,:) ! Temporary vector for analysis 
  REAL, ALLOCATABLE :: svals(:)      ! Singular values of Ainv
  REAL, ALLOCATABLE :: work(:)       ! Work array for syevTYPE
  INTEGER, ALLOCATABLE :: ipiv(:)    ! vector of pivot indices for GESVTYPE
  INTEGER, SAVE :: mythread, nthreads  ! Thread variables for OpenMP

!$OMP THREADPRIVATE(mythread, nthreads, lastdomain, allocflag, screenout)


! *******************
! *** Preparation ***
! *******************

#if defined (_OPENMP)
  nthreads = omp_get_num_threads()
  mythread = omp_get_thread_num()
#else
  nthreads = 1
  mythread = 0
#endif

  ! Control screen output
  IF (lastdomain<domain_p .AND. lastdomain>-1) THEN
     screenout = .false.
  ELSE
     screenout = .true.

     ! In case of OpenMP, let only thread 0 write output to the screen
     IF (mythread>0) screenout = .false.

     ! Output, only in case of OpenMP parallelization
#if defined (_OPENMP)
     IF (screenout .AND. screen > 0) THEN
        WRITE (*,'(a, 5x,a,i5,a)') 'PDAF', '--- Use OpenMP parallelization with ', nthreads, ' threads'
     END IF
#endif
  END IF


! ************************
! *** Compute residual ***
! ***   d = y - H x    ***
! ************************

  CALL PDAF_timeit(12, 'new')

  haveobsB: IF (dim_obs_l > 0) THEN
     ! *** The residual only exists for domains with observations ***

     ALLOCATE(resid_l(dim_obs_l))
     ALLOCATE(obs_l(dim_obs_l))
     ALLOCATE(HXbar_l(dim_obs_l))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3 * dim_obs_l)

     ! Restrict mean obs. state onto local observation space
     CALL U_g2l_obs(domain_p, step, dim_obs_f, dim_obs_l, HXbar_f, HXbar_l)

     ! get local observation vector
     CALL U_init_obs_l(domain_p, step, dim_obs_l, obs_l)

     ! Get residual as difference of observation and observed state
     resid_l = obs_l - HXbar_l

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

  haveobsA: IF (dim_obs_l > 0) THEN
     ! *** The contribution of observation matrix ist only ***
     ! *** computed for domains with observations          ***

     CALL PDAF_timeit(30, 'new')

     ! *** Compute HL = [Hx_1 ... Hx_N] Omega
     ALLOCATE(HL_l(dim_obs_l, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l * dim_ens)

     ENS: DO member = 1, dim_ens
        ! [Hx_1 ... Hx_(r+1)] for local analysis domain
        CALL U_g2l_obs(domain_p, step, dim_obs_f, dim_obs_l, HX_f(:, member), &
             HL_l(:, member))
     END DO ENS

     ! *** Set the value of the forgetting factor  ***
     ! *** Inserted here, because HL_l is required ***
     IF (type_forget == 2) THEN
        CALL PDAF_set_forget_local(domain_p, step, dim_obs_l, dim_ens, HL_l, &
             HXbar_l, resid_l, obs_l, U_init_n_domains_p, U_init_obsvar_l, &
             forget)
     ENDIF
     DEALLOCATE(HXbar_l)

     ! Complete HL = [Hx_1 ... Hx_N] Omega
     CALL PDAF_estkf_AOmega(dim_obs_l, dim_ens, HL_l)

     CALL PDAF_timeit(30, 'old')
     CALL PDAF_timeit(31, 'new')


     ! ***                RiHL = Rinv HL                 ***
     ! *** this is implemented as a subroutine thus that ***
     ! *** Rinv does not need to be allocated explicitly ***
     ALLOCATE(RiHL_l(dim_obs_l, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l * rank)

     CALL U_prodRinvA_l(domain_p, step, dim_obs_l, rank, obs_l, HL_l, RiHL_l)
     DEALLOCATE(obs_l)
 
     ! *** Initialize Ainv = (N-1) I ***
     Ainv_l = 0.0
     DO i = 1, rank
        Ainv_l(i,i) = REAL(dim_ens - 1)
     END DO

     ! ***             T        ***
     ! ***  Compute  HL  RiHL   ***
     ALLOCATE(tmp_Ainv_l(rank, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank**2)

     CALL gemmTYPE('t', 'n', rank, rank, dim_obs_l, &
          1.0, HL_l, dim_obs_l, RiHL_l, dim_obs_l, &
          0.0, tmp_Ainv_l, rank)

     DEALLOCATE(HL_l)

  ELSE haveobsA
     ! *** For domains with dim_obs_l=0 there is no ***
     ! *** direct observation-contribution to Ainv  ***

     CALL PDAF_timeit(31, 'new')

     ! *** Initialize Ainv = (N-1) I ***
     Ainv_l = 0.0
     DO i = 1, rank
        Ainv_l(i,i) = REAL(dim_ens - 1)
     END DO

     ! No observation-contribution to Ainv from this domain
     ALLOCATE(tmp_Ainv_l(rank, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank**2)

     tmp_Ainv_l = 0.0

  END IF haveobsA

  ! *** Complete computation of Ainv  ***
  ! ***   -1                T         ***
  ! ***  A  = forget I  + HL RiHL     ***
  Ainv_l = forget * Ainv_l + tmp_Ainv_l

  CALL PDAF_timeit(31, 'old')
  CALL PDAF_timeit(10, 'old')


! ***********************************************
! *** Compute weight for model state update   ***
! ***                                         ***
! ***              T                    f     ***
! ***    w = A RiHL d  with d = (y - H x )    ***
! ***********************************************

  CALL PDAF_timeit(13, 'new')

  ! *** Compute RiHLd = RiHL^T d ***
  ALLOCATE(RiHLd_l(rank))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank)

  haveobsC: IF (dim_obs_l > 0) THEN
     ! *** RiHLd_p/=0 only with observations ***

     CALL gemvTYPE('t', dim_obs_l, rank, 1.0, RiHL_l, &
          dim_obs_l, resid_l, 1, 0.0, RiHLd_l, 1)

     DEALLOCATE(RiHL_l, resid_l)

  ELSE haveobsC

     RiHLd_l = 0.0

  END IF haveobsC


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

  typeainv1: IF (type_sqrt==1) THEN
     ! *** Variant 1: Solve Ainv w= RiHLd for w

     ALLOCATE(ipiv(rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'i', rank)

     ! save matrix Ainv
     tmp_Ainv_l = Ainv_l

     ! call solver (gesvTYPE - LU solver)
     CALL gesvTYPE(rank, 1, tmp_Ainv_l, rank, ipiv, &
          RiHLd_l, rank, lib_info)
     DEALLOCATE(ipiv)

  ELSE typeainv1
     ! *** Variant 2: Invert Ainv using SVD

     ALLOCATE(svals(rank))
     ALLOCATE(work(3 * rank))
     ldwork = 3 * rank
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 4 * rank)

     ! save matrix Ainv
     tmp_Ainv_l = Ainv_l

     ! Compute SVD of Ainv
     CALL syevTYPE('v', 'l', rank, tmp_Ainv_l, rank, svals, work, ldwork, lib_info)

     DEALLOCATE(work)

     ! Compute product A RiHLd
     IF (lib_info==0) THEN
        ALLOCATE(VRiHLd_l(rank))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank)

        CALL gemvTYPE('t', rank, rank, 1.0, tmp_Ainv_l, &
             rank, RiHLd_l, 1, 0.0, VRiHLd_l, 1)
     
        DO row = 1,rank
           VRiHLd_l(row) = VRiHLd_l(row) / svals(row)
        END DO
  
        CALL gemvTYPE('n', rank, rank, 1.0, tmp_Ainv_l, &
             rank, VRiHLd_l, 1, 0.0, RiHLd_l, 1)

        DEALLOCATE(svals, VRiHLd_l)
     END IF
  END IF typeainv1

  ! *** check if SVD was successful
  IF (lib_info == 0) THEN
     flag = 0
  ELSE
     WRITE (*, '(/5x, a, i10, a/)') &
          'PDAF-ERROR(1): Domain ', domain_p, ' Problem in computation of analysis weights !!!'
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

     ALLOCATE(TRiHLd_l(dim_ens, 1))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)

     CALL PDAF_estkf_OmegaA(rank, 1, RiHLd_l, TRiHLd_l)
     DEALLOCATE(RiHLd_l)

     CALL PDAF_timeit(13, 'old')

     ! *****************************
     ! *** Update state estimate ***
     ! *****************************

     CALL PDAF_timeit(14, 'new')

     CALL gemvTYPE('n', dim_l, dim_ens, 1.0, ens_l, &
          dim_l, TRiHLd_l, 1, 0.0, state_inc_l, 1)
     DEALLOCATE(TRiHLd_l)
     
     ! Shift ensemble
     DO col = 1, dim_ens
        DO row = 1, dim_l
           ens_l(row, col) = ens_l(row, col) + state_inc_l(row)
        END DO
     END DO
     
     IF (incremental == 0) THEN
        ! update state here if incremental updating is not used
        state_l = state_l + state_inc_l
     END IF

     CALL PDAF_timeit(14, 'old')

  ELSE check1

     CALL PDAF_timeit(13, 'old')

  END IF check1

! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(tmp_Ainv_l)

  IF (allocflag == 0) allocflag = 1

  ! Store domain index
  lastdomain = domain_p

END SUBROUTINE PDAF_lestkf_analysis_fixed
