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
!$Id: PDAF-D_lseik_analysis.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_lseik_analysis --- Perform LSEIK analysis step
!
! !INTERFACE:
SUBROUTINE PDAF_lseik_analysis(domain_p, step, dim_l, dim_obs_f, dim_obs_l, &
     dim_ens, rank, state_l, Uinv_l, ens_l, HX_f, &
     HXbar_f, state_inc_l, forget, U_g2l_obs, U_init_obs_l, &
     U_prodRinvA_l, U_init_obsvar_l, U_init_n_domains_p, screen, incremental, &
     type_forget, flag)

! !DESCRIPTION:
! Analysis step of the LSEIK filter
! with adaptive forgetting factor.
!
! Variant for domain decomposed states.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2005-09 - Lars Nerger - Initial code
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
  REAL, INTENT(inout) :: state_l(dim_l)        ! State on local analysis domain
  REAL, INTENT(inout) :: Uinv_l(rank, rank)    ! Inverse of matrix U
  REAL, INTENT(in) :: ens_l(dim_l, dim_ens)    ! Local state ensemble
  REAL, INTENT(in) :: HX_f(dim_obs_f, dim_ens) ! PE-local full observed state ens.
  REAL, INTENT(in) :: HXbar_f(dim_obs_f)       ! PE-local full observed ens. mean
  REAL, INTENT(in) :: state_inc_l(dim_l)       ! Local state increment
  REAL, INTENT(inout) :: forget      ! Forgetting factor
  INTEGER, INTENT(in) :: screen      ! Verbosity flag
  INTEGER, INTENT(in) :: incremental ! Control incremental updating
  INTEGER, INTENT(in) :: type_forget ! Type of forgetting factor
  INTEGER, INTENT(inout) :: flag     ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_g2l_obs, &   ! Restrict full obs. vector to local analysis domain
       U_init_obs_l, &       ! Init. observation vector on local analysis domain
       U_init_obsvar_l, &    ! Initialize local mean observation error variance
       U_init_n_domains_p, & ! Provide PE-local number of local analysis domains
       U_prodRinvA_l         ! Provide product R^-1 A for local analysis domain

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update
! Calls: U_g2l_obs
! Calls: U_init_obs_l
! Calls: U_prodRinvA_l
! Calls: PDAF_set_forget_local
! Calls: PDAF_seik_matrixT
! Calls: PDAF_seik_Uinv
! Calls: PDAF_seik_TtimesA
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: gemvTYPE (BLAS; dgemv or sgemv dependent on precision)
! Calls: gesvTYPE (LAPACK; dgesv or sgesv dependent on precision)
!EOP
       
! *** local variables ***
  INTEGER :: member                  ! counters
  INTEGER, SAVE :: allocflag = 0     ! Flag whether first time allocation is done
  INTEGER, SAVE :: lastdomain = -1   ! store domain index
  LOGICAL :: screenout = .true.      ! Whether to print information to stdout
  REAL, ALLOCATABLE :: HL_l(:,:)     ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: RiHL_l(:,:)   ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: Uinv_inc(:,:) ! local Uinv
  REAL, ALLOCATABLE :: resid_l(:)    ! observation residual
  REAL, ALLOCATABLE :: obs_l(:)      ! local observation vector
  REAL, ALLOCATABLE :: HXbar_l(:)    ! state projected onto obs. space
  REAL, ALLOCATABLE :: RiHLd_l(:)    ! local RiHLd
  REAL, ALLOCATABLE :: TRiHLd_l(:,:) ! Temporary vector for analysis 
  REAL, ALLOCATABLE :: Uinv_l_tmp(:,:) ! Temporary storage of Uinv
  INTEGER, ALLOCATABLE :: ipiv(:)      ! vector of pivot indices for GESV
  INTEGER :: gesv_info                 ! control flag for GESV
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
        WRITE (*,'(a, 5x,a,i5,a)') &
             'PDAF', '--- Use OpenMP parallelization with ', nthreads, ' threads'
     END IF
#endif
  END IF


  ! ************************
  ! *** Compute residual ***
  ! ***   d = y - H x    ***
  ! ************************

  CALL PDAF_timeit(12, 'new')
  ALLOCATE(resid_l(dim_obs_l))
  ALLOCATE(obs_l(dim_obs_l))
  ALLOCATE(HXbar_l(dim_obs_l))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3 * dim_obs_l)

  ! Restrict mean obs. state onto local observation space
  CALL U_g2l_obs(domain_p, step, dim_obs_f, dim_obs_l, HXbar_f, HXbar_l)

  ! get local observation vector
  CALL U_init_obs_l(domain_p, step, dim_obs_l, obs_l)

  ! get residual as difference of observation and
  ! projected state
  resid_l = obs_l - HXbar_l

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
  CALL PDAF_timeit(30, 'new')

  ! HL = [Hx_1 ... Hx_(r+1)] T
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

  ! HL = [Hx_1 ... Hx_(r+1)] T
  CALL PDAF_seik_matrixT(dim_obs_l, dim_ens, HL_l)

  CALL PDAF_timeit(30, 'old')
  CALL PDAF_timeit(31, 'new')


  ! ***                RiHL = Rinv HL                 ***
  ! *** this is implemented as a subroutine thus that ***
  ! *** Rinv does not need to be allocated explicitly ***
  ALLOCATE(RiHL_l(dim_obs_l, rank))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l * rank)

  CALL U_prodRinvA_l(domain_p, step, dim_obs_l, rank, obs_l, HL_l, RiHL_l)
  DEALLOCATE(obs_l)
 
  ! *** Initialize Uinv = (r+1) T^T T ***
  CALL PDAF_seik_Uinv(rank, Uinv_l)

  ! *** Finish computation of Uinv  ***
  ! ***   -1          -1    T       ***
  ! ***  U  = forget U  + HL RiHL   ***

  ALLOCATE(Uinv_inc(rank, rank))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank**2)
  CALL gemmTYPE('t', 'n', rank, rank, dim_obs_l, &
       1.0, HL_l, dim_obs_l, RiHL_l, dim_obs_l, &
       0.0, Uinv_inc, rank)

  DEALLOCATE(HL_l)

  Uinv_l = forget * Uinv_l + Uinv_inc

  DEALLOCATE(Uinv_inc)

  CALL PDAF_timeit(31, 'old')
  CALL PDAF_timeit(10, 'old')


! ************************************
! ***      update model state      ***
! ***                              ***
! ***  a   f          T         f  ***
! *** x = x + L U RiHV  (y - H x ) ***
! ***                              ***
! ************************************

  CALL PDAF_timeit(13, 'new')
  ! ************************
  ! *** RiHLd = RiHV^T d ***
  ! ************************
  ALLOCATE(RiHLd_l(rank))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * rank)

  CALL gemvTYPE('t', dim_obs_l, rank, 1.0, RiHL_l, &
       dim_obs_l, resid_l, 1, 0.0, RiHLd_l, 1)

  DEALLOCATE(RiHL_l, resid_l)


  ! ****************************************
  ! *** Compute  w = U RiHLd  by solving ***
  ! ***           -1                     ***
  ! ***          U  w = RiHLd            ***
  ! *** for w. We use the LAPACK         ***
  ! *** routine GESV.                    ***
  ! ****************************************

  ALLOCATE(Uinv_l_tmp(rank, rank))
  ALLOCATE(ipiv(rank))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank**2)
  IF (allocflag == 0) CALL PDAF_memcount(3, 'i', rank)

  ! save matrix Uinv
  Uinv_l_tmp = Uinv_l

  ! call solver (GESV - LU solver)
  CALL gesvTYPE(rank, 1, Uinv_l_tmp, rank, ipiv, &
       RiHLd_l, rank, gesv_info)
  DEALLOCATE(Uinv_l_tmp, ipiv)


  ! *** check if solve was successful
  update: IF (gesv_info /= 0) THEN
     WRITE (*, '(/5x, a, i10, 39a/)') &
          'PDAF-ERROR(1): Domain ', domain_p, 'Problem in solve for state analysis !!!'
     flag = 1
  ELSE

     ! **************************
     ! *** Compute vector T w ***
     ! **************************

     ALLOCATE(TRiHLd_l(dim_ens, 1))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)

     CALL PDAF_seik_TtimesA(rank, 1, RiHLd_l, TRiHLd_l)
     DEALLOCATE(RiHLd_l)

     CALL PDAF_timeit(13, 'old')


     ! **************************
     ! *** Update model state ***
     ! ***    a   f           ***
     ! ***   x = x + L RiHLd  ***
     ! **************************

     CALL PDAF_timeit(14, 'new')

     CALL gemvTYPE('n', dim_l, dim_ens, 1.0, ens_l, &
          dim_l, TRiHLd_l, 1, 0.0, state_inc_l, 1)
     DEALLOCATE(TRiHLd_l)

     IF (incremental == 0) THEN
        ! update state only if incremental updating is not used
        state_l = state_l + state_inc_l
     END IF

     CALL PDAF_timeit(14, 'old')
    
  END IF update


! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1

  ! Store domain index
  lastdomain = domain_p

END SUBROUTINE PDAF_lseik_analysis
