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
! !ROUTINE: PDAF_seik_analysis_newT --- Perform SEIK analysis step
!
! !INTERFACE:
SUBROUTINE PDAF_seik_analysis_newT(step, dim_p, dim_obs_p, dim_ens, rank, &
     state_p, Uinv, ens_p, state_inc_p, &
     HL_p, HXbar_p, obs_p, forget, U_prodRinvA, &
     screen, incremental, debug, flag)

! !DESCRIPTION:
! Analysis step of the SEIK filter
! with adaptive forgetting factor.
!
! Variant for domain decomposed states
! and with more efficient implementation of XT.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2003-11 - Lars Nerger - Initial code
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
  USE PDAF_mod_filtermpi, &
       ONLY: mype, MPIerr, COMM_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step         ! Current time step
  INTEGER, INTENT(in) :: dim_p        ! PE-local dimension of model state
  INTEGER, INTENT(in) :: dim_obs_p    ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens      ! Size of ensemble
  INTEGER, INTENT(in) :: rank         ! Rank of initial covariance matrix
  REAL, INTENT(inout) :: state_p(dim_p)           ! PE-local model state
  REAL, INTENT(inout) :: Uinv(rank, rank)         ! Inverse of eigenvalue matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)    ! PE-local state ensemble
  REAL, INTENT(inout) :: state_inc_p(dim_p)       ! PE-local state analysis increment
  REAL, INTENT(inout) :: HL_p(dim_obs_p, dim_ens) ! PE-local observed ensemble (perturbations)
  REAL, INTENT(in)    :: HXbar_p(dim_obs_p)       ! PE-local observed state
  REAL, INTENT(in)    :: obs_p(dim_obs_p)         ! PE-local observation vector
  REAL, INTENT(in)    :: forget       ! Forgetting factor
  INTEGER, INTENT(in) :: screen       ! Verbosity flag
  INTEGER, INTENT(in) :: incremental  ! Control incremental updating
  INTEGER, INTENT(in) :: debug        ! Flag for writing debug output
  INTEGER, INTENT(inout) :: flag      ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA             ! Provide product R^-1 A

! !CALLING SEQUENCE:
! Called by: PDAF_seik_update
! Calls: U_prodRinvA
! Calls: PDAF_set_forget
! Calls: PDAF_seik_matrixT
! Calls: PDAF_seik_Uinv
! Calls: PDAF_seik_TtimesA
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: gemvTYPE (BLAS; dgemv or sgemv dependent on precision)
! Calls: gesvTYPE (LAPACK; dgesv or sgesv dependent on precision)
! Calls: MPI_allreduce (MPI)
!EOP

! *** local variables ***
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  REAL, ALLOCATABLE :: RiHL_p(:,:)    ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: Uinv_p(:,:)    ! Uinv for PE-local domain
  REAL, ALLOCATABLE :: Uinv_inc(:,:)  ! Increment for Uinv
  REAL, ALLOCATABLE :: innov_p(:)     ! PE-local observation innovation
  REAL, ALLOCATABLE :: RiHLd(:)       ! Temporary vector for analysis 
  REAL, ALLOCATABLE :: RiHLd_p(:)     ! PE-local RiHLd
  REAL, ALLOCATABLE :: TRiHLd(:,:)    ! Temporary vector for analysis 
  REAL, ALLOCATABLE :: temp_Uinv(:,:) ! Temporary storage of Uinv
  INTEGER, ALLOCATABLE :: ipiv(:)     ! Vector of pivot indices for GESV
  INTEGER :: gesv_info                ! Control flag for GESV

  
! **********************
! *** INITIALIZATION ***
! **********************

  CALL PDAF_timeit(51, 'new')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_seik_analysis -- START'

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, i7, 3x, a)') &
          'PDAF ', step, 'Assimilating observations - SEIK'
  END IF


! **************************
! *** Compute innovation ***
! ***     d = y - H x    ***
! **************************
  
  haveobsB: IF (dim_obs_p > 0) THEN
     ! The innovation only exists for domains with observations

     CALL PDAF_timeit(10, 'new')

     ALLOCATE(innov_p(dim_obs_p))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p)

     innov_p = obs_p - HXbar_p

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_seik_analysis:', debug, &
             'innovation d(1:min(dim_obs_p,10))', innov_p(1:min(dim_obs_p,10))
        WRITE (*,*) '++ PDAF-debug PDAF_seik_analysis:', debug, &
             'MIN/MAX of innovation', MINVAL(innov_p), MAXVAL(innov_p)
     END IF

     CALL PDAF_timeit(10, 'old')
  END IF haveobsB

  CALL PDAF_timeit(51, 'old')


! ******************************************
! ***   Compute analyzed matrix Uinv     ***
! ***                                    ***
! ***  -1                T        T  -1  ***
! *** U  = forget*N T T + (HL)  R  (HL)  ***
! ***  i                      i  i     i ***
! ******************************************

  CALL PDAF_timeit(11, 'new')

  haveobsA: IF (dim_obs_p > 0) THEN
     ! *** The contribution of observation matrix ist only ***
     ! *** computed for domains with observations          ***

     CALL PDAF_timeit(51, 'new')

     ! Project observed ensemble onto error space: HL = [Hx_1 ... Hx_N] T
     CALL PDAF_seik_matrixT(dim_obs_p, dim_ens, HL_p)

     CALL PDAF_timeit(51, 'old')


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

     CALL PDAF_timeit(31, 'new')
     CALL PDAF_timeit(51, 'new')

    ! *** Initialize Uinv = N T^T T ***
     CALL PDAF_seik_Uinv(rank, Uinv)


     ! *** Finish computation of Uinv  ***
     ! ***   -1          -1    T       ***
     ! ***  U  = forget U  + HL RiHL   ***
     ALLOCATE(Uinv_p(rank, rank))
     ALLOCATE(Uinv_inc(rank, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * rank**2)

     CALL gemmTYPE('t', 'n', rank, rank, dim_obs_p, &
          1.0, HL_p, dim_obs_p, RiHL_p, dim_obs_p, &
          0.0, Uinv_p, rank)

  ELSE haveobsA
     ! *** For domains with dim_obs_p=0 there is no ***
     ! *** direct observation-contribution to Uinv  ***

     CALL PDAF_timeit(31, 'new')
     CALL PDAF_timeit(51, 'new')

     ! Initialize Uinv = N T^T T 
     CALL PDAF_seik_Uinv(rank, Uinv)

     ALLOCATE(Uinv_p(rank, rank))
     ALLOCATE(Uinv_inc(rank, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * rank**2)

     ! No observation-contribution to Uinv from this domain
     Uinv_p = 0.0

  END IF haveobsA

  ! get total sum on all filter PEs
  CALL MPI_allreduce(Uinv_p, Uinv_inc, rank * rank, &
       MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

  Uinv = forget * Uinv + Uinv_inc

  DEALLOCATE(Uinv_p, Uinv_inc)

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_seik_analysis:', debug, '  U^-1', Uinv

  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(31, 'old')
  CALL PDAF_timeit(11, 'old')


! ************************************
! ***      update model state      ***
! ***                              ***
! ***  a   f          T         f  ***
! *** x = x + L U RiHV  (y - H x ) ***
! ***                              ***
! ************************************

  CALL PDAF_timeit(51, 'new')
  CALL PDAF_timeit(12, 'new')

  ! ************************
  ! *** RiHLd = RiHV^T d ***
  ! ************************

  ALLOCATE(RiHLd_p(rank))
  ALLOCATE(RiHLd(rank))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * rank)

  haveobsC: IF (dim_obs_p > 0) THEN
     ! *** RiHLd_p/=0 only with observations ***
    
     ! local products (partial sum)
     CALL gemvTYPE('t', dim_obs_p, rank, 1.0, RiHL_p, &
          dim_obs_p, innov_p, 1, 0.0, RiHLd_p, 1)
     DEALLOCATE(RiHL_p, innov_p)
  ELSE haveobsC
     RiHLd_p = 0.0
  END IF haveobsC

  ! get total sum on all filter PEs
  CALL MPI_allreduce(RiHLd_p, RiHLd, rank, &
       MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

  DEALLOCATE(RiHLd_p)


  ! ****************************************
  ! *** Compute  w = U RiHLd  by solving ***
  ! ***           -1                     ***
  ! ***          U  w = RiHLd            ***
  ! *** for w. We use the LAPACK         ***
  ! *** routine GESV.                    ***
  ! ****************************************

  ALLOCATE(temp_Uinv(rank, rank))
  ALLOCATE(ipiv(rank))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank**2)
  IF (allocflag == 0) CALL PDAF_memcount(3, 'i', rank)

  ! save matrix Uinv
  temp_Uinv = Uinv

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_seik_analysis:', debug, &
       '  Invert U^-1 using solver GESV'

  ! call solver (GESV - LU solver)
  CALL gesvTYPE(rank, 1, temp_Uinv, rank, ipiv, &
       RiHLd, rank, gesv_info)
  DEALLOCATE(temp_Uinv, ipiv)

  ! *** check if solve was successful
  update: IF (gesv_info /= 0) THEN
     WRITE (*, '(/5x, a/)') 'PDAF-ERROR(1): Problem in solve for state analysis !!!'
     flag = 1

     CALL PDAF_timeit(12, 'old')
  ELSE

     ! **************************
     ! *** Compute vector T w ***
     ! **************************

     ALLOCATE(TRiHLd(dim_ens, 1))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)
    
     CALL PDAF_seik_TtimesA(rank, 1, RiHLd, TRiHLd)
     DEALLOCATE(RiHLd)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_seik_analysis:', debug, '  wbar', TRiHLd

     CALL PDAF_timeit(12, 'old')


     ! **************************
     ! *** Update model state ***
     ! ***    a   f           ***
     ! ***   x = x + LT Tw    ***
     ! **************************

     CALL PDAF_timeit(13, 'new')

     CALL gemvTYPE('n', dim_p, dim_ens, 1.0, ens_p, &
          dim_p, TRiHLd, 1, 0.0, state_inc_p, 1)
     DEALLOCATE(TRiHLd)
     
     IF (incremental == 0) THEN
        ! update state here if incremental updating is not used
        state_p = state_p + state_inc_p
     END IF

     CALL PDAF_timeit(13, 'old')
     
  END IF update

  CALL PDAF_timeit(51, 'old')


! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_seik_analysis -- END'

END SUBROUTINE PDAF_seik_analysis_newT
