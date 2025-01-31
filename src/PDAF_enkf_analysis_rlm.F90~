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
! !ROUTINE: PDAF_enkf_analysis_rlm --- Perform EnKF analysis step
!
! !INTERFACE:
SUBROUTINE PDAF_enkf_analysis_rlm(step, dim_p, dim_obs_p, dim_ens, rank_ana, &
     state_p, ens_p, HXB, forget, U_init_dim_obs, U_obs_op, &
     U_add_obs_err, U_init_obs, U_init_obs_covar, screen, flag)

! !DESCRIPTION:
! Analysis step of ensemble Kalman filter with 
! representer-type formulation.  This version is 
! for large observation dimension  in which HP 
! is not explicitely computed.  It is optimal if 
! the number of observations is larger than half 
! of the ensemble size.
! The final ensemble update uses a block
! formulation to reduce memory requirements.
!
! Variant for domain decomposition.
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
  USE PDAF_mod_filter, &
       ONLY: obs_member, debug
  USE PDAF_mod_filtermpi, &
       ONLY: mype, npes_filter, MPIerr, COMM_filter
  USE PDAFomi, &
       ONLY: omi_n_obstypes => n_obstypes, omi_omit_obs => omit_obs, &
       PDAFomi_gather_obsdims
  USE PDAF_analysis_utils, &
       ONLY: PDAF_omit_obs_omi

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step      ! Current time step
  INTEGER, INTENT(in)  :: dim_p     ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p ! PE-local dimension of observation vector
  INTEGER, INTENT(in)  :: dim_ens   ! Size of state ensemble
  INTEGER, INTENT(in)  :: rank_ana  ! Rank to be considered for inversion of HPH
  REAL, INTENT(inout)  :: state_p(dim_p)        ! PE-local ensemble mean state
  REAL, INTENT(inout)  :: ens_p(dim_p, dim_ens) ! PE-local state ensemble
  REAL, INTENT(inout)  :: HXB(dim_ens, dim_ens) ! Ensemble tranformation matrix
  REAL, INTENT(in)     :: forget    ! Forgetting factor
  INTEGER, INTENT(in)  :: screen    ! Verbosity flag
  INTEGER, INTENT(inout) :: flag    ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, & ! Initialize dimension of observation vector
       U_obs_op, &              ! Observation operator
       U_init_obs, &            ! Initialize observation vector
       U_init_obs_covar, &      ! Initialize observation error covariance matrix
       U_add_obs_err            ! Add observation error covariance matrix

! !CALLING SEQUENCE:
! Called by: PDAF_enkf_update
! Calls: U_init_dim_obs
! Calls: U_obs_op
! Calls: U_add_obs_err
! Calls: PDAF_enkf_gather_resid
! Calls: PDAF_enkf_obs_ensemble
! Calls: PDAF_timeit
! Calls: PDAF_memcount
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: gesvTYPE (LAPACK; dgesv or sgesv dependent on precision)
! Calls: syevxTYPE (LAPACK; dsyevx or ssyevx dependent on precision)
! Calls: MPI_allreduce (MPI)
!EOP

! *** local variables ***
  INTEGER :: i, j, member ! counters
  INTEGER :: dim_obs      ! global dimension of observation vector
  REAL :: invdim_ens      ! inverse of ensemble size
  REAL :: invdim_ensm1    ! inverse of ensemble size minus 1
  REAL :: sqrtinvforget   ! square root of inverse forgetting factor
  INTEGER, SAVE :: allocflag = 0    ! Flag whether first time allocation is done
  INTEGER, SAVE :: allocflag_b = 0  ! Flag whether first time allocation is done
  REAL, ALLOCATABLE :: HPH(:,:)        ! Temporary matrix for analysis
  REAL, ALLOCATABLE :: XminMean_b(:,:) ! Temporary matrix for analysis
  REAL, ALLOCATABLE :: HX(:,:)         ! H(ensstate)-H(meanstate)
  REAL, ALLOCATABLE :: resid(:,:)      ! ensemble of global residuals
  REAL, ALLOCATABLE :: resid_p(:,:)    ! ensemble of local residuals
  REAL, ALLOCATABLE :: m_state_p(:)    ! PE-local observed state vector
  REAL, ALLOCATABLE :: HXmean_p(:)     ! Temporary matrix for analysis
  INTEGER, ALLOCATABLE :: ipiv(:)      ! vector of pivot indices
  REAL, ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
  INTEGER :: sgesv_info                ! output flag of SGESV

  ! *** Variables for variant using pseudo inverse with eigendecompositon
  REAL, ALLOCATABLE :: eval(:)         ! vector of eigenvalues
  REAL, ALLOCATABLE :: rwork(:)        ! workarray for eigenproblem
  REAL, ALLOCATABLE :: evec(:,:)       ! matrix of eigenvectors
  REAL, ALLOCATABLE :: evec_temp(:,:)  ! matrix of eigenvectors
  REAL, ALLOCATABLE :: repres(:,:)     ! matrix of representer vectors
  INTEGER :: syev_info     ! output flag of eigenproblem routine
  REAL    :: VL, VU         ! temporary variables for SYEVX (never really used)
  INTEGER :: Ilower, Iupper ! variables defining the interval of eigenvalues
  REAL    :: abstol         ! error tolerance for eigenvalue problem
  INTEGER :: nEOF           ! number of EOFs as computed by SYEVX
  INTEGER, ALLOCATABLE :: iwork(:)     ! workarray for SYEVX
  INTEGER, ALLOCATABLE :: ifail(:)     ! workarray for SYEVX
  REAL, EXTERNAL :: DLAMCH             ! function to specify tolerance of SYEVX
  REAL    :: eval_inv                  ! inverse of an eigenvalue
  INTEGER :: maxblksize, blklower, blkupper ! Variables for block formulation


! **********************
! *** INITIALIZATION ***
! **********************

  CALL PDAF_timeit(51, 'new')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_enkf_analysis -- START'

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, i7, 3x, a)') &
          'PDAF ', step, 'Assimilating observations - EnKF large-dim_obs version using transform matrix'
     IF (rank_ana > 0) THEN
        WRITE (*, '(a, 5x, a, i5)') &
             'PDAF', '--- use pseudo inverse of HPH, rank= ', rank_ana
     ELSE
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- use HPH directly'
     END IF
  END IF

  ! init numbers
  invdim_ens = 1.0 / REAL(dim_ens)
  invdim_ensm1 = 1.0 / (REAL(dim_ens - 1))
  sqrtinvforget = SQRT(1.0 / forget)

  CALL PDAF_timeit(51, 'old')


! *********************************
! *** Get observation dimension ***
! *********************************

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_enkf_analysis:', debug, '  dim_p', dim_p
     WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_enkf_analysis -- call init_dim_obs'
  END IF

  CALL PDAF_timeit(43, 'new')
  CALL U_init_dim_obs(step, dim_obs_p)
  CALL PDAF_timeit(43, 'old')

  CALL PDAF_timeit(51, 'new')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_enkf_analysis:', debug, '  dim_obs_p', dim_obs_p

  IF (screen > 2) THEN
     WRITE (*, '(a, 5x, a13, 1x, i6, 1x, a, i10)') &
          'PDAF', '--- PE-domain', mype, 'dimension of observation vector', dim_obs_p
  END IF

  ! Get global dimension of observation vector
  IF (npes_filter>1) THEN
     CALL MPI_allreduce(dim_obs_p, dim_obs, 1, MPI_INTEGER, MPI_SUM, &
          COMM_filter, MPIerr)
  ELSE
     ! This is a work around for working with nullmpi.F90
     dim_obs = dim_obs_p
  END IF


! **********************************
! *** Compute representer vector ***
! ***                            ***
! *** We compute the ensemble of ***
! *** representer vectors b by   ***
! *** solving                    ***
! ***        T                   ***
! *** (H P H  + R) b  = y - H x  ***
! **********************************

  ! *******************************************************
  ! *** Compute mean forecasted state                   ***
  ! *** for normal EnKF using ensemble mean as forecast ***
  ! *******************************************************
  CALL PDAF_timeit(11, 'new')

  state_p = 0.0
  DO member = 1, dim_ens
     DO i = 1, dim_p
        state_p(i) = state_p(i) + invdim_ens * ens_p(i, member)
     END DO
  END DO

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_enkf_analysis:', debug, 'forecast ensemble mean (1:min(dim_p,6)):', &
          state_p(1:min(dim_p,6))
  END IF

  CALL PDAF_timeit(11, 'old')
  CALL PDAF_timeit(10, 'new')


  ! ****************************************************************
  ! *** Omit observations if innovation is too large             ***
  ! *** This step also initializes obs_p, whic his not used here ***
  ! ****************************************************************

  IF (omi_omit_obs) THEN
     ALLOCATE(obs_p(dim_obs_p))

     CALL PDAF_omit_obs_omi(dim_p, dim_obs_p, dim_ens, state_p, ens_p, &
          obs_p, U_init_obs, U_obs_op, 0, screen)

     DEALLOCATE(obs_p)
  END IF


  ! **********************************************
  ! *** We directly compute the matrices       ***
  ! ***                                T       ***
  ! ***   HP = H P     and  HPH = H P H        ***
  ! *** as ensemble means by just projecting   ***
  ! *** the state ensemble onto observation    ***
  ! *** space. The covariance matrix is not    ***
  ! *** explicitly computed.                   ***
  ! **********************************************
  ALLOCATE(resid_p(dim_obs_p, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * dim_ens)


  ! ***                             T ***
  ! *** get HP = H P and HPH = H P H  ***
  ! *** as ensemble means             ***

  ENSa: DO member = 1, dim_ens
     ! spread out state ensemble according to forgetting factor
     IF (forget /= 1.0) THEN
        ens_p(:, member) = state_p(:) &
             + (ens_p(:, member) - state_p(:)) * sqrtinvforget
     END IF
  END DO ENSa

  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(44, 'new')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_enkf_analysis -- call obs_op', dim_ens, 'times'

  ENSb: DO member = 1, dim_ens
     ! Store member index to make it accessible with PDAF_get_obsmemberid
     obs_member = member

     ! project state ensemble onto observation space
     ! stored in RESID_P for compactness
     CALL U_obs_op(step, dim_p, dim_obs_p, ens_p(:, member), resid_p(:, member))
  END DO ENSb

  CALL PDAF_timeit(44, 'old')
  CALL PDAF_timeit(51, 'new')

  ! compute mean of ensemble projected on obseration space
  ALLOCATE(HXmean_p(dim_obs_p))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p)
    
  CALL PDAF_timeit(33, 'new')
  HXmean_p(:) = 0.0
  DO member = 1, dim_ens
     DO i = 1, dim_obs_p
        HXmean_p(i) = HXmean_p(i) + invdim_ens * resid_p(i, member)
     END DO
  END DO
  CALL PDAF_timeit(33, 'old')
     
  CALL PDAF_timeit(30, 'new')
  ! compute difference between ensemble state projected on obs space
  ! and the corresponding ensemble mean
  DO member = 1, dim_ens
     resid_p(:, member) = resid_p(:, member) - HXmean_p(:)
  END DO
  DEALLOCATE(HXmean_p)

  IF (debug>0) THEN
     DO i = 1, dim_ens
        WRITE (*,*) '++ PDAF-debug PDAF_enkf_analysis:', debug, 'process-local observed ensemble pert, member', i, &
             ' values (1:min(dim_obs_p,6)):', resid_p(1:min(dim_obs_p,6),i)
     END DO
  END IF

  ! Allgather residual
  ALLOCATE(HX(dim_obs, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs * dim_ens)

  CALL PDAF_enkf_gather_resid(dim_obs, dim_obs_p, dim_ens, resid_p, HX)
  CALL PDAF_timeit(30, 'old')

  ! Finish computation of HPH
  ALLOCATE(HPH(dim_obs, dim_obs))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs * dim_obs)

  CALL PDAF_timeit(32, 'new')
  CALL gemmTYPE ('n', 't', dim_obs, dim_obs, dim_ens, &
       invdim_ensm1, HX, dim_obs, HX, dim_obs, &
       0.0, HPH, dim_obs)
  CALL PDAF_timeit(32, 'old')

  CALL PDAF_timeit(51, 'old')

  ! For OMI: Gather global observation dimensions
  IF (omi_n_obstypes > 0) CALL PDAFomi_gather_obsdims()

  ! *** Add observation error covariance ***
  ! ***       HPH^T = (HPH + R)          ***
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_enkf_analysis -- call add_obs_err'

  CALL PDAF_timeit(46, 'new')
  CALL U_add_obs_err(step, dim_obs, HPH)
  CALL PDAF_timeit(46, 'old')

  CALL PDAF_timeit(10, 'old')


! *****************************************
! *** generate ensemble of observations ***
! *****************************************

  CALL PDAF_timeit(15, 'new')
  ! observation ensemble is initialized into the residual matrix
  CALL PDAF_enkf_obs_ensemble(step, dim_obs_p, dim_obs, dim_ens, resid_p, &
       U_init_obs, U_init_obs_covar, screen, flag)
  CALL PDAF_timeit(15, 'old')


! ************************
! *** Compute residual ***
! ***   d = y - H x    ***
! ************************
  ALLOCATE(m_state_p(dim_obs_p))
  ALLOCATE(resid(dim_obs, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p + dim_obs * dim_ens)

  CALL PDAF_timeit(12, 'new')
  ! *** Project state onto observation space and    ***
  ! *** compute observation residual (innovation) d ***
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_enkf_analysis -- call obs_op', dim_ens, 'times'

  DO member = 1, dim_ens
     ! Store member index to make it accessible with PDAF_get_obsmemberid
     obs_member = member

     ! Project state onto observation space
     CALL PDAF_timeit(44, 'new')
     CALL U_obs_op(step, dim_p, dim_obs_p, ens_p(:, member), m_state_p)
     CALL PDAF_timeit(44, 'old')

     ! get residual as difference of observation and
     ! projected state
     CALL PDAF_timeit(51, 'new')
     resid_p(:, member) = resid_p(:, member) - m_state_p(:)
     CALL PDAF_timeit(51, 'old')
  END DO
  DEALLOCATE(m_state_p)

  IF (debug>0) THEN
     DO i = 1, dim_ens
        WRITE (*,*) '++ PDAF-debug PDAF_enkf_analysis:', debug, 'process-local innovation member', i, &
             ' values (1:min(dim_obs_p,6)):', resid_p(1:min(dim_obs_p,6),i)
     END DO
  END IF

  ! Allgather residual
  CALL PDAF_timeit(51, 'new')
  CALL PDAF_enkf_gather_resid(dim_obs, dim_obs_p, dim_ens, resid_p, resid)
  CALL PDAF_timeit(51, 'old')

  DEALLOCATE(resid_p)

  CALL PDAF_timeit(12, 'old')
  CALL PDAF_timeit(14, 'new')
  CALL PDAF_timeit(51, 'new')


  whichupdate: IF (rank_ana > 0) THEN
! **************************************************
! *** Update using pseudo inverse of HPH         ***
! *** by performing incomplete eigendecompostion ***
! *** and using Moore-Penrose inverse of this    ***
! *** matrix                                     ***
! **************************************************

     ! *** Initialization ***
     ALLOCATE(repres(dim_obs, dim_ens))
     ALLOCATE(eval(rank_ana))
     ALLOCATE(evec(dim_obs, rank_ana))
     ALLOCATE(evec_temp(dim_obs, rank_ana))
     ALLOCATE(rwork(8 * dim_obs))
     ALLOCATE(iwork(5 * dim_obs))
     ALLOCATE(ifail(dim_obs))

     IF (allocflag_b == 0) THEN
        ! count allocated memory
        CALL PDAF_memcount(3, 'r', dim_obs * dim_ens + rank_ana + &
             2 * dim_obs * rank_ana + 8 * dim_obs)
        CALL PDAF_memcount(3, 'i', 6 * dim_obs)
        allocflag_b = 1
     END IF

     CALL PDAF_timeit(13, 'new')
     CALL PDAF_timeit(36, 'new')

     ! **************************************
     ! *** compute pseudo inverse of HPH  ***
     ! *** using Moore-Penrose inverse    ***
     ! *** o rank reduced matrix          ***
     ! **************************************

     Iupper = dim_obs
     Ilower = dim_obs - rank_ana + 1
     abstol = 2 * DLAMCH('S')

     ! *** Decompose HPH = eigenvec ev eigenvec^T by   ***
     ! *** computing the RANK_ANA largest eigenvalues  ***
     ! *** and the corresponding eigenvectors          ***
     ! *** We use the LAPACK routine SYEVX             ***
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_enkf_analysis:', debug, &
          '  Compute eigenvalue decomposition of HPH^T'

     CALL syevxTYPE('v', 'i', 'u', dim_obs, HPH, &
          dim_obs, VL, VU, Ilower, Iupper, &
          abstol, nEOF, eval, evec, dim_obs, &
          rwork, 8 * dim_obs, iwork, ifail, syev_info)

     ! check if eigendecomposition was successful
     EVPok: IF (syev_info == 0) THEN
        ! Eigendecomposition OK, continue

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_enkf_resample:', debug, '  eigenvalues', eval

        ! *** store V ***
        evec_temp = evec

        ! *** compute  V diag(ev^(-1)) ***
        DO j = 1, rank_ana
           eval_inv = 1.0 / eval(j)
           DO i = 1, dim_obs
              evec(i, j) = eval_inv * evec(i, j)
           END DO
        END DO
      
        ! *** compute HPH^(-1) ~ V evinv V^T ***
        ! *** HPH^(-1) is stored in HPH      ***
        CALL gemmTYPE('n', 't', dim_obs, dim_obs, rank_ana, &
             1.0, evec, dim_obs, evec_temp, dim_obs, &
             0.0, HPH, dim_obs)

        DEALLOCATE(eval, evec, evec_temp, rwork, iwork, ifail)

        CALL PDAF_timeit(36, 'old')


        ! ****************************************
        ! *** Compute ensemble of representer  ***
        ! *** vectors b as the product         ***
        ! ***           b = invHPH d           ***
        ! ****************************************

        CALL PDAF_timeit(37, 'new')
        CALL gemmTYPE('n', 'n', dim_obs, dim_ens, dim_obs, &
             1.0, HPH, dim_obs, resid, dim_obs, &
             0.0, repres, dim_obs)
        CALL PDAF_timeit(37, 'old')

        CALL PDAF_timeit(13, 'old')


        ! **************************************
        ! *** Update model state ensemble    ***
        ! ***    a   f         f     _    T  ***
        ! ***   x = x + K d = x + (X-X) HX B ***
        ! **************************************

        CALL PDAF_timeit(16, 'new')

        ! *** HXB = HX^T B
        CALL PDAF_timeit(31, 'new')
        CALL gemmTYPE('t', 'n', dim_ens, dim_ens, dim_obs, &
             1.0, HX, dim_obs, repres, dim_obs, &
             0.0, HXB, dim_ens)
        CALL PDAF_timeit(31, 'old')

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_enkf_analysis:', debug, '  transform HXB', HXB

        ! *** Blocking loop for ensemble update ***

        ! Initializations
        maxblksize = 200
        IF (mype == 0) &
             WRITE (*, '(a, 5x, a, i5)') 'PDAF', '--- use blocking with size ', maxblksize

        ! *** XminMean
        ALLOCATE(XminMean_b(maxblksize, dim_ens))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', maxblksize * dim_ens)

        blocking1: DO blklower = 1, dim_p, maxblksize
        
           blkupper = MIN(blklower + maxblksize - 1, dim_p)

           CALL PDAF_timeit(35, 'new')
           ENSc: DO member = 1, dim_ens
              ! initialize XminMean
              XminMean_b(1 : blkupper - blklower + 1, member) = &
                   ens_p(blklower : blkupper, member) - &
                   state_p(blklower : blkupper)
           END DO ENSc
           CALL PDAF_timeit(35, 'old')

           ! *** Update ensemble
           CALL PDAF_timeit(38, 'new')
           CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
                invdim_ensm1, XminMean_b, maxblksize, HXB(1, 1), dim_ens, &
                1.0, ens_p(blklower, 1), dim_p)
           CALL PDAF_timeit(38, 'old')
           
        END DO blocking1

        DEALLOCATE(XminMean_b)

        CALL PDAF_timeit(16, 'old')

     ELSE
        ! Error in the EVP
        CALL PDAF_timeit(32, 'old')
        CALL PDAF_timeit(13, 'old')
     END IF EVPok

     ! *** Clean up ***
     DEALLOCATE(repres)


  ELSE whichupdate
! *******************************************
! *** Update using matrix HPH directly to ***      
! *** compute representer amplitudes b by ***
! *** solving HPH b = d for b.            ***
! *******************************************

     ! ****************************************
     ! *** Compute ensemble of representer  ***
     ! *** vectors b by solving             ***
     ! ***              HPH b = d           ***
     ! *** We use the LAPACK routine GESV   ***
     ! ****************************************
     ALLOCATE(ipiv(dim_obs))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'i', dim_obs)

     CALL PDAF_timeit(13, 'new')
     CALL gesvTYPE(dim_obs, dim_ens, HPH, dim_obs, ipiv, &
          resid, dim_obs, sgesv_info)
     CALL PDAF_timeit(13, 'old')


     ! *** check if solve was successful
     update: IF (sgesv_info /= 0) THEN
        WRITE (*, '(/5x, a/)') 'PDAF-ERROR(2): Problem in solve for Kalman gain !!!'
        flag = 2
     ELSE

     ! **************************************
     ! *** Update model state ensemble    ***
     ! ***    a   f         f     _    T  ***
     ! ***   x = x + K d = x + (X-X) HX B ***
     ! **************************************

        
        CALL PDAF_timeit(16, 'new')

        ! *** HXB = HX^T B
        CALL PDAF_timeit(31, 'new')
        CALL gemmTYPE('t', 'n', dim_ens, dim_ens, dim_obs, &
             1.0, HX, dim_obs, resid, dim_obs, &
             0.0, HXB, dim_ens)
        CALL PDAF_timeit(31, 'old')

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_enkf_analysis:', debug, '  transform', invdim_ensm1*HXB

        ! *** Blocking loop for ensemble update ***

        ! Initializations
        maxblksize = 200
        IF (mype == 0) &
             WRITE (*, '(a, 5x, a, i5)') 'PDAF', '--- use blocking with size ', maxblksize

        ! *** XminMean
        ALLOCATE(XMinMean_b(maxblksize, dim_ens))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', maxblksize * dim_ens)

        blocking2: DO blklower = 1, dim_p, maxblksize
      
           blkupper = MIN(blklower + maxblksize - 1, dim_p)

           CALL PDAF_timeit(35, 'new')
           DO member = 1, dim_ens
              ! initialize XminMean
              XminMean_b(1 : blkupper - blklower + 1, member) = &
                   ens_p(blklower : blkupper, member) - &
                   state_p(blklower : blkupper)
           END DO
           CALL PDAF_timeit(35, 'old')

           ! *** Update ensemble
           CALL PDAF_timeit(38, 'new')
           CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
                invdim_ensm1, XminMean_b, maxblksize, HXB(1, 1), dim_ens, &
                1.0, ens_p(blklower, 1), dim_p)
           CALL PDAF_timeit(38, 'old')
           
        END DO blocking2

        DEALLOCATE(XminMean_b)

        CALL PDAF_timeit(16, 'old')

     END IF update
     DEALLOCATE(ipiv)

  END IF whichupdate

  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(14, 'old')


! ********************
! *** Finishing up ***
! ********************

  ! Clean up
  DEALLOCATE(HX)
  DEALLOCATE(resid)
  DEALLOCATE(HPH)

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_enkf_analysis -- END'

END SUBROUTINE PDAF_enkf_analysis_rlm
