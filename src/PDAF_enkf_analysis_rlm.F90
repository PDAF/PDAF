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
!> Perform EnKF analysis step
!!
!! Analysis step of ensemble Kalman filter with 
!! representer-type formulation.  This version is 
!! for large observation dimension  in which HP 
!! is not explicitely computed.  It is optimal if 
!! the number of observations is larger than half 
!! of the ensemble size.
!! The final ensemble update uses a block
!! formulation to reduce memory requirements.
!!
!! Variant for domain decomposition.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2003-11 - Lars Nerger - Initial code
!! * Later revisions - see svn log
!!
MODULE PDAF_enkf_analysis_rlm

CONTAINS
SUBROUTINE PDAF_enkf_ana_rlm(step, dim_p, dim_obs_p, dim_ens, rank_ana, &
     state_p, ens_p, HZB, HX_p, HXbar_p, obs_p, &
     U_add_obs_err, U_init_obs_covar, screen, debug, flag)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE mpi
  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: mype, npes_filter, MPIerr, COMM_filter
  USE PDAFomi, &
       ONLY: omi_n_obstypes => n_obstypes, PDAFomi_gather_obsdims
  USE PDAF_enkf, &
       ONLY: PDAF_enkf_gather_resid, PDAF_enkf_obs_ensemble

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step          !< Current time step
  INTEGER, INTENT(in) :: dim_p         !< PE-local dimension of model state
  INTEGER, INTENT(in) :: dim_obs_p     !< PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens       !< Size of state ensemble
  INTEGER, INTENT(in) :: rank_ana      !< Rank to be considered for inversion of HPH
  REAL, INTENT(inout) :: state_p(dim_p)           !< PE-local ensemble mean state
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)    !< PE-local state ensemble
  REAL, INTENT(inout) :: HZB(dim_ens, dim_ens)    ! Ensemble tranformation matrix
  REAL, INTENT(in)    :: HX_p(dim_obs_p, dim_ens) !< PE-local observed ensemble
  REAL, INTENT(in)    :: HXbar_p(dim_obs_p)       !< PE-local observed state
  REAL, INTENT(in)    :: obs_p(dim_obs_p)         !< PE-local observation vector
  INTEGER, INTENT(in) :: screen        !< Verbosity flag
  INTEGER, INTENT(in) :: debug         !< Flag for writing debug output
  INTEGER, INTENT(inout) :: flag       !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_obs_covar, &      !< Initialize observation error covariance matrix
       U_add_obs_err                   !< Add observation error covariance matrix

! *** local variables ***
  INTEGER :: i, j, member              ! counters
  INTEGER :: dim_obs                   ! global dimension of observation vector
  REAL :: invdim_ens                   ! inverse of ensemble size
  REAL :: invdim_ensm1                 ! inverse of ensemble size minus 1
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  INTEGER, SAVE :: allocflag_b = 0     ! Flag whether first time allocation is done
  REAL, ALLOCATABLE :: HPH(:,:)        ! Temporary matrix for analysis
  REAL, ALLOCATABLE :: XminMean_b(:,:) ! Temporary matrix for analysis
  REAL, ALLOCATABLE :: HZ(:,:)         ! H(ensstate)-H(meanstate)
  REAL, ALLOCATABLE :: resid(:,:)      ! ensemble of global residuals
  REAL, ALLOCATABLE :: resid_p(:,:)    ! ensemble of local residuals
  INTEGER, ALLOCATABLE :: ipiv(:)      ! vector of pivot indices
  INTEGER :: sgesv_info                ! output flag of SGESV

  ! *** Variables for variant using pseudo inverse with eigendecompositon
  REAL, ALLOCATABLE :: eval(:)         ! vector of eigenvalues
  REAL, ALLOCATABLE :: rwork(:)        ! workarray for eigenproblem
  REAL, ALLOCATABLE :: evec(:,:)       ! matrix of eigenvectors
  REAL, ALLOCATABLE :: evec_temp(:,:)  ! matrix of eigenvectors
  REAL, ALLOCATABLE :: repres(:,:)     ! matrix of representer vectors
  INTEGER :: syev_info                 ! output flag of eigenproblem routine
  REAL    :: VL, VU                    ! temporary variables for SYEVX (never really used)
  INTEGER :: Ilower, Iupper            ! variables defining the interval of eigenvalues
  REAL    :: abstol                    ! error tolerance for eigenvalue problem
  INTEGER :: nEOF                      ! number of EOFs as computed by SYEVX
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
          'PDAF ', step, 'EnKF analysis - large-dim_obs version using transform matrix'
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

  ! *** Get global dimension of observation vector ***
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

  CALL PDAF_timeit(10, 'new')


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
     
  CALL PDAF_timeit(30, 'new')

  ! Initialize array of observed ensemble perturbations
  DO member = 1, dim_ens
     resid_p(:, member) = HX_p(:, member) - HXbar_p(:)
  END DO

  IF (debug>0) THEN
     DO i = 1, dim_ens
        WRITE (*,*) '++ PDAF-debug PDAF_enkf_analysis:', debug, &
             'process-local observed ensemble pert, member', i, &
             ' values (1:min(dim_obs_p,6)):', resid_p(1:min(dim_obs_p,6),i)
     END DO
  END IF

  ! Allgather global array of observed ensemble perturbations
  ALLOCATE(HZ(dim_obs, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs * dim_ens)

  CALL PDAF_enkf_gather_resid(dim_obs, dim_obs_p, dim_ens, resid_p, HZ)
  CALL PDAF_timeit(30, 'old')

  ! Finish computation of HPH
  ALLOCATE(HPH(dim_obs, dim_obs))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs * dim_obs)

  CALL PDAF_timeit(32, 'new')
  CALL gemmTYPE ('n', 't', dim_obs, dim_obs, dim_ens, &
       invdim_ensm1, HZ, dim_obs, HZ, dim_obs, &
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

  CALL PDAF_timeit(11, 'new')
  ! observation ensemble is initialized into the residual matrix
  CALL PDAF_enkf_obs_ensemble(step, dim_obs_p, dim_obs, dim_ens, resid_p, &
       obs_p, U_init_obs_covar, screen, flag)
  CALL PDAF_timeit(11, 'old')


! *************************************
! *** Compute matrix of innovations ***
! ***         D = Y - H X           ***
! *************************************

  CALL PDAF_timeit(12, 'new')
  CALL PDAF_timeit(51, 'new')

  ALLOCATE(resid(dim_obs, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs * dim_ens)

  resid_p(:,:) = resid_p(:,:) - HX_p(:,:)

  IF (debug>0) THEN
     DO i = 1, dim_ens
        WRITE (*,*) '++ PDAF-debug PDAF_enkf_analysis:', debug, 'process-local innovation member', i, &
             ' values (1:min(dim_obs_p,6)):', resid_p(1:min(dim_obs_p,6),i)
     END DO
  END IF

  ! Allgather residual
  CALL PDAF_enkf_gather_resid(dim_obs, dim_obs_p, dim_ens, resid_p, resid)

  DEALLOCATE(resid_p)

  CALL PDAF_timeit(12, 'old')


  whichupdate: IF (rank_ana > 0) THEN
! **************************************************
! *** Update using pseudo inverse of HPH         ***
! *** by performing incomplete eigendecompostion ***
! *** and using Moore-Penrose inverse of this    ***
! *** matrix                                     ***
! **************************************************

     CALL PDAF_timeit(13, 'new')

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

     CALL PDAF_timeit(35, 'new')

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

        CALL PDAF_timeit(35, 'old')


        ! ****************************************
        ! *** Compute ensemble of representer  ***
        ! *** vectors b as the product         ***
        ! ***           b = invHPH d           ***
        ! ****************************************

        CALL PDAF_timeit(36, 'new')
        CALL gemmTYPE('n', 'n', dim_obs, dim_ens, dim_obs, &
             1.0, HPH, dim_obs, resid, dim_obs, &
             0.0, repres, dim_obs)
        CALL PDAF_timeit(36, 'old')


        ! **************************************
        ! *** Update model state ensemble    ***
        ! ***    a   f         f     _    T  ***
        ! ***   x = x + K d = x + (X-X) HZ B ***
        ! **************************************

        ! *** HZB = HZ^T B
        CALL PDAF_timeit(37, 'new')
        CALL gemmTYPE('t', 'n', dim_ens, dim_ens, dim_obs, &
             1.0, HZ, dim_obs, repres, dim_obs, &
             0.0, HZB, dim_ens)
        CALL PDAF_timeit(37, 'old')

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_enkf_analysis:', debug, '  transform HZB', HZB

        CALL PDAF_timeit(13, 'old')
        CALL PDAF_timeit(14, 'new')

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

           ENSc: DO member = 1, dim_ens
              ! initialize XminMean
              XminMean_b(1 : blkupper - blklower + 1, member) = &
                   ens_p(blklower : blkupper, member) - &
                   state_p(blklower : blkupper)
           END DO ENSc

           ! *** Update ensemble
           CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
                invdim_ensm1, XminMean_b, maxblksize, HZB(1, 1), dim_ens, &
                1.0, ens_p(blklower, 1), dim_p)
           
        END DO blocking1

        CALL PDAF_timeit(14, 'old')

        DEALLOCATE(XminMean_b)

     ELSE
        ! Error in the EVP
        CALL PDAF_timeit(32, 'old')

     END IF EVPok

     ! *** Clean up ***
     DEALLOCATE(repres)

  ELSE whichupdate
! *******************************************
! *** Update using matrix HPH directly to ***      
! *** compute representer amplitudes b by ***
! *** solving HPH b = d for b.            ***
! *******************************************

     CALL PDAF_timeit(13, 'new')

     ! ****************************************
     ! *** Compute ensemble of representer  ***
     ! *** vectors b by solving             ***
     ! ***              HPH b = d           ***
     ! *** We use the LAPACK routine GESV   ***
     ! ****************************************
     ALLOCATE(ipiv(dim_obs))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'i', dim_obs)

     CALL PDAF_timeit(33, 'new')
     CALL gesvTYPE(dim_obs, dim_ens, HPH, dim_obs, ipiv, &
          resid, dim_obs, sgesv_info)
     CALL PDAF_timeit(33, 'old')


     ! *** check if solve was successful
     update: IF (sgesv_info /= 0) THEN
        WRITE (*, '(/5x, a/)') 'PDAF-ERROR(2): Problem in solve for Kalman gain !!!'
        flag = 2

        CALL PDAF_timeit(13, 'old')
     ELSE

     ! **************************************
     ! *** Update model state ensemble    ***
     ! ***    a   f         f     _    T  ***
     ! ***   x = x + K d = x + (X-X) HZ B ***
     ! **************************************


        ! *** HZB = HZ^T B
        CALL PDAF_timeit(34, 'new')
        CALL gemmTYPE('t', 'n', dim_ens, dim_ens, dim_obs, &
             1.0, HZ, dim_obs, resid, dim_obs, &
             0.0, HZB, dim_ens)
        CALL PDAF_timeit(34, 'old')

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_enkf_analysis:', debug, '  transform', invdim_ensm1*HZB

        CALL PDAF_timeit(13, 'old')
        CALL PDAF_timeit(14, 'new')

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

           DO member = 1, dim_ens
              ! initialize XminMean
              XminMean_b(1 : blkupper - blklower + 1, member) = &
                   ens_p(blklower : blkupper, member) - &
                   state_p(blklower : blkupper)
           END DO

           ! *** Update ensemble
           CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
                invdim_ensm1, XminMean_b, maxblksize, HZB(1, 1), dim_ens, &
                1.0, ens_p(blklower, 1), dim_p)
           
        END DO blocking2

        CALL PDAF_timeit(14, 'old')

        DEALLOCATE(XminMean_b)

     END IF update
     DEALLOCATE(ipiv)

  END IF whichupdate

  CALL PDAF_timeit(51, 'old')


! ********************
! *** Finishing up ***
! ********************

  ! Clean up
  DEALLOCATE(HZ)
  DEALLOCATE(resid)
  DEALLOCATE(HPH)

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_enkf_analysis -- END'

END SUBROUTINE PDAF_enkf_ana_rlm

!> Smoother extension for EnKF
!!
!! Smoother extension for the ensemble Kalman filter (EnKF).
!! The routine uses the matrix Ainv computed by the filter analysis
!! to perform the smoothing on past ensembles.
!!
!! Variant for domain decomposed states.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2013-04 - Lars Nerger - Initial code
!! * 2025-02 - Lars Nerger - merged into module PDAF_enkf_analysis_rlm
!! * Later revisions - see svn log
!!
SUBROUTINE PDAF_smoother_enkf(dim_p, dim_ens, dim_lag, Ainv, sens_p, &
     cnt_maxlag, forget, screen)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: mype
  USE PDAF_analysis_utils, &
       ONLY: PDAF_subtract_colmean

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p        !< PE-local dimension of model state
  INTEGER, INTENT(in) :: dim_ens      !< Size of ensemble
  INTEGER, INTENT(in) :: dim_lag      !< Number of past time instances for smoother
  REAL, INTENT(in)   :: Ainv(dim_ens, dim_ens)           !< Weight matrix for ensemble transformation
  REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag) !< PE-local smoother ensemble
  INTEGER, INTENT(inout) :: cnt_maxlag                   !< Count available number of time steps for smoothing
  REAL, INTENT(in)    :: forget       !< Forgetting factor
  INTEGER, INTENT(in) :: screen       !< Verbosity flag

! *** local variables ***
  INTEGER :: member, lagcol           ! Counters
  INTEGER :: n_lags                   ! Available number of time instances for smoothing
  INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  REAL :: fact
  REAL, ALLOCATABLE :: ens_blk(:,:)   ! Temporary block of state ensemble
  REAL, ALLOCATABLE :: W_smooth(:,:)  ! Weight matrix for smoothing

  
! **********************
! *** INITIALIZATION ***
! **********************

  ! Determine number of time instances for smoothing
  IF (cnt_maxlag >= dim_lag) THEN
     ! Already performed enough analysis to smooth over full lag
     n_lags = dim_lag
  ELSE
     ! Not yet enough analysis steps to smoother over full lag
     n_lags = cnt_maxlag
  END IF

  IF (mype == 0 .AND. screen > 0 .AND. n_lags > 0) THEN
     WRITE (*, '(a, 5x, a, i8)') 'PDAF', 'Perform smoothing up to lag', n_lags
  END IF

  ! init scale factor for weight matrix
  fact = sqrt(forget) / (Real(Dim_ens - 1))


! **********************************************
! *** Compute transform matrix for smoothing ***
! ***                                        ***
! *** W_smooth = (1 - 1_N) A_filter          ***
! **********************************************

  havelag: IF (n_lags > 0) THEN

     ALLOCATE(W_smooth(dim_ens, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

     W_smooth = Ainv

     ! Part 4: T W
     CALL PDAF_subtract_colmean(dim_ens, dim_ens, W_smooth)
  

! **********************************************
! *** Perform smoothing                      ***
! *** Transform state ensemble at past times ***
! ***          a    f     _                  ***
! ***         X  = X + (X-X) W_smooth        ***
! **********************************************

     ! Use block formulation for transformation
     maxblksize = 200
     ALLOCATE(ens_blk(maxblksize, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', maxblksize * dim_ens)
     ens_blk = 0.0 ! This is not really necessary, but there was a case having a problem without it. 

     ! *** Smooth for all available lags ***
     smoothing: DO lagcol = 1, n_lags

        ! Use block formulation for transformation
        blocking: DO blklower = 1, dim_p, maxblksize
           
           blkupper = MIN(blklower + maxblksize - 1, dim_p)

           ! Store former analysis ensemble
           DO member = 1, dim_ens
              ens_blk(1 : blkupper-blklower+1, member) &
                   = sens_p(blklower : blkupper, member, lagcol)
           END DO

           !                        a(i)   a(i-1)    a(i-1)
           ! Transform ensemble:   X    = X       + X       W_smooth
           CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
                fact, ens_blk(1, 1), maxblksize, W_smooth, dim_ens, &
                1.0, sens_p(blklower, 1, lagcol), dim_p)

        END DO blocking

     END DO smoothing

     DEALLOCATE(ens_blk, W_smooth)

  END IF havelag


! ********************
! *** Finishing up ***
! ********************
  
  ! Increment maxlag counter
  cnt_maxlag = cnt_maxlag + 1

  ! Set flag for memory counting
  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_smoother_enkf

END MODULE PDAF_enkf_analysis_rlm
