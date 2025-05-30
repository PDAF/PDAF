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
!> Perform LEnKF analysis step
!!
!! Analysis step of ensemble Kalman filter with 
!! representer-type formulation.  In this version 
!! HP is explicitly computed.  This variant is 
!! optimal if the number of observations is 
!! smaller than or equal to half of the ensemble 
!! size.
!! The final ensemble update uses a block
!! formulation to reduce memory requirements.
!!
!! Variant for domain decomposition.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2003-10 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAF_lenkf_analysis_rsm

CONTAINS
SUBROUTINE PDAF_lenkf_ana_rsm(step, dim_p, dim_obs_p, dim_ens, rank_ana, &
     state_p, ens_p, HX_p, HXbar_p, obs_p, &
     U_add_obs_err, U_init_obs_covar, U_localize, screen, debug, flag)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE mpi
  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_parallel, &
       ONLY: mype, npes_filter, MPIerr, COMM_filter
  USE PDAFomi_obs_f, &
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
  REAL, INTENT(in)    :: HX_p(dim_obs_p, dim_ens) !< PE-local observed ensemble
  REAL, INTENT(in)    :: HXbar_p(dim_obs_p)       !< PE-local observed state
  REAL, INTENT(in)    :: obs_p(dim_obs_p)         !< PE-local observation vector
  INTEGER, INTENT(in) :: screen        !< Verbosity flag
  INTEGER, INTENT(in) :: debug         !< Flag for writing debug output
  INTEGER, INTENT(inout) :: flag       !< Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_obs_covar, &      !< Initialize observation error covariance matrix
       U_add_obs_err, &                !< Add observation error covariance matrix
       U_localize                      !< Apply localization to HP and HPH^T

! *** local variables ***
  INTEGER :: i, j, member              ! counters
  INTEGER :: dim_obs                   ! global dimension of observation vector
  REAL :: invdim_ens                   ! inverse of ensemble size
  REAL :: invdim_ensm1                 ! inverse of ensemble size minus 1
  INTEGER, SAVE :: allocflag = 0       ! Flag for first-time allocation
  INTEGER, SAVE :: allocflag_b = 0     ! Flag for first-time allocation
  REAL, ALLOCATABLE :: HP_p(:,:)       ! Temporary matrix for analysis
  REAL, ALLOCATABLE :: HPH(:,:)        ! Temporary matrix for analysis
  REAL, ALLOCATABLE :: XminMean_p(:,:) ! Temporary matrix for analysis
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
  REAL,EXTERNAL :: DLAMCH              ! function to specify tolerance of SYEVX
  REAL    :: eval_inv                  ! inverse of an eigenvalue


! **********************
! *** INITIALIZATION ***
! **********************

  CALL PDAF_timeit(51, 'new')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lenkf_analysis -- START'

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, i7, 3x, a)') &
          'PDAF ', step, 'Localized EnKF analysis - small-dim_obs version'
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
  ALLOCATE(XminMean_p(dim_p, dim_ens))
  IF (allocflag == 0) &
       CALL PDAF_memcount(3, 'r', dim_obs_p * dim_ens + dim_p * dim_ens)
  
  ! ***                             T ***
  ! *** get HP = H P and HPH = H P H  ***
  ! *** as ensemble means             ***

  ! Initialize array of ensemble perturbations
  ENSa: DO member = 1, dim_ens
     XminMean_p(:, member) = ens_p(:, member) - state_p(:)
  END DO ENSa

  CALL PDAF_timeit(30, 'new')

  ! Initialize array of observed ensemble perturbations
  DO member = 1, dim_ens
     resid_p(:, member) = HX_p(:, member) - HXbar_p(:)
  END DO

  IF (debug>0) THEN
     DO i = 1, dim_ens
        WRITE (*,*) '++ PDAF-debug PDAF_lenkf_analysis:', debug, &
             'process-local observed ensemble pert, member', i, &
             ' values (1:min(dim_obs_p,6)):', resid_p(1:min(dim_obs_p,6),i)
     END DO
  END IF

  ! Allgather global residual
  ALLOCATE(resid(dim_obs, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs * dim_ens)

  CALL PDAF_enkf_gather_resid(dim_obs, dim_obs_p, dim_ens, resid_p, resid)
  CALL PDAF_timeit(30, 'old')

  ! Finish computation of HP and HPH
  ALLOCATE(HP_p(dim_obs, dim_p))
  ALLOCATE(HPH(dim_obs, dim_obs))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs * dim_p + dim_obs * dim_obs)

  CALL PDAF_timeit(31, 'new')
  CALL gemmTYPE('n', 't', dim_obs, dim_p, dim_ens, &
       invdim_ensm1, resid, dim_obs, XminMean_p, dim_p, &
       0.0, HP_p, dim_obs)
  CALL PDAF_timeit(31, 'old')

  CALL PDAF_timeit(32, 'new')
  CALL gemmTYPE('n', 't', dim_obs, dim_obs, dim_ens, &
       invdim_ensm1, resid, dim_obs, resid, dim_obs, &
       0.0, HPH, dim_obs)
  CALL PDAF_timeit(32, 'old')

  CALL PDAF_timeit(51, 'old')

  DEALLOCATE(XminMean_p)

  ! For OMI: Gather global observation dimensions
  IF (omi_n_obstypes > 0) CALL PDAFomi_gather_obsdims()

  ! Apply localization
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lenkf_analysis -- call localize_covar'

  CALL PDAF_timeit(45, 'new')
  CALL U_localize(dim_p, dim_obs, HP_p, HPH)
  CALL PDAF_timeit(45, 'old')

  ! *** Add observation error covariance ***
  ! ***       HPH^T = (HPH + R)          ***
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lenkf_analysis -- call add_obs_err'

  CALL PDAF_timeit(46, 'new')
  CALL U_add_obs_err(step, dim_obs, HPH)
  CALL PDAF_timeit(46, 'old')

  CALL PDAF_timeit(10, 'old')
  CALL PDAF_timeit(51, 'new')


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

  ! *** Project state onto observation space and    ***
  ! *** compute observation residual (innovation) d ***
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lenkf_analysis -- call obs_op', dim_ens, 'times'

  resid_p(:,:) = resid_p(:,:) - HX_p(:,:)

  IF (debug>0) THEN
     DO i = 1, dim_ens
        WRITE (*,*) '++ PDAF-debug PDAF_lenkf_analysis:', debug, 'process-local innovation member', i, &
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
        CALL PDAF_memcount(3, 'r', dim_obs * dim_ens + rank_ana &
             + 2 * dim_obs * rank_ana + 8 * dim_obs)
        CALL PDAF_memcount(3, 'i', 6 * dim_obs)
        allocflag_b = 1
     END IF
    
     CALL PDAF_timeit(35,'new')

     ! **************************************
     ! *** compute pseudo inverse of HPH  ***
     ! *** using Moore-Penrose inverse    ***
     ! *** of rank reduced matrix         ***
     ! **************************************

     Iupper = dim_obs
     Ilower = dim_obs - rank_ana + 1
     abstol = 2 * DLAMCH('S')

     ! *** Decompose HPH = eigenvec ev eigenvec^T by   ***
     ! *** computing the RANK_ANA largest eigenvalues  ***
     ! *** and the corresponding eigenvectors          ***
     ! *** We use the LAPACK routine SYEVX             ***
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lenkf_analysis:', debug, &
          '  Compute eigenvalue decomposition of HPH^T'

     CALL syevxTYPE('v', 'i', 'u', dim_obs, HPH, &
          dim_obs, VL, VU, Ilower, Iupper, &
          abstol, nEOF, eval, evec, dim_obs, &
          rwork, 8 * dim_obs, iwork, ifail, syev_info)

     ! check if eigendecomposition was successful
     EVPok: IF (syev_info == 0) THEN
        ! Eigendecomposition OK, continue

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_lenkf_resample:', debug, '  eigenvalues', eval

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

        IF (debug>0) THEN
           DO i = 1, dim_ens
              WRITE (*,*) '++ PDAF-debug PDAF_lenkf_analysis:', debug, 'representer member', i, &
                   ' values (1:min(dim_p,6)):', repres(1:min(dim_p,6),i)
           END DO
        END IF

        CALL PDAF_timeit(36, 'old')
        CALL PDAF_timeit(13, 'old')


        ! ***********************************
        ! *** Update model state ensemble ***
        ! ***          a   f              ***
        ! ***         x = x + K d         ***
        ! ***********************************

        CALL PDAF_timeit(14, 'new')

        CALL gemmTYPE('t', 'n', dim_p, dim_ens, dim_obs, &
             1.0, HP_p, dim_obs, repres, dim_obs, &
             1.0, ens_p, dim_p)

        CALL PDAF_timeit(14, 'old')

     ELSE
        CALL PDAF_timeit(35, 'old')
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
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lenkf_analysis:', debug, &
          '  Compute representers using solver GESV'

     CALL PDAF_timeit(33, 'new')
     CALL gesvTYPE(dim_obs, dim_ens, HPH, dim_obs, ipiv, &
          resid, dim_obs, sgesv_info)
     CALL PDAF_timeit(33, 'old')

     IF (debug>0) THEN
        DO i = 1, dim_ens
           WRITE (*,*) '++ PDAF-debug PDAF_lenkf_analysis:', debug, 'representer member', i, &
                ' values (1:min(dim_obs,6)):', resid(1:min(dim_obs,6),i)
        END DO
     END IF
     CALL PDAF_timeit(13, 'old')

     ! *** check if solve was successful
     update: IF (sgesv_info /= 0) THEN
        WRITE (*, '(/a, 3x, a/)') 'PDAF', '!!! Problem in solve for Kalman gain !!!'
        flag = 2
     ELSE

        ! ***********************************
        ! *** Update model state ensemble ***
        ! ***    a   f         f    T     ***
        ! ***   x = x + K d = x + HP b    ***
        ! ***********************************

        CALL PDAF_timeit(14, 'new')
        CALL gemmTYPE('t', 'n', dim_p, dim_ens, dim_obs, &
             1.0, HP_p, dim_obs, resid, dim_obs, &
             1.0, ens_p, dim_p)
        CALL PDAF_timeit(14, 'old')

        DEALLOCATE(resid)
     END IF update
      
     DEALLOCATE(ipiv)
     
  END IF whichupdate

  CALL PDAF_timeit(51, 'old')


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(HP_p)
  DEALLOCATE(HPH)

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lenkf_analysis -- END'

END SUBROUTINE PDAF_lenkf_ana_rsm

END MODULE PDAF_lenkf_analysis_rsm
