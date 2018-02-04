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
!$Id: PDAF-D_enkf_analysis_rsm.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_enkf_analysis_rsm --- Perform EnKF analysis step
!
! !INTERFACE:
SUBROUTINE PDAF_enkf_analysis_rsm(step, dim_p, dim_obs_p, dim_ens, rank_ana, &
     state_p, ens_p, forget, U_init_dim_obs, U_obs_op, &
     U_add_obs_err, U_init_obs, U_init_obs_covar, screen, flag)


! !DESCRIPTION:
! Analysis step of ensemble Kalman filter with 
! representer-type formulation.  In this version 
! HP is explicitly computed.  This variant is 
! optimal if the number of observations is 
! smaller than or equal to half of the ensemble 
! size.
! The final ensemble update uses a block
! formulation to reduce memory requirements.
!
! Variant for domain decomposition.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2003-10 - Lars Nerger - Initial code
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
  USE PDAF_mod_filter, &
       ONLY: obs_member
  USE PDAF_mod_filtermpi, &
       ONLY: mype, npes_filter, MPIerr, COMM_filter, MPI_SUM, MPI_INTEGER

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step      ! Current time step
  INTEGER, INTENT(in)  :: dim_p     ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p ! PE-local dimension of observation vector
  INTEGER, INTENT(in)  :: dim_ens   ! Size of state ensemble
  INTEGER, INTENT(in)  :: rank_ana  ! Rank to be considered for inversion of HPH
  REAL, INTENT(inout)  :: state_p(dim_p)        ! PE-local ensemble mean state
  REAL, INTENT(inout)  :: ens_p(dim_p, dim_ens) ! PE-local state ensemble
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
  INTEGER :: i, j, member  ! counters
  INTEGER :: dim_obs       ! global dimension of observation vector
  REAL :: invdim_ens       ! inverse of ensemble size
  REAL :: invdim_ensm1     ! inverse of ensemble size minus 1
  REAL :: sqrtinvforget    ! square root of inverse forgetting factor
  INTEGER, SAVE :: allocflag = 0       ! Flag for first-time allocation
  INTEGER, SAVE :: allocflag_b = 0     ! Flag for first-time allocation
  REAL, ALLOCATABLE :: HP_p(:,:)       ! Temporary matrix for analysis
  REAL, ALLOCATABLE :: HPH(:,:)        ! Temporary matrix for analysis
  REAL, ALLOCATABLE :: XminMean_p(:,:) ! Temporary matrix for analysis
  REAL, ALLOCATABLE :: resid(:,:)      ! ensemble of global residuals
  REAL, ALLOCATABLE :: resid_p(:,:)    ! ensemble of local residuals
  REAL, ALLOCATABLE :: m_state_p(:)    ! PE-local observed state
  REAL, ALLOCATABLE :: HXmean_p(:)     ! Temporary vector for analysis
  INTEGER, ALLOCATABLE :: ipiv(:)      ! vector of pivot indices
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
  REAL,EXTERNAL :: DLAMCH   ! function to specify tolerance of SYEVX
  REAL    :: eval_inv       ! inverse of an eigenvalue


! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, i7, 3x, a)') &
          'PDAF ', step, 'Assimilating observations - EnKF small-m version'
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


! *********************************
! *** Get observation dimension ***
! *********************************

  CALL U_init_dim_obs(step, dim_obs_p)

  IF (screen > 0) THEN
     WRITE (*, '(a, 5x, a13, 1x, i6, 1x, a, i10)') &
          'PDAF','--- PE-domain', mype, 'dimension of observation vector', dim_obs_p
  END IF

  ! Get global dimension of observation vector
  ! Using reals is a work around such that this also works with nullmpi.F90
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

  CALL PDAF_timeit(11, 'old')
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
  ENSa: DO member = 1, dim_ens

     ! spread out state ensemble according to forgetting factor
     IF (forget /= 1.0) THEN
        ens_p(:, member) = state_p(:) &
             + (ens_p(:, member) - state_p(:)) * sqrtinvforget
     END IF

     ! initialize XMINMEAN
     XminMean_p(:, member) = ens_p(:, member) - state_p(:)

  END DO ENSa

  ENSb: DO member = 1, dim_ens
     ! Store member index to make it accessible with PDAF_get_obsmemberid
     obs_member = member

     ! project state ensemble onto observation space
     ! stored in RESID_P for compactness
     CALL U_obs_op(step, dim_p, dim_obs_p, ens_p(:, member), resid_p(:, member))
  END DO ENSb

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

  ! Allgather residual
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

  DEALLOCATE(XminMean_p)


  ! *** Add observation error covariance ***
  ! ***       HPH^T = (HPH + R)          ***
  CALL U_add_obs_err(step, dim_obs, HPH)

  CALL PDAF_timeit(10, 'old')


! *****************************************
! *** generate ensemble of observations ***
! *****************************************

  CALL PDAF_timeit(15, 'new')
  ! observation ensemble is initialized into the residual matrix
  CALL PDAF_enkf_obs_ensemble(step, dim_obs_p, dim_obs, dim_ens, resid_p, &
       U_init_obs, U_init_obs_covar, screen, flag)
  CALL PDAF_timeit(15, 'old')


! ***********************************
! *** Compute matrix of residuals ***
! ***         D = Y - H X         ***
! ***********************************
  ALLOCATE(m_state_p(dim_obs_p))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p)

  CALL PDAF_timeit(12, 'new')
  ! *** Project state onto observation space and    ***
  ! *** compute observation residual (innovation) d ***
  DO member = 1, dim_ens
     ! Store member index to make it accessible with PDAF_get_obsmemberid
     obs_member = member

     ! Project state onto observation space
     CALL U_obs_op(step, dim_p, dim_obs_p, ens_p(:, member), m_state_p)

     ! get residual as difference of observation and
     ! projected state
     resid_p(:, member) = resid_p(:, member) - m_state_p(:)
  END DO
  DEALLOCATE(m_state_p)

  ! Allgather residual
  CALL PDAF_enkf_gather_resid(dim_obs, dim_obs_p, dim_ens, resid_p, resid)

  DEALLOCATE(resid_p)

  CALL PDAF_timeit(12, 'old')
  CALL PDAF_timeit(14, 'new')


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
        CALL PDAF_memcount(3, 'r', dim_obs * dim_ens + rank_ana &
             + 2 * dim_obs * rank_ana + 8 * dim_obs)
        CALL PDAF_memcount(3, 'i', 6 * dim_obs)
        allocflag_b = 1
     END IF
    
     CALL PDAF_timeit(13, 'new')
     CALL PDAF_timeit(36,'new')

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
     CALL syevxTYPE('v', 'i', 'u', dim_obs, HPH, &
          dim_obs, VL, VU, Ilower, Iupper, &
          abstol, nEOF, eval, evec, dim_obs, &
          rwork, 8 * dim_obs, iwork, ifail, syev_info)

     ! check if eigendecomposition was successful
     EVPok: IF (syev_info == 0) THEN
        ! Eigendecomposition OK, continue

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


        ! ***********************************
        ! *** Update model state ensemble ***
        ! ***          a   f              ***
        ! ***         x = x + K d         ***
        ! ***********************************

        CALL PDAF_timeit(16, 'new')

        CALL gemmTYPE('t', 'n', dim_p, dim_ens, dim_obs, &
             1.0, HP_p, dim_obs, repres, dim_obs, &
             1.0, ens_p, dim_p)
 
        CALL PDAF_timeit(16, 'old')

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
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs)

     CALL PDAF_timeit(13, 'new')
     CALL gesvTYPE(dim_obs, dim_ens, HPH, dim_obs, ipiv, &
          resid, dim_obs, sgesv_info)
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

        CALL PDAF_timeit(16, 'new')
        CALL gemmTYPE('t', 'n', dim_p, dim_ens, dim_obs, &
             1.0, HP_p, dim_obs, resid, dim_obs, &
             1.0, ens_p, dim_p)
        CALL PDAF_timeit(16, 'old')

        DEALLOCATE(resid)
     END IF update
      
     DEALLOCATE(ipiv)
     
  END IF whichupdate

  CALL PDAF_timeit(14, 'old')


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(HP_p)
  DEALLOCATE(HPH)

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_enkf_analysis_rsm
