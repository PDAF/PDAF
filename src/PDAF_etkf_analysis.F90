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
!> ETKF analysis cf. Hunt et al. (2007)
!!
!! Analysis step of the ETKF following Hunt et al., Efficient data 
!! assimilation for spatiotemporal chaos: A local ensemble transform 
!! Kalman filter. Physica D 230 (2007) 112-126.
!!
!! The implementation also supports an adaptive forgetting factor.
!!
!! Variant for domain decomposed states. 
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2009-07 - Lars Nerger - Initial code
!! * Later revisions - see svn log
!!
SUBROUTINE PDAF_etkf_analysis(step, dim_p, dim_obs_p, dim_ens, &
     state_p, Ainv, ens_p, state_inc_p, &
     HZ_p, HXbar_p, obs_p, forget, U_prodRinvA, &
     screen, incremental, type_trans, debug, flag)

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
  USE PDAF_analysis_utils, &
       ONLY: PDAF_subtract_rowmean, PDAF_subtract_colmean

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step         !< Current time step
  INTEGER, INTENT(in) :: dim_p        !< PE-local dimension of model state
  INTEGER, INTENT(in ) :: dim_obs_p   !< PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens      !< Size of ensemble
  REAL, INTENT(out)   :: state_p(dim_p)           !< on exit: PE-local forecast state
  REAL, INTENT(out)   :: Ainv(dim_ens, dim_ens)   !< on exit: weight matrix for ensemble transformation
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)    !< PE-local state ensemble
  REAL, INTENT(inout) :: state_inc_p(dim_p)       !< PE-local state analysis increment
  REAL, INTENT(inout) :: HZ_p(dim_obs_p, dim_ens) !< PE-local observed ensemble
  REAL, INTENT(in) :: HXbar_p(dim_obs_p)          !< PE-local observed state
  REAL, INTENT(in) :: obs_p(dim_obs_p)            !< PE-local observation vector
  REAL, INTENT(in)    :: forget       !< Forgetting factor
  INTEGER, INTENT(in) :: screen       !< Verbosity flag
  INTEGER, INTENT(in) :: incremental  !< Control incremental updating
  INTEGER, INTENT(in) :: type_trans   !< Type of ensemble transformation
  INTEGER, INTENT(in) :: debug        !< Flag for writing debug output
  INTEGER, INTENT(inout) :: flag      !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA             !< Provide product R^-1 A
       
! *** local variables ***
  INTEGER :: i, col, row              ! counters
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  INTEGER :: syev_info                ! Status flag for SYEV
  INTEGER :: ldwork                   ! Size of work array for SYEV
  INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
  REAL :: sqrtNm1                     ! Temporary variable: sqrt(dim_ens-1)
  REAL, ALLOCATABLE :: RiHZ_p(:,:)    ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: innov_p(:)     ! PE-local observation innovation
  REAL, ALLOCATABLE :: RiHZd(:)       ! Temporary vector for analysis 
  REAL, ALLOCATABLE :: RiHZd_p(:)     ! PE-local RiHZd
  REAL, ALLOCATABLE :: VRiHZd(:)      ! Temporary vector for analysis
  REAL, ALLOCATABLE :: tmp_Ainv(:,:)  ! Temporary storage of Ainv
  REAL, ALLOCATABLE :: Usqrt(:, :)    ! Square-root of matrix U
  REAL, ALLOCATABLE :: rndmat(:,:)    ! Temporary random matrix
  REAL, ALLOCATABLE :: ens_blk(:,:)   ! Temporary block of state ensemble
  REAL, ALLOCATABLE :: svals(:)       ! Singular values of Ainv
  REAL, ALLOCATABLE :: work(:)        ! Work array for SYEV
  INTEGER :: incremental_dummy        ! Dummy variable to avoid compiler warning
  REAL :: state_inc_p_dummy(1)        ! Dummy variable to avoid compiler warning


! **********************
! *** INITIALIZATION ***
! **********************

  CALL PDAF_timeit(51, 'new')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_etkf_analysis -- START'

  ! Initialize variable to prevent compiler warning
  incremental_dummy = incremental
  state_inc_p_dummy(1) = state_inc_p(1)


! **************************
! *** Compute innovation ***
! ***     d = y - H x    ***
! **************************

  IF (dim_obs_p > 0) THEN
     ! The innovatiopn only exists for domains with observations
     
     CALL PDAF_timeit(10, 'new')

     ALLOCATE(innov_p(dim_obs_p))

     innov_p = obs_p - HXbar_p

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_etkf_analysis:', debug, &
             'innovation d(1:min(dim_obs_p,10))', innov_p(1:min(dim_obs_p,10))
        WRITE (*,*) '++ PDAF-debug PDAF_etkf_analysis:', debug, &
             'MIN/MAX of innovation', MINVAL(innov_p), MAXVAL(innov_p)
     END IF
     CALL PDAF_timeit(10, 'old')
  END IF

  CALL PDAF_timeit(51, 'old')


! **********************************************
! ***   Compute analyzed matrix Ainv         ***
! ***                                        ***
! ***     -1                 T  -1           ***
! ***    U  = forget I + (HZ)  R   HZ        ***
! **********************************************

  CALL PDAF_timeit(11, 'new')

  ALLOCATE(Usqrt(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

  haveobsA: IF (dim_obs_p > 0) THEN
     ! *** The contribution of observation matrix ist only ***
     ! *** computed for domains with observations          ***

     CALL PDAF_timeit(51, 'new')

     ! Subtract mean from observed ensemble: HZ = [Hx_1 ... Hx_N] T
     CALL PDAF_subtract_rowmean(dim_obs_p, dim_ens, HZ_p)

     CALL PDAF_timeit(51, 'old')

     ! ***                RiHZ = Rinv HZ                
     ! *** This is implemented as a subroutine thus that
     ! *** Rinv does not need to be allocated explicitly.
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_etkf_analysis -- call prodRinvA_l'

     ALLOCATE(RiHZ_p(dim_obs_p, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * dim_ens)

     CALL PDAF_timeit(48, 'new')
     CALL U_prodRinvA(step, dim_obs_p, dim_ens, obs_p, HZ_p, RiHZ_p)
     CALL PDAF_timeit(48, 'old')

     CALL PDAF_timeit(51, 'new')

     ! *** Initialize Ainv = (N-1) I ***
     Ainv = 0.0
     DO i = 1, dim_ens
        Ainv(i, i) = REAL(dim_ens - 1)
     END DO

     ! ***             T        ***
     ! ***  Compute  HZ  RiHZ   ***

     CALL gemmTYPE('t', 'n', dim_ens, dim_ens, dim_obs_p, &
          1.0, HZ_p, dim_obs_p, RiHZ_p, dim_obs_p, &
          0.0, Usqrt, dim_ens)

  ELSE haveobsA
     ! *** For domains with dim_obs_p=0 there is no ***
     ! *** direct observation-contribution to Ainv  ***
 
     CALL PDAF_timeit(51, 'new')
    
     ! *** Initialize Ainv = (N-1) I ***
     Ainv = 0.0
     DO i = 1, dim_ens
        Ainv(i, i) = REAL(dim_ens - 1)
     END DO

     ! No observation-contribution to Ainv from this domain
     Usqrt = 0.0

  END IF haveobsA

  ! get total sum on all filter PEs
  ALLOCATE(tmp_Ainv(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

  CALL MPI_allreduce(Usqrt, tmp_Ainv, dim_ens**2, &
       MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

  ! *** Complete computation of Ainv ***
  ! ***   -1          -1    T        ***
  ! ***  U  = forget U  + HZ RiHZ    ***
  Ainv = forget * Ainv + tmp_Ainv

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_etkf_analysis:', debug, '  A^-1', Ainv

  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(11, 'old')


! ***********************************************
! *** Compute weight for model state update   ***
! ***                                         ***
! ***              T                    f     ***
! ***    w = U RiHZ d  with d = (y - H x )    ***
! ***                                         ***
! ***********************************************

  CALL PDAF_timeit(51, 'new')
  CALL PDAF_timeit(12, 'new')

  ! *** Subtract ensemble mean from ensemble matrix ***
  ! ***          Z = [x_1, ..., x_N] T              ***
  CALL PDAF_subtract_rowmean(dim_p, dim_ens, ens_p)
  
  ! *** Compute RiHZd = RiHZ^T d ***
  ALLOCATE(RiHZd_p(dim_ens))
  ALLOCATE(RiHZd(dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 2 * dim_ens)

  haveobsC: IF (dim_obs_p > 0) THEN
     ! *** RiHZd_p/=0 only with observations ***

     ! local products (partial sum)
     CALL gemvTYPE('t', dim_obs_p, dim_ens, 1.0, RiHZ_p, &
          dim_obs_p, innov_p, 1, 0.0, RiHZd_p, 1)

     DEALLOCATE(RiHZ_p, innov_p)

  ELSE haveobsC

     RiHZd_p = 0.0

  END IF haveobsC

  ! get total sum on all filter PEs
  CALL MPI_allreduce(RiHZd_p, RiHZd, dim_ens, &
       MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_etkf_analysis:', debug, '  (HXT R^-1)^T d', RiHZd

  DEALLOCATE(RiHZd_p)

 
  ! *** Compute weight vector for state analysis:        ***
  ! ***          w = U RiHZd                             ***
  ! *** Use singular value decomposition of Ainv         ***
  ! ***        Ainv = ASB^T                              ***
  ! *** Then: U = A S^(-1) B                             ***
  ! *** The decomposition is also used for the symmetric ***
  ! *** square-root for the ensemble transformation.     ***

  ! Invert Ainv using SVD
  ALLOCATE(svals(dim_ens))
  ALLOCATE(work(3 * dim_ens))
  ldwork = 3 * dim_ens
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3 * dim_ens)
    
  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_etkf_analysis:', debug, &
       '  Compute eigenvalue decomposition of A^-1'
    
  ! Compute SVD of Ainv
  CALL syevTYPE('v', 'l', dim_ens, Ainv, dim_ens, svals, work, ldwork, syev_info)

  DEALLOCATE(work)

  ! Check if SVD was successful
  IF (syev_info == 0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_etkf_resample:', debug, '  eigenvalues', svals

     flag = 0
  ELSE
     WRITE (*, '(/5x, a/)') 'PDAF-ERROR(1): Problem in SVD of inverse of U !!!'
     flag = 1
  END IF

  ! *** Compute w = U RiHZd stored in RiHZd
  check0: IF (flag == 0) THEN

     ALLOCATE(VRiHZd(dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)

     CALL gemvTYPE('t', dim_ens, dim_ens, 1.0, Ainv, &
          dim_ens, RiHZd, 1, 0.0, VRiHZd, 1)
     
     DO row = 1, dim_ens
        VRiHZd(row) = VRiHZd(row) / svals(row)
     END DO
  
     CALL gemvTYPE('n', dim_ens, dim_ens, 1.0, Ainv, &
          dim_ens, VRiHZd, 1, 0.0, RiHZd, 1)

     DEALLOCATE(VRiHZd)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_etkf_analysis:', debug, '  A(HXT R^-1)^T d', RiHZd

  END IF check0
     
  CALL PDAF_timeit(12, 'old')


! ************************************************
! ***     Transform state ensemble             ***
! ***              a   _f   f                  ***
! ***             X  = X + X  W                ***
! *** The weight matrix W is stored in Ainv.   ***
! ************************************************

! *** Prepare weight matrix for ensemble transformation ***

  check1: IF (flag == 0) THEN

     IF (mype == 0 .AND. screen > 0) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', 'Perform ensemble transformation'
     END IF

     CALL PDAF_timeit(20, 'new')

     ! Part 1: square-root of U
     DO col = 1, dim_ens
        DO row = 1, dim_ens
           tmp_Ainv(row, col) = Ainv(row, col) / SQRT(svals(col))
        END DO
     END DO

     sqrtNm1 = SQRT(REAL(dim_ens-1))
     CALL gemmTYPE('n', 't', dim_ens, dim_ens, dim_ens, &
          sqrtNm1, tmp_Ainv, dim_ens, Ainv, dim_ens, &
          0.0, Usqrt, dim_ens)

     ! Part 2 - Optional 
     ! Multiply by orthogonal random matrix with eigenvector (1,...,1)^T
     multrnd: IF (type_trans == 2) THEN
        WRITE (*,'(a, 5x,a)') 'PDAF', '--- Apply random rotation to ensemble'

        ALLOCATE(rndmat(dim_ens, dim_ens))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)
        
        ! Initialize random matrix
        CALL PDAF_generate_rndmat(dim_ens, rndmat, 2)

        CALL gemmTYPE('n', 'n', dim_ens, dim_ens, dim_ens, &
             1.0, Usqrt, dim_ens, rndmat, dim_ens, &
             0.0, tmp_Ainv, dim_ens)

        DEALLOCATE(rndmat)
     ELSE
        ! Non-random case
        tmp_Ainv = Usqrt
     END IF multrnd
     

     ! Part 3: W = sqrt(U) + w
     DO col = 1, dim_ens
        DO row = 1, dim_ens
           Ainv(row, col) = tmp_Ainv(row, col) + RiHZd(row)
        END DO
     END DO

     DEALLOCATE(tmp_Ainv, svals)
     DEALLOCATE(RiHZd, Usqrt)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_etkf_resample:', debug, '  transform', Ainv

     CALL PDAF_timeit(20, 'old')


! *** Perform ensemble transformation ***

     CALL PDAF_timeit(21, 'new')

     ! Use block formulation for transformation
     maxblksize = 200
     IF (mype == 0 .AND. screen > 0) &
          WRITE (*, '(a, 5x, a, i5)') &
          'PDAF', '--- use blocking with size ', maxblksize
        
     ALLOCATE(ens_blk(maxblksize, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', maxblksize * dim_ens)

     blocking: DO blklower = 1, dim_p, maxblksize
           
        blkupper = MIN(blklower + maxblksize - 1, dim_p)

        ! Store forecast ensemble
        DO col = 1, dim_ens
           ens_blk(1 : blkupper - blklower + 1, col) &
                = ens_p(blklower : blkupper, col)
        END DO
           
        DO col = 1,dim_ens
           ens_p(blklower : blkupper, col) = state_p(blklower : blkupper)
        END DO

        !                        a  _f   f
        ! Transform ensemble:   X = X + X  W

        CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
             1.0, ens_blk(1, 1), maxblksize, Ainv(1, 1), dim_ens, &
             1.0, ens_p(blklower, 1), dim_p)

     END DO blocking

     CALL PDAF_timeit(21, 'old')

     DEALLOCATE(ens_blk)

  END IF check1


! ********************
! *** Finishing up ***
! ********************

  ! Apply T from left side to allow for smoothing
  CALL PDAF_subtract_colmean(dim_ens, dim_ens, Ainv)

  CALL PDAF_timeit(51, 'old')

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_etkf_analysis -- END'

END SUBROUTINE PDAF_etkf_analysis
