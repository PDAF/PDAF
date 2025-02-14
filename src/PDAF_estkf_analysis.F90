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
!> ESTKF analysis/ensemble transformation
!!
!! Analysis step of the ESTKF with direct
!! transformation of the forecast into the 
!! analysis ensemble. This variant does not
!! compute the analysis state, but only the
!! analysis ensemble, whose mean is the analysis
!! state.
!! Supported is also the adaptive forgetting factor.
!!
!! Variant for domain decomposed states.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2011-09 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAF_estkf_analysis

CONTAINS
SUBROUTINE PDAF_estkf_ana(step, dim_p, dim_obs_p, dim_ens, rank, &
     state_p, Ainv, ens_p, state_inc_p, &
     HL_p, HXbar_p, obs_p, forget, U_prodRinvA, &
     screen, incremental, type_sqrt, type_trans, TA, debug, flag)

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
       ONLY: PDAF_seik_Omega, PDAF_estkf_AOmega, PDAF_estkf_OmegaA

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step         !< Current time step
  INTEGER, INTENT(in) :: dim_p        !< PE-local dimension of model state
  INTEGER, INTENT(inout) :: dim_obs_p !< PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens      !< Size of ensemble
  INTEGER, INTENT(in) :: rank         !< Rank of initial covariance matrix
  REAL, INTENT(inout) :: state_p(dim_p)           !< on exit: PE-local forecast mean state
  REAL, INTENT(inout) :: Ainv(rank, rank)         !< Inverse of matrix A - temporary use only
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)    !< PE-local state ensemble
  REAL, INTENT(inout) :: state_inc_p(dim_p)       !< PE-local state analysis increment
  REAL, INTENT(inout) :: HL_p(dim_obs_p, dim_ens) !< PE-local observed ensemble
  REAL, INTENT(in)    :: HXbar_p(dim_obs_p)       !< PE-local observed state
  REAL, INTENT(in)    :: obs_p(dim_obs_p)         !< PE-local observation vector
  REAL, INTENT(in)    :: forget       !< Forgetting factor
  INTEGER, INTENT(in) :: screen       !< Verbosity flag
  INTEGER, INTENT(in) :: incremental  !< Control incremental updating
  INTEGER, INTENT(in) :: type_sqrt    !< Type of square-root of A
                                      !< (0): symmetric sqrt; (1): Cholesky decomposition
  INTEGER, INTENT(in) :: type_trans   !< Type of ensemble transformation
  REAL, INTENT(inout) :: TA(dim_ens, dim_ens)     !< Ensemble transformation matrix
  INTEGER, INTENT(in) :: debug        !< Flag for writing debug output
  INTEGER, INTENT(inout) :: flag      !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA             !< Provide product R^-1 with some matrix


! *** local variables ***
  INTEGER :: i, j, col, row           ! counters
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  INTEGER :: lib_info                 ! Status flag for LAPACK calls
  INTEGER :: ldwork                   ! Size of work array for syevTYPE
  INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
  REAL    :: fac                      ! Temporary variable
  LOGICAL :: storeOmega = .FALSE.     ! Store matrix Omega instead of recomputing it
  LOGICAL, SAVE :: firsttime = .TRUE. ! Indicates first call to resampling
  REAL, ALLOCATABLE :: RiHL_p(:,:)    ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: innov_p(:)     ! PE-local observation residual
  REAL, ALLOCATABLE :: RiHLd(:)       ! Temporary vector for analysis 
  REAL, ALLOCATABLE :: RiHLd_p(:)     ! PE-local RiHLd
  REAL, ALLOCATABLE :: VRiHLd(:)      ! Temporary vector for analysis
  REAL, ALLOCATABLE :: Asqrt(:, :)    ! Square-root of matrix A
  REAL, ALLOCATABLE :: Ainv_p(:,:)    ! Ainv for PE-local domain
  REAL, ALLOCATABLE :: tmp_Ainv(:,:)  ! Temporary storage of Ainv
  REAL, ALLOCATABLE :: Omega(:,:)     ! Orthogonal matrix Omega
  REAL, ALLOCATABLE :: OmegaT(:,:)    ! Transpose of Omega
  REAL, SAVE, ALLOCATABLE :: OmegaTsave(:,:) ! Saved transpose of Omega
  REAL, ALLOCATABLE :: ens_blk(:,:)   ! Temporary blocked state ensemble
  REAL, ALLOCATABLE :: svals(:)       ! Singular values of Ainv
  REAL, ALLOCATABLE :: work(:)        ! Work array for SYEVTYPE
  INTEGER, ALLOCATABLE :: ipiv(:)     ! vector of pivot indices for GESVTYPE
  REAL :: state_inc_p_dummy(1)        ! Dummy variable to avoid compiler warning


! **********************
! *** INITIALIZATION ***
! **********************

  CALL PDAF_timeit(51, 'new')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_estkf_analysis -- START'

  ! Initialize variable to prevent compiler warning
  state_inc_p_dummy(1) = state_inc_p(1)


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
        WRITE (*,*) '++ PDAF-debug PDAF_estkf_analysis:', debug, &
             'innovation d(1:min(dim_obs_p,10))', innov_p(1:min(dim_obs_p,10))
        WRITE (*,*) '++ PDAF-debug PDAF_estkf_analysis:', debug, &
             'MIN/MAX of innovation', MINVAL(innov_p), MAXVAL(innov_p)
     END IF

     CALL PDAF_timeit(10, 'old')
  END IF haveobsB

  CALL PDAF_timeit(51, 'old')


! *************************************************
! ***   Compute analyzed matrix Ainv            ***
! ***                                           ***
! ***  -1                       T  -1           ***
! *** A  = forget*(N-1) I + (HL)  R  (HL)       ***
! ***                                           ***
! *************************************************

  CALL PDAF_timeit(11, 'new')

  haveobsA: IF (dim_obs_p > 0) THEN
     ! *** The contribution of observation matrix ist only ***
     ! *** computed for domains with observations          ***

     CALL PDAF_timeit(51, 'new')

     ! Compute HL = [Hx_1 ... Hx_N] T
     CALL PDAF_estkf_AOmega(dim_obs_p, dim_ens, HL_p)

     CALL PDAF_timeit(51, 'old')
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

  ELSE haveobsA
     ! *** For domains with dim_obs_p=0 there is no ***
     ! *** direct observation-contribution to Ainv  ***

     CALL PDAF_timeit(31, 'new')
     CALL PDAF_timeit(51, 'new')

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
     
  END IF haveobsA

  ! get total sum on all filter PEs
  CALL MPI_allreduce(Ainv_p, tmp_Ainv, rank * rank, &
       MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

  ! *** Complete computation of Ainv  ***
  ! ***   -1                T         ***
  ! ***  A  = forget I  + HL RiHL     ***

  Ainv = forget * Ainv + tmp_Ainv

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_estkf_analysis:', debug, '  A^-1', Ainv

  DEALLOCATE(Ainv_p)

  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(31, 'old')
  CALL PDAF_timeit(11, 'old')


! ***********************************************
! *** Compute weight for model state update   ***
! ***                                         ***
! ***              T                    f     ***
! ***    w = A RiHL d  with d = (y - H x )    ***
! ***********************************************

  CALL PDAF_timeit(51, 'new')
  CALL PDAF_timeit(12, 'new')

  ! *** RiHLd = RiHL^T d ***
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

     ! Compute SVD of Ainv
     CALL syevTYPE('v', 'l', rank, Ainv, rank, svals, work, ldwork, lib_info)

     DEALLOCATE(work)

     ! Compute product RiHLd A
     IF (lib_info==0) THEN
        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_estkf_analysis:', debug, '  eigenvalues', svals

        ALLOCATE(VRiHLd(rank))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank)

        CALL gemvTYPE('t', rank, rank, 1.0, Ainv, &
             rank, RiHLd, 1, 0.0, VRiHLd, 1)
     
        DO row = 1,rank
           VRiHLd(row) = VRiHLd(row) / svals(row)
        END DO
  
        CALL gemvTYPE('n', rank, rank, 1.0, Ainv, &
             rank, VRiHLd, 1, 0.0, RiHLd, 1)

        DEALLOCATE(VRiHLd)

     END IF
  END IF typeainv1

  CALL PDAF_timeit(12, 'old')

  ! *** check whether solve was successful
  IF (lib_info == 0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_estkf_analysis:', debug, '  A(HXT R^-1)^T d', RiHLd

     flag = 0
  ELSE
     WRITE (*, '(/5x, a/)') 'PDAF-ERROR(1): Problem in computation of analysis weights!!!'
     flag = 1
  END IF


! ****************************************************************
! ***     Transform state ensemble                             ***
! ***              a   _f   f                                  ***
! ***             X  = X + X  W                                ***
! *** with                                   -T      T         ***
! ***          W = Omega (RiHLd + sqrt(N-1) C   Omega )        ***
! ****************************************************************

! *** Prepare weight matrix for ensemble transformation ***

  check1: IF (flag == 0) THEN

     IF (mype == 0 .AND. screen > 0) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', 'Perform ensemble transformation'
        IF (type_sqrt == 1) THEN
           WRITE (*, '(a, 5x, a)') 'PDAF', '--- use Cholesky square-root of A'
        ELSE
           WRITE (*, '(a, 5x, a)') 'PDAF', '--- use symmetric square-root of A'
        END IF
     END IF

     CALL PDAF_timeit(20, 'new')
     CALL PDAF_timeit(32, 'new')

     ! Part 1: square-root of A

     ! Asqrt is allocated with dim_ens cols, because this is 
     ! required further below. Now only rank columns are used
     ALLOCATE(Asqrt(rank, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens * rank)

     typeainv2: IF (type_sqrt == 1) THEN
        ! Variant, if Ainv has been inverted above by solving

        ! *** Store Ainv for temporary use
        tmp_Ainv(:,:) = Ainv(:,:)

        ! Cholesky decomposition of tmp_Ainv = C C^T
        CALL potrfTYPE('l', rank, tmp_Ainv, rank, lib_info)

        ! check if Cholesky decomposition was successful
        CholeskyOK: IF (lib_info == 0) THEN
           ! Decomposition OK, continue
           flag = 0
        ELSE
           ! Decomposition failed
           WRITE (*, '(/5x, a/)') &
                'PDAF-ERROR(3): Problem with Cholesky decomposition of Ainv !!!'
           flag = 3
        ENDIF CholeskyOK

     ELSE typeainv2
        ! Variant, if SVD inversion of Ainv has been performed

        DO col = 1, rank
           DO row = 1, rank
              Asqrt(row, col) = Ainv(row, col) / SQRT(svals(col))
           END DO
        END DO

        CALL gemmTYPE('n', 't', rank, rank, rank, &
             1.0, Asqrt, rank, Ainv, rank, &
             0.0, tmp_Ainv, rank)
        DEALLOCATE(svals)

        ! Set flag
        flag = 0

     END IF typeainv2

     CALL PDAF_timeit(32, 'old')

  END IF check1

  
  ! *** Part 2: Initialize the transformation Omega ***

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
           ! *** At first call to ESTKF_ANALYSIS initialize  ***
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
              Asqrt = TRANSPOSE(Omega)
              ! store transposed Omega
              OmegaTsave = Asqrt
           END IF

           firsttime  = .FALSE.
      
           DEALLOCATE(Omega)

        ELSE first
           IF (mype == 0 .AND. screen > 0) &
                WRITE (*, '(a, 5x, a)') 'PDAF', '--- use stored Omega'
           IF (type_sqrt == 1) THEN
              OmegaT = OmegaTsave
           ELSE
              Asqrt = OmegaTsave
           END IF
        END IF first

     ELSE Omega_store
     IF (debug>0) THEN
        IF (type_sqrt == 1) THEN
           WRITE (*,*) '++ PDAF-debug PDAF_estkf_analysis:', debug, '  Omega^T', OmegaT
        ELSE
           WRITE (*,*) '++ PDAF-debug PDAF_estkf_analysis:', debug, '  Omega^T', Asqrt
        END IF
     END IF

        ! *** Initialize the matrix Omega in SEIK_Omega ***
        ! *** each time ESTKF_ANALYSIS is called        ***

        ALLOCATE(Omega(dim_ens, rank))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens * rank)
      
        ! *** Generate uniform orthogonal matrix OMEGA ***
        CALL PDAF_seik_Omega(rank, Omega, type_trans, screen)
        
        ! transpose Omega
        IF (type_sqrt == 1) THEN
           OmegaT = TRANSPOSE(Omega)
        ELSE
           Asqrt = TRANSPOSE(Omega)
        END IF

        DEALLOCATE(Omega)

     END IF Omega_store
     CALL PDAF_timeit(33, 'old')


     ! *** Part 3: Product A^(1/2) Omega ***
    
     CALL PDAF_timeit(34, 'new')
     IF (type_sqrt == 1) THEN
        ! A = (Omega C^(-1)) by solving Ct A = OmegaT for A
        CALL trtrsTYPE('L', 'T', 'N', rank, dim_ens, &
             tmp_Ainv, rank, OmegaT, rank, lib_info)
     ELSE
        ! TMP_AINV already contains matrix C (no more inversion)

        CALL gemmTYPE('n', 'n', rank, dim_ens, rank, &
             1.0, tmp_Ainv, rank, Asqrt, rank, &
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

     DEALLOCATE(Asqrt)
  END IF check2

  check3: IF (flag == 0) THEN

     ! *** Part 4: Add RiHLd and multiply by scaling factor

     CALL PDAF_timeit(35, 'new')

     ! *** Add RiHLd to At
     fac = SQRT(REAL(dim_ens - 1))

     IF (incremental /= 2) THEN
        DO j = 1, dim_ens
           DO i = 1, rank
              OmegaT(i,j) = fac * OmegaT(i,j) + RiHLd(i)
           END DO
        END DO
     ELSE
        ! For ensemble 3D-Var update only ensemble perturbations
        DO j = 1, dim_ens
           DO i = 1, rank
              OmegaT(i,j) = fac * OmegaT(i,j)
           END DO
        END DO
     END IF

     DEALLOCATE(RiHLd)

     ! *** Omega A^T (A^T stored in OmegaT) ***
     CALL PDAF_estkf_OmegaA(rank, dim_ens, OmegaT, TA)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_estkf_analysis:', debug, '  transform', TA

     CALL PDAF_timeit(35, 'old')
     CALL PDAF_timeit(20, 'old')

! *** Perform ensemble transformation ***

     CALL PDAF_timeit(21, 'new')

     ! Use block formulation for transformation
     maxblksize = 200
     IF (mype == 0 .AND. screen > 0) &
          WRITE (*, '(a, 5x, a, i5)') &
          'PDAF', '--- Ensemble update: use blocking with size ', maxblksize
        
     ALLOCATE(ens_blk(maxblksize, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', maxblksize * dim_ens)

     blocking: DO blklower = 1, dim_p, maxblksize
           
        blkupper = MIN(blklower + maxblksize - 1, dim_p)

        ! Store old state ensemble
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
             1.0, ens_blk(1, 1), maxblksize, TA(1, 1), dim_ens, &
             1.0, ens_p(blklower, 1), dim_p)

     END DO blocking

     CALL PDAF_timeit(21, 'old')

     DEALLOCATE(ens_blk)
     DEALLOCATE(OmegaT)

  ELSE

     CALL PDAF_timeit(20, 'old')

  END IF check3

  CALL PDAF_timeit(51, 'old')


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(tmp_Ainv)

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_estkf_analysis -- END'

END SUBROUTINE PDAF_estkf_ana

END MODULE PDAF_estkf_analysis
