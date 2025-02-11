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
!> Perform SEIK analysis and resampling - old formulation
!!
!! Analysis step of the SEIK filter
!! with adaptive forgetting factor.
!!
!! Variant for domain decomposed states. 
!! Old formulation regarding application of T.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! 2003-10 - Lars Nerger - Initial code
!! Later revisions - see svn log
!!
MODULE PDAF_seik_analysis

CONTAINS
SUBROUTINE PDAFseik_analysis(step, dim_p, dim_obs_p, dim_ens, rank, &
     state_p, Uinv, ens_p, state_inc_p, &
     HL_p, HXbar_p, obs_p, forget, U_prodRinvA, &
     screen, incremental, debug, flag)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE mpi
  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: MPIerr, COMM_filter
  USE PDAF_analysis_utils, &
       ONLY: PDAF_seik_matrixT, PDAF_seik_Uinv

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step         !< Current time step
  INTEGER, INTENT(in) :: dim_p        !< PE-local dimension of model state
  INTEGER, INTENT(in) :: dim_obs_p    !< PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens      !< Size of ensemble
  INTEGER, INTENT(in) :: rank         !< Rank of initial covariance matrix
  REAL, INTENT(inout) :: state_p(dim_p) !< PE-local model state
  REAL, INTENT(inout) :: Uinv(rank, rank)         !< Inverse of eigenvalue matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)    !< PE-local state ensemble
  REAL, INTENT(inout) :: state_inc_p(dim_p)       !< PE-local state analysis increment
  REAL, INTENT(inout) :: HL_p(dim_obs_p, dim_ens) !< PE-local observed ensemble (perturbations)
  REAL, INTENT(in)    :: HXbar_p(dim_obs_p)       !< PE-local observed state
  REAL, INTENT(in)    :: obs_p(dim_obs_p)         !< PE-local observation vector
  REAL, INTENT(in)    :: forget       !< Forgetting factor
  INTEGER, INTENT(in) :: screen       !< Verbosity flag
  INTEGER, INTENT(in) :: incremental  !< Control incremental updating
  INTEGER, INTENT(in) :: debug        !< Flag for writing debug output
  INTEGER, INTENT(inout) :: flag      !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA             !< Provide product R^-1 A
       
! *** local variables ***
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  REAL, ALLOCATABLE :: RiHL_p(:,:)    ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: Uinv_p(:,:)    ! local Uinv
  REAL, ALLOCATABLE :: Uinv_inc(:,:)  ! increment for Uinv
  REAL, ALLOCATABLE :: innov_p(:)     ! observation residual
  REAL, ALLOCATABLE :: RiHLd(:)       ! Temporary vector for analysis 
  REAL, ALLOCATABLE :: RiHLd_p(:)     ! local RiHLd
  REAL, ALLOCATABLE :: temp_Uinv(:,:) ! Temporary storage of Uinv
  INTEGER, ALLOCATABLE :: ipiv(:)     ! vector of pivot indices for GESV
  INTEGER :: gesv_info                ! control flag for GESV


! **********************
! *** INITIALIZATION ***
! **********************

  CALL PDAF_timeit(51, 'new')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_seik_analysis -- START'


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


! **********************************************
! ***   Compute analyzed matrix Uinv         ***
! ***                                        ***
! ***  -1                T        T  -1      ***
! *** U  = forget*(r+1) T T + (HL)  R  (HL)  ***
! ***  i                          i  i     i ***
! **********************************************

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

     ! *** Initialize Uinv = (r+1) T^T T ***
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
  CALL MPI_allreduce(Uinv_p, Uinv_inc, rank**2, &
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

  ! ************************************************
  ! *** Compute matrix L = [x_1, ..., x_(r+1)] T ***
  ! ************************************************

  CALL PDAF_seik_matrixT(dim_p, dim_ens, ens_p)

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
  ! *** routine GESV.                   ***
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

  CALL PDAF_timeit(12, 'old')

  ! *** check if solve was successful
  update: IF (gesv_info /= 0) THEN
     WRITE (*, '(/5x,a/)') 'PDAF-ERROR(1): Problem in solve for state analysis !!!'
     flag = 1
  ELSE
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_seik_analysis:', debug, '  wbar', RiHLd

     CALL PDAF_timeit(13, 'new')

     ! **************************
     ! *** Update model state ***
     ! ***    a   f           ***
     ! ***   x = x + L RiHLd  ***
     ! **************************

     CALL gemvTYPE('n', dim_p, rank, 1.0, ens_p, &
          dim_p, RiHLd, 1, 0.0, state_inc_p, 1)
     DEALLOCATE(RiHLd)
    
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

END SUBROUTINE PDAFseik_analysis



!> Perform ensemble transformation in SEIK
!!
!! Routine for ensemble transformation in the SEIK filter.
!! The routine transforms a forecast ensemble to represent
!! the analysis state und the analysis covariance
!! matrix given in factored form P = L U L$^T$.
!!
!! Variant for domain decomposition.
!! Old formulation regarding application of T.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2003-10 - Lars Nerger - Initial code
!! * Later revisions - see svn log
!!
SUBROUTINE PDAFseik_resample(subtype, dim_p, dim_ens, rank, Uinv, &
     state_p, ensT_p, type_sqrt, type_trans, Nm1vsN, screen, flag)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: mype
  USE PDAF_mod_filter, &
       ONLY: debug
  USE PDAF_analysis_utils, &
       ONLY: PDAF_seik_Omega

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: subtype      !< Filter subtype
  INTEGER, INTENT(in) :: dim_p        !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens      !< Size of ensemble
  INTEGER, INTENT(in) :: rank         !< Rank of initial covariance matrix
  REAL, INTENT(inout) :: Uinv(rank, rank)       !< Inverse of matrix U
  REAL, INTENT(inout) :: state_p(dim_p)         !< PE-local model state
  REAL, INTENT(inout) :: ensT_p(dim_p, dim_ens) !< PE-local ensemble times T
  INTEGER, INTENT(in) :: type_sqrt    !< Type of square-root of A
                                      !< (0): symmetric sqrt; (1): Cholesky decomposition
  INTEGER, INTENT(in) :: type_trans   !< Type of ensemble transformation
  INTEGER, INTENT(in) :: Nm1vsN       !< Flag which normalization of P ist used in SEIK
  INTEGER, INTENT(in) :: screen       !< Verbosity flag
  INTEGER, INTENT(inout) :: flag      !< Status flag
       
! *** local variables ***
  INTEGER :: i, j, row, col           ! Counters
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  INTEGER :: lib_info                 ! Status flags for library calls
  INTEGER :: ldwork                   ! Size of work array for SYEV
  INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
  REAL :: fac                         ! Temporary variable: sqrt(dim_ens) or sqrt(rank)
  REAL :: rdim_ens                    ! Size of ensemble in real format
  LOGICAL :: storeOmega = .FALSE.     ! Store matrix Omega or recompute it
  LOGICAL, SAVE :: firsttime = .TRUE. ! Indicates first call to resampling
  REAL, ALLOCATABLE :: Omega(:,:)     ! Orthogonal matrix Omega
  REAL, ALLOCATABLE :: OmegaT(:,:)    ! Transpose of Omega
  REAL, SAVE, ALLOCATABLE :: OmegaTsave(:,:) ! Saved transpose of Omega
  REAL, ALLOCATABLE :: ens_blk(:,:)          ! Temporary blocked state ensemble
  REAL, ALLOCATABLE :: tempUinv(:,:)         ! Temporary matrix Uinv
  REAL, ALLOCATABLE :: Ttrans(:,:)           ! Temporary matrix T^T
  REAL, ALLOCATABLE :: svals(:)   ! Singular values of Uinv
  REAL, ALLOCATABLE :: work(:)    ! Work array for SYEV
  REAL, ALLOCATABLE :: Usqrt(:,:) ! Temporary for square-root of U


! **********************
! *** INITIALIZATION ***
! **********************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_seik_resample -- START'

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a)') 'PDAF', 'Transform state ensemble'
     IF (type_sqrt == 1) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- use Cholesky square-root of U'
     ELSE
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- use symmetric square-root of U'
     END IF
  END IF

  CALL PDAF_timeit(20, 'new')
  CALL PDAF_timeit(32, 'new')

! ************************************
! *** Compute square-root of U     ***
! ************************************

  ! initialize Uinv for internal use
  ALLOCATE(tempUinv(rank, rank))
  IF (allocflag == 0) CALL PDAF_memcount(4, 'r', rank ** 2)
  IF (subtype /= 3) THEN
     tempUinv(:,:) = Uinv(:,:)
  ELSE
     rdim_ens = REAL(dim_ens)

     ! Initialize matrix T^T
     ALLOCATE(Ttrans(rank, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(4, 'r', dim_ens * rank)
     DO i = 1, rank
        DO j = 1, dim_ens
           Ttrans(i, j) = -1.0 / rdim_ens
        END DO
     END DO
     DO i = 1, rank
        Ttrans(i, i) = Ttrans(i, i) + 1.0
     END DO
     CALL gemmTYPE('n', 't', rank, rank, dim_ens, &
          rdim_ens, Ttrans, rank, Ttrans, rank, &
          0.0, tempUinv, rank)
     DEALLOCATE(Ttrans)
  END IF

  ! Usqrt is allocated with dim_ens cols, because this is 
  ! required further below. Now only rank columns are used
  ALLOCATE(Usqrt(rank, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens * rank)

  typesqrtU: IF (type_sqrt == 1) THEN
     ! Compute square-root by Cholesky-decomposition

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_seik_resample:', debug, &
          '  Compute Cholesky decomposition of U^-1'

     CALL potrfTYPE('l', rank, tempUinv, rank, lib_info)

  ELSE
     ! Compute symmetric square-root by SVD of Uinv

     ALLOCATE(svals(rank))
     ALLOCATE(work(3 * rank))
     ldwork = 3 * rank
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3 * rank)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_seik_resample:', debug, &
          '  Compute eigenvalue decomposition of U^-1'

     ! Compute SVD of Uinv
     CALL syevTYPE('v', 'l', rank, Uinv, rank, svals, work, ldwork, lib_info)

     DEALLOCATE(work)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_seik_resample:', debug, '  eigenvalues', svals

     DO col = 1, rank
        DO row = 1, rank
           Usqrt(row, col) = Uinv(row, col) / SQRT(svals(col))
        END DO
     END DO

     CALL gemmTYPE('n', 't', rank, rank, rank, &
          1.0, Usqrt, rank, Uinv, rank, &
          0.0, tempUinv, rank)
     DEALLOCATE(svals)

  END IF typesqrtU

  CALL PDAF_timeit(32, 'old')


  ! check if computation of square-root was successful
  CholeskyOK: IF (lib_info == 0) THEN
     ! Decomposition OK, continue

  
! *************************************************
! *** Generate ensemble of interpolating states ***
! *************************************************

     ! allocate fields
     ALLOCATE(OmegaT(rank, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(4, 'r', dim_ens * rank)
     
     IF (storeOmega .AND. allocflag == 0) THEN
        ALLOCATE(OmegaTsave(rank, dim_ens))
        CALL PDAF_memcount(4, 'r', dim_ens * rank)
     END IF

     CALL PDAF_timeit(33, 'new')
     Omega_store: IF (storeOmega) THEN

        first: IF (firsttime) THEN
           ! *** At first call to SEIK_RESAMPLE initialize   ***
           ! *** the matrix Omega in SEIK_Omega and store it ***

           ALLOCATE(Omega(dim_ens, rank))
           IF (allocflag == 0) CALL PDAF_memcount(4, 'r', dim_ens * rank)

           ! *** Generate uniform orthogonal matrix OMEGA ***
           CALL PDAF_seik_Omega(rank, Omega, type_trans, screen)
        
           ! transpose Omega
           IF (type_sqrt == 1) THEN
              OmegaT = TRANSPOSE(Omega)
              ! store transposed Omega
              OmegaTsave = OmegaT
           ELSE
              Usqrt = TRANSPOSE(Omega)
              ! store transposed Omega
              OmegaTsave = Usqrt
           END IF

           firsttime = .FALSE.
      
           DEALLOCATE(Omega)

        ELSE first

           IF (mype == 0 .AND. screen > 0) &
                WRITE (*, '(a, 5x, a)') 'PDAF', '--- use stored Omega'
           IF (type_sqrt == 1) THEN
              OmegaT = OmegaTsave
           ELSE
              Usqrt = OmegaTsave
           END IF

        END IF first

     ELSE Omega_store

        ! *** Initialize the matrix Omega in SEIK_Omega ***
        ! *** each time SEIK_RESAMPLE is called         ***

        ALLOCATE(Omega(dim_ens, rank))
        IF (allocflag == 0) CALL PDAF_memcount(4, 'r', dim_ens * rank)

        ! *** Generate uniform orthogonal matrix OMEGA ***
        CALL PDAF_seik_Omega(rank, Omega, type_trans, screen)
        
        ! transpose Omega
        IF (type_sqrt == 1) THEN
           OmegaT = TRANSPOSE(Omega)
        ELSE
           Usqrt = TRANSPOSE(Omega)
        END IF

        DEALLOCATE(Omega)

     END IF Omega_store
     IF (debug>0) THEN
        IF (type_sqrt == 1) THEN
           WRITE (*,*) '++ PDAF-debug PDAF_seik_update:', debug, '  Omega^T', OmegaT
        ELSE
           WRITE (*,*) '++ PDAF-debug PDAF_seik_update:', debug, '  Omega^T', Usqrt
        END IF
     END IF

     CALL PDAF_timeit(33, 'old')


     ! ***     Generate ensemble of states                             ***
     ! *** x_i = x + sqrt(FAC) L (Omega C^(-1))t                       ***
     ! *** Here FAC depends on the use definition of the covariance    ***
     ! *** matrix P using a factor N^-1 or (N-1)^-1.                   ***
    
     CALL PDAF_timeit(34, 'new')
     IF (type_sqrt == 1) THEN
        ! A = (Omega C^(-1)) by solving Ct A = OmegaT for A
        CALL trtrsTYPE('L', 'T', 'N', rank, dim_ens, &
          tempUinv, rank, OmegaT, rank, lib_info)
     ELSE
        ! TMP_UINV already contains matrix C (no more inversion)

        CALL gemmTYPE('n', 'n', rank, dim_ens, rank, &
             1.0, tempUinv, rank, Usqrt, rank, &
             0.0, OmegaT, rank)

        lib_info = 0

     END IF
     CALL PDAF_timeit(34, 'old')
     CALL PDAF_timeit(20, 'old')
     
     ! check if solve was successful
     solveOK: IF (lib_info == 0) THEN
        ! Solve for A OK, continue

        CALL PDAF_timeit(21, 'new')

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_seik_resample:', debug, '  transform', OmegaT

        ! *** Block formulation for resampling
        maxblksize = 200
        IF (mype == 0 .AND. screen > 0) &
             WRITE (*, '(a, 5x, a, i5)') 'PDAF', '--- use blocking with size ', maxblksize
        
        ALLOCATE(ens_blk(maxblksize, dim_ens))
        IF (allocflag == 0) CALL PDAF_memcount(4, 'r', maxblksize * dim_ens)

        blocking: DO blklower = 1, dim_p, maxblksize
           
           blkupper = MIN(blklower + maxblksize - 1, dim_p)

           ! Store old state ensemble
           DO col = 1, dim_ens
              ens_blk(1 : blkupper - blklower + 1, col) &
                   = ensT_p(blklower : blkupper, col)
           END DO
           
           DO col = 1,dim_ens
              ensT_p(blklower : blkupper, col) = state_p(blklower : blkupper)
           END DO

           ! *** X = state+ sqrt(FAC) state_ens T A^T (A^T stored in OmegaT) ***
           IF (Nm1vsN == 1) THEN
              ! Use factor (N-1)^-1
              fac = SQRT(REAL(dim_ens - 1))
           ELSE
              ! Use factor N^-1
              fac = SQRT(REAL(dim_ens))
           END IF

           CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, rank, &
                fac, ens_blk(1, 1), maxblksize, OmegaT(1, 1), rank, &
                1.0, ensT_p(blklower, 1), dim_p)

        END DO blocking

        CALL PDAF_timeit(21, 'old')

        DEALLOCATE(ens_blk)

     ELSE SolveOK

        ! Solve for A failed
        WRITE (*, '(/5x, a/)') &
             'PDAF-ERROR(2): Problem with solve for A in SEIK_RESAMPLE !!!'
        flag = 2

     ENDIF SolveOK

     DEALLOCATE(OmegaT)

  ELSE CholeskyOK

     ! eigendecomposition failed
     WRITE (*, '(/5x, a/)') &
          'PDAF-ERROR(1): Problem with Cholesky decomposition of Uinv !!!'
     flag = 1
     
  ENDIF CholeskyOK


! ****************
! *** clean up ***
! ****************

  DEALLOCATE(Usqrt)

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_seik_resample -- END'

END SUBROUTINE PDAFseik_resample

END MODULE PDAF_seik_analysis
