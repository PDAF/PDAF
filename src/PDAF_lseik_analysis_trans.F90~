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
!$Id: PDAF_lseik_analysis_trans.F90 1147 2023-03-12 16:14:34Z lnerger $
!BOP
!
! !ROUTINE: PDAF_lseik_analysis_trans --- LSEIK analysis/transformation
!
! !INTERFACE:
SUBROUTINE PDAF_lseik_analysis_trans(domain_p, step, dim_l, dim_obs_f, dim_obs_l, &
     dim_ens, rank, state_l, Uinv_l, ens_l, HX_f, &
     HXbar_f, state_inc_l, OmegaT_in, forget, U_g2l_obs, &
     U_init_obs_l, U_prodRinvA_l, U_init_obsvar_l, U_init_n_domains_p, &
     screen, incremental, type_forget, type_sqrt, flag)

! !DESCRIPTION:
! Analysis step of the LSEIK filter with direct
! transformation of the forecast into the 
! analysis ensemble. This variant does not
! compute the analysis state, but only the
! analysis ensemble, whose mean is the analysis
! state.
! Supported is also the adaptive forgetting factor.
!
! Variant for domain decomposed states.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2009-07 - Lars Nerger - Initial code
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
       ONLY: Nm1vsN, obs_member, debug
  USE PDAF_mod_filtermpi, &
       ONLY: mype
  USE PDAFomi, &
       ONLY: omi_omit_obs => omit_obs
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
  REAL, INTENT(inout) :: Uinv_l(rank, rank)    ! Inverse of matrix U - temporary use only
  REAL, INTENT(inout) :: ens_l(dim_l, dim_ens) ! Local state ensemble
  REAL, INTENT(in) :: HX_f(dim_obs_f, dim_ens) ! PE-local full observed state ens.
  REAL, INTENT(in) :: HXbar_f(dim_obs_f)       ! PE-local full observed ens. mean
  REAL, INTENT(in) :: state_inc_l(dim_l)       ! Local state increment
  REAL, INTENT(inout) :: OmegaT_in(rank, dim_ens) ! Matrix Omega
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
! Calls: syevTYPE (LAPACK; dsyev or ssyev dependent on precision)
! Calls: potrfTYPE (LAPACK; dpotrf or spotrf dependent on precision)
! Calls: trtrsTYPE (LAPACK; dtrtrs or strtrs dependent on precision)
!EOP
       
! *** local variables ***
  INTEGER :: i, j, member, col, row    ! Counters
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  INTEGER :: lib_info                  ! Status flag for LAPACK calls
  INTEGER :: ldwork                    ! Size of work array for SYEV
  INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
  REAL    :: fac                       ! Temporary variable sqrt(dim_ens) or sqrt(rank)
  INTEGER, SAVE :: lastdomain = -1     ! store domain index
  LOGICAL, SAVE :: screenout = .true.  ! Whether to print information to stdout
  REAL, ALLOCATABLE :: HL_l(:,:)       ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: RiHL_l(:,:)     ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: resid_l(:)      ! observation residual
  REAL, ALLOCATABLE :: obs_l(:)        ! local observation vector
  REAL, ALLOCATABLE :: HXbar_l(:)      ! state projected onto obs. space
  REAL, ALLOCATABLE :: RiHLd_l(:)      ! local RiHLd
  REAL, ALLOCATABLE :: VRiHLd_l(:)     ! Temporary vector for analysis
  REAL, ALLOCATABLE :: tmp_Uinv_l(:,:) ! Temporary storage of Uinv
  REAL, ALLOCATABLE :: OmegaT(:,:)     ! Transpose of Omega
  REAL, ALLOCATABLE :: TA(:,:)         ! Temporary matrix
  REAL, ALLOCATABLE :: ens_blk(:,:)    ! Temporary blocked state ensemble
  REAL, ALLOCATABLE :: svals(:)        ! Singular values of Uinv
  REAL, ALLOCATABLE :: work(:)         ! Work array for SYEV
  INTEGER, ALLOCATABLE :: ipiv(:)      ! vector of pivot indices for GESV
  INTEGER, SAVE :: mythread, nthreads  ! Thread variables for OpenMP
  INTEGER :: incremental_dummy         ! Dummy variable to avoid compiler warning
  REAL :: state_inc_l_dummy(1)         ! Dummy variable to avoid compiler warning

!$OMP THREADPRIVATE(mythread, nthreads, lastdomain, allocflag, screenout)


! *******************
! *** Preparation ***
! *******************

  CALL PDAF_timeit(51, 'new')

  ! Initialize variable to prevent compiler warning
  incremental_dummy = incremental
  state_inc_l_dummy(1) = state_inc_l(1)

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
     IF (screenout) THEN
        WRITE (*,'(a, 5x,a,i5,a)') &
             'PDAF', '--- Use OpenMP parallelization with ', nthreads, ' threads'
     END IF
#endif
  END IF

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lseik_analysis -- START'

  CALL PDAF_timeit(51, 'old')


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
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lseik_analysis -- call g2l_obs for mean'

     CALL PDAF_timeit(46, 'new')
     obs_member = 0
     CALL U_g2l_obs(domain_p, step, dim_obs_f, dim_obs_l, HXbar_f, HXbar_l)
     CALL PDAF_timeit(46, 'old')

     ! get local observation vector
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lseik_analysis -- call init_obs_l'

     CALL PDAF_timeit(47, 'new')
     CALL U_init_obs_l(domain_p, step, dim_obs_l, obs_l)
     CALL PDAF_timeit(47, 'old')

     ! Get residual as difference of observation and observed state
     CALL PDAF_timeit(51, 'new')
     resid_l = obs_l - HXbar_l
     CALL PDAF_timeit(51, 'old')

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lseik_analysis:', debug, '  innovation d_l', resid_l

     ! Omit observations with too high innovation
     IF (omi_omit_obs) CALL PDAFomi_omit_by_inno_l_cb(domain_p, dim_obs_l, resid_l, obs_l)

  END IF haveobsB

  CALL PDAF_timeit(12, 'old')


! *************************************************
! ***   Compute analyzed matrix Uinv            ***
! ***                                           ***
! ***  -1              T        T  -1           ***
! *** U  = forget fac T T + (HL)  R  (HL)       ***
! ***  i                        i  i     i      ***
! ***                                           ***
! *** Here FAC is a scaling factor according    ***
! *** to the definition of the ensemble         ***
! *** covariance scaled by N^-1 (original SEIK) ***
! *** (N-1)^-1 (SEIK as ensemble KF)            ***
! *************************************************

  CALL PDAF_timeit(10, 'new')

  haveobsA: IF (dim_obs_l > 0) THEN
     ! *** The contribution of observation matrix ist only ***
     ! *** computed for domains with observations          ***

     CALL PDAF_timeit(30, 'new')

     ! *** Compute HL = [Hx_1 ... Hx_N] T
     ALLOCATE(HL_l(dim_obs_l, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l * dim_ens)

     CALL PDAF_timeit(46, 'new')

     IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lseik_analysis -- call g2l_obs', dim_ens, 'times'

     ENS: DO member = 1, dim_ens
        ! Store member index
        obs_member = member

        ! [Hx_1 ... Hx_(r+1)] for local analysis domain
        CALL U_g2l_obs(domain_p, step, dim_obs_f, dim_obs_l, HX_f(:, member), &
             HL_l(:, member))
     END DO ENS

     CALL PDAF_timeit(46, 'old')

     ! *** Set the value of the forgetting factor  ***
     ! *** Inserted here, because HL_l is required ***
     CALL PDAF_timeit(51, 'new')
     IF (type_forget == 2) THEN
        CALL PDAF_set_forget_local(domain_p, step, dim_obs_l, dim_ens, HL_l, &
             HXbar_l, resid_l, obs_l, U_init_n_domains_p, U_init_obsvar_l, &
             forget)
     ENDIF
     DEALLOCATE(HXbar_l)

     ! Complete HL = [Hx_1 ... Hx_(r+1)] T
     CALL PDAF_seik_matrixT(dim_obs_l, dim_ens, HL_l)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lseik_analysis:', debug, '  HXT_l', HL_l(:, 1:dim_ens-1)

     CALL PDAF_timeit(51, 'old')
     CALL PDAF_timeit(30, 'old')
     CALL PDAF_timeit(31, 'new')


     ! ***                RiHL = Rinv HL                 ***
     ! *** this is implemented as a subroutine thus that ***
     ! *** Rinv does not need to be allocated explicitly ***
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lseik_analysis -- call prodRinvA_l'

     ALLOCATE(RiHL_l(dim_obs_l, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l * rank)

     CALL PDAF_timeit(48, 'new')
     CALL U_prodRinvA_l(domain_p, step, dim_obs_l, rank, obs_l, HL_l, RiHL_l)
     CALL PDAF_timeit(48, 'old')

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lseik_analysis:', debug, '  R^-1(HXT_l)', RiHL_l

     DEALLOCATE(obs_l)

     CALL PDAF_timeit(51, 'new')

     ! *** Initialize Uinv = fac T^T T ***
     CALL PDAF_seik_Uinv(rank, Uinv_l)

     ! ***             T        ***
     ! ***  Compute  HL  RiHL   ***
     ALLOCATE(tmp_Uinv_l(rank, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank**2)

     CALL gemmTYPE('t', 'n', rank, rank, dim_obs_l, &
          1.0, HL_l, dim_obs_l, RiHL_l, dim_obs_l, &
          0.0, tmp_Uinv_l, rank)

     DEALLOCATE(HL_l)

     CALL PDAF_timeit(51, 'old')

  ELSE haveobsA
     ! *** For domains with dim_obs_l=0 there is no ***
     ! *** direct observation-contribution to Uinv  ***

     CALL PDAF_timeit(31, 'new')
     CALL PDAF_timeit(51, 'new')

     ! Initialize Uinv = fac T^T T 
     CALL PDAF_seik_Uinv(rank, Uinv_l)

     ! No observation-contribution to Uinv from this domain
     ALLOCATE(tmp_Uinv_l(rank, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank**2)

     tmp_Uinv_l = 0.0

     CALL PDAF_timeit(51, 'old')

  END IF haveobsA

  ! *** Complete computation of Uinv  ***
  ! ***   -1          -1    T         ***
  ! ***  U  = forget U  + HL RiHL     ***
  CALL PDAF_timeit(51, 'new')
  Uinv_l = forget * Uinv_l + tmp_Uinv_l
  CALL PDAF_timeit(51, 'old')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lseik_analysis:', debug, '  U^-1_l', Uinv_l

  CALL PDAF_timeit(31, 'old')
  CALL PDAF_timeit(10, 'old')


! ***********************************************
! *** Compute weight for model state update   ***
! ***                                         ***
! ***              T                    f     ***
! ***    w = U RiHL d  with d = (y - H x )    ***
! ***********************************************

  CALL PDAF_timeit(13, 'new')
  CALL PDAF_timeit(51, 'new')

  ! *** Computer RiHLd = RiHV^T d ***
  ALLOCATE(RiHLd_l(rank))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank)

  haveobsC: IF (dim_obs_l > 0) THEN
     ! *** RiHLd_p/=0 only with observations ***

     CALL gemvTYPE('t', dim_obs_l, rank, 1.0, RiHL_l, &
          dim_obs_l, resid_l, 1, 0.0, RiHLd_l, 1)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lseik_analysis:', debug, '  (HXT_l R^-1)^T d_l', RiHLd_l

     DEALLOCATE(RiHL_l, resid_l)

  ELSE haveobsC

     RiHLd_l = 0.0

  END IF haveobsC


  ! *** Compute weight vector for state analysis:     ***
  ! ***          w = U RiHLd                          ***
  ! *** For this, two variants are implemented:       ***
  ! *** 1. solve for w in:                            ***
  ! ***           -1                                  ***
  ! ***          U  w = RiHLd                         ***
  ! ***   We use the LAPACK routine GESV              ***
  ! *** 2. Compute singular value decomposition       ***
  ! ***   of Uinv: Uinv = ASB^T                       ***
  ! ***   Then: U = A S^(-1) B                        ***
  ! ***   This is combined with a symmetric           ***
  ! ***   square-root for the ensemble transformation ***

  typeuinv1: IF (type_sqrt==1) THEN
     ! *** Variant 1: Solve Uinv w= RiHLd for w

     ALLOCATE(ipiv(rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'i', rank)

     ! save matrix Uinv
     tmp_Uinv_l = Uinv_l

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lseik_analysis:', debug, &
          '  Invert U^-1_l using solver GESV'

     ! call solver (GESV - LU solver)
     CALL gesvTYPE(rank, 1, tmp_Uinv_l, rank, ipiv, &
          RiHLd_l, rank, lib_info)
     DEALLOCATE(ipiv)

  ELSE typeuinv1
     ! *** Variant 2: Invert Uinv using SVD

     ALLOCATE(svals(rank))
     ALLOCATE(work(3 * rank))
     ldwork = 3 * rank
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 4 * rank)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lseik_analysis:', debug, &
          '  Compute eigenvalue decomposition of U^-1_l'

     ! Compute SVD of Uinv
     CALL syevTYPE('v', 'l', rank, Uinv_l, rank, svals, work, ldwork, lib_info)

     DEALLOCATE(work)

     ! Compute product U RiHLd
     IF (lib_info==0) THEN
        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_lseik_analysis:', debug, '  eigenvalues', svals

        ALLOCATE(VRiHLd_l(rank))
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank)

        CALL gemvTYPE('t', rank, rank, 1.0, Uinv_l, &
             rank, RiHLd_l, 1, 0.0, VRiHLd_l, 1)
     
        DO row = 1,rank
           VRiHLd_l(row) = VRiHLd_l(row) / svals(row)
        END DO
  
        CALL gemvTYPE('n', rank, rank, 1.0, Uinv_l, &
             rank, VRiHLd_l, 1, 0.0, RiHLd_l, 1)

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_lseik_analysis:', debug, '  U(HXT_l R^-1)^T d_l', RiHLd_l

        DEALLOCATE(VRiHLd_l)
     END IF
  END IF typeuinv1

  CALL PDAF_timeit(13, 'old')

  ! *** check whether SVD was successful
  IF (lib_info == 0) THEN
     flag = 0
  ELSE
     WRITE (*, '(/5x, a, i10, a/)') &
          'PDAF-ERROR(1): Domain ', domain_p, ' Problem in computation of analysis weights !!!'
     flag = 1
  END IF


! ****************************************************************
! ***     Transform state ensemble                             ***
! ***              a   _f   f                                  ***
! ***             X  = X + X  W                                ***
! *** with                               -T      T             ***
! ***          W = T (RiHLd + sqrt(FAC) C   Omega )            ***
! *** Here FAC depends on the use definition of the covariance ***
! *** matrix P using a factor (r+1)^-1 or r^-1.                ***
! ****************************************************************

! *** Prepare weight matrix for ensemble transformation ***

  check1: IF (flag == 0) THEN

     ! allocate fields
     ALLOCATE(OmegaT(rank, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank * dim_ens)

     IF (mype == 0 .AND. screen > 0 .AND. screenout) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', 'Transform state ensemble'
        IF (type_sqrt == 1) THEN
           WRITE (*, '(a, 5x, a)') 'PDAF', '--- use Cholesky square-root of U'
        ELSE
           WRITE (*, '(a, 5x, a)') 'PDAF', '--- use symmetric square-root of U'
        END IF
     END IF

     CALL PDAF_timeit(20, 'new')
     CALL PDAF_timeit(32, 'new')

     ! Part 1: square-root of U
     typeuinv2: IF (type_sqrt == 1) THEN
        ! Variant, if Uinv has been inverted above by solving

        ! Store Uinv for temporary use
        tmp_Uinv_l(:, :) = Uinv_l(:, :)

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_lseik_analysis:', debug, &
             '  Compute Cholesky decomposition of U^-1_l'

        ! Cholesky decomposition of tmp_Uinv_l = C C^T
        CALL potrfTYPE('l', rank, tmp_Uinv_l, rank, lib_info)

        ! check if Cholesky decomposition was successful
        CholeskyOK: IF (lib_info == 0) THEN
           ! Decomposition OK, continue
           flag = 0
        ELSE
           ! Decomposition failed
           WRITE (*, '(/5x, a, i10, a/)') &
                'PDAF-ERROR(3): Problem with Cholesky decomposition of Uinv - domain ', &
                domain_p, ' !!!'
           flag = 3
        ENDIF CholeskyOK

     ELSE typeuinv2
        ! Variant, if SVD inversion of Uinv has been performed

        ! Use OmegaT as temporary array
        DO col = 1, rank
           DO row = 1, rank
              OmegaT(row, col) = Uinv_l(row, col) / SQRT(svals(col))
           END DO
        END DO

        CALL gemmTYPE('n', 't', rank, rank, rank, &
             1.0, OmegaT, rank, Uinv_l, rank, &
             0.0, tmp_Uinv_l, rank)
        DEALLOCATE(svals)

        ! Set flag
        flag = 0

     END IF typeuinv2

     CALL PDAF_timeit(32, 'old')

  END IF check1

  check2: IF (flag == 0) THEN

     ! *** Part 2: Product U^(1/2) Omega ***
    
     CALL PDAF_timeit(34, 'new')
     IF (type_sqrt == 1) THEN
        ! Initialize the matrix Omega from argument OmegaT_in
        OmegaT = OmegaT_in

       ! A = (Omega C^(-1)) by solving Ct A = OmegaT for A
        CALL trtrsTYPE('L', 'T', 'N', rank, dim_ens, &
             tmp_Uinv_l, rank, OmegaT, rank, lib_info)
     ELSE
        ! TEMP_UINV already contains matrix C (no more inversion)

        CALL gemmTYPE('n', 'n', rank, dim_ens, rank, &
             1.0, tmp_Uinv_l, rank, OmegaT_in, rank, &
             0.0, OmegaT, rank)

        lib_info = 0

     END IF

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lseik_analysis:', debug, &
          '  transform for perturbations', OmegaT

     CALL PDAF_timeit(34, 'old')

     ! check whether solve was successful
     solveOK: IF (lib_info == 0) THEN
        ! Solve for A OK, continue
        flag = 0
     ELSE
        ! Solve for A failed
        WRITE (*, '(/5x, a, i10, a/)') &
             'PDAF-ERROR(2): Problem in computation of transformation matrix - domain ', &
             domain_p, ' !!!'
        flag = 2

        CALL PDAF_timeit(20, 'old')
     END IF solveOK

  END IF check2

  check3: IF (flag == 0) THEN

     CALL PDAF_timeit(35, 'new')

     ! *** Part 4: Add RiHLd and multiply by scaling factor

     IF (Nm1vsN == 1) THEN
        ! Use factor (N-1)^-1
        fac = SQRT(REAL(dim_ens - 1))
     ELSE
        ! Use factor N^-1
        fac = SQRT(REAL(dim_ens))
     END IF

     ! *** Add RiHLd to A^T stored in OmegaT
     DO j = 1, dim_ens
        DO i = 1, rank
           OmegaT(i, j) = fac * OmegaT(i, j) + RiHLd_l(i)
        END DO
     END DO
     DEALLOCATE(RiHLd_l)
      
     ! *** T A^T (A^T stored in OmegaT_l) ***
     ALLOCATE(TA(dim_ens, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

     CALL PDAF_seik_TtimesA(rank, dim_ens, OmegaT, TA)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lseik_analysis:', debug, '  transform', TA

     CALL PDAF_timeit(35, 'old')
     CALL PDAF_timeit(20, 'old')


! *** Perform ensemble transformation ***

     ! Use block formulation for transformation
     maxblksize = 200
     IF (mype == 0 .AND. screen > 0 .AND. screenout) &
          WRITE (*, '(a, 5x, a, i5)') &
          'PDAF', '--- use blocking with size ', maxblksize

     ALLOCATE(ens_blk(maxblksize, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', maxblksize * dim_ens)

     blocking: DO blklower = 1, dim_l, maxblksize

        blkupper = MIN(blklower + maxblksize - 1, dim_l)

        ! Store old state ensemble
        CALL PDAF_timeit(21, 'new')
        DO col = 1, dim_ens
           ens_blk(1 : blkupper - blklower + 1, col) &
                = ens_l(blklower : blkupper, col)
        END DO

        DO col = 1, dim_ens
           ens_l(blklower : blkupper, col) = state_l(blklower : blkupper)
        END DO
        CALL PDAF_timeit(21, 'old')

        !                        a  _f   f    T
        ! Transform ensemble:   X = X + X  T(A )
        CALL PDAF_timeit(22, 'new')

        CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
             1.0, ens_blk, maxblksize, TA, dim_ens, &
             1.0, ens_l(blklower, 1), dim_l)

        CALL PDAF_timeit(22, 'old')

     END DO blocking
        
     DEALLOCATE(ens_blk, TA)
     DEALLOCATE(OmegaT)

  END IF check3

  CALL PDAF_timeit(51, 'old')


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(tmp_Uinv_l)

  IF (allocflag == 0) allocflag = 1

  ! Store domain index
  lastdomain = domain_p

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lseik_analysis -- END'

END SUBROUTINE PDAF_lseik_analysis_trans
