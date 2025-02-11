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
!$Id: PDAF_letkf_analysis_fixed.F90 1147 2023-03-12 16:14:34Z lnerger $
!BOP
!
! !ROUTINE: PDAF_letkf_analysis_fixed --- LETKF with state update/no transform
!
! !INTERFACE:
SUBROUTINE PDAF_letkf_analysis_fixed(domain_p, step, dim_l, dim_obs_l, &
     dim_ens, state_l, Ainv_l, ens_l, HZ_l, &
     HXbar_l, obs_l, state_inc_l, forget, &
     U_prodRinvA_l, &
     incremental, screen, debug, flag)

! !DESCRIPTION:
! Analysis step of the ESTKF with direct update
! of the state estimate, but no ensemble 
! transformation. The ensemble is only shifted
! to represent the analysis state. This variant
! is used for the filter variant with a fixed
! covariance matrix. This variant bases on the
! LETKF implementation suing the T-matrix.
!
! The implementation also supports an adaptive forgetting factor.
!
! Variant for domain decomposed states.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! __Revision history:__
! 2012-03 - Lars Nerger - Initial code
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
  USE PDAF_mod_filtermpi, &
       ONLY: mype
  USE PDAF_analysis_utils, &
       ONLY: PDAF_subtract_rowmean, PDAF_subtract_colmean
#if defined (_OPENMP)
  USE omp_lib, &
       ONLY: omp_get_num_threads, omp_get_thread_num
#endif

  IMPLICIT NONE

! *** Arguments ***
!  Variable naming scheme:
!    suffix _p: Denotes a full variable on the PE-local domain
!    suffix _l: Denotes a local variable on the current analysis domain
  INTEGER, INTENT(in) :: domain_p      !< Current local analysis domain
  INTEGER, INTENT(in) :: step          !< Current time step
  INTEGER, INTENT(in) :: dim_l         !< State dimension on local analysis domain
  INTEGER, INTENT(in) :: dim_obs_l     !< Size of obs. vector on local ana. domain
  INTEGER, INTENT(in) :: dim_ens       !< Size of ensemble 
  REAL, INTENT(inout) :: state_l(dim_l)           !< Local forecast state
  REAL, INTENT(out)   :: Ainv_l(dim_ens, dim_ens) !< on exit: local weight matrix for ensemble transformation
  REAL, INTENT(inout) :: ens_l(dim_l, dim_ens)    !< Local state ensemble
  REAL, INTENT(inout) :: HZ_l(dim_obs_l, dim_ens) !< Local observed state ensemble (perturbation)
  REAL, INTENT(in)    :: HXbar_l(dim_obs_l)       !< Local observed ensemble mean
  REAL, INTENT(in)    :: obs_l(dim_obs_l)         !< Local observation vector
  REAL, INTENT(in)    :: state_inc_l(dim_l)       !< Local state increment
  REAL, INTENT(inout) :: forget        !< Forgetting factor
  INTEGER, INTENT(in) :: incremental   !< Control incremental updating
  INTEGER, INTENT(in) :: screen        !< Verbosity flag
  INTEGER, INTENT(in) :: debug         !< Flag for writing debug output
  INTEGER, INTENT(inout) :: flag       !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA_l            !< Provide product R^-1 A for local analysis domain
       
! *** local variables ***
  INTEGER :: i, col, row               ! Counters
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  INTEGER :: syev_info                 ! Status flag for SYEV
  INTEGER :: ldwork                    ! Size of work array for SYEV
  INTEGER, SAVE :: lastdomain = -1     ! store domain index
  LOGICAL, SAVE :: screenout = .true.  ! Whether to print information to stdout
  REAL, ALLOCATABLE :: RiHZ_l(:,:)     ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: innov_l(:)      ! local observation innovation
  REAL, ALLOCATABLE :: RiHZd_l(:)      ! local RiHZd
  REAL, ALLOCATABLE :: VRiHZd_l(:)     ! Temporary vector for analysis
  REAL, ALLOCATABLE :: tmp_Ainv_l(:,:) ! Temporary storage of Ainv
  REAL, ALLOCATABLE :: svals(:)        ! Singular values of Ainv
  REAL, ALLOCATABLE :: work(:)         ! Work array for SYEV
  INTEGER, SAVE :: mythread, nthreads  ! Thread variables for OpenMP
  INTEGER :: screen_dummy              ! Dummy variable to avoid compiler warning

!$OMP THREADPRIVATE(mythread, nthreads, lastdomain, allocflag, screenout)


! *******************
! *** Preparation ***
! *******************

  CALL PDAF_timeit(51, 'new')

  ! Initialize variable to prevent compiler warning
  screen_dummy = screen

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
     IF (mype == 0 .AND. screen > 0 .AND. screenout) THEN
        WRITE (*,'(a, 5x, a, i5, a)') &
             'PDAF', '--- Use OpenMP parallelization with ', nthreads, ' threads'
     END IF
#endif
  END IF

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_analysis -- START'

  CALL PDAF_timeit(51, 'old')


! **************************
! *** Compute innovation ***
! ***   d = y - H x      ***
! **************************

  CALL PDAF_timeit(16, 'new')
  CALL PDAF_timeit(20, 'new')

  haveobsB: IF (dim_obs_l > 0) THEN
     ! *** The residual only exists for domains with observations ***

     CALL PDAF_timeit(51, 'new')

     ALLOCATE(innov_l(dim_obs_l))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l)

     innov_l = obs_l - HXbar_l

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_letkf_analysis:', debug, '  innovation d_l', innov_l

     CALL PDAF_timeit(51, 'old')

  END IF haveobsB

  CALL PDAF_timeit(20, 'old')


! **********************************************
! ***   Compute analyses matrix Ainv         ***
! ***                                        ***
! ***     -1                 T  -1           ***
! ***    A  = forget I + (HZ)  R   HZ        ***
! **********************************************

  haveobsA: IF (dim_obs_l > 0) THEN
     ! *** The contribution of observation matrix ist only ***
     ! *** computed for domains with observations          ***

     CALL PDAF_timeit(21, 'new')
     CALL PDAF_timeit(51, 'new')

     ! Subtract ensemble mean: HZ = [Hx_1 ... Hx_N] T
     CALL PDAF_subtract_rowmean(dim_obs_l, dim_ens, HZ_l)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_letkf_analysis:', debug, '  HXT_l', HZ_l(:, 1:dim_ens-1)

     CALL PDAF_timeit(51, 'old')
     CALL PDAF_timeit(21, 'old')


     ! ***                RiHZ = Rinv HZ                
     ! *** This is implemented as a subroutine thus that
     ! *** Rinv does not need to be allocated explicitly.
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_analysis -- call prodRinvA_l'

     ALLOCATE(RiHZ_l(dim_obs_l, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l * dim_ens)

     CALL PDAF_timeit(48, 'new')
     CALL U_prodRinvA_l(domain_p, step, dim_obs_l, dim_ens, obs_l, HZ_l, RiHZ_l)
     CALL PDAF_timeit(48, 'old')

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_letkf_analysis:', debug, '  R^-1(HXT_l)', RiHZ_l

     CALL PDAF_timeit(21, 'new')
     CALL PDAF_timeit(51, 'new')

     ! *** Initialize Ainv = (N-1) I ***
     Ainv_l = 0.0
     DO i = 1, dim_ens
        Ainv_l(i, i) = REAL(dim_ens - 1)
     END DO

     ! ***             T        ***
     ! ***  Compute  HZ  RiHZ   ***

     ALLOCATE(tmp_Ainv_l(dim_ens, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

     CALL gemmTYPE('t', 'n', dim_ens, dim_ens, dim_obs_l, &
          1.0, HZ_l, dim_obs_l, RiHZ_l, dim_obs_l, &
          0.0, tmp_Ainv_l, dim_ens)

     CALL PDAF_timeit(51, 'old')

  ELSE haveobsA
     ! *** For domains with dim_obs_l=0 there is no ***
     ! *** direct observation-contribution to Ainv  ***

     CALL PDAF_timeit(21, 'new')
     CALL PDAF_timeit(51, 'new')

     ! *** Initialize Ainv = (N-1) I ***
     Ainv_l = 0.0
     DO i = 1, dim_ens
        Ainv_l(i, i) = REAL(dim_ens - 1)
     END DO

     ! No observation-contribution to Ainv from this domain
     ALLOCATE(tmp_Ainv_l(dim_ens, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

     tmp_Ainv_l = 0.0

     CALL PDAF_timeit(51, 'old')

  END IF haveobsA

  ! *** Complete computation of Ainv  ***
  ! ***   -1          -1    T         ***
  ! ***  A  = forget A  + HZ RiHZ     ***
  CALL PDAF_timeit(51, 'new')
  Ainv_l = forget * Ainv_l + tmp_Ainv_l
  CALL PDAF_timeit(51, 'old')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_letkf_analysis:', debug, '  A^-1_l', Ainv_l

  CALL PDAF_timeit(21, 'old')
  CALL PDAF_timeit(16, 'old')


! ***********************************************
! *** Compute weight for model state update   ***
! ***                                         ***
! ***              T                    f     ***
! ***    w = A RiHZ d  with d = (y - H x )    ***
! ***                                         ***
! ***********************************************

  CALL PDAF_timeit(17, 'new')
  CALL PDAF_timeit(51, 'new')

  ! *** Compute RiHZd = RiHZ^T d ***
  ALLOCATE(RiHZd_l(dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)

  haveobsC: IF (dim_obs_l > 0) THEN
     ! *** RiHLd_p/=0 only with observations ***

     CALL gemvTYPE('t', dim_obs_l, dim_ens, 1.0, RiHZ_l, &
          dim_obs_l, innov_l, 1, 0.0, RiHZd_l, 1)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_letkf_analysis:', debug, '  (HXT_l R^-1)^T d_l', RiHZd_l

     DEALLOCATE(RiHZ_l, innov_l)

  ELSE haveobsC

     RiHZd_l = 0.0

  END IF haveobsC


  ! *** Compute weight vector for state analysis:        ***
  ! ***          w = A RiHZd                             ***
  ! *** Use singular value decomposition of Ainv         ***
  ! ***        Ainv = USV^T                              ***
  ! *** Then: A = U S^(-1) V                             ***
  ! *** The decomposition is also used for the symmetric ***
  ! *** square-root for the ensemble transformation.     ***

  ! *** Invert Ainv using SVD
  ALLOCATE(svals(dim_ens))
  ALLOCATE(work(3 * dim_ens))
  ldwork = 3 * dim_ens
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3 * dim_ens)

  ! save matrix Ainv
  tmp_Ainv_l = Ainv_l

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_letkf_analysis:', debug, &
       '  Compute eigenvalue decomposition of A^-1_l'

  ! Compute SVD of Ainv
  CALL syevTYPE('v', 'l', dim_ens, tmp_Ainv_l, dim_ens, svals, work, ldwork, syev_info)

  DEALLOCATE(work)

  ! *** check if SVD was successful
  IF (syev_info == 0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_letkf_analysis:', debug, '  eigenvalues', svals

     flag = 0
  ELSE
     WRITE (*, '(/5x, a, i10, a/)') &
          'PDAF-ERROR(1): Domain ', domain_p, ' Problem in EVD of inverse of A !!!'
     flag = 1
  END IF

  ! *** Compute w = A RiHZd stored in RiHZd
  check0: IF (flag == 0) THEN

     ALLOCATE(VRiHZd_l(dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)

     CALL gemvTYPE('t', dim_ens, dim_ens, 1.0, tmp_Ainv_l, &
          dim_ens, RiHZd_l, 1, 0.0, VRiHZd_l, 1)
     
     DO row = 1, dim_ens
        VRiHZd_l(row) = VRiHZd_l(row) / svals(row)
     END DO
  
     CALL gemvTYPE('n', dim_ens, dim_ens, 1.0, tmp_Ainv_l, &
          dim_ens, VRiHZd_l, 1, 0.0, RiHZd_l, 1)

     DEALLOCATE(svals, tmp_Ainv_l, VRiHZd_l)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_letkf_analysis:', debug, '  A(HXT_l R^-1)^T d_l', RiHZd_l

  END IF check0


! ************************************
! ***      update model state      ***
! ***                              ***
! ***     a   f   f                ***
! ***    x = x + X  Omega RiHLd    ***
! ***                              ***
! ************************************

  check1: IF (flag == 0) THEN

     ! **************************
     ! *** Compute vector T w ***
     ! **************************

     CALL PDAF_subtract_colmean(dim_ens, 1, RiHZd_l)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_letkf_analysis:', debug, '  transform wbar_l', RiHZd_l

     CALL PDAF_timeit(17, 'old')


     ! *****************************
     ! *** Update state estimate ***
     ! *****************************

     CALL PDAF_timeit(18, 'new')

     CALL gemvTYPE('n', dim_l, dim_ens, 1.0, ens_l, &
          dim_l, RiHZd_l, 1, 0.0, state_inc_l, 1)
     DEALLOCATE(RiHZd_l)
     
     ! Shift ensemble
     DO col = 1, dim_ens
        DO row = 1, dim_l
           ens_l(row, col) = ens_l(row, col) + state_inc_l(row)
        END DO
     END DO

     IF (incremental == 0) THEN
        ! update state here if incremental updating is not used
        state_l = state_l + state_inc_l
     END IF

     CALL PDAF_timeit(18, 'old')

  ELSE check1

     CALL PDAF_timeit(17, 'old')

  END IF check1

  CALL PDAF_timeit(51, 'old')


! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1

  ! Store domain index
  lastdomain = domain_p

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_analysis -- END'

END SUBROUTINE PDAF_letkf_analysis_fixed
