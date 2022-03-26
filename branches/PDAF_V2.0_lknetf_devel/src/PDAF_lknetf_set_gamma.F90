! Copyright (c) 2014-2021 Lars Nerger
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
!BOP
!
! !ROUTINE: PDAF_lknetf_set_gamma --- compute hybrid weight gamma for hybrid LKNETF
!
! !INTERFACE:
SUBROUTINE PDAF_lknetf_set_gamma(domain_p, dim_obs_l, dim_ens, &
     HX_l, HXbar_l, weights, type_hyb, hybrid_a_x, hnorm, &
     alpha_X, n_eff_out, skew_mabs, kurt_mabs, &
     screen, flag)

! !DESCRIPTION:
! Compute hybrid weight for LKNETF
!
! Variant for domain decomposed states.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2018-01 - Lars Nerger - 
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
  INTEGER, INTENT(in) :: dim_obs_l   ! Size of obs. vector on local ana. domain
  INTEGER, INTENT(in) :: dim_ens     ! Size of ensemble 
  REAL, INTENT(in) :: HX_l(dim_obs_l, dim_ens)  ! local observed state ens.
  REAL, INTENT(in) :: HXbar_l(dim_obs_l)        ! local mean observed ensemble
  REAL, INTENT(in) :: weights(dim_ens) ! Weight vector
  INTEGER, INTENT(in) :: type_hyb      ! Type of hybrid weight
  REAL, INTENT(in) :: hybrid_a_x       ! Prescribed hybrid weight for state transformation
  REAL, INTENT(in) :: hnorm            ! Hybrid weight for covariance transformation
  REAL, INTENT(inout) :: alpha_X(1)    ! Hybrid weight for state transformation
  REAL, INTENT(inout) :: n_eff_out(1)  ! Effective ensemble size
  REAL, INTENT(inout) :: skew_mabs(1)  ! Mean absolute skewness
  REAL, INTENT(inout) :: kurt_mabs(1)  ! Mean absolute kurtosis
  INTEGER, INTENT(in) :: screen        ! Verbosity flag
  INTEGER, INTENT(inout) :: flag       ! Status flag


! !CALLING SEQUENCE:
! Called by: PDAF_lknetf_compute_alpha
! Calls: PDAF_timeit
! Calls: PDAF_memcount
!EOP
       
! *** local variables ***
  INTEGER :: i                       ! Counter
  INTEGER, SAVE :: allocflag=0       ! Flag whether first time allocation is done
  LOGICAL :: screenout = .TRUE.      ! Whether to print information to stdout
  REAL :: n_eff                      ! Effective sample size
  INTEGER , SAVE :: lastdomain=-1    !save domain index
  INTEGER, SAVE :: mythread, nthreads  ! Thread variables for OpenMP
  REAL :: hfac
  REAL, PARAMETER :: pi=3.14159265358979
  REAL :: skew(1), kurt(1)
  REAL :: alpha_kurt, kurt_limit
  REAL :: alpha_skew, alpha_adap, alpha_stat

!$OMP THREADPRIVATE(mythread, nthreads, lastdomain, allocflag, screenout)


! *******************
! *** Preparation ***
! *******************

#if defined (_OPENMP)
  nthreads = omp_get_num_threads()
  mythread = omp_get_thread_num()
#else
  nthreads = 1
  mythread = 0
#endif

  ! Control screen output
  IF (lastdomain<domain_p .AND. lastdomain>-1) THEN
     screenout = .FALSE.
  ELSE
     screenout = .TRUE.

     ! In case of OpenMP, let only thread 0 write output to the screen
     IF (mythread>0) screenout = .FALSE.

     ! Output, only in case of OpenMP parallelization
  END IF


  ! **********************************************
  ! *** Compute particle weights as likelihood ***
  ! **********************************************

  CALL PDAF_timeit(22, 'new')
  ! Get residual as difference of observation and observed state for 
  ! each ensemble member only on domains where observations are availible
  
  ! Compute adaptive hybrid weights
  IF (type_hyb==2 .OR. type_hyb==12 .OR. type_hyb==112 .OR. type_hyb==212 &
       .OR. type_hyb==312 .OR. type_hyb==712) THEN
     CALL PDAF_lknetf_alpha_neff(dim_ens, weights, hybrid_a_x, alpha_x(1))
     alpha_adap = alpha_x(1)
  END IF

  CALL PDAF_timeit(22, 'old')


  ! Compute effective ensemble size
  CALL PDAF_diag_effsample(dim_ens, weights, n_eff)
  n_eff_out(1) = n_eff

  ! Compute ensemble statistics for observed ensemble
  DO i = 1, dim_obs_l
     CALL PDAF_diag_ensstats(dim_obs_l, dim_ens, i, &
          HXbar_l, HX_l, skew, kurt, flag)
     skew_mabs = skew_mabs + ABS(skew)
     kurt_mabs = kurt_mabs + ABS(kurt)
  END DO
  skew_mabs = skew_mabs / dim_obs_l
  kurt_mabs = kurt_mabs / dim_obs_l


  ! *********************
  ! *** Screen output ***
  ! *********************

  IF (mype == 0 .AND. screen > 0 .AND. screenout) THEN
     IF (type_hyb==0) THEN
        WRITE (*, '(a, 5x, a, f12.3)') 'PDAF', 'Set hybrid weight in LETKF to ', hybrid_a_x
     ELSEIF (type_hyb==1) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- linear dependence on N_eff/N'
     ELSEIF (type_hyb==2) THEN
        WRITE (*, '(a, 6x, a, es10.2)') 'PDAF', '--- Hybrid weight from N_eff/N>=', hybrid_a_x
     ELSEIF (type_hyb==3) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- quadratic dependence on 1 - N_eff/N'
     ELSEIF (type_hyb==4) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- square-root dependence on 1 - N_eff/N'
     ELSEIF (type_hyb==5) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- square-root dependence on N_eff/N; 1 for N_eff=N'
     ELSEIF (type_hyb==6) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- sine dependence on N_eff/N with minimum constraint'
     ELSEIF (type_hyb==7) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- square-root dependence on N_eff/N, minimum limit'
     ELSEIF (type_hyb==8) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- linear dependence 1 - 0.5 N_eff/N'
     ELSEIF (type_hyb==9) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- quadratic dependence 1 - 0.5 (N_eff/N)^2'
     ELSEIF (type_hyb==10) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - (1-limit)* skewness/sqrt(N)'
     ELSEIF (type_hyb==11) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - skewness/sqrt(N) with minimum limit'
     ELSEIF (type_hyb==12) THEN
        WRITE (*, '(a, 5x, a, f12.3)') 'PDAF', '--- dependence 1 - skewness/sqrt(N) with min. from N_eff/N>=', hybrid_a_x
     ELSEIF (type_hyb==13) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - skewness/sqrt(N) with min. from 1-N_eff/N'
     ELSEIF (type_hyb==14) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- linear dependence on (N_eff-1)/(N-1)'
     ELSEIF (type_hyb==15) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- linear dependence on (N_eff-1)/(N-1), limit 0.95'
     ELSEIF (type_hyb==110) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - (1-limit)* kurtosis/sqrt(N)'
     ELSEIF (type_hyb==111) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - kurtosis/sqrt(N) with minimum limit'
     ELSEIF (type_hyb==112) THEN
        WRITE (*, '(a, 5x, a, f12.3)') 'PDAF', '--- dependence 1 - kurtosis/sqrt(N) with min. from N_eff/N>=', hybrid_a_x
     ELSEIF (type_hyb==113) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - kurtosis/sqrt(N) with min. from 1-N_eff/N'
     ELSEIF (type_hyb==212) THEN
        WRITE (*, '(a, 5x, a, f12.3)') 'PDAF', '--- dependence 1 - kurtnorm/sqrt(N) with min. from N_eff/N>=', hybrid_a_x
     ELSEIF (type_hyb==213) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - kurtnorm/sqrt(N) with min. from 1-N_eff/N'
     ELSEIF (type_hyb==312) THEN
        WRITE (*, '(a, 5x, a, f12.3)') 'PDAF', '--- dependence 1 - ensstat/sqrt(N) with min. from N_eff/N>=', hybrid_a_x
     ELSEIF (type_hyb==313) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - ensstat/sqrt(N) with min. from 1-N_eff/N'
     ELSEIF (type_hyb==413) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- min of 1 - ensstat/sqrt(N) and 1-N_eff/N'
     ELSEIF (type_hyb==513) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - ensstat/sqrt(N) with min. from 1-(N_eff-1)/(N-1)'
     ELSEIF (type_hyb==613) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - ensstat/sqrt(N) with min. from 1-(N_eff-1)/(N-1) B'
     ELSEIF (type_hyb==712) THEN
        WRITE (*, '(a, 5x, a, f12.3, a, f12.1)') 'PDAF', &
             '--- dependence 1 - ensstat/sqrt(hnorm) with min. from N_eff/N>=', hybrid_a_x, ', hnorm=', hnorm
     ELSEIF (type_hyb==713) THEN
        WRITE (*, '(a, 5x, a, a, f12.1)') 'PDAF', &
             '--- dependence 1 - ensstat/sqrt(hnorm) with min. from 1-N_eff/N', ', hnorm', hnorm
     END IF
  END IF


! *********************************************
! *** Compute adaptive hybridization weight ***
! *********************************************

  IF (type_hyb>0) THEN
     ! Compute adaptive hybrid weight for state update

     IF (type_hyb==1) THEN
        alpha_x(1) = (1.0 - n_eff/REAL(dim_ens))*(hybrid_a_x)
     elseIF (type_hyb==2) THEN
        alpha_x(1) = alpha_adap
     ELSEIF (type_hyb==3) THEN
        alpha_x(1) = ((1.0 - n_eff/REAL(dim_ens))**2)*(hybrid_a_x)
     ELSEIF (type_hyb==4 .OR. type_hyb==5) THEN
        hfac = 1.0 - n_eff/REAL(dim_ens)
        IF (hfac>0.0) THEN
           alpha_x(1) = SQRT(1.0 - n_eff/REAL(dim_ens))*(hybrid_a_x)
        ELSE
           IF (type_hyb==4) THEN
              alpha_x(1) = 0.0
           ELSEIF (type_hyb==5) THEN
              alpha_x(1) = 1.0
           ENDIF
        ENDIF
     ELSEIF (type_hyb==6) THEN
        alpha_x(1) = n_eff/REAL(dim_ens)
        alpha_x(1) = 1.0 - SIN(alpha_x(1)*pi)*(1.0-hybrid_a_x)
        IF (alpha_x(1)<hybrid_a_x) alpha_x(1) = hybrid_a_x
        IF (alpha_x(1)>1.0) alpha_x(1) = 1.0
     ELSEIF (type_hyb==7) THEN
        hfac = 1.0 - n_eff/REAL(dim_ens)
        IF (hfac>0.0) THEN
           alpha_x(1) = SQRT(1.0 - n_eff/REAL(dim_ens))
        ELSE
           alpha_x(1) = 0.0
        ENDIF
        IF (alpha_x(1)<hybrid_a_x) alpha_x(1) = hybrid_a_x
     ELSEIF (type_hyb==8) THEN
        alpha_x(1) = (1.0 - 0.5*n_eff/REAL(dim_ens))*(hybrid_a_x)
     ELSEIF (type_hyb==9) THEN
        alpha_x(1) = (1.0 - 0.5*(n_eff/REAL(dim_ens))**2)*(hybrid_a_x)
     ELSEIF (type_hyb==10) THEN
        alpha_x(1) = 1.0 - skew_mabs(1)/SQRT(REAL(dim_ens))*(1.0-hybrid_a_x)
        IF (alpha_x(1) < (hybrid_a_x)) alpha_x(1) = hybrid_a_x
     ELSEIF (type_hyb==11) THEN
        alpha_x(1) = 1.0 - skew_mabs(1)/SQRT(REAL(dim_ens))
        IF (alpha_x(1) < (hybrid_a_x)) alpha_x(1) = hybrid_a_x
     ELSEIF (type_hyb==12) THEN
        alpha_skew = 1.0 - skew_mabs(1)/SQRT(REAL(dim_ens))
        IF (alpha_skew<0.0) alpha_skew = 0.0
        alpha_x(1) = MAX(alpha_skew, alpha_adap)
     ELSEIF (type_hyb==13) THEN
        alpha_skew = 1.0 - skew_mabs(1)/SQRT(REAL(dim_ens))
        IF (alpha_skew<0.0) alpha_skew = 0.0
        alpha_adap = 1.0 - n_eff/REAL(dim_ens)
        alpha_x(1) = MAX(alpha_skew, alpha_adap)
     ELSEIF (type_hyb==14) THEN
        alpha_x(1) = (1.0 - (n_eff-1.0)/REAL(dim_ens-1))*(hybrid_a_x)
        IF (alpha_x(1)<0.0) alpha_x(1) = 0.0
     ELSEIF (type_hyb==15) THEN
        alpha_x(1) = (1.0 - (n_eff-1.0)/REAL(dim_ens-1))*(hybrid_a_x)
        IF (alpha_x(1)<0.0) alpha_x(1) = 0.0
        IF (alpha_x(1)>0.95) alpha_x(1)=0.95
     ELSEIF (type_hyb==110) THEN
        kurt_limit = 0.5*(REAL(dim_ens-3))/(REAL(dim_ens-2))*skew_mabs(1)*skew_mabs(1) - REAL(dim_ens)/0.5-3.0
        alpha_x(1) = 1.0 - kurt_mabs(1)/kurt_limit*(1.0-hybrid_a_x)
        IF (alpha_x(1) < (hybrid_a_x)) alpha_x(1) = hybrid_a_x
     ELSEIF (type_hyb==111) THEN
        kurt_limit = 0.5*(REAL(dim_ens-3))/(REAL(dim_ens-2))*skew_mabs(1)*skew_mabs(1) - REAL(dim_ens)/0.5-3.0
        alpha_x(1) = 1.0 - kurt_mabs(1)/kurt_limit
        IF (alpha_x(1) < (hybrid_a_x)) alpha_x(1) = hybrid_a_x
     ELSEIF (type_hyb==112) THEN
        kurt_limit = 0.5*(REAL(dim_ens-3))/(REAL(dim_ens-2))*skew_mabs(1)*skew_mabs(1) - REAL(dim_ens)/0.5-3.0
        alpha_kurt = 1.0 - kurt_mabs(1)/kurt_limit
        IF (alpha_kurt<0.0) alpha_kurt = 0.0
        alpha_x(1) = MAX(alpha_kurt, alpha_adap)
     ELSEIF (type_hyb==113) THEN
        kurt_limit = 0.5*(REAL(dim_ens-3))/(REAL(dim_ens-2))*skew_mabs(1)*skew_mabs(1) - REAL(dim_ens)/0.5-3.0
        alpha_kurt = 1.0 - kurt_mabs(1)/kurt_limit
        IF (alpha_kurt<0.0) alpha_kurt = 0.0
        alpha_adap = 1.0 - n_eff/REAL(dim_ens)
        alpha_x(1) = MAX(alpha_kurt, alpha_adap)
     ELSEIF (type_hyb==212) THEN
        kurt_limit = REAL(dim_ens)
        alpha_kurt = 1.0 - kurt_mabs(1)/kurt_limit
        IF (alpha_kurt<0.0) alpha_kurt = 0.0
        alpha_x(1) = MAX(alpha_kurt, alpha_adap)
     ELSEIF (type_hyb==213) THEN
        kurt_limit = REAL(dim_ens)
        alpha_kurt = 1.0 - kurt_mabs(1)/kurt_limit
        IF (alpha_kurt<0.0) alpha_kurt = 0.0
        alpha_adap = 1.0 - n_eff/REAL(dim_ens)
        alpha_x(1) = MAX(alpha_kurt, alpha_adap)
     ELSEIF (type_hyb==312) THEN
        alpha_skew = 1.0 - skew_mabs(1)/SQRT(REAL(dim_ens))
        IF (alpha_skew<0.0) alpha_skew = 0.0
        alpha_kurt = 1.0 - kurt_mabs(1)/REAL(dim_ens)
        IF (alpha_kurt<0.0) alpha_kurt = 0.0
        alpha_stat = MIN(alpha_skew, alpha_kurt)
        alpha_x(1) = MAX(alpha_stat, alpha_adap)
     ELSEIF (type_hyb==313) THEN
        alpha_skew = 1.0 - skew_mabs(1)/SQRT(REAL(dim_ens))
        IF (alpha_skew<0.0) alpha_skew = 0.0
        alpha_kurt = 1.0 - kurt_mabs(1)/REAL(dim_ens)
        IF (alpha_kurt<0.0) alpha_kurt = 0.0
        alpha_stat = MIN(alpha_skew, alpha_kurt)
        alpha_adap = 1.0 - n_eff/REAL(dim_ens)
        alpha_x(1) = MAX(alpha_stat, alpha_adap)
     ELSEIF (type_hyb==413) THEN
        alpha_skew = 1.0 - skew_mabs(1)/SQRT(REAL(dim_ens))
        IF (alpha_skew<0.0) alpha_skew = 0.0
        alpha_kurt = 1.0 - kurt_mabs(1)/REAL(dim_ens)
        IF (alpha_kurt<0.0) alpha_kurt = 0.0
        alpha_stat = MIN(alpha_skew, alpha_kurt)
        alpha_adap = 1.0 - n_eff/REAL(dim_ens)
        IF (alpha_adap<0.0) alpha_adap=0.0
        alpha_x(1) = MIN(alpha_stat, alpha_adap)
     ELSEIF (type_hyb==513) THEN
        alpha_skew = 1.0 - skew_mabs(1)/SQRT(REAL(dim_ens))
        IF (alpha_skew<0.0) alpha_skew = 0.0
        alpha_kurt = 1.0 - kurt_mabs(1)/REAL(dim_ens)
        IF (alpha_kurt<0.0) alpha_kurt = 0.0
        alpha_stat = MIN(alpha_skew, alpha_kurt)
        alpha_adap = 1.0 - (n_eff-1.0)/REAL(dim_ens-1)
        alpha_x(1) = MAX(alpha_stat, alpha_adap)
     ELSEIF (type_hyb==613) THEN
        alpha_skew = 1.0 - skew_mabs(1)/SQRT(REAL(dim_ens))
        IF (alpha_skew<0.0) alpha_skew = 0.0
        alpha_kurt = 1.0 - kurt_mabs(1)/REAL(dim_ens)
        IF (alpha_kurt<0.0) alpha_kurt = 0.0
        alpha_stat = MIN(alpha_skew, alpha_kurt)
        alpha_adap = 1.0 - (n_eff-1.0)/REAL(dim_ens-1)
        alpha_x(1) = MAX(alpha_stat, alpha_adap)
        IF (alpha_x(1)>0.95) alpha_x(1)=0.95
     ELSEIF (type_hyb==712) THEN
        alpha_skew = 1.0 - skew_mabs(1)/SQRT(hnorm)
        IF (alpha_skew<0.0) alpha_skew = 0.0
        alpha_kurt = 1.0 - kurt_mabs(1)/hnorm
        IF (alpha_kurt<0.0) alpha_kurt = 0.0
        alpha_stat = MIN(alpha_skew, alpha_kurt)
        alpha_x(1) = MAX(alpha_stat, alpha_adap)
     ELSEIF (type_hyb==713) THEN
        alpha_skew = 1.0 - skew_mabs(1)/SQRT(hnorm)
        IF (alpha_skew<0.0) alpha_skew = 0.0
        alpha_kurt = 1.0 - kurt_mabs(1)/hnorm
        IF (alpha_kurt<0.0) alpha_kurt = 0.0
        alpha_stat = MIN(alpha_skew, alpha_kurt)
        alpha_adap = 1.0 - n_eff/REAL(dim_ens)
        alpha_x(1) = MAX(alpha_stat, alpha_adap)
     END IF

  ELSE

     ! fixed hybrid weight
     alpha_x(1) = hybrid_a_x

  END IF


! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1

  lastdomain = domain_p

END SUBROUTINE PDAF_lknetf_set_gamma
