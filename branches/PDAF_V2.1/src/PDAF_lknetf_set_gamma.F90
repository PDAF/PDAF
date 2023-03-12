! Copyright (c) 2014-2023 Lars Nerger
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
     HX_l, HXbar_l, weights, type_hyb, hyb_g, hyb_k, &
     gamma, n_eff_out, maSkew, maKurt, &
     screen, flag)

! !DESCRIPTION:
! Compute hybrid weight for LKNETF
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
  USE PDAF_mod_filter, &
       ONLY: debug
#if defined (_OPENMP)
  USE omp_lib, &
       ONLY: omp_get_num_threads, omp_get_thread_num
#endif

  IMPLICIT NONE

! !ARGUMENTS:
! ! Variable naming scheme:
! !   suffix _p: Denotes a full variable on the PE-local domain
! !   suffix _l: Denotes a local variable on the current analysis domain
  INTEGER, INTENT(in) :: domain_p      ! Current local analysis domain
  INTEGER, INTENT(in) :: dim_obs_l     ! Size of obs. vector on local ana. domain
  INTEGER, INTENT(in) :: dim_ens       ! Size of ensemble 
  REAL, INTENT(in) :: HX_l(dim_obs_l, dim_ens)  ! local observed state ens.
  REAL, INTENT(in) :: HXbar_l(dim_obs_l)        ! local mean observed ensemble
  REAL, INTENT(in) :: weights(dim_ens) ! Weight vector
  INTEGER, INTENT(in) :: type_hyb      ! Type of hybrid weight
  REAL, INTENT(in) :: hyb_g            ! Prescribed hybrid weight for state transformation
  REAL, INTENT(in) :: hyb_k            ! Scale factor kappa (for type_hyb 3 and 4)
  REAL, INTENT(inout) :: gamma(1)      ! Hybrid weight for state transformation
  REAL, INTENT(inout) :: n_eff_out(1)  ! Effective ensemble size
  REAL, INTENT(inout) :: maSkew(1)     ! Mean absolute skewness
  REAL, INTENT(inout) :: maKurt(1)     ! Mean absolute kurtosis
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
  REAL, PARAMETER :: pi=3.14159265358979    ! Pi
  REAL :: skew(1), kurt(1)           ! Skewness and kurtosis of observed ensemble
  REAL :: kurt_limit                 ! Asymptotic value of kurtosis
  REAL :: gamma_kurt                 ! Gamma from kurtosis
  REAL :: gamma_skew                 ! Gamma from Skewness
  REAL :: gamma_Neff                 ! Gamma from effective sample size
  REAL :: gamma_stat                 ! Gamma from combining kurtosis and skewness

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

  ! Dummy initialization to prevent compiler warning
  gamma_Neff = hyb_g


  ! *********************
  ! *** Screen output ***
  ! *********************

  IF (mype == 0 .AND. screen > 0 .AND. screenout) THEN
     IF (type_hyb==0) THEN
        WRITE (*, '(a, 5x, a, f8.3)') 'PDAF', 'Set hybrid weight in LETKF to ', hyb_g
     ELSE
        WRITE (*, '(a, 5x, a)') 'PDAF', 'Compute adaptive hybrid weight'
     END IF

     ! First four methods are discussed in Nerger (2022)
     IF (type_hyb==1) THEN
        WRITE (*, '(a, 5x, a, f8.3)') 'PDAF', '--- gamma_lin: (1 - N_eff/N)*scale, scale=', hyb_g
     ELSEIF (type_hyb==2) THEN
        WRITE (*, '(a, 5x, a, f8.3)') 'PDAF', '--- gamma_alpha: hybrid weight from N_eff/N>=', hyb_g
     ELSEIF (type_hyb==3) THEN
        WRITE (*, '(a, 5x, a, f8.3, a, f8.2)') 'PDAF', &
             '--- gamma_ska: 1 - min(s,k)/sqrt(kappa) with N_eff/N>=', hyb_g, ', kappa=', hyb_k
     ELSEIF (type_hyb==4) THEN
        WRITE (*, '(a, 5x, a, a, f8.3)') 'PDAF', &
             '--- gamma_sklin: 1 - min(s,k)/sqrt(kappa) >= 1-N_eff/N', ', kappa', hyb_k

     ! Additional methods
     ELSEIF (type_hyb==23) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- quadratic dependence on 1 - N_eff/N'
     ELSEIF (type_hyb==24) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- square-root dependence on 1 - N_eff/N; 0 for N_eff=N'
     ELSEIF (type_hyb==25) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- square-root dependence on 1 - N_eff/N; 1 for N_eff=N'
     ELSEIF (type_hyb==26) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- sine dependence on N_eff/N with minimum constraint'
     ELSEIF (type_hyb==27) THEN
        WRITE (*, '(a, 5x, a,f8.3)') 'PDAF', '--- square-root dependence on N_eff/N, minimum limit', hyb_g
     ELSEIF (type_hyb==28) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- linear dependence 1 - 0.5 N_eff/N'
     ELSEIF (type_hyb==29) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- quadratic dependence 1 - 0.5 (N_eff/N)^2'
     ELSEIF (type_hyb==10) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - (1-limit)* skewness/sqrt(N)'
     ELSEIF (type_hyb==11) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - skewness/sqrt(N) with minimum limit'
     ELSEIF (type_hyb==12) THEN
        WRITE (*, '(a, 5x, a, f12.3)') 'PDAF', '--- dependence 1 - skewness/sqrt(N) with min. from N_eff/N>=', hyb_g
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
        WRITE (*, '(a, 5x, a, f12.3)') 'PDAF', '--- dependence 1 - kurtosis/sqrt(N) with min. from N_eff/N>=', hyb_g
     ELSEIF (type_hyb==113) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - kurtosis/sqrt(N) with min. from 1-N_eff/N'
     ELSEIF (type_hyb==212) THEN
        WRITE (*, '(a, 5x, a, f12.3)') 'PDAF', '--- dependence 1 - kurtnorm/sqrt(N) with min. from N_eff/N>=', hyb_g
     ELSEIF (type_hyb==213) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - kurtnorm/sqrt(N) with min. from 1-N_eff/N'
     ELSEIF (type_hyb==312) THEN
        WRITE (*, '(a, 5x, a, f12.3)') 'PDAF', '--- dependence 1 - ensstat/sqrt(N) with min. from N_eff/N>=', hyb_g
     ELSEIF (type_hyb==313) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - ensstat/sqrt(N) with min. from 1-N_eff/N'
     ELSEIF (type_hyb==413) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- min of 1 - ensstat/sqrt(N) and 1-N_eff/N'
     ELSEIF (type_hyb==513) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - ensstat/sqrt(N) with min. from 1-(N_eff-1)/(N-1)'
     ELSEIF (type_hyb==613) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - ensstat/sqrt(N) with min. from 1-(N_eff-1)/(N-1) B'
     END IF
  END IF


  ! **********************************************
  ! *** Compute particle weights as likelihood ***
  ! **********************************************

  CALL PDAF_timeit(54, 'new')
  ! Get residual as difference of observation and observed state for 
  ! each ensemble member only on domains where observations are availible
  
  ! Compute adaptive hybrid weights
  IF (type_hyb==2 .OR. type_hyb==3 .OR. type_hyb==12 .OR. type_hyb==112 &
       .OR. type_hyb==212 .OR. type_hyb==312) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, &
          'PDAF_lknetf_compute_gamma -- Determine gamma for N_eff/N>=limit iteratively'
     CALL PDAF_lknetf_alpha_neff(dim_ens, weights, hyb_g, gamma(1))
     gamma_Neff = gamma(1)

     IF (debug>0 .AND. type_hyb/=2) &
          WRITE (*,*) '++ PDAF-debug PDAF_lknetf_compute_gamma:', debug, '  gamma for N_eff/N>=limit', gamma_Neff
  END IF

  CALL PDAF_timeit(54, 'old')


  ! Compute effective ensemble size
  CALL PDAF_diag_effsample(dim_ens, weights, n_eff)
  n_eff_out(1) = n_eff

  ! Compute ensemble statistics for observed ensemble
  DO i = 1, dim_obs_l
     CALL PDAF_diag_ensstats(dim_obs_l, dim_ens, i, &
          HXbar_l, HX_l, skew, kurt, flag)
     maSkew = maSkew + ABS(skew)
     maKurt = maKurt + ABS(kurt)
  END DO
  maSkew = maSkew / dim_obs_l
  maKurt = maKurt / dim_obs_l

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lknetf_compute_gamma:', debug, '  MAS, MAK', maSkew, maKurt


! *********************************************
! *** Compute adaptive hybridization weight ***
! *********************************************

  IF (type_hyb>0) THEN
     ! Compute adaptive hybrid weight for state update

     ! First four methods are discussed in Nerger (2022)
     IF (type_hyb==1) THEN
        gamma(1) = (1.0 - n_eff/REAL(dim_ens))*(hyb_g)
     ELSEIF (type_hyb==2) THEN
        gamma(1) = gamma_Neff
     ELSEIF (type_hyb==3) THEN
        gamma_skew = 1.0 - maSkew(1)/SQRT(hyb_k)
        IF (gamma_skew<0.0) gamma_skew = 0.0
        gamma_kurt = 1.0 - maKurt(1)/hyb_k
        IF (gamma_kurt<0.0) gamma_kurt = 0.0
        gamma_stat = MIN(gamma_skew, gamma_kurt)
        gamma(1) = MAX(gamma_stat, gamma_Neff)
     ELSEIF (type_hyb==4) THEN
        gamma_skew = 1.0 - maSkew(1)/SQRT(hyb_k)
        IF (gamma_skew<0.0) gamma_skew = 0.0
        gamma_kurt = 1.0 - maKurt(1)/hyb_k
        IF (gamma_kurt<0.0) gamma_kurt = 0.0
        gamma_stat = MIN(gamma_skew, gamma_kurt)
        gamma_Neff = 1.0 - n_eff/REAL(dim_ens)
        gamma(1) = MAX(gamma_stat, gamma_Neff)

     ! Additional methods - experimental
     ELSEIF (type_hyb==23) THEN
        gamma(1) = ((1.0 - n_eff/REAL(dim_ens))**2)*(hyb_g)
     ELSEIF (type_hyb==24 .OR. type_hyb==25) THEN
        IF (1.0 - n_eff/REAL(dim_ens) > 0.0) THEN
           gamma(1) = SQRT(1.0 - n_eff/REAL(dim_ens))*(hyb_g)
        ELSE
           IF (type_hyb==24) THEN
              gamma(1) = 0.0
           ELSEIF (type_hyb==25) THEN
              gamma(1) = 1.0
           ENDIF
        ENDIF
     ELSEIF (type_hyb==26) THEN
        gamma(1) = n_eff/REAL(dim_ens)
        gamma(1) = 1.0 - SIN(gamma(1)*pi)*(1.0-hyb_g)
        IF (gamma(1)<hyb_g) gamma(1) = hyb_g
        IF (gamma(1)>1.0) gamma(1) = 1.0
     ELSEIF (type_hyb==27) THEN
        IF (1.0 - n_eff/REAL(dim_ens) > 0.0) THEN
           gamma(1) = SQRT(1.0 - n_eff/REAL(dim_ens))
        ELSE
           gamma(1) = 0.0
        ENDIF
        IF (gamma(1)<hyb_g) gamma(1) = hyb_g
     ELSEIF (type_hyb==28) THEN
        gamma(1) = (1.0 - 0.5*n_eff/REAL(dim_ens))*(hyb_g)
     ELSEIF (type_hyb==29) THEN
        gamma(1) = (1.0 - 0.5*(n_eff/REAL(dim_ens))**2)*(hyb_g)
     ELSEIF (type_hyb==10) THEN
        gamma(1) = 1.0 - maSkew(1)/SQRT(REAL(dim_ens))*(1.0-hyb_g)
        IF (gamma(1) < (hyb_g)) gamma(1) = hyb_g
     ELSEIF (type_hyb==11) THEN
        gamma(1) = 1.0 - maSkew(1)/SQRT(REAL(dim_ens))
        IF (gamma(1) < (hyb_g)) gamma(1) = hyb_g
     ELSEIF (type_hyb==12) THEN
        gamma_skew = 1.0 - maSkew(1)/SQRT(REAL(dim_ens))
        IF (gamma_skew<0.0) gamma_skew = 0.0
        gamma(1) = MAX(gamma_skew, gamma_Neff)
     ELSEIF (type_hyb==13) THEN
        gamma_skew = 1.0 - maSkew(1)/SQRT(REAL(dim_ens))
        IF (gamma_skew<0.0) gamma_skew = 0.0
        gamma_Neff = 1.0 - n_eff/REAL(dim_ens)
        gamma(1) = MAX(gamma_skew, gamma_Neff)
     ELSEIF (type_hyb==14) THEN
        gamma(1) = (1.0 - (n_eff-1.0)/REAL(dim_ens-1))*(hyb_g)
        IF (gamma(1)<0.0) gamma(1) = 0.0
     ELSEIF (type_hyb==15) THEN
        gamma(1) = (1.0 - (n_eff-1.0)/REAL(dim_ens-1))*(hyb_g)
        IF (gamma(1)<0.0) gamma(1) = 0.0
        IF (gamma(1)>0.95) gamma(1)=0.95
     ELSEIF (type_hyb==110) THEN
        kurt_limit = 0.5*(REAL(dim_ens-3))/(REAL(dim_ens-2))*maSkew(1)*maSkew(1) - REAL(dim_ens)/0.5-3.0
        gamma(1) = 1.0 - maKurt(1)/kurt_limit*(1.0-hyb_g)
        IF (gamma(1) < (hyb_g)) gamma(1) = hyb_g
     ELSEIF (type_hyb==111) THEN
        kurt_limit = 0.5*(REAL(dim_ens-3))/(REAL(dim_ens-2))*maSkew(1)*maSkew(1) - REAL(dim_ens)/0.5-3.0
        gamma(1) = 1.0 - maKurt(1)/kurt_limit
        IF (gamma(1) < (hyb_g)) gamma(1) = hyb_g
     ELSEIF (type_hyb==112) THEN
        kurt_limit = 0.5*(REAL(dim_ens-3))/(REAL(dim_ens-2))*maSkew(1)*maSkew(1) - REAL(dim_ens)/0.5-3.0
        gamma_kurt = 1.0 - maKurt(1)/kurt_limit
        IF (gamma_kurt<0.0) gamma_kurt = 0.0
        gamma(1) = MAX(gamma_kurt, gamma_Neff)
     ELSEIF (type_hyb==113) THEN
        kurt_limit = 0.5*(REAL(dim_ens-3))/(REAL(dim_ens-2))*maSkew(1)*maSkew(1) - REAL(dim_ens)/0.5-3.0
        gamma_kurt = 1.0 - maKurt(1)/kurt_limit
        IF (gamma_kurt<0.0) gamma_kurt = 0.0
        gamma_Neff = 1.0 - n_eff/REAL(dim_ens)
        gamma(1) = MAX(gamma_kurt, gamma_Neff)
     ELSEIF (type_hyb==212) THEN
        kurt_limit = REAL(dim_ens)
        gamma_kurt = 1.0 - maKurt(1)/kurt_limit
        IF (gamma_kurt<0.0) gamma_kurt = 0.0
        gamma(1) = MAX(gamma_kurt, gamma_Neff)
     ELSEIF (type_hyb==213) THEN
        kurt_limit = REAL(dim_ens)
        gamma_kurt = 1.0 - maKurt(1)/kurt_limit
        IF (gamma_kurt<0.0) gamma_kurt = 0.0
        gamma_Neff = 1.0 - n_eff/REAL(dim_ens)
        gamma(1) = MAX(gamma_kurt, gamma_Neff)
     ELSEIF (type_hyb==312) THEN
        gamma_skew = 1.0 - maSkew(1)/SQRT(REAL(dim_ens))
        IF (gamma_skew<0.0) gamma_skew = 0.0
        gamma_kurt = 1.0 - maKurt(1)/REAL(dim_ens)
        IF (gamma_kurt<0.0) gamma_kurt = 0.0
        gamma_stat = MIN(gamma_skew, gamma_kurt)
        gamma(1) = MAX(gamma_stat, gamma_Neff)
     ELSEIF (type_hyb==313) THEN
        gamma_skew = 1.0 - maSkew(1)/SQRT(REAL(dim_ens))
        IF (gamma_skew<0.0) gamma_skew = 0.0
        gamma_kurt = 1.0 - maKurt(1)/REAL(dim_ens)
        IF (gamma_kurt<0.0) gamma_kurt = 0.0
        gamma_stat = MIN(gamma_skew, gamma_kurt)
        gamma_Neff = 1.0 - n_eff/REAL(dim_ens)
        gamma(1) = MAX(gamma_stat, gamma_Neff)
     ELSEIF (type_hyb==413) THEN
        gamma_skew = 1.0 - maSkew(1)/SQRT(REAL(dim_ens))
        IF (gamma_skew<0.0) gamma_skew = 0.0
        gamma_kurt = 1.0 - maKurt(1)/REAL(dim_ens)
        IF (gamma_kurt<0.0) gamma_kurt = 0.0
        gamma_stat = MIN(gamma_skew, gamma_kurt)
        gamma_Neff = 1.0 - n_eff/REAL(dim_ens)
        IF (gamma_Neff<0.0) gamma_Neff=0.0
        gamma(1) = MIN(gamma_stat, gamma_Neff)
     ELSEIF (type_hyb==513) THEN
        gamma_skew = 1.0 - maSkew(1)/SQRT(REAL(dim_ens))
        IF (gamma_skew<0.0) gamma_skew = 0.0
        gamma_kurt = 1.0 - maKurt(1)/REAL(dim_ens)
        IF (gamma_kurt<0.0) gamma_kurt = 0.0
        gamma_stat = MIN(gamma_skew, gamma_kurt)
        gamma_Neff = 1.0 - (n_eff-1.0)/REAL(dim_ens-1)
        gamma(1) = MAX(gamma_stat, gamma_Neff)
     ELSEIF (type_hyb==613) THEN
        gamma_skew = 1.0 - maSkew(1)/SQRT(REAL(dim_ens))
        IF (gamma_skew<0.0) gamma_skew = 0.0
        gamma_kurt = 1.0 - maKurt(1)/REAL(dim_ens)
        IF (gamma_kurt<0.0) gamma_kurt = 0.0
        gamma_stat = MIN(gamma_skew, gamma_kurt)
        gamma_Neff = 1.0 - (n_eff-1.0)/REAL(dim_ens-1)
        gamma(1) = MAX(gamma_stat, gamma_Neff)
        IF (gamma(1)>0.95) gamma(1)=0.95
     END IF

  ELSE

     ! fixed hybrid weight
     gamma(1) = hyb_g

  END IF


! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1

  lastdomain = domain_p

END SUBROUTINE PDAF_lknetf_set_gamma
