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
!$Id$
!BOP
!
! !ROUTINE: PDAF_set_forget_local - Set local adaptive forgetting factor
!
! !INTERFACE:
SUBROUTINE PDAF_set_forget_local(domain, step, dim_obs_l, dim_ens, &
     HX_l, HXbar_l, obs_l, U_init_obsvar_l, forget, aforget)

! !DESCRIPTION:
! Dynamically set the global forgetting factor individually for
! each local analysis domain in LSEIK. 
! This is a typical implementation that tries to ensure
! statistical consistency by enforcing the condition\\
! var\_resid = 1/forget var\_ens + var\_obs\\
! where var\_res is the variance of the innovation residual,
! var\_ens is the ensemble-estimated variance, and
! var\_obs is the observation error variance.\\
! This variant is not proven to improve the estimates!
!
! __Revision history:__
! 2006-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_mod_filtermpi, &
       ONLY: mype
#if defined (_OPENMP)
  USE omp_lib, &
       ONLY: omp_get_num_threads, omp_get_thread_num
#endif

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain           ! Current local analysis domain
  INTEGER, INTENT(in) :: step             ! Current time step
  INTEGER, INTENT(in) :: dim_obs_l        ! Dimension of local observation vector
  INTEGER, INTENT(in) :: dim_ens          ! Ensemble size
  REAL, INTENT(in) :: HX_l(dim_obs_l, dim_ens) ! Local observed ensemble
  REAL, INTENT(in) :: HXbar_l(dim_obs_l)  ! Local observed state estimate
  REAL, INTENT(in) :: obs_l(dim_obs_l)    ! Local observation vector
  REAL, INTENT(in) :: forget              ! Prescribed forgetting factor
  REAL, INTENT(out) :: aforget            ! Adaptive forgetting factor

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_obsvar_l             ! Initialize local mean obs. error variance

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_analysis
! Calls: U_init_obsvar_l
!EOP
  
! *** local variables ***
  INTEGER :: i, j                     ! Counters
  REAL :: var_ens, var_resid, var_obs ! Variances
  INTEGER, SAVE :: first = 1          ! Flag for very first call to routine
  INTEGER, SAVE :: lastdomain = -1    ! store domain index
  LOGICAL, SAVE :: screenout = .true. ! Whether to print information to stdout
  INTEGER, SAVE :: mythread, nthreads ! Thread variables for OpenMP
  REAL :: forget_neg, forget_max, forget_min ! limiting values of forgetting factor


! **********************
! *** INITIALIZATION ***
! **********************

  ! Define limiting values of forgetting factor
  ! These are set very arbitrarily for now
  forget_neg = forget
  forget_max = 100.0
  forget_min = 0.01

#if defined (_OPENMP)
  nthreads = omp_get_num_threads()
  mythread = omp_get_thread_num()
#else
  nthreads = 1
  mythread = 0
#endif


  ! Control screen output
  IF (lastdomain<domain .AND. lastdomain>-1) THEN
     screenout = .false.
  ELSE
     screenout = .true.

     ! In case of OpenMP, let only thread 0 write output to the screen
     IF (mythread>0) screenout = .false.

     ! Output, only in case of OpenMP parallelization
  END IF


! ****************************************************
! *** Initialize adaptive local forgetting factors ***
! ****************************************************

  IF (screenout) THEN
     ! At first call during each forecast phase
     IF (mype == 0) THEN
        WRITE (*, '(a, 5x, a)') &
             'PDAF', '--- Apply dynamically estimated local forgetting factors'
        WRITE (*, '(a, 9x, a, es10.2)') &
             'PDAF', 'Maximum limit for forgetting factor', forget_max
        WRITE (*, '(a, 9x, a, es10.2)') &
             'PDAF', 'Minimum limit for forgetting factor', forget_min
        WRITE (*, '(a, 9x, a, es10.2)') &
             'PDAF', 'Forgetting factor if var(obs)>var(resid)', forget_neg
     ENDIF
          
     ! Set flag
     first = 0
  ENDIF
  lastdomain = domain


  ! ************************************************************
  ! *** Compute optimal forgetting factor for current domain ***
  ! ************************************************************

  ! *** Compute mean ensemble variance ***

  ! local
  var_ens = 0.0
  DO i = 1, dim_obs_l
     DO j = 1, dim_ens
        var_ens = var_ens + (HXbar_l(i) - HX_l(i, j)) ** 2
     ENDDO
  ENDDO
  var_ens = var_ens / REAL(dim_ens - 1) / REAL(dim_obs_l)


  ! *** Compute mean of observation-minus-forecast residual ***
   
  var_resid = 0.0
  DO i = 1, dim_obs_l
     var_resid = var_resid + (obs_l(i) - HXbar_l(i)) ** 2
  ENDDO
  var_resid = var_resid / REAL(dim_obs_l)


  ! *** Compute mean observation variance ***

  ! Get mean observation error variance
  CALL PDAF_timeit(52, 'new')
  CALL U_init_obsvar_l(domain, step, dim_obs_l, obs_l, var_obs)
  CALL PDAF_timeit(52, 'old')

  ! *** Compute optimal forgetting factor ***
  aforget = var_ens / (var_resid - var_obs)

  ! Apply special condition if observation variance is larger than residual variance
  IF (aforget < 0.0) aforget = forget_neg

  ! Impose upper limit for forgetting factor
  ! - the value for this is quite arbitary
  IF (aforget > forget_max) aforget = forget_max

  ! Impose lower limit for forgetting factor
  ! - the value for this is quite arbitary
  IF (aforget < forget_min) aforget = forget_min

END SUBROUTINE PDAF_set_forget_local
