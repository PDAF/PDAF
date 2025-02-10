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
!> Set adaptive forgetting factor
!!
!! Dynamically set the global forgetting factor.
!! This is a typical implementation that tries to ensure
!! statistical consistency by enforcing the condition\\
!! var\_resid = 1/forget var\_ens + var\_obs\\
!! where var\_res is the variance of the innovation residual,
!! var\_ens is the ensemble-estimated variance, and
!! var\_obs is the observation error variance.\\
!! This routine is used in SEIK. It can also be used in LSEIK. 
!! In this case a forgetting factor for the PE-local domain is 
!! computed. An alternative for LSEIK is PDAF\_set\_forget\_local, 
!! which computes a forgetting factor for each local analysis 
!! domain. The implementation used in both routines is 
!! experimental and not proven to improve the estimates.
!!
!! __Revision history:__
!! * 2006-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE PDAF_set_forget(step, localfilter, dim_obs_p, dim_ens, mens_p, &
     mstate_p, obs_p, U_init_obsvar, forget_in, forget_out, &
     screen)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE mpi
  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_mod_filtermpi, &
       ONLY: mype, npes_filter, MPIerr, COMM_filter

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                      !< Current time step
  INTEGER, INTENT(in) :: dim_obs_p                 !< Dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens                   !< Ensemble size
  REAL, INTENT(in) :: mens_p(dim_obs_p, dim_ens)   !< Observed PE-local ensemble
  REAL, INTENT(in) :: mstate_p(dim_obs_p)          !< Observed PE-local mean state
  REAL, INTENT(in) :: obs_p(dim_obs_p)             !< Observation vector
  REAL, INTENT(in) :: forget_in                    !< Prescribed forgetting factor
  REAL, INTENT(out) :: forget_out                  !< Adaptively estimated forgetting factor
  INTEGER, INTENT(in) :: localfilter               !< Whether filter is domain-local
  INTEGER, INTENT(in) :: screen                    !< Verbosity flag

! *** External subroutine ***
!  (PDAF-internal name, real name is defined in the call to PDAF)
  EXTERNAL :: U_init_obsvar                        !< Initialize mean obs. error variance
  
! *** local variables ***
  INTEGER :: i, j                            ! Counters
  INTEGER :: dim_obs                         ! Global observation dimension for non-local filters
                                             ! PE-local dimension for local filters
  INTEGER :: dim_obs_g                       ! Global observation dimension
  REAL :: var_ens_p, var_ens                 ! Variance of ensemble
  REAL :: var_resid_p, var_resid             ! Variance of residual
  REAL :: var_obs                            ! Variance of observation errors
  REAL :: forget_neg, forget_max, forget_min ! Limiting values of forgetting factor


! **********************
! *** INITIALIZATION ***
! **********************

  ! Define limiting values of forgetting factor
  ! These are set very arbitrarily for now
  forget_neg = forget_in
  forget_max = 100.0
  forget_min = 0.01

  IF (mype == 0) THEN
     WRITE (*, '(a, 5x, a)') &
          'PDAF', '--- Apply global adaptive forgetting factor'
     WRITE (*, '(a, 9x, a, es10.2)') &
          'PDAF', 'Maximum limit for forgetting factor', forget_max
     WRITE (*, '(a, 9x, a, es10.2)') &
          'PDAF', 'Minimum limit for forgetting factor', forget_min
     WRITE (*, '(a, 9x, a, es10.2)') &
          'PDAF', 'Forgetting factor if var(obs) > var(resid)', forget_neg
  ENDIF


! ******************************************
! *** Compute adaptive forgetting factor ***
! ******************************************

  ! *** Compute mean ensemble variance ***

  CALL PDAF_timeit(51, 'new')

  IF (npes_filter>1) THEN
     CALL MPI_allreduce(dim_obs_p, dim_obs_g, 1, &
          MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
  ELSE
     dim_obs_g = dim_obs_p
  END IF

  IF (localfilter==0) THEN
     ! global fitlers 
     dim_obs = dim_obs_g
  ELSE
     ! domain-local filters
     dim_obs = dim_obs_p
  ENDIF

  IF (dim_obs_g > 0) THEN

     ! local
     var_ens_p = 0.0
     IF (dim_obs > 0) THEN
        DO i = 1, dim_obs_p
           DO j = 1, dim_ens
              var_ens_p = var_ens_p + (mstate_p(i) - mens_p(i, j)) ** 2
           ENDDO
        ENDDO
        var_ens_p = var_ens_p / REAL(dim_ens - 1) / REAL(dim_obs)
     END IF

     IF (localfilter==0) THEN
        ! global 
        CALL MPI_allreduce(var_ens_p, var_ens, 1, MPI_REALTYPE, MPI_SUM, &
             COMM_filter, MPIerr)
     ELSE
        ! For domain-local filters use only PE-local variance
        var_ens = var_ens_p
     ENDIF

     ! *** Compute mean of innovation ***
   
     ! Compute variance
     var_resid_p = 0.0
     IF (dim_obs > 0) THEN
        DO i = 1, dim_obs_p
           var_resid_p = var_resid_p + (obs_p(i) - mstate_p(i)) ** 2
        ENDDO
        var_resid_p = var_resid_p / REAL(dim_obs)
     END IF
     
     IF (localfilter==0) THEN
        ! global 
        CALL MPI_allreduce(var_resid_p, var_resid, 1, &
             MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)
     ELSE
        ! For domain-local filters use only PE-local variance
        var_resid = var_resid_p
     ENDIF

     CALL PDAF_timeit(51, 'old')

     ! *** Compute mean observation variance ***

     ! Get mean observation error variance
     CALL PDAF_timeit(49, 'new')
     IF (dim_obs_p>0) THEN
        CALL U_init_obsvar(step, dim_obs_p, obs_p, var_obs)
     ELSE
        var_obs=0.0
     END IF
     CALL PDAF_timeit(49, 'old')

     CALL PDAF_timeit(51, 'new')

     ! *** Compute optimal forgetting factor ***
     IF (var_resid>0.0 .AND. var_obs>0.0) THEN
        forget_out = var_ens / (var_resid - var_obs)
     ELSE
        forget_out = forget_in
     END IF

     ! Apply special condition if observation variance is larger than residual variance
     IF (forget_out < 0.0) forget_out = forget_neg

     ! Impose upper limit for forgetting factor
     IF (forget_out > forget_max) forget_out = forget_max

     ! Impose lower limit for forgetting factor
     IF (forget_out < forget_min) forget_out = forget_min

  ELSE
     ! No observations available in full model domain
     forget_out = forget_in
  END IF


! ********************
! *** FINISHING UP ***
! ********************

  IF (mype == 0 .AND. screen>0) THEN
     WRITE (*, '(a, 9x, a, es10.2)') &
          'PDAF', '--> Computed forgetting factor', forget_out
  ENDIF

  CALL PDAF_timeit(51, 'old')
   
END SUBROUTINE PDAF_set_forget
