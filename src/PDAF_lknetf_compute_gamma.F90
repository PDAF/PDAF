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
! !ROUTINE: PDAF_lknetf_compute_gamma --- compute hybrid weight gamma for hybrid LKNETF
!
! !INTERFACE:
SUBROUTINE PDAF_lknetf_compute_gamma(domain_p, step, dim_obs_l, dim_ens, &
     HX_l, HXbar_l, obs_l, type_hyb, hyb_g, hyb_k, &
     gamma, n_eff_out, skew_mabs, kurt_mabs, &
     U_likelihood_l, screen, flag)

! !DESCRIPTION:
! This routine computes the hybrid weight gamma for the
! two-step variants HNK and HKN. 
! For this, it first computes the particle weights and
! then it calls PDAF_lknetf_set_gamma to get gamma
! according to the chosen hybridication weight type. 
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
  USE PDAF_mod_filter, &
       ONLY: obs_member

  IMPLICIT NONE

! !ARGUMENTS:
! ! Variable naming scheme:
! !   suffix _p: Denotes a full variable on the PE-local domain
! !   suffix _l: Denotes a local variable on the current analysis domain
  INTEGER, INTENT(in) :: domain_p    ! Current local analysis domain
  INTEGER, INTENT(in) :: step        ! Current time step
  INTEGER, INTENT(in) :: dim_obs_l   ! Size of obs. vector on local ana. domain
  INTEGER, INTENT(in) :: dim_ens     ! Size of ensemble 
  REAL, INTENT(in) :: HX_l(dim_obs_l, dim_ens)  ! local observed state ens.
  REAL, INTENT(in) :: HXbar_l(dim_obs_l)        ! local mean observed ensemble
  REAL, INTENT(in) :: obs_l(dim_obs_l) ! Local observation vector
  INTEGER, INTENT(in) :: type_hyb      ! Type of hybrid weight
  REAL, INTENT(in) :: hyb_g       ! Prescribed hybrid weight for state transformation
  REAL, INTENT(in) :: hyb_k            ! Hybrid weight for covariance transformation
  REAL, INTENT(inout) :: gamma(1)    ! Hybrid weight for state transformation
  REAL, INTENT(inout) :: n_eff_out(1)  ! Effective ensemble size
  REAL, INTENT(inout) :: skew_mabs(1)  ! Mean absolute skewness
  REAL, INTENT(inout) :: kurt_mabs(1)  ! Mean absolute kurtosis
  INTEGER, INTENT(in) :: screen        ! Verbosity flag
  INTEGER, INTENT(inout) :: flag       ! Status flag


! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_likelihood_l           ! Compute observation likelihood for an ensemble member

! !CALLING SEQUENCE:
! Called by: PDAF_lknetf_step_update
! Calls: PDAF_timeit
! Calls: PDAF_memcount
!EOP
       
! *** local variables ***
  INTEGER :: i, member               ! Counters
  INTEGER, SAVE :: allocflag=0       ! Flag whether first time allocation is done
  REAL :: weight                     ! Ensemble weight (likelihood)
  REAL, ALLOCATABLE :: resid_i(:)    ! Residual
  REAL, ALLOCATABLE :: weights(:)    ! weight vector
  REAL :: total_weight               ! Sum of weight

!$OMP THREADPRIVATE(allocflag)


! **********************************************
! *** Compute particle weights as likelihood ***
! **********************************************

  ! Allocate weight vector
  ALLOCATE(weights(dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)

  ! Allocate temporal array for obs-ens_i
  ALLOCATE(resid_i(dim_obs_l))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l)

  CALL PDAF_timeit(22, 'new')
  ! Get residual as difference of observation and observed state for 
  ! each ensemble member only on domains where observations are availible

  CALC_w: DO member = 1,dim_ens

     ! Store member index to make it accessible with PDAF_get_obsmemberid
     obs_member = member
     
     ! Calculate local residual  
     resid_i = obs_l - HX_l(:,member)

     ! Compute likelihood
     CALL PDAF_timeit(49, 'new')
     CALL U_likelihood_l(domain_p, step, dim_obs_l, obs_l, resid_i, weight)
     CALL PDAF_timeit(49, 'old')
     weights(member) = weight

  END DO CALC_w

  ! Normalize weights
  total_weight = 0.0
  DO i = 1, dim_ens
     total_weight = total_weight + weights(i)
  END DO

  IF (total_weight /= 0.0) THEN
     weights = weights / total_weight
  ELSE
     ! ERROR: weights are zero
     flag = 1
     WRITE(*,'(/5x,a/)') 'PDAF-ERROR (1): Zero weights in LKNETF analysis step'
  END IF

  CALL PDAF_timeit(22, 'old')


! *******************************
! *** Set hybrid weight gamma ***
! *******************************

 CALL PDAF_lknetf_set_gamma(domain_p, dim_obs_l, dim_ens, &
    HX_l, HXbar_l, weights, type_hyb, hyb_g, hyb_k, &
    gamma, n_eff_out, skew_mabs, kurt_mabs, &
    screen, flag)


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(weights, resid_i)

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_lknetf_compute_gamma
