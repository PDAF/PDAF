! Copyright (c) 2014-2018 Paul Kirchgessner
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
! !ROUTINE: PDAF_ewpf_ew_step --- equal-weights step of EWPF
!
! !INTERFACE:
SUBROUTINE PDAF_ewpf_ew_step(step, dim_p, dim_obs, dim_ens,&
     ens, weights, U_obs_op, U_init_dim_obs,&
     U_init_obs, U_solve_invHQHTpR, &
     U_adjoint_obs_op, U_prodqA, U_randVec, U_prodRinvA)

! !DESCRIPTION:
! Compute equivalent weight step of EWPF.
!
! Variant for ETKF with domain decomposition.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2014-05 - Paul Kirchgessner - Initial code
! Later revisions - see svn log
!
! !USES:
#include "typedefs.h"

!  use random  !module with functions to generate random numbers. (licence free)
  USE qsort_c_module 
  USE PDAF_mod_ewpf, &
       ONLY: keep

  IMPLICIT NONE

  INTEGER, INTENT(in) :: step
  INTEGER, INTENT(in) :: dim_p
  INTEGER, INTENT(in) :: dim_obs
  INTEGER, INTENT(in) :: dim_ens
  REAL, INTENT(inout), DIMENSION(dim_p,dim_ens) :: ens
  REAL, INTENT(inout), DIMENSION(dim_ens) :: weights

  EXTERNAL ::  U_init_dim_obs,    &  ! Get observation dimension
       U_init_obs,        &  ! Initilize observations
       U_obs_op,          &  ! Observation operator
       U_solve_invHQHTpR, &  ! Compute inverse of (HQH'+R)
       U_adjoint_obs_op,  &  ! Compute adjoint of observation operator
       U_prodqA,          &  ! Compute product of Q and some matrix A
       U_randVec,         &  ! Get a random vector with covariance Q
       U_prodRinvA           ! Compute the product of R inverse with some matrix A

  ! Local Variables 
  INTEGER :: particle 
  INTEGER, ALLOCATABLE :: index_out(:)
  REAL, ALLOCATABLE :: Hens(:,:)
  REAL, ALLOCATABLE :: Kymx(:,:)
  REAL :: w, alpha 
  REAL, ALLOCATABLE :: c(:), csorted(:)
  REAL :: cmax
  REAL, ALLOCATABLE :: xT(:)
  REAL, ALLOCATABLE :: a(:),b(:)
  REAL, ALLOCATABLE :: HtKymx(:), QHtKymx(:)
  REAL, ALLOCATABLE :: HQHtKymx(:),Rm1HQHTKymx(:)
  REAL, ALLOCATABLE :: betan(:)
  REAL  :: dRm1dt
  REAL, EXTERNAL :: ddot
  INTEGER :: i
  INTEGER :: cnt
  REAL, ALLOCATABLE :: ens_tmp(:,:)
  REAL, ALLOCATABLE :: observations(:)

  WRITE (*, '(a, 1x, 64a)') 'PDAF',('-', i = 1, 64)
  WRITE (*, '(a, 5x, a)') 'PDAF','+++++ Equal weights step +++++'

! Normalize weights: NOT NEEDED?!
! weights = exp(-weights + maxval(weights))
! weights = exp(-weights)
! weights = weights/sum(weights)
! weights = -log(weights)

  ALLOCATE(HtKymx(dim_p))
  ALLOCATE(QHtKymx(dim_p))
  ALLOCATE(betan(dim_p))
  ALLOCATE(index_out(dim_ens))
  ALLOCATE(a(dim_ens))
  ALLOCATE(b(dim_ens))
  ALLOCATE(c(dim_ens))
  ALLOCATE(csorted(dim_ens))

  a = 0 
  b = 0

  ! Initilize observations
  CALL U_init_dim_obs(step, dim_obs)
  ALLOCATE(observations(dim_obs)) ! Allocate observations
  ALLOCATE(Hens(dim_obs, dim_ens))
  ALLOCATE(Kymx(dim_obs, dim_ens))
  ALLOCATE(xT(dim_obs))
  ALLOCATE(HQHtKymx(dim_obs))
  ALLOCATE(Rm1HQHTKymx(dim_obs))

  CALL U_init_obs(step, dim_obs, observations)

  CALC_max_weight: DO particle = 1,dim_ens

     CALL U_obs_op(step,dim_p, dim_obs, ens(:,particle), Hens(:,particle))
     Hens(:,particle) = observations - Hens(:,particle)

     ! calc w1 = (HQH^T+R)^-1*(y-Hx)
     CALL U_solve_invHQHTpR(step, dim_obs, dim_p, dim_ens, Hens(:,particle), &
          Kymx(:,particle), ens(:,particle)) 
        ! the current model state needs to be in the interface for possible linearization
     c(particle)=0.5*DOT_PRODUCT(Kymx(:,particle), Hens(:,particle))
     c(particle) =weights(particle)+c(particle)    
  ENDDO CALC_max_weight

  ! sort weights
  csorted = c
  CALL qsortC(csorted)
  cmax = csorted(NINT((keep)*dim_ens)) !here cmin, since its the -log(w)

  particle_loop: DO particle = 1, dim_ens
     compute_correction: IF (c(particle) <= cmax) THEN
        ! Since R is symetric R^(-1)x instead of (x^TR^-1)^T is calculated
        CALL U_prodRinvA(step, dim_obs, 1, observations, Hens(:,particle), xT)
        dRm1dt = DOT_PRODUCT(xT,Hens(:,particle)) 

        b(particle) = 0.5*dRm1dt-cmax+weights(particle)
        ! Finished calculating b - a below
        CALL U_adjoint_obs_op(step, dim_obs, dim_p, Kymx(:,particle), HtKymx)
        CALL U_prodQA(dim_p, HtKymx, ens(:,particle), QHtKymx)
        CALL U_obs_op(step, dim_p, dim_obs, QHtKymx, HQHtKymx)
        CALL U_prodRinvA(step, dim_obs, 1, observations, HQHtKymx, Rm1HQHtKymx)

        a(particle) =  0.5*DOT_PRODUCT(Hens(:,particle), Rm1HQHtKymx)

        alpha = 1.0+SQRT(1.0-b(particle)/a(particle)+ 10.0D-6)

        ! Update determistic part of the state
        ens(:,particle) = ens(:,particle)+alpha*QHtKymx  

        weights(particle) = weights(particle) + &
             (alpha**2.0 - 2.0*alpha)*a(particle) +0.5*dRm1dt

        ! skip the Gaussian distribution, since its numerically zero! 
     ELSE
        ! these particles are resampled.
        weights(particle) = HUGE(REAL(1.0))
     ENDIF compute_correction !if(uniform)

     ! Add random forcing to the particles
     betan = 0.0
     CALL U_randvec(dim_p, betan, 1) 
     ens(:,particle) = ens(:,particle) + betan

  END DO particle_loop

  ! Resample particles:
  CALL PDAF_ewpf_resample(dim_ens, weights, index_out)

  ALLOCATE(ens_tmp(dim_p, dim_ens))

  ens_tmp = ens ! > should be possible to omit storing the whole ensemble again. 

  DO i = 1,dim_ens
     IF( i /= index_out(i)) THEN
        ens(:,i) = ens_tmp(:, index_out(i))
     ENDIF
  ENDDO


! ********************
! *** finishing up ***
! ********************

  DEALLOCATE(observations)
  DEALLOCATE(ens_tmp)

  DEALLOCATE(Hens, Kymx, xT, HQHtKymx, Rm1HQHTKymx)
  DEALLOCATE(HtKymx, QHtKymx, betan)
  DEALLOCATE(index_out, a, b, c, csorted)

END SUBROUTINE PDAF_ewpf_ew_step















