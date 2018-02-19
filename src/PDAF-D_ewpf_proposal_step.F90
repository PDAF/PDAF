! Copyright (c) 2004-2018 Lars Nerger and Paul Kirchgessner
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
! !ROUTINE: PDAF_ewpf_proposal_step --- Proposal step of EWPF
!
! !INTERFACE:
SUBROUTINE PDAF_ewpf_proposal_step(cnt_steps, nsteps, dim_p,&
     U_collect_state, U_prodRinvA, U_prodQA, U_obs_op, U_init_dim_obs,&
     U_init_obs, U_randvec, U_adjoint_obs_op, U_distribute_state, &
     outflag)

! !DESCRIPTION:
! Implementation of the proposal step for the EWPF.
!
! Variant for EWPF with domain decomposition.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2014-05 - Paul Kirchgessner - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_mod_ewpf, &
       ONLY: dim_ens, subtype_filter,step, &
       weights, step_obs, screen, eofV, member ,&
       local_dim_ens, HX_last, type_nudging, observation, dim_obs, &
       state_last, weight
  USE PDAF_mod_filtermpi, &
       ONLY: filterpe,mype_world, mype_couple 
 
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: cnt_steps ! Steps until next observation comes
  INTEGER, INTENT(in) :: nsteps    ! Steps between two observations
  INTEGER , INTENT(in) :: dim_p
  INTEGER, INTENT(out) :: outflag  ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_collect_state, &
       U_prodRinvA,&
       U_prodQA, &
       U_obs_op, & 
       U_init_dim_obs, &
       U_init_obs, &
       U_randvec, & 
       U_adjoint_obs_op, &
       U_distibute_state

! Local variables
  REAL, ALLOCATABLE :: state_p(:)
  REAL, ALLOCATABLE :: betan(:)           !Random error
  REAL, ALLOCATABLE :: Hymx(:)            !y-H(x^(n-1))
  REAL, ALLOCATABLE :: RiHx(:)            !R^-1 (x^(n-1))
  REAL, ALLOCATABLE :: QHTRiHx(:)         !QH^T(HQH^T+R)^(-1)(y-H(psi^(n-1)))
  REAL, ALLOCATABLE :: HTRiHx(:)      
  REAL :: weight_loc
  REAL :: nudge_strength
  INTEGER :: i
  REAL :: weights_gather(dim_ens)
  REAL :: inno1,inno2
  REAL :: weight_tmp(1)


 ! ****************************************
 ! *** Perform proposal step on all PEs ***
 ! ****************************************

  ! Allocate arrays
  ALLOCATE(state_p(dim_p))
  ALLOCATE(betan(dim_p))
  ALLOCATE(QHTRiHx(dim_p))
  ALLOCATE(HTRiHx(dim_p))
  ALLOCATE(Hymx(dim_obs))
  ALLOCATE(Rihx(dim_obs))
  IF( cnt_steps == 1) THEN
     ALLOCATE(state_last(dim_p))
  ENDIF

  ! Initialize arrays
  weight_loc = 0
  QHTRiHx = 0
  HTRiHx = 0

  ! Collect the current state
  CALL U_collect_state(dim_p, state_p)

  ! Initilize forcing strength
  CALL PDAF_ewpf_Btau(cnt_steps, 0, nsteps, nudge_strength, type_nudging)

  ! Initilize random forcing  
  CALL U_randvec(dim_p, betan, 0) 

  ! Initialize observations and weights at first timestep after the analysis
  doproposal: IF (cnt_steps == 1) THEN

     ! Allocate future observations for the use in proposal
     CALL U_init_dim_obs(step_obs, dim_obs)
     IF ( dim_obs > 0 ) THEN
        ALLOCATE(observation(dim_obs)) ! Allocate observations
        ALLOCATE(HX_last(dim_obs))     ! Prepare vector for future observations 
        CALL U_init_obs(step_obs, dim_obs, observation)
     ENDIF
     ! Initialize weights at first time step
     weights = -LOG(REAL(1./dim_ens))
      
  ELSEIF (cnt_steps > 1 ) THEN doproposal

     ! Compute innovation
     Hymx = observation - Hx_last
     inno1 = SQRT(SUM(HYmx**2)) 
 
     IF (nudge_strength > 0) THEN
        CALL U_prodRinvA(step_obs, dim_obs, 1, observation, Hymx,RiHx)
        CALL U_adjoint_obs_op(dim_obs, dim_p, dim_ens, RiHx, HTRiHx)
        CALL U_ProdQA(dim_p, HTRiHx, state_p, QHTRiHx)
     ENDIF

     state_p = state_p + nudge_strength*QHTRiHx + betan

     QHTRiHx = (nudge_strength**2)*QHTRiHx+ 2*nudge_strength*betan

     CALL dgemv('T', dim_p, 1, 1.0, HTRiHx, dim_p, &
          QHTRiHx, 1, 0, weight_loc, 1)
  
  ENDIF doproposal

  weight_loc = 0.5*weight_loc
  CALL U_distribute_state(dim_p,state_p)

  IF (.NOT.filterpe) THEN
     weight_tmp(1) = weight_loc
     weight = weight + weight_loc ! for new ewpf formulation
     CALL PDAF_ewpf_gather_weights(1, weight_tmp(1), screen)
  ELSE
     weights_gather = 0
     weights(1) = weights(1)+ weight_loc 
     weight = weight + weight_loc ! For new formulation of ewpf
     CALL PDAF_ewpf_gather_weights(dim_ens, weights_gather, screen)
     weights = weights + weights_gather
  END IF

! Normalize the weights on the filterpe
  IF ( filterpe ) THEN
   ! weights = weights-minval(weights)
   ! weights = exp(-weights)
   ! weights = weights/real(sum(weights)) 
   ! weights = -log(weights)
  ENDIF

  ! Compute Hx for the use in the next cycle. Use current state vector
  CALL U_obs_op(step_obs,dim_p,dim_obs, state_p,Hx_last)

  state_last = state_p
  

  IF (cnt_steps == nsteps-1) THEN
     DEALLOCATE(HX_last)
     DEALLOCATE(observation)
     DEALLOCATE(state_last)
  ENDIF

  outflag = 0

  DEALLOCATE(RiHX)
  DEALLOCATE(Hymx)
  DEALLOCATE(state_p, betan, QHTRiHx, HTRiHx)

END SUBROUTINE PDAF_ewpf_proposal_step

