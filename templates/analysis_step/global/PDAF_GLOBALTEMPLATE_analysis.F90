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
!> Analysis routine for DA method GLOBALTEMPLATE
!!
!! This routine computes the analysis of the DA method.
!! Thus, for ensemble DA, it transforms the forecast ensemble
!! into the analysis ensemble.
!!
!! ADAPTING THE TEMPLATE:
!! This template contains a few typical steps of ensemble filters
!! On this basis one can implement another DA method. Below we 
!! describe the steps that are included in this code template.
!!
!! __Revision history:__
!! * 2024-12 - Lars Nerger - Initial code for template based on ETKF
!! * Later revisions - see repository log
!!
SUBROUTINE PDAF_GLOBALTEMPLATE_analysis(step, dim_p, dim_obs_p, dim_ens, &
     state_p, Ainv, ens_p, forget, &
     U_init_dim_obs, U_obs_op, U_init_obs, U_prodRinvA, &
     screen, flag)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE mpi
  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: mype, MPIerr, COMM_filter
  USE PDAF_mod_filter, &
       ONLY: type_trans, obs_member, debug
  USE PDAFomi, &
       ONLY: omi_n_obstypes => n_obstypes, omi_omit_obs => omit_obs

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step         ! Current time step
  INTEGER, INTENT(in) :: dim_p        ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_p   ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens      ! Size of ensemble
  REAL, INTENT(out)   :: state_p(dim_p)          ! on exit: PE-local forecast state
  REAL, INTENT(out)   :: Ainv(dim_ens, dim_ens)  ! on entry: uninitialized
                                      ! on exit: weight matrix for ensemble transformation
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  REAL, INTENT(in)    :: forget       ! Forgetting factor
  INTEGER, INTENT(in) :: screen       ! Verbosity flag
  INTEGER, INTENT(inout) :: flag      ! Status flag

! *** External subroutines ***
! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_dim_obs, &  ! Initialize dimension of observation vector
       U_obs_op, &               ! Observation operator
       U_init_obsvar, &          ! Initialize mean observation error variance
       U_init_obs, &             ! Initialize observation vector
       U_prodRinvA               ! Provide product R^-1 A

! *** Local variables ***
  INTEGER :: i, member, col, row      ! Counters
  INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
  REAL :: invdimens                   ! Inverse overall ensemble size
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  REAL, ALLOCATABLE :: HX_p(:,:)      ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: HXmean_p(:)    ! PE-local observed state
  REAL, ALLOCATABLE :: innov_p(:)     ! PE-local observation innovation
  REAL, ALLOCATABLE :: obs_p(:)       ! PE-local observation vector
  REAL, ALLOCATABLE :: RiHZ_p(:,:)    ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: tmp_Ainv(:,:)  ! Temporary storage of Ainv
  REAL, ALLOCATABLE :: rndmat(:,:)    ! Temporary random matrix
  REAL, ALLOCATABLE :: ens_blk(:,:)   ! Temporary block of state ensemble
  

! **********************
! *** INITIALIZATION ***
! **********************

  CALL PDAF_timeit(51, 'new')

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 1x, i7, 3x, a)') &
          'PDAF', step, 'Assimilating observations - GLOBALTEMPLATE'
  END IF


! ***********************************
! *** Compute mean forecast state ***
! ***********************************

  CALL PDAF_timeit(11, 'new')

  state_p = 0.0
  invdimens = 1.0 / REAL(dim_ens)
  DO member = 1, dim_ens
     DO row = 1, dim_p
        state_p(row) = state_p(row) + invdimens * ens_p(row, member)
     END DO
  END DO
  
  CALL PDAF_timeit(11, 'old')
  CALL PDAF_timeit(51, 'new')


! *********************************
! *** Get observation dimension ***
! *********************************

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++ TEMPLATE:                                                 +++
! +++ The number of observations is initialized here. With OMI  +++
! +++ this is the number of observation for the sub-domain of   +++
! +++ an MPI process.                                           +++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  CALL PDAF_timeit(15, 'new')
  CALL U_init_dim_obs(step, dim_obs_p)
  CALL PDAF_timeit(15, 'old')

  IF (screen > 2) THEN
     WRITE (*, '(a, 5x, a13, 1x, i6, 1x, a, i10)') &
          'PDAF', '--- PE-domain', mype, 'dimension of observation vector', dim_obs_p
  END IF


! **************************
! *** Compute innovation ***
! ***    d = y - H x     ***
! **************************

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++ TEMPLATE:                                                  +++
! +++ We include this part in the template because it is generic +++
! +++ and will likely be needed for most ensemble DA methods     +++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  CALL PDAF_timeit(12, 'new')
  
  haveobsA: IF (dim_obs_p > 0) THEN
     ! *** The innovation only exists for domains with observations ***

     ! Allocate observation-related arrays
     ALLOCATE(innov_p(dim_obs_p))
     ALLOCATE(obs_p(dim_obs_p))
     ALLOCATE(HX_p(dim_obs_p, dim_ens))
     ALLOCATE(HXmean_p(dim_obs_p))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3 * dim_obs_p + dim_obs_p*dim_ens)

     ! Compute observed ensemble states
     ! apply H to each ensemble state; then average (correct approach for nonlinear H)
     CALL PDAF_timeit(44, 'new')
     ENS1: DO member = 1, dim_ens
        ! Store member index to make it accessible with PDAF_get_obsmemberid
        obs_member = member

        ! Apply observation operator to ensemble member
        CALL U_obs_op(step, dim_p, dim_obs_p, ens_p(:, member), HX_p(:, member))
     END DO ENS1
     CALL PDAF_timeit(44, 'old')

     ! Compute mean of observed ensemble
     CALL PDAF_timeit(51, 'new')

     HXmean_p = 0.0
     DO member = 1, dim_ens
        DO row = 1, dim_obs_p
           HXmean_p(row) = HXmean_p(row) + invdimens * HX_p(row, member)
        END DO
     END DO

     CALL PDAF_timeit(51, 'old')

     ! get observation vector
     CALL PDAF_timeit(50, 'new')
     CALL U_init_obs(step, dim_obs_p, obs_p)
     CALL PDAF_timeit(50, 'old')

     ! Get innovation as difference of observation and mean observed state
     CALL PDAF_timeit(51, 'new')
     innov_p = obs_p - HXmean_p
     CALL PDAF_timeit(51, 'old')

     ! Omit observations with too high innovation
     IF (omi_omit_obs)  THEN
        CALL PDAF_timeit(51, 'new')
        CALL PDAFomi_omit_by_inno_cb(dim_obs_p, innov_p, obs_p)
        CALL PDAF_timeit(51, 'old')
     END IF

  ELSE IF (dim_obs_p == 0) THEN

     ! We need to call the observation operator also for dim_obs_p=0.
     ! This allows for observations in which the obervation operator
     ! include a global MPI communication
     ALLOCATE(HX_p(1,1))
     DO member = 1, dim_ens
        ! Store member index to make it accessible with PDAF_get_obsmemberid
        obs_member = member

        ! Call observation operator (it will only set a pointer for OMI)
        CALL U_obs_op(step, dim_p, dim_obs_p, ens_p(:, member), HX_p(:, member))
     END DO

  END IF haveobsA

  CALL PDAF_timeit(12, 'old')


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++ TEMPLATE:                                                    +++
! +++ At this point the observations, the observed ensemble states +++
! +++ and the innovation are initialized. Now one can use these    +++
! +++ to e.g. calculate a transformation matrix for the ensemble,  +++
! +++ which is specific for each DA method.                        +++
! +++ Below we just describe some routines that can be used in     +++
! +++ the calculations.                                            +++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! **********************************************
! ***  Calculate transformation matrix Ainv  ***
! **********************************************

  CALL PDAF_timeit(10, 'new')

  haveobsB: IF (dim_obs_p > 0) THEN
     ! *** The contribution of observation matrix ist only ***
     ! *** computed for domains with observations          ***

! +++ TEMPLATE: To compute the array of ensemble perturbations
! +++ one can use PDAF_subtract_rowmean which replaces the values of 
! +++ the input array HX_p by the perturbations.

     ! Subtract ensemble mean: HX = [Hx_1 ... Hx_N] T
     CALL PDAF_timeit(51, 'new')
     CALL PDAF_subtract_rowmean(dim_obs_p, dim_ens, HX_p)
     CALL PDAF_timeit(51, 'old')

     ! ***                RiHZ = Rinv HZ                
     ! *** This is implemented as a subroutine thus that
     ! *** Rinv does not need to be allocated explicitly.

! +++ TEMPLATE: All DA methods take observation error into account.
! +++ One variant, e.g. used in ESTKF and ETKF is to multiply with
! +++ the inverse observation error covariance matrix. The routine
! +++ U_prodRinvA compute this product. This is a call-back routine.
! +++ If PDAF-OMI is used this routine is provided by OMI.

     ALLOCATE(RiHZ_p(dim_obs_p, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_p * dim_ens)

     CALL PDAF_timeit(48, 'new')
     CALL U_prodRinvA(step, dim_obs_p, dim_ens, obs_p, HX_p, RiHZ_p)
     CALL PDAF_timeit(48, 'old')

! +++ TEMPLATE: Note that the template does not compute the full
! +++ matrix Ainv here because the exact operation might depend on 
! +++ the DA method. 

  ELSE haveobsB
     ! *** For domains with dim_obs_p=0 there is no ***
     ! *** direct observation-contribution to Ainv  ***

     CALL PDAF_timeit(51, 'new')

! +++ TEMPLATE: Note that this implementation is what is done
! +++ in the ETKF. It can be specific for the DA method.

     ! *** Initialize Ainv = (N-1) I ***
     Ainv = 0.0
     DO i = 1, dim_ens
        Ainv(i, i) = REAL(dim_ens - 1)
     END DO

     CALL PDAF_timeit(51, 'old')

  END IF haveobsB

  CALL PDAF_timeit(51, 'new')

! +++ TEMPLATE: The calculation of Ainv before is local for each
! +++ MPI process. One needs a global sum to get the overall value

  ! get total sum on all filter PEs
  ALLOCATE(tmp_Ainv(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

  CALL MPI_allreduce(Ainv, tmp_Ainv, dim_ens**2, &
       MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)

! +++ TEMPLATE: Note that the calculation here at first would give
! +++ the actual matrix A and one would need to invert it. In the
! +++ ETKF this is done using an eigenvalue decomposition, which
! +++ is not included here.

  ! *** Complete computation of Ainv ***
  Ainv = tmp_Ainv

  CALL PDAF_timeit(51, 'old')

  CALL PDAF_timeit(10, 'old')


  ! Optional 
  ! Multiply by orthogonal random matrix with eigenvector (1,...,1)^T
  multrnd: IF (type_trans == 2) THEN
     CALL PDAF_timeit(51, 'new')

     WRITE (*,'(a, 5x, a)') 'PDAF', '--- Apply random rotation to ensemble'

! +++ TEMPLATE: One might want to apply a random rotation to Ainv
! +++ PDAF_general_rndmat provides a matrix with such random rotation
! +++ which is used in this example

     ALLOCATE(rndmat(dim_ens, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)
        
     ! Initialize random matrix
     CALL PDAF_generate_rndmat(dim_ens, rndmat, 2)

     CALL gemmTYPE('n', 'n', dim_ens, dim_ens, dim_ens, &
          1.0, Ainv, dim_ens, rndmat, dim_ens, &
          0.0, tmp_Ainv, dim_ens)

     Ainv = tmp_Ainv

     DEALLOCATE(rndmat)

     CALL PDAF_timeit(51, 'old')
  END IF multrnd


! ************************************************
! ***     Transform state ensemble             ***
! ***              a   _f   f                  ***
! ***             X  = X + X  W                ***
! *** The weight matrix W is stored in Ainv.   ***
! ************************************************

  CALL PDAF_timeit(51, 'new')

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a)') 'PDAF', 'Perform ensemble transformation'
  END IF

! +++ TEMPLATE: Finally one multiplies the ensemble (or ensemble perturbations)
! +++ with the transformation matrix Ainv. Since this overwrites the ensemble 
! +++ matrix ens_p we use he a blocked variant. A block of 'blocksize' rows
! +++ of ens_p is updates at a time so that only a small temporary matrix
! +++ (ens_blk) is required. The value of 'maxblksize' is harcoded, but could
! +++ be changed.
! +++ The example below computes the product and adds the forecast mean state.
! +++ Thus, the update of the ensemble mean and perturbations is done together.

  ! Use block formulation for transformation
  maxblksize = 200
  IF (mype == 0 .AND. screen > 0) &
       WRITE (*, '(a, 5x, a, i5)') &
       'PDAF', '--- use blocking with size ', maxblksize
        
  ALLOCATE(ens_blk(maxblksize, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', maxblksize * dim_ens)

  blocking: DO blklower = 1, dim_p, maxblksize
           
     blkupper = MIN(blklower + maxblksize - 1, dim_p)

     ! Store forecast ensemble
     CALL PDAF_timeit(21, 'new')
     DO col = 1, dim_ens
        ens_blk(1 : blkupper - blklower + 1, col) &
             = ens_p(blklower : blkupper, col)
     END DO

     ! Store mean forecast in ensemble matrix
     DO col = 1,dim_ens
        ens_p(blklower : blkupper, col) = state_p(blklower : blkupper)
     END DO
     CALL PDAF_timeit(21, 'old')

     !                        a  _f   f
     ! Transform ensemble:   X = X + X  Ainv
     CALL PDAF_timeit(22, 'new')

     CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
          1.0, ens_blk(1, 1), maxblksize, Ainv(1, 1), dim_ens, &
          1.0, ens_p(blklower, 1), dim_p)

     CALL PDAF_timeit(22, 'old')

  END DO blocking

  CALL PDAF_timeit(51, 'old')


! ********************
! *** Finishing up ***
! ********************

! +++ TEMPLATE: One should be careful to deallocate all allocated arrays
! +++ We collect most deallocates here

  IF (dim_obs_p > 0) THEN
     DEALLOCATE(obs_p)
     DEALLOCATE(HXmean_p)
     DEALLOCATE(innov_p)
     DEALLOCATE(RiHZ_p)
  END IF
  DEALLOCATE(HX_p)
  DEALLOCATE(tmp_Ainv)
  DEALLOCATE(ens_blk)

! +++ TEMPLATE: Below is generic operation that is required
! +++ memory counting work

  ! Set flag that allocation was already done once (used for memory counting)
  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_GLOBALTEMPLATE_analysis
