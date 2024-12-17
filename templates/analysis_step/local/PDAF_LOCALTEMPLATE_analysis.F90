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
!> Analysis routine for DA method LOCALTEMPLATE
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
SUBROUTINE PDAF_LOCALTEMPLATE_analysis(domain_p, step, dim_l, dim_obs_f, dim_obs_l, &
     dim_ens, state_l, Ainv_l, ens_l, HX_f, &
     HXmean_f, rndmat, forget, &
     U_g2l_obs, U_init_obs_l, U_prodRinvA_l, U_init_n_domains_p, &
     screen, flag)

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
       ONLY: type_trans, obs_member, debug
  USE PDAFomi, &
       ONLY: omi_omit_obs => omit_obs
#if defined (_OPENMP)
  USE omp_lib, &
       ONLY: omp_get_num_threads, omp_get_thread_num
#endif

  IMPLICIT NONE

! *** Arguments ***
! Variable naming scheme:
!    suffix _p: Denotes a full variable on the PE-local domain
!    suffix _l: Denotes a local variable on the current analysis domain
!    suffix _f: Denotes a full variable of all observations required for the
!               analysis loop on the PE-local domain
  INTEGER, INTENT(in) :: domain_p    ! Current local analysis domain
  INTEGER, INTENT(in) :: step        ! Current time step
  INTEGER, INTENT(in) :: dim_l       ! State dimension on local analysis domain
  INTEGER, INTENT(in) :: dim_obs_f   ! PE-local dimension of full observation vector
  INTEGER, INTENT(in) :: dim_obs_l   ! Size of obs. vector on local ana. domain
  INTEGER, INTENT(in) :: dim_ens     ! Size of ensemble 
  REAL, INTENT(inout) :: state_l(dim_l)           ! local forecast state
  REAL, INTENT(out)   :: Ainv_l(dim_ens, dim_ens) ! on entry: uninitialized
                               ! on exit: local weight matrix for ensemble transformation
  REAL, INTENT(inout) :: ens_l(dim_l, dim_ens)    ! Local state ensemble
  REAL, INTENT(in) :: HX_f(dim_obs_f, dim_ens) ! PE-local full observed state ens.
  REAL, INTENT(in) :: HXmean_f(dim_obs_f)       ! PE-local full observed ens. mean
  REAL, INTENT(inout) :: rndmat(dim_ens, dim_ens) ! Global random rotation matrix
  REAL, INTENT(inout) :: forget      ! Forgetting factor
  INTEGER, INTENT(in) :: screen      ! Verbosity flag
  INTEGER, INTENT(inout) :: flag     ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_g2l_obs, &   ! Restrict full obs. vector to local analysis domain
       U_init_obs_l, &       ! Init. observation vector on local analysis domain
       U_init_n_domains_p, & ! Provide PE-local number of local analysis domains
       U_prodRinvA_l         ! Provide product R^-1 A for local analysis domain
       
! *** Local variables ***
  INTEGER :: i, member, col, row       ! Counters
  INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
  INTEGER :: syev_info                 ! Status flag for SYEV
  INTEGER :: ldwork                    ! Size of work array for SYEV
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  INTEGER, SAVE :: lastdomain = -1     ! store domain index
  INTEGER, SAVE :: mythread, nthreads  ! Thread variables for OpenMP for controlling screen output
  LOGICAL, SAVE :: screenout = .true.  ! Whether to print information to stdout
  REAL, ALLOCATABLE :: HX_l(:,:)       ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: HXmean_l(:)     ! state projected onto obs. space
  REAL, ALLOCATABLE :: RiHZ_l(:,:)     ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: innov_l(:)      ! local observation innovation
  REAL, ALLOCATABLE :: obs_l(:)        ! local observation vector
  REAL, ALLOCATABLE :: tmp_Ainv_l(:,:) ! Temporary storage of Ainv
  REAL, ALLOCATABLE :: ens_blk(:,:)    ! Temporary block of state ensemble

!$OMP THREADPRIVATE(mythread, nthreads, lastdomain, allocflag, screenout)


! *******************
! *** Preparation ***
! *******************

  CALL PDAF_timeit(51, 'new')

! +++ TEMPLATE:
! +++ Here we determine whether information is printed to the screen
! +++ in case of using OpenMP parallelization only thread 0 is writing

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
        WRITE (*,'(a, 5x, a, i5, a)') &
             'PDAF', '--- Use OpenMP parallelization with ', nthreads, ' threads'
     END IF
#endif
  END IF

  CALL PDAF_timeit(51, 'old')


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

  haveobsB: IF (dim_obs_l > 0) THEN
     ! *** The innovation only exists for domains with observations ***

     ALLOCATE(innov_l(dim_obs_l))
     ALLOCATE(obs_l(dim_obs_l))
     ALLOCATE(HX_l(dim_obs_l, dim_ens))
     ALLOCATE(HXmean_l(dim_obs_l))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3 * dim_obs_l + dim_obs_l * dim_ens)

     ! Compute observed ensemble states
     ! Restrict ensemble of observed states onto local observation space
     CALL PDAF_timeit(46, 'new')
     CALL PDAF_timeit(30, 'new')
     ENS: DO member = 1, dim_ens
        ! Store member index
        obs_member = member

        ! [Hx_1 ... Hx_N] for local analysis domain
        CALL U_g2l_obs(domain_p, step, dim_obs_f, dim_obs_l, HX_f(:, member), &
             HX_l(:, member))
     END DO ENS
     CALL PDAF_timeit(30, 'old')

     ! Restrict mean observed state onto local observation space
     obs_member = 0
     CALL U_g2l_obs(domain_p, step, dim_obs_f, dim_obs_l, HXmean_f, HXmean_l)
     CALL PDAF_timeit(46, 'old')

     ! get local observation vector
     CALL PDAF_timeit(47, 'new')
     CALL U_init_obs_l(domain_p, step, dim_obs_l, obs_l)
     CALL PDAF_timeit(47, 'old')

     ! Get innovation as difference of observation and observed state
     CALL PDAF_timeit(51, 'new')
     innov_l = obs_l - HXmean_l
     CALL PDAF_timeit(51, 'old')

     ! Omit observations with too high innovation
     IF (omi_omit_obs) CALL PDAFomi_omit_by_inno_l_cb(domain_p, dim_obs_l, innov_l, obs_l)

  END IF haveobsB

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

! ***************************************************************
! ***  Calculate transformation matrix Ainv for local domain  ***
! ***************************************************************

  CALL PDAF_timeit(10, 'new')

  haveobsA: IF (dim_obs_l > 0) THEN
     ! *** The contribution of observation matrix ist only ***
     ! *** computed for domains with observations          ***

     CALL PDAF_timeit(30, 'new')
     CALL PDAF_timeit(51, 'new')

! +++ TEMPLATE: To compute the array of ensemble perturbations
! +++ one can use PDAF_subtract_rowmean which replaces the values of 
! +++ the input array HX_l by the perturbations.

     ! Subtract ensemble mean: HZ = [Hx_1 ... Hx_N] T
     CALL PDAF_subtract_rowmean(dim_obs_l, dim_ens, HX_l)

     CALL PDAF_timeit(51, 'old')
     CALL PDAF_timeit(30, 'old')


     ! ***                RiHZ = Rinv HZ                
     ! *** This is implemented as a subroutine thus that
     ! *** Rinv does not need to be allocated explicitly.

! +++ TEMPLATE: All DA methods take observation error into account.
! +++ One variant, e.g. used in LESTKF and LETKF is to multiply with
! +++ the inverse observation error covariance matrix. The routine
! +++ U_prodRinvA_l compute this product and in addition applied
! +++ localization weights to the observation error covariance matirx. 
! +++ This is a call-back routine.
! +++ If PDAF-OMI is used this routine is provided by OMI.

     ALLOCATE(RiHZ_l(dim_obs_l, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l * dim_ens)

     CALL PDAF_timeit(48, 'new')
     CALL U_prodRinvA_l(domain_p, step, dim_obs_l, dim_ens, obs_l, HX_l, RiHZ_l)
     CALL PDAF_timeit(48, 'old')

! +++ TEMPLATE: Note that the template does not compute the full
! +++ matrix Ainv here because the exact operation might depend on 
! +++ the DA method. 

  ELSE haveobsA
     ! *** For domains with dim_obs_l=0 there is no ***
     ! *** direct observation-contribution to Ainv  ***

     CALL PDAF_timeit(51, 'new')

     ! *** Initialize Ainv = (N-1) I ***
     Ainv_l = 0.0
     DO i = 1, dim_ens
        Ainv_l(i, i) = REAL(dim_ens - 1)
     END DO
     CALL PDAF_timeit(51, 'old')

  END IF haveobsA

  CALL PDAF_timeit(10, 'old')


  ! Optional 
  ! Multiply by orthogonal random matrix with eigenvector (1,...,1)^T
  multrnd: IF (type_trans == 2) THEN
     CALL PDAF_timeit(51, 'new')

! +++ TEMPLATE: One can applly here the random rotation matrix
! +++ that was initialized in PDAF_LOCALTEMPLATE_update.

     ALLOCATE(tmp_Ainv_l(dim_ens, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

     CALL gemmTYPE('n', 'n', dim_ens, dim_ens, dim_ens, &
          1.0, Ainv_l, dim_ens, rndmat, dim_ens, &
          0.0, tmp_Ainv_l, dim_ens)

     Ainv_l = tmp_Ainv_l

     DEALLOCATE(tmp_Ainv_l)

     CALL PDAF_timeit(51, 'old')
  END IF multrnd


! ************************************************
! ***     Transform state ensemble             ***
! ***              a   _f   f                  ***
! ***             X  = X + X  W                ***
! *** The weight matrix W is stored in Ainv_l. ***
! ************************************************

  CALL PDAF_timeit(51, 'new')

  IF (mype == 0 .AND. screen > 0 .AND. screenout) THEN
     WRITE (*, '(a, 5x, a)') 'PDAF', 'Perform ensemble transformation'
  END IF


! +++ TEMPLATE: Finally one multiplies the ensemble (or ensemble perturbations)
! +++ with the transformation matrix Ainv. Since this overwrites the ensemble 
! +++ matrix ens_p we use he a blocked variant. A block of 'blocksize' rows
! +++ of ens_p is updates at a time so that only a small temporary matrix
! +++ (ens_blk) is required. The value of 'maxblksize' is harcoded, but could
! +++ be changed. Note, that the blocking might not be required for the local
! +++ ensemble if it is small enough.
! +++ The example below computes the product and adds the forecast mean state.
! +++ Thus, the update of the ensemble mean and perturbations is done together.

  ! Use block formulation for transformation
  maxblksize = 200
  IF (mype == 0 .AND. screen > 0 .AND. screenout) &
       WRITE (*, '(a, 5x, a, i5)') &
       'PDAF', '--- use blocking with size ', maxblksize

  ALLOCATE(ens_blk(maxblksize, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(4, 'r', maxblksize * dim_ens)

  blocking: DO blklower = 1, dim_l, maxblksize

     blkupper = MIN(blklower + maxblksize - 1, dim_l)

     ! Store forecast ensemble
     CALL PDAF_timeit(21, 'new')
     DO col = 1, dim_ens
        ens_blk(1 : blkupper - blklower + 1, col) &
             = ens_l(blklower : blkupper, col)
     END DO

     ! Store mean forecast in ensemble matrix
     DO col = 1, dim_ens
        ens_l(blklower : blkupper, col) = state_l(blklower : blkupper)
     END DO
     CALL PDAF_timeit(21, 'old')

     !                        a  _f   f
     ! Transform ensemble:   X = X + X  TW
     CALL PDAF_timeit(22, 'new')

     CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
          1.0, ens_blk, maxblksize, Ainv_l, dim_ens, &
          1.0, ens_l(blklower, 1), dim_l)

     CALL PDAF_timeit(22, 'old')

  END DO blocking
        

  CALL PDAF_timeit(51, 'old')


! ********************
! *** Finishing up ***
! ********************

! +++ TEMPLATE: One should be careful to deallocate all allocated arrays
! +++ We collect most deallocates here

  IF (dim_obs_l > 0) THEN
     DEALLOCATE(obs_l)
     DEALLOCATE(HXmean_l)
     DEALLOCATE(HX_l)
     DEALLOCATE(innov_l)
     DEALLOCATE(RiHZ_l)
  END IF
  DEALLOCATE(ens_blk)

! +++ TEMPLATE: Below are generic operations that are required
! +++ to make screen output and memory counting work

  ! Set flag that allocation was already done once (used for memory counting)
  IF (allocflag == 0) allocflag = 1

  ! Store domain index to control screen output
  lastdomain = domain_p

END SUBROUTINE PDAF_LOCALTEMPLATE_analysis
