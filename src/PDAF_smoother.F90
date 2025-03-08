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
!> Smoother routines for PDAF framework and for square-root filters
!!
!! This module contains the routines related to smoothing
!! in PDAF and the smoother updates for global and local 
!! square-root filters.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! 2012-05 - Lars Nerger - Initial code
!! Other revisions - see repository log
!!
MODULE PDAF_smoother

CONTAINS
!> Smoother extension for square-root filters
!!
!! Smoother extension for the ensemble square-root filters (ETKF, ESTKF). 
!! The routine uses the matrix Ainv computed by the filter analysis
!! to perform the smoothing on past ensembles.
!!
!! Variant for domain decomposed states.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! 2012-05 - Lars Nerger - Initial code
!! Other revisions - see repository log
!!
SUBROUTINE PDAF_smoothing(dim_p, dim_ens, dim_lag, Ainv, sens_p, &
     cnt_maxlag, forget, screen)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_parallel, &
       ONLY: mype

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p         !< PE-local dimension of model state
  INTEGER, INTENT(in) :: dim_ens       !< Size of ensemble
  INTEGER, INTENT(in) :: dim_lag       !< Number of past time instances for smoother
  REAL, INTENT(in)   :: Ainv(dim_ens, dim_ens)           !< Weight matrix for ensemble transformation
  REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag) !< PE-local smoother ensemble
  INTEGER, INTENT(inout) :: cnt_maxlag !< Count available number of time steps for smoothing
  REAL, INTENT(in)    :: forget        !< Forgetting factor
  INTEGER, INTENT(in) :: screen        !< Verbosity flag

! *** local variables ***
  INTEGER :: member, col, row, lagcol  !< Counters
  INTEGER :: n_lags                    !< Available number of time instances for smoothing
  INTEGER :: maxblksize, blkupper, blklower  !< Variables for blocked ensemble update
  INTEGER, SAVE :: allocflag = 0       !< Flag whether first time allocation is done
  REAL :: invdimens                    !< Inverse of global ensemble size
  REAL, ALLOCATABLE :: ens_blk(:,:)    !< Temporary block of state ensemble
  REAL, ALLOCATABLE :: W_smooth(:,:)   !< Weight matrix for smoothing

  
! **********************
! *** INITIALIZATION ***
! **********************

  ! Determine number of time instances for smoothing
  IF (cnt_maxlag >= dim_lag) THEN
     ! Already performed enough analysis to smooth over full lag
     n_lags = dim_lag
  ELSE
     ! Not yet enough analysis steps to smoother over full lag
     n_lags = cnt_maxlag
  END IF

  IF (mype == 0 .AND. screen > 0 .AND. n_lags > 0) THEN
     WRITE (*, '(a, 5x, a, i8)') 'PDAF', 'Perform smoothing up to lag ', n_lags
  END IF


! **********************************************
! *** Compute transform matrix for smoothing ***
! ***                                        ***
! *** W_smooth = A_filter + 1_N              ***
! **********************************************

  havelag: IF (n_lags > 0) THEN

     invdimens = 1.0 / REAL(dim_ens)

     ALLOCATE(W_smooth(dim_ens, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

     DO col = 1, dim_ens
        DO row = 1, dim_ens
            W_smooth(row, col) = forget * Ainv(row, col) + invdimens
        END DO
     END DO
  

! **********************************************
! *** Perform smoothing                      ***
! *** Transform state ensemble at past times ***
! ***              a    f                    ***
! ***             X  = X  W_smooth           ***
! **********************************************

     ! Use block formulation for transformation
     maxblksize = 200
     ALLOCATE(ens_blk(maxblksize, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', maxblksize * dim_ens)
     lagcol=1

     ! *** Smooth for all available lags ***
     smoothing: DO lagcol = 1, n_lags

        ! Use block formulation for transformation
        blocking: DO blklower = 1, dim_p, maxblksize
           
           blkupper = MIN(blklower + maxblksize - 1, dim_p)

           ! Store former analysis ensemble
           DO member = 1, dim_ens
              ens_blk(1 : blkupper-blklower+1, member) &
                   = sens_p(blklower : blkupper, member, lagcol)
           END DO

           !                        a   f
           ! Transform ensemble:   X = X  W_smooth
           CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
                1.0, ens_blk(1, 1), maxblksize, W_smooth, dim_ens, &
                0.0, sens_p(blklower, 1, lagcol), dim_p)

        END DO blocking

     END DO smoothing

     DEALLOCATE(ens_blk, W_smooth)

  END IF havelag


! ********************
! *** Finishing up ***
! ********************
  
  ! Increment maxlag counter
  cnt_maxlag = cnt_maxlag + 1

  ! Set flag for memory counting
  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_smoothing

!-------------------------------------------------------------------------------
!> Smoother extension for local square-root filters
!!
!! Smoother extension for the ensemble square-root filters (ETKF, ESTKF). 
!! The routine uses the matrix Ainv computed by the filter analysis
!! to perform the smoothing on past ensembles.
!!
!! Variant for domain decomposed states.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2012-05 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_smoothing_local(domain_p, step, dim_p, dim_l, dim_ens, &
     dim_lag, Ainv, ens_l, sens_p, cnt_maxlag, &
     U_g2l_state, U_l2g_state, forget, screen)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_parallel, &
       ONLY: mype
#if defined (_OPENMP)
  USE omp_lib, &
       ONLY: omp_get_num_threads, omp_get_thread_num
#endif

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p      !< Current local analysis domain
  INTEGER, INTENT(in) :: step          !< Current time step
  INTEGER, INTENT(in) :: dim_p         !< PE-local dimension of model state
  INTEGER, INTENT(in) :: dim_l         !< State dimension on local analysis domain
  INTEGER, INTENT(in) :: dim_ens       !< Size of ensemble
  INTEGER, INTENT(in) :: dim_lag       !< Number of past time instances for smoother
  REAL, INTENT(in)   :: Ainv(dim_ens, dim_ens)  !< Weight matrix for ensemble transformation
  REAL, INTENT(inout) :: ens_l(dim_l, dim_ens)  !< local past ensemble (temporary)
  REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag)   !< PE-local smoother ensemble
  INTEGER, INTENT(inout) :: cnt_maxlag !< Count available number of time steps for smoothing
  REAL, INTENT(in)    :: forget        !< Forgetting factor
  INTEGER, INTENT(in) :: screen        !< Verbosity flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_g2l_state, & !< Get state on local ana. domain from global state
       U_l2g_state           !< Init full state from state on local analysis domain

! *** local variables ***
  INTEGER :: member, col, row, lagcol  ! Counters
  INTEGER :: n_lags                    ! Available number of time instances for smoothing
  INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  INTEGER, SAVE :: first = 1           ! Flag for very first call to routine
  INTEGER, SAVE :: domain_save = 1     ! Index of domain from last call to routine
  REAL :: invdimens                    ! Inverse of global ensemble size
  REAL, ALLOCATABLE :: ens_blk(:,:)    ! Temporary block of state ensemble
  REAL, ALLOCATABLE :: W_smooth(:,:)   ! Weight matrix for smoothing
  INTEGER, SAVE :: mythread, nthreads  ! Thread variables for OpenMP

!$OMP THREADPRIVATE(mythread, nthreads, allocflag, first, domain_save)

  
! **********************
! *** INITIALIZATION ***
! **********************

#if defined (_OPENMP)
  nthreads = omp_get_num_threads()
  mythread = omp_get_thread_num()
#else
  nthreads = 1
  mythread = 0
#endif

  ! Determine number of time instances for smoothing
  IF (cnt_maxlag >= dim_lag) THEN
     ! Already performed enough analysis to smooth over full lag
     n_lags = dim_lag
  ELSE
     ! Not yet enough analysis steps to smooth over full lag
     n_lags = cnt_maxlag
  END IF

  IF ((domain_p <= domain_save) .OR. (first == 1)) THEN
     IF (mype == 0 .AND. screen > 0 .AND. n_lags > 0 .AND. mythread==0) THEN
        WRITE (*, '(a, 5x, a, i8)') 'PDAF', 'Perform smoothing up to lag ', n_lags
     END IF
  END IF


! **********************************************
! *** Compute transform matrix for smoothing ***
! ***                                        ***
! *** W_smooth = A_filter + 1_N              ***
! **********************************************

  havelag: IF (n_lags > 0) THEN

     invdimens = 1.0 / REAL(dim_ens)

     ALLOCATE(W_smooth(dim_ens, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

     DO col = 1, dim_ens
        DO row = 1, dim_ens
           W_smooth(row, col) = forget*Ainv(row, col) + invdimens
        END DO
     END DO


! **********************************************
! *** Perform smoothing                      ***
! *** Transform state ensemble at past times ***
! ***              a    f                    ***
! ***             X  = X  W_smooth           ***
! **********************************************

     ! Use block formulation for transformation
     maxblksize = 200
     ALLOCATE(ens_blk(maxblksize, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', maxblksize * dim_ens)
     lagcol=1

     ! *** Smooth for all available lags ***
     smoothing: DO lagcol = 1, n_lags

        ! *** Get local ensemble ***
        CALL PDAF_timeit(15, 'new')
        DO member = 1, dim_ens
           CALL U_g2l_state(step, domain_p, dim_p, sens_p(:, member, lagcol), dim_l, &
                ens_l(:, member))
        END DO
        CALL PDAF_timeit(15, 'old')

        ! Use block formulation for transformation
        blocking: DO blklower = 1, dim_l, maxblksize
           
           blkupper = MIN(blklower + maxblksize - 1, dim_l)

           ! Store former analysis ensemble
           DO member = 1, dim_ens
              ens_blk(1 : blkupper-blklower+1, member) &
                   = ens_l(blklower : blkupper, member)
           END DO

           !                        a   f
           ! Transform ensemble:   X = X  W_smooth
           CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
                1.0, ens_blk(1, 1), maxblksize, W_smooth, dim_ens, &
                0.0, ens_l(blklower, 1), dim_l)

        END DO blocking

        ! *** Initialize global ensemble ***
        CALL PDAF_timeit(16, 'new')
        DO member = 1, dim_ens
           CALL U_l2g_state(step, domain_p, dim_l, ens_l(:, member), dim_p, &
                sens_p(:, member, lagcol))
        END DO
        CALL PDAF_timeit(16, 'old')

     END DO smoothing
     
     DEALLOCATE(ens_blk, W_smooth)

  END IF havelag


! ********************
! *** Finishing up ***
! ********************
  
  ! Increment maxlag counter
  IF ((domain_p <= domain_save) .OR. (first == 1)) THEN

     cnt_maxlag = cnt_maxlag + 1
          
     ! Set flag
     first = 0
  END IF
  domain_save = domain_p

  ! Set flag for memory counting
  IF (allocflag == 0) allocflag = 1

END SUBROUTINE PDAF_smoothing_local

!-------------------------------------------------------------------------------
!> Shift ensemble states in ensemble array for smoothing
!!
!! Routine to store a previous analysis ensemble for smoothing.
!! The storage is performed by shifting all previous ensembles
!! by one. Thus the previous ensembles are stored consecutively.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! 2012-05 - Lars Nerger - Initial code
!! Other revisions - see repository log
!!
SUBROUTINE PDAF_smoother_shift(dim_p, dim_ens, dim_lag, ens_p, sens_p, cnt_maxlag, screen)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_mod_parallel, &
       ONLY: mype

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p         !< PE-local dimension of model state
  INTEGER, INTENT(in) :: dim_ens       !< Size of ensemble
  INTEGER, INTENT(in) :: dim_lag       !< Number of past time instances for smoother
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens, 1)         !< PE-local state ensemble
  REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag)  !< PE-local smoother ensemble
  INTEGER, INTENT(inout) :: cnt_maxlag !< Count available number of time steps for smoothing
  INTEGER, INTENT(in) :: screen        !< Verbosity flag

! *** local variables ***
  INTEGER :: lag             ! Counters
  INTEGER :: n_lags          ! Available number of tiem instances for smoothing
  
  
! **********************
! *** INITIALIZATION ***
! **********************

  ! Determine number of time instances for smoothing
  IF (cnt_maxlag >= dim_lag) THEN
     ! Already performed enough analysis to smooth over full lag
     n_lags = dim_lag
  ELSE
     ! Not yet enough analysis steps to smoothe over full lag
     n_lags = cnt_maxlag
  END IF

  IF (mype == 0 .AND. screen > 0 .AND. n_lags > 0) THEN
     WRITE (*, '(a, 5x, a)') 'PDAF', 'Store previous analysis for smoother'
  END IF


! ***********************
! *** Shift ensembles ***
! ***********************

  ! Shift past ensembles
  DO lag = n_lags-1, 1, -1

     IF (mype == 0 .AND. screen > 2) &
          write (*,*) 'PDAF: smoother: shift column', lag, 'to ',lag+1

     sens_p(:,:,lag+1) = sens_p(:,:,lag)

  END DO

  ! Store current ensemble
  IF (n_lags > 0) THEN

     IF (mype == 0 .AND. screen > 2) &
          write (*,*) 'PDAF: smoother: store current ensemble in smoother array'

     sens_p(:,:,1) = ens_p(:,:,1)
  END IF


! ********************
! *** Finishing up ***
! ********************

END SUBROUTINE PDAF_smoother_shift

END MODULE PDAF_smoother
