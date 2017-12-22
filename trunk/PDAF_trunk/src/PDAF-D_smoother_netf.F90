! Copyright (c) 2004-1016 Lars Nerger
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
!$Id: PDAF-D_smoother_netf.F90 1672 2016-12-02 09:09:38Z lnerger $
!BOP
!
! !ROUTINE: PDAF_smoother_netf --- Smoother extension for NETF
!
! !INTERFACE:
SUBROUTINE PDAF_smoother_netf(dim_p, dim_ens, dim_lag, Ainv, sens_p, &
     cnt_maxlag, forget, screen)

! !DESCRIPTION:
! Smoother extension for the ensemble square-root filters (ETKF, ESTKF). 
! The routine uses the matrix Ainv computed by the filter analysis
! to perform the smoothing on past ensembles.
!
! Variant for domain decomposed states.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2016-08 - Lars Nerger - Initial code adating PDAF_smoother
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

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p        ! PE-local dimension of model state
  INTEGER, INTENT(in) :: dim_ens      ! Size of ensemble
  INTEGER, INTENT(in) :: dim_lag      ! Number of past time instances for smoother
  REAL, INTENT(in)   :: Ainv(dim_ens, dim_ens)  ! Weight matrix for ensemble transformation
  REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag)   ! PE-local smoother ensemble
  INTEGER, INTENT(inout) :: cnt_maxlag ! Count available number of time steps for smoothing
  REAL, INTENT(in)    :: forget       ! Forgetting factor
  INTEGER, INTENT(in) :: screen       ! Verbosity flag

! !CALLING SEQUENCE:
! Called by: PDAF_etkf_update, PDAF_estkf_update
! Calls: PDAF_timeit
! Calls: PDAF_memcount
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
!EOP

! *** local variables ***
  INTEGER :: member, col, row, lagcol ! Counters
  INTEGER :: n_lags                   ! Available number of time instances for smoothing
  INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
  INTEGER, SAVE :: allocflag = 0      ! Flag whether first time allocation is done
  REAL :: invdimens                   ! Inverse of global ensemble size
  REAL, ALLOCATABLE :: ens_blk(:,:)   ! Temporary block of state ensemble
  REAL, ALLOCATABLE :: W_smooth(:,:) ! Weight matrix for smoothing

  
! **********************
! *** INITIALIZATION ***
! **********************

  CALL PDAF_timeit(17, 'new')

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
            W_smooth(row, col) = Ainv(row, col)
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

  CALL PDAF_timeit(17, 'old')

END SUBROUTINE PDAF_smoother_netf
