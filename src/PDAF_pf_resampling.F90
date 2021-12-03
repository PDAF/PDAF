! Copyright (c) 2019-2021 Lars Nerger
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
! !ROUTINE: PDAF_pf_resampling --- Get particle indices for resampling
!
! !INTERFACE:
SUBROUTINE PDAF_pf_resampling(method, Nin, Nout, weights, IDs, screen)

! !DESCRIPTION:
! Determine particle indices for resampling. Implemented
! are three sampling schemes:
! (1) probabilistic resampling
! (2) stochastic unversal resampling
! (3) residual resampling
! The schemes follow the algorithms in 
! Vetra-Carvalho et al., State-of-the-art stochastic data
! assimilation methods for high-dimensional non-Gaussian problems.
! Tellus A, 70:1, 1445364, doi:10.1080/16000870.2018.1445364
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2019-05 - Lars Nerger initial code
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
  INTEGER, INTENT(in) :: method       ! Choose resampling method
                                      ! (1) probabilistic resampling
  INTEGER, INTENT(in) :: Nin          ! number of particles
  INTEGER, INTENT(in) :: Nout         ! number of particles to be resampled
  REAL, INTENT(in)    :: weights(Nin) ! Weights
  INTEGER, INTENT(out) :: IDs(Nout)   ! Indices of resampled ensmeble states
  INTEGER, INTENT(in) :: screen       ! Verbosity flag

! !CALLING SEQUENCE:
!EOP
       
! *** local variables ***
  INTEGER :: i, j                    ! Loop counters
  INTEGER :: c, Nr                   ! Counter for resampling
  INTEGER, SAVE :: first = 1         ! flag for init of random number seed
  INTEGER, SAVE :: iseed(4)          ! seed array for random number routine
  REAL :: rndval                     ! Random value
  REAL, ALLOCATABLE :: w_acc(:)      ! accumulated weights
  INTEGER, ALLOCATABLE :: w_i(:)     ! Integer weights
  REAL, ALLOCATABLE :: w_r(:)        ! residual weights


! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a)') &
          'PDAF', 'Resample particles'
     IF (method == 1) THEN
        WRITE (*, '(a, 5x, a)') &
             'PDAF', '--- Probabilistic resampling'
     ELSE IF (method == 2) THEN
        WRITE (*, '(a, 5x, a)') &
             'PDAF', '--- Stachastic universal resampling'
     ELSE IF (method == 3) THEN
        WRITE (*, '(a, 5x, a)') &
             'PDAF', '--- Residual resampling'
     END IF
  END IF

  IF (first==1) THEN
     ! Set random seed
     iseed(1)=15
     iseed(2)=25
     iseed(3)=47
     iseed(4)=17
     first = 0
  END IF

  ! Initialize accumulated weights
  ALLOCATE(w_acc(Nin))
  w_acc = 0.0

  ! Initialize rasampling IDs
  IDs = 0

   
  Rtype: IF (method == 1) THEN

! ********************************
! *** Probabilistic resampling ***
! ********************************

     ! Get accumulated weights
     w_acc(1) = weights(1)

     DO i = 2, Nin
        W_acc(i) = w_acc(i-1) + weights(i)
     END DO

     c = 1
     DO j = 1, Nout

        ! Init random number
        CALL larnvTYPE(1, iseed, 1, rndval)

        ! Determine index
        checkacc: DO i = 1, Nin
           IF (rndval > w_acc(i)) THEN
              c = c + 1
           ELSE
              EXIT checkacc
           END IF
        END DO checkacc

        IDs(j) = c
        c = 1

     END DO


  ELSE IF (method == 2) THEN

! ***************************************
! *** Stochastic universal resampling ***
! ***************************************

     ! Get accumulated weights
     w_acc(1) = weights(1)

     DO i = 2, Nin
        w_acc(i) = w_acc(i-1) + weights(i)
     END DO

     ! Init random number
     CALL larnvTYPE(1, iseed, 1, rndval)
     rndval = rndval / Nin


     c = 1

     DO j = 1, Nout

        ! Determine index
        checkaccB: DO i = 1, Nin
           IF (rndval > w_acc(i)) THEN
              c = c + 1
           ELSE
              EXIT checkaccB
           END if
        END DO checkaccB

        IDs(j) = c
        
        rndval = rndval + 1.0/REAL(Nin)

        c = 1

     END DO

  ELSE IF (method == 3) THEN

! ***************************
! *** Residual resampling ***
! ***************************

     ALLOCATE(w_i(Nin))
     ALLOCATE(w_r(Nin))

     DO i = 1, Nin
        w_i(i) = FLOOR(weights(i) * Nout)
     END DO

     ! First round setting resamping indices
     Nr = Nout
     c = 1
     IDs = 0
     DO j = 1, Nout
        IF (w_i(j) > 0) THEN
           IDs(c:c+w_i(j)-1) = j
           c = c + w_i(j)
           Nr = Nr - w_i(j)
        END IF
     END DO

     ! Now perform probabilistic resampling for remaining Nr indices

     ! Determine residual weights 

     IF (Nr > 0) THEN
        DO i = 1, Nin
           w_r(i) = (weights(i)*Nout - REAL(w_i(i))) / REAL (Nr)
        END DO
     END IF

     ! Get accumulated weights
     w_acc(1) = w_r(1)

     DO i = 2, Nin
        w_acc(i) = w_acc(i-1) + w_r(i)
     END DO

     c = 1
     DO j = Nout-Nr+1, Nout

        ! Init random number
        CALL larnvTYPE(1, iseed, 1, rndval)
        rndval = rndval * REAL(Nout) / REAL(Nin)

        ! Determine index
        checkaccC: DO i = 1, Nin
           IF (rndval > w_acc(i)) THEN
              c = c + 1
           ELSE
              EXIT checkaccC
           END IF
        END DO checkaccC

        IDs(j) = c
        c = 1

     END DO

     DEALLOCATE(w_i, w_r)

    END IF Rtype


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(w_acc)

END SUBROUTINE PDAF_pf_resampling
