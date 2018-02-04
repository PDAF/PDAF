! Copyright (c) 2004-2016 Lars Nerger
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
!$Id: PDAF_enkf_omega.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_enkf_omega --- Generate random matrix with special properties
!
! !INTERFACE:
SUBROUTINE PDAF_enkf_omega(seed, r, dim_ens, omega, norm, &
     otype, screen)

! !DESCRIPTION:
! Generate a random matrix OMEGA for the initilization
! of ensembles for the EnKF.
!
! The following properties can be set:\\
! 1. Simply fill the matrix with random numbers from a 
!   Gaussian distribution with mean zero and unit 
!   variance. (This corresponds to the simple random 
!   initialization of EnKF.)\\
! 2. Constrain the columns of OMEGA to be of unit norm
!   (This corrects error in the variance estimates caused
!   by the finite number of random numbers.)\\
! 3. Constrain the columns of OMEGA to be of norm dim\_ens$^{-1/2}$
!   (This corrects variance errors as in 2)\\
! 4. Project columns of OMEGA to be orthogonal to the vector
!   $(1,....,1)^T$ by Householder reflections. (This assures
!   that the mean of the generated ensemble equals the
!   prescribed state estimate.)\\
! Property 4 can be combined with either property 2 or 3.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2005-04 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: seed(4)  ! Seed for random number generation
  INTEGER, INTENT(in) :: r        ! Approximated rank of covar matrix
  INTEGER, INTENT(in) :: dim_ens  ! Ensemble size
  REAL, INTENT(inout) :: omega(dim_ens,r)  ! Random matrix
  REAL, INTENT(inout) :: norm     ! Norm for ensemble transformation
  INTEGER, INTENT(in) :: otype    ! Type of omega:
                                  ! (1) Simple Gaussian random matrix
                                  ! (2) Columns of unit norm
                                  ! (3) Columns of norm dim_ens^(-1/2)
                                  ! (4) Projection orthogonal (1,..,1)^T
                                  ! (6) Combination of 2 and 4
                                  ! (7) Combination of 3 and 4
                                  ! (8) Rows of sum 0 and variance 1
  INTEGER, INTENT(in) :: screen  ! Control verbosity

! !CALLING SEQUENCE:
! Called by: Used-provided ensemble initialization routine for EnKF
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: larnvTYPE (BLAS; dlarnv or slarnv dependent on precision)
!EOP

!  *** local variables ***
  INTEGER :: i, j     ! Counters
  INTEGER, SAVE :: iseed(4)            ! seed array for random number routine
  REAL, ALLOCATABLE :: rndvec(:)       ! Vector of random numbers
  REAL, ALLOCATABLE :: house(:,:)      ! Householder matrix
  REAL, ALLOCATABLE :: Omega_tmp(:,:)  ! Temporary OMEGA for projection
  REAL :: colnorm     ! Norm of matrix column
  REAL :: rownorm     ! Norm of matrix row


! **********************
! *** INITIALIZATION ***
! **********************

  IF (seed(1) >= 0) THEN
     ! Use given seed
     iseed=seed
  END IF

  ALLOCATE(rndvec(dim_ens-1))

! *** Initialize random matrix ***

 DO i = 1, dim_ens
     ! Fill row-wise to be consistent with old sampling formulation
     CALL larnvTYPE(3, iseed, r, rndvec)
      Omega(i, :) = rndvec
  END DO
  
! *** Normalize columns of Omega ***
  normcols: IF (otype == 2 .OR. otype == 3 .OR. otype == 6 .OR. otype == 7) THEN
     IF ((screen > 0) .AND. (otype == 2 .OR. otype == 6)) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', '--- EnKF_omega: Normalize columns of random matrix'
     ELSE IF (screen > 0) THEN
        WRITE (*, '(a, 5x,a)') &
             'PDAF', '--- EnKF_omega: Normalize columns of random matrix to dim_ens^(-1/2)'
     END IF

     DO j = 1, r
        ! Compute norm
        colnorm = 0.0
        DO i = 1, dim_ens
           colnorm = colnorm + Omega(i, j)**2
        END DO
        IF (otype == 3 .OR. otype == 7) THEN
           ! Set column norm to 1/sqrt(dim_ens)
           colnorm = colnorm / REAL(dim_ens)
        END IF
        colnorm = SQRT(colnorm)

        ! Perform normalization
        DO i = 1, dim_ens
           Omega(i, j) = Omega(i, j) / colnorm
        END DO
     END DO
  END IF normcols


! *** Project columns orthogonal to (1,1,...,1)^T ***
  doproject: IF (otype == 4 .OR. otype == 6 .OR. otype == 7) THEN
     IF (screen > 0) &
          WRITE (*, '(a, 5x, a)') &
          'PDAF', '--- EnKF_omega: Project columns orthogonal to (1,...,1)^T'

     ALLOCATE(Omega_tmp(dim_ens, r))
     ALLOCATE(house(dim_ens, dim_ens))

     ! Store Omega
     Omega_tmp = Omega

     ! Initialize Householder matrix
     housecolb: DO j = 1, dim_ens
        houserowb: DO i = 1, dim_ens
           house(i, j) = -1.0 / REAL(dim_ens)
        END DO houserowb
     END DO housecolb
     DO j = 1, dim_ens
        house(j, j) = house(j, j) + 1.0
     END DO
     
     ! Perform reflection
     CALL gemmTYPE ('n', 'n', dim_ens, r, dim_ens, &
          1.0, house, dim_ens, omega_tmp, dim_ens, &
          0.0, omega, dim_ens)

     DEALLOCATE(Omega_tmp, house)

  END IF doproject

  rowzero: IF (otype == 8) THEN
     IF (screen > 0) &
          WRITE (*, '(a, 5x, a)') &
          'PDAF', '--- EnKF_omega: Ensure that row sums are zero'

     DO i = 1, dim_ens

        rownorm = 0.0

        DO j = 1, r
           rownorm = rownorm + Omega(i,j)
        ENDDO
        rownorm = rownorm / REAL(r)

        DO j = 1, r
           Omega(i,j) = Omega(i,j) - rownorm
        ENDDO
        
     END DO

     IF (screen > 0) &
          WRITE (*, '(a, 5x, a)') &
          'PDAF', '--- EnKF_omega: Ensure that variance in rows is one'

     DO i = 1, dim_ens

        rownorm = 0.0

        DO j = 1, r
           rownorm = rownorm + Omega(i,j)*Omega(i,j)
        ENDDO
        rownorm = rownorm / REAL(r-1)
        rownorm = SQRT(rownorm)

        DO j = 1, r
           Omega(i,j) = Omega(i,j) / rownorm
        ENDDO
        
     END DO
  END IF rowzero


! *** Initialize norm for ensemble transformation ***
  IF (otype == 1) THEN
     norm = 1.0
  ELSEIF (otype == 2) THEN
     norm = SQRT(REAL(dim_ens - 1))
  ELSEIF (otype == 3) THEN
     norm = 1.0
  ELSEIF (otype == 4) THEN
     norm = 1.0
  ELSEIF (otype == 6) THEN
     norm = SQRT(REAL(dim_ens - 1))
  ELSEIF (otype == 7) THEN
     norm = 1.0
  END IF


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(rndvec)

END SUBROUTINE PDAF_enkf_omega
