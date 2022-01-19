! Copyright (c) 2004-2021 Lars Nerger, lars.nerger@awi.de
!
! This routine is free software: you can redistribute it and/or modify
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
! License along with this software.  If not, see <http://www.gnu.org/licenses/>.
!
!$Id$
!BOP
!
! !ROUTINE: PDAF_SampleEns --- Sample an ensemble from EOF modes
!
! !INTERFACE:
SUBROUTINE PDAF_SampleEns(dim, dim_ens, modes, svals, state, &
     ens, verbose, flag)

! !DESCRIPTION:
! This routine generates an ensemble of model states from a provided
! mean state and EOF modes (singular vectors of a peturbation matrix)
! and singular values. The resulting ensemble is the 2nd order-exact
! sample covariance matrix and mean state. 
! the ensemble state vectors are computed as
!   $ens_i = state + sqrt(dim_ens-1) modes (\Omega C)^T$
! where $C$ holds in its diagonal the singular values ($svals$). $\Omega$
! is an orthogonal transformation matrix that preserves the mean state.
! The generated ensemble fulfills the condition for the state error
! covariance matrix
!   $P = 1/(sqrt(dim_ens-1)  \sum_{i=1}^{dim\_ens} (ens_i - state)(ens_i - state)^T$
!
! !REVISION HISTORY:
! 2014-05 - Lars Nerger - Initial code for SANGOMA based on PDAF example code.
! 2016-05 - Lars Nerger - Back-porting to PDAF
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_mod_filter, &
       ONLY: Nm1vsN

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim                   ! Size of state vector
  INTEGER, INTENT(in) :: dim_ens               ! Size of ensemble
  REAL, INTENT(inout) :: modes(dim, dim_ens-1) ! Array of EOF modes
  REAL, INTENT(in)    :: svals(dim_ens-1)      ! Vector of singular values
  REAL, INTENT(inout) :: state(dim)            ! PE-local model state
  REAL, INTENT(out)   :: ens(dim, dim_ens)     ! State ensemble
  INTEGER, INTENT(in) :: verbose               ! Verbosity flag
  INTEGER, INTENT(inout) :: flag               ! Status flag
!EOP

! *** local variables ***
  INTEGER :: row, col                 ! counters
  REAL, ALLOCATABLE :: omega(:,:)     ! Transformation matrix Omega
  REAL :: fac                         ! Square-root of dim_ens-1


! **********************
! *** INITIALIZATION ***
! **********************

  IF (verbose>0) THEN
     WRITE (*,'(a, 5x,a)') 'PDAF', '******************************************************'
     WRITE (*,'(a, 5x,a)') 'PDAF', '*                  PDAF_SampleEns                    *'
     WRITE (*,'(a, 5x,a)') 'PDAF', '*                                                    *'
     WRITE (*,'(a, 5x,a)') 'PDAF', '*  Sample an ensemble with 2nd-order exact sampling  *'
     WRITE (*,'(a, 5x,a)') 'PDAF', '*    based on the Sangoma tool Sangoma_SampleEns.    *'
     WRITE (*,'(a, 5x,a)') 'PDAF', '******************************************************'

     ! *** Generate full ensemble on filter-PE 0 ***
     WRITE (*, '(/a, 5x, a)') 'PDAF', 'Sample state ensemble from covariance matrix'
     WRITE (*, '(a, 5x, a)') 'PDAF', 'given as EOF vectors and singular values'
     WRITE (*, '(a, 5x, a, i5)') 'PDAF', '--- Ensemble size:  ', dim_ens
     WRITE (*, '(a, 5x, a, i5)') 'PDAF', '--- number of EOFs: ', dim_ens-1
  END IF

  ! allocate memory for temporary fields
  ALLOCATE(omega(dim_ens, dim_ens-1))


! ********************************************************
! *** Generate ensemble by transformation of EOF modes ***
! ********************************************************

  ! *** Generate uniform orthogonal matrix OMEGA ***
  CALL PDAF_seik_omega(dim_ens-1, Omega, 1, 1)

  ! ***      Generate ensemble of states                  ***
  ! *** ens_i = state + sqrt(dim_ens-1) modes (Omega C)^T ***

  ! A = Omega C
  DO col = 1, dim_ens-1
     DO row = 1, dim_ens
        Omega(row, col) = Omega(row, col) * svals(col)
     END DO
  END DO
      
 ! ens = state + sqrt(dim_ens-1) modes A^T
  DO col = 1, dim_ens
     ens(:, col) = state(:)
  END DO

  IF (Nm1vsN == 1) THEN
     fac = SQRT(REAL(dim_ens-1))
  ELSE
     fac = SQRT(REAL(dim_ens))
  END IF
  CALL gemmTYPE('n', 't', dim, dim_ens, dim_ens-1, &
       fac, modes, dim, Omega, dim_ens, &
       1.0, ens, dim)


! ****************
! *** clean up ***
! ****************

  DEALLOCATE(omega)

  flag = 0

END SUBROUTINE PDAF_SampleEns
