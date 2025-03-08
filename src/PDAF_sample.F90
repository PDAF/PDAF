! Copyright (c) 2004-2025 Lars Nerger, lars.nerger@awi.de
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
!
!> Module provideing routines for ensmeble sampling
!!
MODULE PDAF_sample

CONTAINS

!> Perform EOF decomposition of state trajectory
!!
!! This routine performs an EOF analysis by singular value decomposition. It is
!! used to prepare a covariance matrix for initializing an ensemble.  For
!! the decomposition a multivariate scaling can be performed by 
!! 'PDAF\_MVNormalize' to ensure that all fields in the state vectors have
!! unit variance. 
!!
!! To use this routine, one has to initialize the array 'states' holding in
!! each column a perturbation vector (state - mean) from a state trajectory. 
!! Outputs are the arrays of singular values (svals) and left singular vectors
!! (svec). The singular values are scaled by sqrt(1/(nstates-1)). With this,
!! $svec * svals^2 * svec^T$ is the covariance matrix. In addition, the standard
!! deviation of the field variations (stddev) is an output array.
!! To use the multivariate normalization one has to define the number of
!! different fields in the state (nfields), the dimension of each fields and
!! the offset of field from the start of each state vector.
!!
!! The routine uses the LAPACK routine 'dgesvd' to compute the singular value
!! decomposition.
!!
!! __Revision history:__
!! * 2012-09 - L. Nerger - Initial code for SANGOMA based on PDAF
!! * 2013-11 - L. Nerger - Adaption to SANGOMA data model
!! * 2016-05 - L. Nerger - Back-porting to PDAF
!! * 2019-11 - L. Nerger - Clarification that 'states' is destroyed by the SVD
!!
SUBROUTINE PDAF_eofcovar(dim, nstates, nfields, dim_fields, offsets, &
     remove_mstate, do_mv, states, stddev, svals, &
     svec, meanstate, verbose, status)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_mod_core, &
       ONLY: debug

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim                 !< Dimension of state vector
  INTEGER, INTENT(in) :: nstates             !< Number of state vectors
  INTEGER, INTENT(in) :: nfields             !< Number of fields in state vector
  INTEGER, INTENT(in) :: dim_fields(nfields) !< Size of each field
  INTEGER, INTENT(in) :: offsets(nfields)    !< Start position of each field
  INTEGER, INTENT(in) :: do_mv               !< 1: Do multivariate scaling; 0: no scaling
     !< nfields, dim_fields and offsets are only used if do_mv=1
  INTEGER, INTENT(in) :: remove_mstate       !< 1: subtract mean state from states
     !< before computing EOFs; 0: don't remove and don't touch meanstate
  REAL, INTENT(inout)  :: states(dim, nstates)  !< State perturbations
     !< the array content will be destroyed by the singular value decomposition
  REAL, INTENT(out) :: stddev(nfields)       !< Standard deviation of field variability
     !< Without multivariate scaling (do_mv=0), it is stddev = 1.0
  REAL, INTENT(out) :: svals(nstates)        !< Singular values divided by sqrt(nstates-1)
  REAL, INTENT(out) :: svec(dim, nstates)    !< Singular vectors
  REAL, INTENT(inout) :: meanstate(dim)      !< Mean state (only changed if remove_mstate=1)
  INTEGER, INTENT(in) :: verbose             !< Verbosity flag
  INTEGER, INTENT(out) :: status             !< Status flag

! *** local variables ***
  INTEGER :: i                      ! Counter
  INTEGER :: stat                   ! internal status flag
  INTEGER :: ldwork                 ! variable for SVD routine 
  REAL, ALLOCATABLE :: work(:)      ! work array for SVD
  REAL :: svdV                      ! right singular values (not referenced)

  
! **********************
! *** INITIALIZATION ***
! **********************

  IF (verbose>0) THEN
     WRITE (*,'(a, 5x,a)') 'PDAF', '***********************************************'
     WRITE (*,'(a, 5x,a)') 'PDAF', '*                PDAF_EOFCovar                *'
     WRITE (*,'(a, 5x,a)') 'PDAF', '*                                             *'
     WRITE (*,'(a, 5x,a)') 'PDAF', '*    Compute EOF decomposition of a matrix    *'
     WRITE (*,'(a, 5x,a)') 'PDAF', '* based on the Sangoma tool Sangoma_EOFCovar. *'
     WRITE (*,'(a, 5x,a)') 'PDAF', '***********************************************'
  END IF

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_eofcovar -- START'
     WRITE (*,*) '++ PDAF-debug PDAF_eofcovar:', debug, '  dim', dim
     WRITE (*,*) '++ PDAF-debug PDAF_eofcovar:', debug, '  nstates', nstates
     WRITE (*,*) '++ PDAF-debug PDAF_eofcovar:', debug, '  nfields', nfields
     WRITE (*,*) '++ PDAF-debug PDAF_eofcovar:', debug, '  dim_fields', dim_fields
     WRITE (*,*) '++ PDAF-debug PDAF_eofcovar:', debug, '  offsets', offsets
     WRITE (*,*) '++ PDAF-debug PDAF_eofcovar:', debug, '  states(1,:)', states(1,:)
     WRITE (*,*) '++ PDAF-debug PDAF_eofcovar:', debug, &
          '  Note: If REAL values appear incorrect, please check if you provide them with the correct precision'
  END IF


! *************************
! *** Remove mean state ***
! *************************

  removemean: IF (remove_mstate == 1) THEN

     IF (verbose>0) &
          WRITE (*,'(/a,5x,a)') 'PDAF', 'EOFCOVAR: Compute and subtract mean state ----------------------'

     ! *** compute mean state ***
     meanstate = 0.0
     DO i = 1, nstates
        meanstate(:) = meanstate(:) + states(:, i) / REAL(nstates)
     END DO

     ! *** get peturbation matrix ***
     DO i = 1, nstates
        states(:,i) = states(:,i) - meanstate(:)
     END DO

  END IF removemean


! ******************************************
! *** Perform multivariate normalization ***
! ******************************************

  stat = 0

  multivar: IF (do_mv == 1) THEN

     IF (verbose>0) &
          WRITE (*,'(/a,5x,a)') 'PDAF', 'EOFCOVAR: Perform multivariate normalization -------------------'

     DO i = 1, nfields
        CALL PDAF_MVNormalize(1, dim, dim_fields(i), offsets(i), nstates, &
             states, stddev(i), status)

        if (verbose>0) &
             WRITE (*,'(a,5x,a,i5,a,es12.4)') 'PDAF', 'Field', i, ': standard deviation ', stddev(i)

        stat = stat + status
     END DO

  ELSE
     stddev = 1.0
  END IF multivar


! *********************************************************
! *** Singular value decomposition of covariance matrix ***
! ***                                                   ***
! *** The covariance matrix is given by the state       ***
! *** sequences X of k states as                        ***
! ***          -1    _     _ T        T                 ***
! *** P = (k-1)   (X-X) (X-X)  = U L U      (EVP)       ***
! ***                                                   ***
! *** We compute the singular value decomposition       ***
! ***     _        T            -1    2  T              ***
! ***   X-X = U S V ;  P = (k-1)   U S  U               ***
! ***                                                   ***
! ***                         -1/2                      ***
! *** and we return U and (k-1)    S.                   ***
! *********************************************************

  IF (stat==0) THEN

     IF (verbose>0) &
          WRITE (*,'(/a,5x,a)') 'PDAF', 'EOFCOVAR: Compute SVD ------------------------------------------'

     ! Allocate work array and work space size
     ALLOCATE(work(MAX(3 * MIN(dim, nstates) + &
          MAX(dim, nstates), 5 * MIN(dim, nstates))))
     ldwork = MAX(3 * MIN(dim, nstates) + &
          MAX(dim, nstates), 5 * MIN(dim, nstates))
    
     ! Do decomposition
     CALL gesvdTYPE('s', 'n', dim, nstates, states, &
          dim, svals, svec, dim, svdV, &
          dim, work, ldwork, status)
  
     ! *** scale singular values ***
     DO i = 1, nstates
        svals(i) = svals(i) / SQRT(REAL(nstates - 1))
     END DO

     stat = stat+status

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_eofcovar:', debug, '  svals', svals
     END IF

  END IF

! ********************************
! *** Rescale singular vectors ***
! ********************************

  do_rescale: IF (do_mv == 1 .AND. stat == 0) THEN

     IF (verbose>0) &
          WRITE (*,'(/a, 5x,a)') 'PDAF', 'EOFCOVAR: Re-scale singular vectors according to stddev --------'

     DO i = 1, nfields
        CALL PDAF_MVNormalize(2, dim, dim_fields(i), offsets(i), nstates-1, &
             svec, stddev(i), status)
     END DO

     stat = stat+status

  END IF do_rescale


! ****************
! *** Clean up ***
! ****************

  DEALLOCATE(work)

  status = stat

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_eofcovar -- END'
  END IF

END SUBROUTINE PDAF_eofcovar

!--------------------------------------------------------------------------
!> Perform multivariate normalization
!!
!! This routine performs multivariate normalization and re-scaling.
!! It has two modes:
!!
!! mode=1: 
!! In this case, the routine computes the standard deviation of a field
!! inside the array 'states' holding in each column a state vector. The standard
!! deviation is computed over all columns if the state vector array. Then, the
!! field is normalized for unit standard deviation by dividing the values by
!! the standard deviation. The standard deviation is provided on output
!! together with the scaled array 'states'
!!
!! mode=2:
!! In this case the input variable 'stddev' is used to rescale the
!! corresponding part of the array 'states'. Usually 'stddev' is obtained by
!! a call with mode=1 before. 
!!
!! __Revision history:__
!! * 2012-09 - Lars Nerger - Initial code for SANGOMA based on PDAF
!! * 2013-11 - L. Nerger - Adaption to SANGOMA data model
!! * 2016-05 - Lars Nerger - Back-porting to PDAF
!!
SUBROUTINE PDAF_mvnormalize(mode, dim_state, dim_field, offset, &
     ncol, states, stddev, status)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: mode       !< Mode: (1) normalize, (2) re-scale
  INTEGER, INTENT(in) :: dim_state  !< Dimension of state vector
  INTEGER, INTENT(in) :: dim_field  !< Dimension of a field in state vector
  INTEGER, INTENT(in) :: offset     !< Offset of field in state vector
  INTEGER, intent(in) :: ncol       !< Number of columns in array states
  REAL, INTENT(inout) :: states(dim_state, ncol)  !< State vector array
  REAL, INTENT(inout) :: stddev     !< Standard deviation of field
                                    !< stddev is output for mode=1 and input for mode=2
  INTEGER, INTENT(out) :: status    !< Status flag (0=success)


! *** local variables ***
  INTEGER :: i, j       ! Counters


  modes: IF (mode == 1) THEN

! *****************************
! *** Perform normalization ***
! *****************************

     ! *** Compute STDDEV of field ***
     stddev = 0.0

     DO j = 1, ncol
        DO i = 1, dim_field
           stddev = stddev + states(i+offset, j)**2
        END DO
     END DO

     stddev = SQRT(stddev / REAL(ncol * dim_field))

     ! *** Normalize field ***
     DO j = 1, ncol
        DO i = 1, dim_field
           states(i+offset, j) = states(i+offset, j) / stddev
        END DO
     END DO

     ! Set status flag for success
     status = 0

  ELSE IF (mode == 2) THEN modes

! **************************
! *** Perform re-scaling ***
! **************************

     DO j = 1, ncol
        DO i = 1, dim_field
           states(i+offset, j) = states(i+offset, j) * stddev
        END DO
     END DO

     ! Set status flag for success
     status = 0

  ELSE modes

     ! invalid 'mode'
     status = 1

  END IF modes


END SUBROUTINE PDAF_mvnormalize

!--------------------------------------------------------------------------
!> Sample an ensemble from EOF modes
!!
!! This routine generates an ensemble of model states from a provided
!! mean state and EOF modes (singular vectors of a peturbation matrix)
!! and singular values. The resulting ensemble is the 2nd order-exact
!! sample covariance matrix and mean state. 
!! the ensemble state vectors are computed as
!!   $ens_i = state + sqrt(dim_ens-1) modes (\Omega C)^T$
!! where $C$ holds in its diagonal the singular values ($svals$). $\Omega$
!! is an orthogonal transformation matrix that preserves the mean state.
!! The generated ensemble fulfills the condition for the state error
!! covariance matrix
!!   $P = 1/(sqrt(dim_ens-1)  \sum_{i=1}^{dim\_ens} (ens_i - state)(ens_i - state)^T$
!!
!! __Revision history:__
!! * 2014-05 - Lars Nerger - Initial code for SANGOMA based on PDAF example code.
!! * 2016-05 - Lars Nerger - Back-porting to PDAF
!!
SUBROUTINE PDAF_SampleEns(dim, dim_ens, modes, svals, state, &
     ens, verbose, flag)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_mod_core, &
       ONLY: debug
  USE PDAF_analysis_utils, &
       ONLY: PDAF_seik_Omega

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim                   !< Size of state vector
  INTEGER, INTENT(in) :: dim_ens               !< Size of ensemble
  REAL, INTENT(inout) :: modes(dim, dim_ens-1) !< Array of EOF modes
  REAL, INTENT(in)    :: svals(dim_ens-1)      !< Vector of singular values
  REAL, INTENT(inout) :: state(dim)            !< PE-local model state
  REAL, INTENT(out)   :: ens(dim, dim_ens)     !< State ensemble
  INTEGER, INTENT(in) :: verbose               !< Verbosity flag
  INTEGER, INTENT(inout) :: flag               !< Status flag

! *** local variables ***
  INTEGER :: row, col                 ! counters
  REAL, ALLOCATABLE :: Omega(:,:)     ! Transformation matrix Omega
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

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_sampleens -- START'
     WRITE (*,*) '++ PDAF-debug PDAF_samplens:', debug, '  dim', dim
     WRITE (*,*) '++ PDAF-debug PDAF_samplens:', debug, '  modes(1,:)', modes(1,:)
     WRITE (*,*) '++ PDAF-debug PDAF_samplens:', debug, '  svals', svals(:)
     WRITE (*,*) '++ PDAF-debug PDAF_samplens:', debug, &
          '  Note: If REAL values appear incorrect, please check if you provide them with the correct precision'
  END IF

  ! allocate memory for temporary fields
  ALLOCATE(Omega(dim_ens, dim_ens-1))


! ********************************************************
! *** Generate ensemble by transformation of EOF modes ***
! ********************************************************

  ! *** Generate uniform orthogonal matrix OMEGA ***
  CALL PDAF_seik_Omega(dim_ens-1, Omega, 1, verbose)

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

  fac = SQRT(REAL(dim_ens-1))
  CALL gemmTYPE('n', 't', dim, dim_ens, dim_ens-1, &
       fac, modes, dim, Omega, dim_ens, &
       1.0, ens, dim)


! ****************
! *** clean up ***
! ****************

  DEALLOCATE(Omega)

  flag = 0

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_sampleens -- END'
  END IF

END SUBROUTINE PDAF_SampleEns

END MODULE PDAF_sample
