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
! !ROUTINE: PDAF_EOFCovar --- Perform EOF decomposition of state trajectory
!
! !INTERFACE:
SUBROUTINE PDAF_eofcovar(dim, nstates, nfields, dim_fields, offsets, &
     remove_mstate, do_mv, states, stddev, svals, &
     svec, meanstate, verbose, status)

! !DESCRIPTION:
! This routine performs an EOF analysis by singular value decomposition. It is
! used to prepare a covariance matrix for initializing an ensemble.  For
! the decomposition a multivariate scaling can be performed by 
! 'PDAF\_MVNormalize' to ensure that all fields in the state vectors have
! unit variance. 
!
! To use this routine, one has to initialize the array 'states' holding in
! each column a perturbation vector (state - mean) from a state trajectory. 
! Outputs are the arrays of singular values (svals) and left singular vectors
! (svec). The singular values are scaled by sqrt(1/(nstates-1)). With this,
! $svec * svals^2 * svec^T$ is the covariance matrix. In addition, the standard
! deviation of the field variations (stddev) is an output array.
! To use the multivariate normalization one has to define the number of
! different fields in the state (nfields), the dimension of each fields and
! the offset of field from the start of each state vector.
!
! The routine uses the LAPACK routine 'dgesvd' to compute the singular value
! decomposition.
!
! !REVISION HISTORY:
! 2012-09 - L. Nerger - Initial code for SANGOMA based on PDAF
! 2013-11 - L. Nerger - Adaption to SANGOMA data model
! 2016-05 - L. Nerger - Back-porting to PDAF
! 2019-11 - L. Nerger - Clarification that 'states' is destroyed by the SVD
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim                 ! Dimension of state vector
  INTEGER, INTENT(in) :: nstates             ! Number of state vectors
  INTEGER, INTENT(in) :: nfields             ! Number of fields in state vector
  INTEGER, INTENT(in) :: dim_fields(nfields) ! Size of each field
  INTEGER, INTENT(in) :: offsets(nfields)    ! Start position of each field
  INTEGER, INTENT(in) :: do_mv               ! 1: Do multivariate scaling; 0: no scaling
     ! nfields, dim_fields and offsets are only used if do_mv=1
  INTEGER, INTENT(in) :: remove_mstate       ! 1: subtract mean state from states
     ! before computing EOFs; 0: don't remove and don't touch meanstate
  REAL, INTENT(inout)  :: states(dim, nstates)  ! State perturbations
     ! the array content will be destroyed by the singular value decomposition
  REAL, INTENT(out) :: stddev(nfields)       ! Standard deviation of field variability
     ! Without multivariate scaling (do_mv=0), it is stddev = 1.0
  REAL, INTENT(out) :: svals(nstates)        ! Singular values divided by sqrt(nstates-1)
  REAL, INTENT(out) :: svec(dim, nstates)    ! Singular vectors
  REAL, INTENT(inout) :: meanstate(dim)      ! Mean state (only changed if remove_mstate=1)
  INTEGER, INTENT(in) :: verbose             ! Verbosity flag
  INTEGER, INTENT(out) :: status             ! Status flag
!EOP


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

END SUBROUTINE PDAF_eofcovar
