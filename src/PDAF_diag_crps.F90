! Copyright (c) 2012-2023 Lars Nerger, lars.nerger@awi.de
!
! This routine is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version
! 3 of the License, or (at your option) any later version.
!
! This code is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this software.  If not, see <http://www.gnu.org/licenses/>.
!
!$Id$

!> Computation of CRPS
!!
!! This routine computes the continuous ranked probability
!! score (CRPS) and its decomposition into uncertainty and
!! potential CRPS: CRPS = RELI + pot_CRPS. In addition the uncertainty
!! is computed.
!! Resolution can be computed by RESOL = UNCERT - pot_CRPS.
!! A perfectly reliable system gives RELI=0.
!! An informative system gives RESOL ~ UNCERT or pot_CRPS << UNCERT.
!!
!! The computation follows H. Hersbach, Weather and Forecasting
!! 15(2000) 599-570.
!!
!! __Revision history:__
!! * 2021-05 - Lars Nerger - Initial code based on sangoma_ComputeCRPS
!! * Later revisions - see repository log
!! * 2024-04 - Yumeng Chen - refactor; add domain decomposition support
!!

!---------------------------------------------------------------------------
!> CRPS diagnostic routine with original interface
!!
SUBROUTINE PDAF_diag_crps(dim_p, dim_ens, element, oens, obs, &
    CRPS, reli, pot_CRPS, uncert, status)!
#include "typedefs.h"

  USE mpi
  USE PDAF_mod_filtermpi, & 
       ONLY: COMM_filter, mype_filter, npes_filter
  USE SANGOMA_quicksort

  IMPLICIT NONE

  ! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p                !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens              !< Ensemble size
  INTEGER, INTENT(in) :: element              !< index of element in full state vector
       !< If element=0, mean values over dim_p grid points/cases are computed
  REAL, INTENT(in)    :: oens(dim_p, dim_ens) !< State ensemble
  REAL, INTENT(in)    :: obs(dim_p)           !< Observation / truth
  REAL, INTENT(out)   :: CRPS                 !< CRPS
  REAL, INTENT(out)   :: reli                 !< Reliability
  REAL, INTENT(out)   :: pot_CRPS             !< potential CRPS
  REAL, INTENT(out)   :: uncert               !< uncertainty
  INTEGER, INTENT(out) :: status              !< Status flag (0=success)

  CALL PDAF_diag_crps_mpi(dim_p, dim_ens, element, oens, obs, &
       COMM_filter, mype_filter, npes_filter, &
       CRPS, reli, pot_CRPS, uncert, status)!

END SUBROUTINE PDAF_diag_crps

!---------------------------------------------------------------------------
!> CRPS diagnostic routine including MPI-settings in interface
!!
SUBROUTINE PDAF_diag_crps_mpi(dim_p, dim_ens, element, oens, obs, &
    COMM_filter, mype_filter, npes_filter, &
    CRPS, reli, pot_CRPS, uncert, status)
#include "typedefs.h"

  USE mpi
  USE SANGOMA_quicksort

  IMPLICIT NONE

  ! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p                !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens              !< Ensemble size
  INTEGER, INTENT(in) :: element              !< index of element in full state vector
       !< If element=0, mean values over dim_p grid points/cases are computed
  INTEGER, INTENT(in) :: COMM_filter          !< MPI communicator for filter
  INTEGER, INTENT(in) :: mype_filter          !< rank of MPI communicator
  INTEGER, INTENT(in) :: npes_filter          !< size of MPI communicator
  REAL, INTENT(in)    :: oens(dim_p, dim_ens) !< State ensemble
  REAL, INTENT(in)    :: obs(dim_p)           !< Observation / truth
  REAL, INTENT(out)   :: CRPS                 !< CRPS
  REAL, INTENT(out)   :: reli                 !< Reliability
  REAL, INTENT(out)   :: pot_CRPS             !< potential CRPS
  REAL, INTENT(out)   :: uncert               !< uncertainty
  INTEGER, INTENT(out) :: status              !< Status flag (0=success)

  ! *** local variables ***
  INTEGER :: i, k    ! counter
  INTEGER :: istart  ! starting index of grid point/case
  INTEGER :: imax    ! end index of grid point/case
  INTEGER :: dim     ! dimension of the full state vector
  INTEGER :: MPIerr 
  REAL :: wk   ! weight for each grid point/case
  REAL :: x_a  ! truth / verifying analysis (observation)
  REAL :: gi   ! The difference between i+1-th ensemble member and i-th member
  REAL :: oi
  REAL :: o0, one_minus_oN
  REAL :: o0_p, one_minus_oN_p
  REAL :: pval
  INTEGER, ALLOCATABLE :: all_dim_p(:) ! dimensions of the local state vector
  INTEGER, ALLOCATABLE :: dis_dim_p(:) ! dimensions of the local state vector
  REAL, ALLOCATABLE :: one_case(:)   ! ensemble in each case, this is variable x in H. Hersbach (2000)
  REAL, ALLOCATABLE :: allobs(:)
  REAL, ALLOCATABLE :: alpha_p(:), beta_p(:)
  REAL, ALLOCATABLE :: alpha(:), beta(:)

  ! initialise the status flag
  status = 0

  ! initialise crps output
  crps = 0.0
  reli = 0.0
  pot_CRPS = 0.0
  uncert = 0.0

  ! allocate arrays for MPI communication
  ALLOCATE( all_dim_p(npes_filter), dis_dim_p(npes_filter) )
  ! gather the dimension of the local state vector to all_dim_p
  CALL MPI_Allgather(dim_p, 1, MPI_INTEGER, all_dim_p, 1, MPI_INTEGER, COMM_filter, MPIerr)
  ! displacement of the received array used for gatherv
  dis_dim_p(1) = 0
  DO i = 2, npes_filter
      dis_dim_p(i) = dis_dim_p(i - 1) + all_dim_p(i - 1)
  END DO

  ! dimension fo the full state vector
  dim = SUM(all_dim_p)

  ! Set number of element over which CPRS is computed
  IF (element==0) THEN
    istart = 1
    imax = dim_p
    ! weight for each grid point/case
    wk = 1.0/REAL(dim)
  ELSEIF (element<=dim) THEN
    ! error handling
    IF (element < 0) THEN
      status = 100
      WRITE(*, '(a, 5x, a, I4, a)') 'PDAF warning:', &
        'PDAF_diag_crps: element(', element, ') argument must be >= 0.'
      RETURN
    ENDIF
    !
    IF (element <= dis_dim_p(mype_filter + 1) .OR. element > dis_dim_p(mype_filter + 1) + dim_p) THEN
      istart = 1
      imax = 0
    ELSE
      ! index for 
      istart = element - dis_dim_p(mype_filter + 1)
      imax = istart
      wk = 1.0
    END IF
  ELSE
    istart = 1
    imax = 0
    status = 100
    wk = 1.0
    WRITE(*, '(a, 5x, a, I4, a, I4, a)') 'PDAF warning:', &
      'PDAF_diag_crps: element (', element, ') argument must be <= dim_p (', dim_p, ').'
    RETURN
  END IF

  ! Calculate uncertainty based on Eq 20 in H. Hersbach (2000)
  ! Uncertainty is only meaningful when multiple verifying analysis exists
  ! because it is calculated based on the distribution of observations
  IF (element == 0) THEN
    ALLOCATE( allobs(dim), source=0. )
    ! get observations across PEs
    CALL MPI_Allgatherv(obs, dim_p, MPI_REALTYPE, allobs, all_dim_p, dis_dim_p, MPI_REALTYPE, COMM_filter, MPIerr)
    CALL quicksort_d(allobs, dim)
    pval = 0.
    DO k = 1, dim - 1
      pval = pval + wk
      uncert = uncert + (allobs(k+1) - allobs(k)) * pval*(1.0-pval)
    END DO
    DEALLOCATE(allobs)
  END IF

  ! allocate arrays for CRPS calculation
  ALLOCATE(one_case(dim_ens))
  ALLOCATE(alpha_p(0:dim_ens), source=0.)
  ALLOCATE(beta_p(0:dim_ens), source=0.)

  ! initialise values used for summation in the loop
  one_minus_oN_p = 0.
  o0_p = 0.

  ! Loop over grid points/cases
  DO k = istart, imax
    ! Get observation for current case
    x_a = obs(k)

    ! Get sorted ensemble for current case
    DO i = 1, dim_ens
      one_case(i) = oens(k, i)
    END DO
    CALL quicksort_d(one_case, dim_ens)

    ! Outlier cases
    IF (x_a < one_case(1)) THEN
      ! Case 1: obs < all ensemble members
      beta_p(0) = beta_p(0) + wk*(one_case(1) - x_a)
      o0_p = o0_p + wk
    ELSEIF (x_a > one_case(dim_ens)) THEN
      ! Case 2: obs > all ensemble members
      alpha_p(dim_ens) = alpha_p(dim_ens) + wk*(x_a - one_case(dim_ens))
      one_minus_oN_p = one_minus_oN_p + wk
    END IF

    ! Eq. 29 and Eq. 26 in H. Hersbach (2000)
    DO i = 1, dim_ens-1
      alpha_p(i) = alpha_p(i) + wk*MAX( MIN(x_a, one_case(i+1)) - one_case(i), 0.0)
      beta_p(i) = beta_p(i) + wk*MAX(one_case(i+1) - MAX(x_a, one_case(i)), 0.0)
    END DO
  END DO

  ALLOCATE(alpha(0:dim_ens), source=0.)
  ALLOCATE(beta(0:dim_ens), source=0.)
  ! todo: get full alpha and beta
  CALL MPI_Allreduce(alpha_p, alpha, dim_ens, MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)
  CALL MPI_Allreduce(beta_p, beta, dim_ens, MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)
  CALL MPI_Allreduce(one_minus_oN_p, one_minus_oN, 1, MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)
  CALL MPI_Allreduce(o0_p, o0, 1, MPI_REALTYPE, MPI_SUM, COMM_filter, MPIerr)
  DEALLOCATE(one_case, alpha_p, beta_p)

  ! Complete computation of CPRS, reliability and potential CPRS
  ! modify alpha(0) and beta(dim_ens) to accomodate outliers calculation
  ! This is equivalent to Eq. 33 in H. Hersbach (2000)
  IF (alpha(0) /= 0.0) alpha(0) = beta(0) * (1.0/o0 - 1.0)
  IF (beta(dim_ens) /= 0.0) beta(dim_ens) = alpha(dim_ens) * (1.0/one_minus_oN - 1.0)

  DO i = 0, dim_ens
    ! The difference between i+1-th ensemble member and i-th member
    gi = alpha(i) + beta(i)
    IF (gi /= 0.0) THEN
      oi = beta(i) / gi
    ELSE
      oi = 0.0
    END IF
    pval = REAL(i) / REAL(dim_ens)

    crps = crps + alpha(i)*pval*pval + beta(i)*(1.0-pval)*(1.0-pval)
    reli = reli + gi * (oi - pval)*(oi - pval)
    pot_CRPS = pot_CRPS + gi * oi * (1.0-oi)
  END DO

  ! ****************
  ! *** Clean up ***
  ! ****************
  DEALLOCATE(alpha, beta)

END SUBROUTINE PDAF_diag_crps_mpi

!--------------------------------------------------------
!> CRPS routine from PDAF until V2.2.1 without parallelization
!!
SUBROUTINE PDAF_diag_CRPS_nompi(dim, dim_ens, element, oens, obs, &
     CRPS, reli, resol, uncert, status)!

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: dim                !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens            !< Ensemble size
  INTEGER, INTENT(in) :: element            !< ID of element to be used
       !< If element=0, mean values over all elements are computed
  REAL, INTENT(in)    :: oens(dim, dim_ens) !< State ensemble
  REAL, INTENT(in)    :: obs(dim)           !< State ensemble
  REAL, INTENT(out)   :: CRPS               !< CRPS
  REAL, INTENT(out)   :: reli               !< Reliability
  REAL, INTENT(out)   :: resol              !< resolution
  REAL, INTENT(out)   :: uncert             !< uncertainty
  INTEGER, INTENT(out) :: status            !< Status flag (0=success)

! *** local variables ***
  INTEGER :: i, elem, istart                ! Counters
  INTEGER :: imax
  REAL :: dinv
  REAL :: oneobs
  REAL :: gi, oi, pval
  REAL, ALLOCATABLE :: oneens(:)
  REAL, ALLOCATABLE :: allobs(:)
  REAL, ALLOCATABLE :: c_a(:), c_b(:)


! ******************
! *** Initialize ***
! ******************

  ! Initialize status flag
  status = 0


! ********************
! *** Compute CRPS ***
! ********************

  ! Set number of element over which CPRS is computed
  IF (element==0) THEN
     imax = dim
     dinv = 1.0/REAL(dim)
     istart = 1
  ELSEIF (element<=dim) THEN
     imax = element
     istart = element
     dinv = 1.0
  ELSE
     imax = 0
     status = 100
     dinv = 1.0
     istart = 1
  END IF

  ALLOCATE(oneens(dim_ens))
  ALLOCATE(c_a(0:dim_ens))
  ALLOCATE(c_b(0:dim_ens))

  c_a = 0.0
  c_b = 0.0

  ! Loop over elements
  DO elem = istart, imax

     ! Get observation for current element
     oneobs = obs(elem)

     ! Get sorted ensemble for current element
     DO i = 1, dim_ens
        oneens(i) = oens(elem, i)
     END DO
     CALL PDAF_sisort(dim_ens, oneens)

     IF (oneobs < oneens(1)) THEN
        ! Case 1: obs < all ensemble members
        c_a(0) = c_a(0) + dinv*(oneens(1) - oneobs)
        c_b(0) = c_b(0) + dinv
     ELSEIF (oneobs > oneens(dim_ens)) THEN
        ! Case 2: obs > all ensemble members
        c_a(dim_ens) = c_a(dim_ens) + dinv
        c_b(dim_ens) = c_b(dim_ens) + dinv*(oneobs - oneens(dim_ens))
     END IF
     DO i = 1, dim_ens-1
        c_a(i) = c_a(i) + dinv*MAX(oneens(i+1) - MAX(oneobs, oneens(i)), 0.0)
        c_b(i) = c_b(i) + dinv*MAX( MIN(oneobs, oneens(i+1)) - oneens(i), 0.0)
     END DO
  END DO

  haveobs: IF (imax>0) THEN
     
     ! Reinterpretation of c_a(dim_ens) and c_b(0)
     IF (c_b(0) /= 0.0) c_b(0) = c_a(0) * (1.0/c_b(0) - 1.0)
     IF (c_a(dim_ens) /= 0.0) c_a(dim_ens) = c_b(dim_ens) * (1.0/c_a(dim_ens) - 1.0)

     ! Calculate untertainty
     ALLOCATE(allobs(dim))
     allobs = obs

     CALL PDAF_sisort(dim, allobs)
     uncert = 0
     pval = 0.0

     DO i = 1, dim-1
        pval = pval + dinv
        uncert = uncert + (allobs(i+1) - allobs(i)) * pval*(1.0-pval)
     END DO

     ! Complete computation of CPRS, reliability and resolution
     crps = 0.0
     reli = 0.0
     resol = 0.0

     oi = 0.0
     DO i = 0, dim_ens
        gi = c_a(i) + c_b(i)
        IF (gi /= 0.0) THEN
           oi = c_a(i) / gi
        ELSE
           oi = 0.0
        END IF
        pval = REAL(i) / REAL(dim_ens)

        crps = crps + c_b(i)*pval*pval + c_a(i)*(1.0-pval)*(1.0-pval)
        reli = reli + gi * (oi-pval)*(oi-pval)
        resol = resol + gi * oi * (1.0-oi)

     END DO

  ELSE
     crps = 0.0
     reli = 0.0
     resol = 0.0
     uncert = 0.0
  END IF haveobs


! ****************
! *** Clean up ***
! ****************

  DEALLOCATE(oneens, c_a, c_b)

END SUBROUTINE PDAF_diag_CRPS_nompi


!--------------------------------------------------------
!> Sorting routine 
!!
!! Sorts a vector veca(1:n) into ascending numerical order, by
!! straight insertion.
!!
!! For large vectors, this routine will not be efficient.
!!
SUBROUTINE PDAF_sisort(n, veca)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: n 
  REAL, INTENT(inout) :: veca(n)

  INTEGER :: i, j, k 
  REAL :: tmpa
  LOGICAL :: eflag

  DO j = 2, n

    eflag = .FALSE.

    tmpa = veca(j) 

    sortloop: DO i = j-1, 1, -1 
       k = i

       IF(veca(i) <= tmpa) THEN
          eflag = .TRUE.
          EXIT sortloop
       END IF

       veca(i+1) = veca(i) 
    ENDDO sortloop

    IF (.NOT.eflag) k=0

    veca(k+1) = tmpa 

  ENDDO 

END SUBROUTINE PDAF_sisort
