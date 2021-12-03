! Copyright (c) 2012-2021 Lars Nerger, lars.nerger@awi.de
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
!! resolution: CRPS = RELI + RESOL. In addition the uncertainty
!! is computed.
!! A perfectly reliable system gives RELI=0.
!! An informative system gives RESOL << UNCERT.
!!
!! The computation follows H. Hersbach, Weather and Forecasting
!! 15(2000) 599-570. Here, RESOL is equivalent to CPRS_pot.
!!
!! __Revision history:__
!! * 2021-05 - Lars Nerger - Initial code based on sangoma_ComputeCRPS
!! * Later revisions - see repository log
!!
SUBROUTINE PDAF_diag_CRPS(dim, dim_ens, element, oens, obs, &
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
     DO i = 1, dim_ens
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

END SUBROUTINE PDAF_diag_CRPS


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
