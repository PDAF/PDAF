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
!$Id$
!BOP
!
! !ROUTINE: PDAF_local_weight --- Compute weight for localization
!
! !INTERFACE:
SUBROUTINE PDAF_local_weight(wtype, rtype, cradius, sradius, distance, &
     nrows, ncols, A, var_obs, weight, verbose)

! !DESCRIPTION:
! This routine initializates a single weight based on the given
! distance and localization radii for the specified weighting
! type.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2010-06 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: wtype      ! Type of weight function
          ! (0): unit weight (=1 up to distance=cradius)
          ! (1): exponential decrease (1/e at distance=sradius; 0 for distance>cradius)
          ! (2): 5th order polynomial (Gaspari&Cohn 1999; 0 for distance>cradius)
  INTEGER, INTENT(in) :: rtype      ! Type of regulated weighting
          ! (0): no regulation
          ! (1): regulation using mean variance
          ! (2): regulation using variance of single observation point
  REAL, INTENT(in)    :: cradius    ! Cut-off radius
  REAL, INTENT(in)    :: sradius    ! Support radius 
  REAL, INTENT(in)    :: distance   ! Distance to observation
  INTEGER, INTENT(in) :: nrows      ! Number of rows in matrix A
  INTEGER, INTENT(in) :: ncols      ! Number of columns in matrix A
  REAL, INTENT(in) :: A(nrows, ncols) ! Input matrix
  REAL, INTENT(in)    :: var_obs    ! Observation variance
  REAL, INTENT(out)   :: weight     ! Weights
  INTEGER, INTENT(in) :: verbose    ! Verbosity flag
  

! *** Local variables ***
  INTEGER :: i, j                   ! Counters
  REAL    :: cfaci                  ! parameter for initialization of 5th-order polynomial
  REAL    :: meanvar                ! Mean variance in observation domain
  REAL    :: svarpovar              ! Mean state plus observation variance
  REAL    :: var                    ! variance for Gaussian
  REAL, PARAMETER :: pi=3.141592653589793   !Pi


! ********************************
! *** Print screen information ***
! ********************************

  ptype: IF (verbose == 1) THEN
     WRITE (*,'(a, 5x, a)') 'PDAF','Set localization weights'
     IF (wtype == 0) THEN
        WRITE (*, '(a, 5x, a)') &
             'PDAF', '--- Initialize unit weights'
        WRITE (*, '(a, 5x, a, es12.4)') &
             'PDAF', '--- Support radius ', sradius
        IF (cradius < sradius) THEN
           WRITE (*, '(a, 5x, a, es12.4)') &
                'PDAF', '--- Use cut-off radius ', cradius
        END IF
     ELSE IF (wtype == 1) THEN
       WRITE (*, '(a, 5x, a)') &
             'PDAF', '--- Initialize exponential weight function'
        WRITE (*, '(a, 5x, a, es12.4)') &
             'PDAF', '--- Distance for 1/e   ', sradius
        WRITE (*, '(a, 5x, a, es12.4)') &
             'PDAF', '--- Cut-off radius ', cradius
     ELSE IF (wtype == 2) THEN
        WRITE (*, '(a, 5x, a)') &
             'PDAF', '--- Initialize weights by 5th-order polynomial'
        WRITE (*, '(a, 5x, a, es12.4)') &
             'PDAF', '--- Support radius ', sradius
        IF (cradius < sradius) THEN
          WRITE (*, '(a, 5x, a, es12.4)') &
                'PDAF', '--- Use cut-off radius ', cradius
        END IF
     ELSE IF (wtype == 3) THEN
        WRITE (*, '(a, 5x, a)') &
             'PDAF', '--- Initialize weights by scaled Gaussian function'
        WRITE (*, '(a, 5x, a, es12.4)') &
             'PDAF', '--- Standard deviation ', sradius
        WRITE (*, '(a, 5x, a, es12.4)') &
             'PDAF', '--- Cut-off radius ', cradius
     END IF
  END IF ptype


! **************************
! *** Initialize weights ***
! **************************

  t_weight: IF (wtype == 0) THEN
     ! Unit weights

     IF (distance <= cradius) THEN
        weight = 1.0
     ELSE
        weight = 0.0
     END IF

  ELSE IF (wtype == 1) THEN t_weight
     ! Weighting by exponential decrease

     IF (cradius > 0.0 .AND. sradius > 0.0) THEN

        IF (distance <= cradius) THEN
           weight = EXP(-distance / sradius)
        ELSE
           weight = 0.0
        END IF

     ELSE

        IF (distance > 0.0) THEN
           weight = 0.0
        ELSE
           weight = 1.0
        END IF

     END IF

  ELSE IF (wtype == 2) THEN t_weight
     ! Weighting by the square-root of a 5th-order function;
     ! equation (4.10) of Gaspari&Cohn, QJRMS125, 723 (1999)

     cfaci = REAL(sradius) / 2.0

     ! Compute weight
     cradnull: IF (cradius > 0.0 .and. sradius > 0.0) THEN

        cutoff: IF (distance <= cradius) THEN
           IF (distance <= sradius / 2.0) THEN
              weight = -0.25 * (distance / cfaci)**5 &
                   + 0.5 * (distance / cfaci)**4 &
                   + 5.0 / 8.0 * (distance / cfaci)**3 &
                   - 5.0 / 3.0 * (distance / cfaci)**2 + 1.0
           ELSEIF (distance > sradius / 2.0 .AND. distance < sradius * 0.9) THEN
              weight = 1.0 / 12.0 * (distance / cfaci)**5 &
                   - 0.5 * (distance / cfaci)**4 &
                   + 5.0 / 8.0 * (distance / cfaci)**3 &
                   + 5.0 / 3.0 * (distance / cfaci)**2 &
                   - 5.0 * (distance / cfaci) &
                   + 4.0 - 2.0 / 3.0 * cfaci / distance
           ELSEIF (distance >= sradius * 0.9 .AND. distance < sradius) THEN
              ! Ensure that weight is non-negative
              weight = MAX(1.0 / 12.0 * (distance / cfaci)**5 &
                   - 0.5 * (distance / cfaci)**4 &
                   + 5.0 / 8.0 * (distance / cfaci)**3 &
                   + 5.0 / 3.0 * (distance / cfaci)**2 &
                   - 5.0 * (distance / cfaci) &
                   + 4.0 - 2.0 / 3.0 * cfaci / distance, 0.0)
           ELSE
              weight = 0.0
           ENDIF
        ELSE cutoff
           weight = 0.0
        END IF cutoff

     ELSE cradnull

        IF (distance > 0.0) THEN
           weight = 0.0
        ELSE
           weight = 1.0
        END IF

     END IF cradnull

  ELSE IF (wtype == 3) THEN t_weight
     ! Weighting by Gaussian function scaled for w(0)=1.0

     ! Compute weight
     IF (cradius > 0.0 .and. sradius > 0.0) THEN

        IF (distance <= cradius) THEN
           var = sradius*sradius
           weight = exp(-distance*distance/ (2.0*var))
        ELSE
           weight = 0.0
        END IF

     ELSE

        IF (distance > 0.0) THEN
           weight = 0.0
        ELSE
           weight = 1.0
        END IF

     END IF



  END IF t_weight


! *********************************
! *** Compute weight regulation ***
! *********************************

  regweight: IF (rtype == 1) THEN
     ! Regulated weight with a function based on the fixed weight
     ! function, the variance of the observations and the local mean of
     ! the estimated state error variance for the local analysis domain.

     IF (verbose == 1) THEN
        WRITE (*, '(a, 5x, a)') &
             'PDAF', '--- Compute regulated weight'
     END IF

     ! Compute mean variance from ensemble perturbations
     meanvar = 0.0
     DO i = 1, nrows
        DO j = 1, ncols
           meanvar = meanvar + A(i,j)**2
        END DO
     END DO
     meanvar = meanvar / REAL(ncols) / REAL(nrows)

     svarpovar = meanvar + var_obs

     ! Compute regulated weight
     weight = weight * var_obs / svarpovar &
          / (1.0 - (weight * meanvar / svarpovar))

  END IF regweight


END SUBROUTINE PDAF_local_weight
