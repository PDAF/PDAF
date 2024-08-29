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
!$Id: PDAF_set_comm_pdaf.F90 918 2021-12-03 07:42:19Z lnerger $


!> Get value of a correlation function
!!
!! This routine returns the value of the chosen correlation
!! function according to the specified length scale.
!!
!!  This is a core routine of PDAF and
!!  should not be changed by the user   !
!!
!! __Revision history:__
!! * 2024-08 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE PDAF_correlation_function(ctype, length, distance, value)

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: ctype          !< Type of correlation function
                                        !< (1): Gaussian with f(0)=1.0
                                        !< (2): 5th-order polynomial (Gaspace/Cohn, 1999)
  REAL, INTENT(in) :: length            !< Length scale of function
                                        !< (1): standard deviation
                                        !< (2): support length f=0 for distance>length 
  REAL, INTENT(in) :: distance          !< Distance at which the function is evaluated
  REAL, INTENT(out) :: value            !< Value of the function
  

! *** Local variables ***
  REAL :: shalf                         ! Half of support distance for Gaspari-Cohn


! **************************
! *** Set function value ***
! **************************

  IF (ctype == 1) THEN
     ! *********************************************
     ! *** Gaussian function scaled for f(0)=1.0 ***
     ! *********************************************

     ! Compute weight
     IF (length > 0.0) THEN

        value = exp(-distance*distance/ (2.0*length*length))

     ELSE

        IF (distance > 0.0) THEN
           value = 0.0
        ELSE
           value = 1.0
        END IF

     END IF

  ELSEIF (ctype == 2) THEN
     ! ************************************************************************
     ! *** 5th-order polynomial mimicking Gaussian but with compact support ***
     ! *** Equation (4.10) of Gaspari&Cohn, QJRMS125, 723 (1999)            ***
     ! ************************************************************************

     shalf = REAL(length) / 2.0

     ! Evaluate function
     cradnull: IF (length > 0.0) THEN

        cutoff: IF (distance <= length) THEN
           IF (distance <= length / 2.0) THEN
              value = -0.25 * (distance / shalf)**5 &
                   + 0.5 * (distance / shalf)**4 &
                   + 5.0 / 8.0 * (distance / shalf)**3 &
                   - 5.0 / 3.0 * (distance / shalf)**2 + 1.0
           ELSEIF (distance > length / 2.0 .AND. distance < length * 0.9) THEN
              value = 1.0 / 12.0 * (distance / shalf)**5 &
                   - 0.5 * (distance / shalf)**4 &
                   + 5.0 / 8.0 * (distance / shalf)**3 &
                   + 5.0 / 3.0 * (distance / shalf)**2 &
                   - 5.0 * (distance / shalf) &
                   + 4.0 - 2.0 / 3.0 * shalf / distance
           ELSEIF (distance >= length * 0.9 .AND. distance < length) THEN
              value = MAX(1.0 / 12.0 * (distance / shalf)**5 &
                   - 0.5 * (distance / shalf)**4 &
                   + 5.0 / 8.0 * (distance / shalf)**3 &
                   + 5.0 / 3.0 * (distance / shalf)**2 &
                   - 5.0 * (distance / shalf) &
                   + 4.0 - 2.0 / 3.0 * shalf / distance, 0.0)
           ELSE
              value = 0.0
           ENDIF
        ELSE cutoff
           value = 0.0
        END IF cutoff

     ELSE cradnull

        IF (distance > 0.0) THEN
           value = 0.0
        ELSE
           value = 1.0
        END IF

     END IF cradnull

  END IF


END SUBROUTINE PDAF_correlation_function
