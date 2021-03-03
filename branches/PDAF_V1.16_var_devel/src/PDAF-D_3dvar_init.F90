! Copyright (c) 2004-2020 Lars Nerger
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
! !ROUTINE: PDAF_3dvar_init --- PDAF-internal initialization of 3DVAR
!
! !INTERFACE:
SUBROUTINE PDAF_3dvar_init(subtype, param_int, dim_pint, param_real, dim_preal, &
     ensemblefilter, fixedbasis, verbose, outflag)

! !DESCRIPTION:
! Initialization of 3DVAR within PDAF. Performed are:\\
!   - initialize filter-specific parameters\\
!   - print screen information on filter configuration.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2021-03 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filter, &
       ONLY: incremental, dim_ens, forget, type_forget, dim_bias_p, &
       type_trans, dim_lag, observe_ens

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: subtype                ! Sub-type of filter
  INTEGER, INTENT(in) :: dim_pint               ! Number of integer parameters
  INTEGER, INTENT(inout) :: param_int(dim_pint) ! Integer parameter array
  INTEGER, INTENT(in) :: dim_preal              ! Number of real parameters 
  REAL, INTENT(inout) :: param_real(dim_preal)  ! Real parameter array
  LOGICAL, INTENT(out) :: ensemblefilter ! Is the chosen filter ensemble-based?
  LOGICAL, INTENT(out) :: fixedbasis     ! Does the filter run with fixed error-space basis?
  INTEGER, INTENT(in) :: verbose                ! Control screen output
  INTEGER, INTENT(inout):: outflag              ! Status flag

! !CALLING SEQUENCE:
! Called by: PDAF_init_filters
!EOP

! *** local variables ***
  REAL :: param_real_dummy    ! Dummy variable to avoid compiler warning


! ****************************
! *** INITIALIZE VARIABLES ***
! ****************************

  ! Initialize variable to prevent compiler warning
  param_real_dummy = param_real(1)

  ! Size of lag considered for smoother
  IF (dim_pint>=3) THEN
     IF (param_int(3)>0) WRITE (*,'(/5x, a/)') &
          'PDAF-ERROR(10): 3DVAR does not support smoothing!'
     dim_lag = 0
  END IF

  ! Whether incremental updating is performed
  ! We do not have incremental updating for 3DVAR!
  IF (dim_pint>=4) THEN
     incremental = param_int(4)
     IF (param_int(4) /= 0 .AND. param_int(4) /= 1) THEN
        WRITE (*,'(/5x, a/)') &
             'PDAF-ERROR(10): 3DVAR does not yet support incremental updating!'
        outflag = 10
     END IF
  END IF

  ! Store dimension of bias vector
  IF (dim_pint >= 7) THEN
!      dim_bias_p = param_int(7)
     dim_bias_p = 0 ! Temporary - bias correction not yet implemented for 3DVAR
  END IF

  ! Define whether filter is mode-based or ensemble-based
  ensemblefilter = .TRUE.

  ! Initialize flag for fixed-basis filters
!   IF (subtype == 2 .OR. subtype == 3) THEN
!      fixedbasis = .TRUE.
!   ELSE
!      fixedbasis = .FALSE.
!   END IF


! *********************
! *** Screen output ***
! *********************

  filter_pe2: IF (verbose > 0) THEN
  
     WRITE(*, '(/a, 4x, a)') 'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++                      3D-Var                     +++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     ! *** General output ***
     WRITE (*, '(/a, 4x, a)') 'PDAF', '3DVAR configuration'
     WRITE (*, '(a, 9x, a, i1)') 'PDAF', 'filter sub-type = ', subtype
     IF (subtype == 0) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> 3DVAR standard mode using LBFGS'
     ELSE
        WRITE (*, '(/5x, a/)') 'PDAF-ERROR(2): No valid sub type!'
        outflag = 2
     END IF
     IF (incremental == 1) &
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> Perform incremental updating'
!      IF (type_forget == 0) THEN
!         WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> Use fixed forgetting factor:', forget
     WRITE (*, '(a, 12x, a, i5)') 'PDAF', '--> ensemble size:', dim_ens

  END IF filter_pe2

END SUBROUTINE PDAF_3dvar_init
