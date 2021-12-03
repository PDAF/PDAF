! Copyright (c) 2014-2021 Lars Nerger
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
!BOP
!
! !ROUTINE: PDAF_PF_init --- PDAF-internal initialization of PF
!
! !INTERFACE:
SUBROUTINE PDAF_PF_init(subtype, param_int, dim_pint, param_real, dim_preal, &
     ensemblefilter, fixedbasis, verbose, outflag)

! !DESCRIPTION:
! Initialization of PF within PDAF. Performed are:\\
!   - initialize filter-specific parameters\\
!   - print screen information on filter configuration.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2019-05 - Lars Nerger - Initial code based on code for NETF
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filter, &
       ONLY: restype, noise_type, pf_noise_amp, type_forget, forget, &
       type_winf, limit_winf

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: subtype                ! Sub-type of filter
  INTEGER, INTENT(in) :: dim_pint               ! Number of integer parameters
  INTEGER, INTENT(inout) :: param_int(dim_pint) ! Integer parameter array
  INTEGER, INTENT(in) :: dim_preal              ! Number of real parameters 
  REAL, INTENT(inout) :: param_real(dim_preal)  ! Real parameter array
  LOGICAL, INTENT(out) :: ensemblefilter        ! Is the chosen filter ensemble-based?
  LOGICAL, INTENT(out) :: fixedbasis            ! Does the filter run with fixed error-space basis?
  INTEGER, INTENT(in) :: verbose                ! Control screen output
  INTEGER, INTENT(inout):: outflag              ! Status flag

! !CALLING SEQUENCE:
! Called by: PDAF_init_filters
!EOP

! *** local variables ***
  REAL :: param_real_dummy               ! Dummy variable to avoid compiler warning


! ****************************
! *** INITIALIZE VARIABLES ***
! ****************************

  ! Initialize variable to prevent compiler warning
  param_real_dummy = param_real(1)

  ! Set type of resampling
  IF (dim_pint>=3) THEN
     IF (param_int(3) > 0 .AND. param_int(3) < 4) THEN
        restype = param_int(3)
     END IF
  END IF

  ! Set type of resampling
  IF (dim_pint>=4) THEN
     IF (param_int(4) > 0 .AND. param_int(4) < 3) THEN
        noise_type = param_int(4)
     END IF
  END IF

  ! Store type for forgetting factor
  IF (dim_pint >= 5) THEN
     type_forget = param_int(5)
  END IF

  ! Type of weights inflation
  IF (dim_pint >= 6) THEN     
     type_winf = param_int(6)
  END IF

  ! Store value of forgetting factor variable which is noise amplitude here
  pf_noise_amp = forget

  ! forgetting factor
  IF (dim_preal >= 2) THEN
     forget = param_real(2)
  ELSE
     forget = 1.0
  END IF

  ! Strength of weights inflation
  IF (dim_preal >= 3) THEN
     IF (param_real(3) < 0.0) THEN
        WRITE (*,'(/5x,a/)') &
             'PDAF-ERROR(10): Invalid limit for weight inflation!'
        outflag = 10
     END IF
     limit_winf = param_real(3)
  END IF


  ! Define whether filter is mode-based or ensemble-based
  ensemblefilter = .TRUE.

  ! Initialize flag for fixed-basis filters
  fixedbasis = .FALSE.


! *********************
! *** Screen output ***
! *********************

  filter_pe: IF (verbose == 1) THEN
  
     WRITE(*, '(/a)') 'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     WRITE(*, '(a)')  'PDAF    +++           Particle Filter with resampling             +++'
     WRITE(*, '(a)')  'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     ! *** General output ***
     WRITE (*, '(/a, 4x, a)') 'PDAF', 'PF configuration'
     WRITE (*, '(a, 11x, a, i1)') 'PDAF', 'filter sub-type = ', subtype
     IF (subtype == 0) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> PF '
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> offline mode'
     ELSE
        WRITE (*, '(/5x, a/)') 'PDAF-ERROR(2): No valid sub type!'
        outflag = 3
     END IF
     IF (restype == 1) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> Resample using probabilistic resampling'
     ELSE IF (restype == 2) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> Resample using stochastic universal resampling'
     ELSE IF (restype == 3) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> Resample using residual resampling'
     END IF

  END IF filter_pe

END SUBROUTINE PDAF_PF_init
