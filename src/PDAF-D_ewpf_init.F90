! Copyright (c) 2014-2018 Paul Kirchgessner
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
! !ROUTINE: PDAF_ewpf_init --- PDAF-internal initialization of EWPF
!
! !INTERFACE:
SUBROUTINE PDAF_ewpf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
     ensemblefilter, fixedbasis, verbose, outflag)

! !DESCRIPTION:
! Initialization of Equivalent Weights Particle Filter within PDAF.
! Performed are:
!   - initialize filter-specific parameters
!   - print screen information on filter configuration.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2014-05 - Paul Kirchgessner - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_ewpf, &
       ONLY: incremental, dim_ens, dim_p, &
       dim_bias_p, type_nudging, &
       bt, start_nudging, use_model_error, keep
  USE PDAF_mod_filtermpi, &
       ONLY: mype, filterpe

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: subtype                ! Sub-type of filter
  INTEGER, INTENT(in) :: dim_pint               ! Number of integer parameters
  INTEGER, INTENT(inout) :: param_int(dim_pint) ! Integer parameter array
  INTEGER, INTENT(in) :: dim_preal              ! Number of real parameters 
  REAL, INTENT(inout) :: param_real(dim_preal)  ! Real parameter array
  LOGICAL, INTENT(out) :: ensemblefilter ! Is the chosen filter ensemble-based?
  INTEGER, INTENT(in) :: verbose                ! Control screen output
  INTEGER, INTENT(inout):: outflag              ! Status flag
  LOGICAL, INTENT(out) :: fixedbasis     ! Does the filter run with fixed error-space

! !CALLING SEQUENCE:
! Called by: PDAF_init_filters
!EOP


! ****************************
! *** INITIALIZE VARIABLES ***
! ****************************

  ! Set type of proposal step nudging
  IF (dim_pint>=3) THEN
     IF (param_int(3) > 0) THEN
        type_nudging = param_int(3)
     ELSE
        type_nudging = 0 !no nudging
     END IF
  END IF

  ! Nudging strength
  IF (dim_preal>=2) THEN
     bt = param_real(2)
  END IF

  ! Time interval at which nuding is started
  IF (dim_preal>=3) THEN
     start_nudging = param_real(3)
  END IF

  ! Fraction of particles to keep in equivalent-weights step
  IF (dim_preal>=4) THEN
     keep = param_real(4)
  END IF
   
  ! Define whether filter is mode-based or ensemble-based
  ensemblefilter = .TRUE.

  ! Initialize flag for fixed-basis filters
  fixedbasis = .FALSE.
 
  ! Set model error flag as true.
  use_model_error = .TRUE. 


! *********************
! *** Screen output ***
! *********************

  filter_pe2: IF (verbose == 1) THEN
  
     WRITE(*, '(/8x, a)') '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     WRITE(*, '(8x, a)')  '+++   Equivilent weights particle filter (EWPF)     +++'
     WRITE(*, '(8x, a)')  '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     ! *** General output ***
     WRITE (*, '(/9x, a)') 'EWPF configuration'
     WRITE (*, '(15x, a, i1)') 'filter sub-type = ', subtype
     WRITE (*, '(15x, a, f5.1)') 'Start nudging  = ' , start_nudging
     WRITE (*, '(15x, a, f5.1)') 'Nudging strength  = ' , bt
     IF (subtype == 0) THEN
        WRITE (*, '(17x, a)') 'PDAF: Use standard EWPF'
     ELSE
        WRITE (*, '(/5x, a/)') 'PDAF-ERROR(2): No valid sub type!'
        outflag = 2
     END IF

  END IF filter_pe2

END SUBROUTINE PDAF_ewpf_init
