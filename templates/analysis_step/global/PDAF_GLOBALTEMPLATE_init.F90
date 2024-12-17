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
!
!> PDAF-internal initialization of parameters for DA method GLOBALTEMPLATE
!!
!! Initialization of the DA method within PDAF. Performed are:\\
!!   - initialize filter-specific parameters\\
!!   - print screen information on filter configuration.
!! Integer parameters are provided in 'param\_int' and
!! floating point parameters (reals) in 'param\_real'. 
!! The value 'dim\_pint' specifies the index up to which integer
!! parameters are to be considered and the value 'dim_preal' specifies
!! the index up to which real valued parameters are to be considered.
!!
!! In this routine only the parameters are considered that are 
!! specific to a DA method. The generic parameters are taken into
!! account in the general routine PDAF\_init. Specifically these
!! parameters are:\\
!! - param_int(1) - state dimension (dim_p)\\
!! - param_int(2) - ensemble size (dim_ens)\\
!! - param_real(1) - forgetting factor (forget) 
!!
!! ADAPTING THE TEMPLATE:
!! For parameters that are specific for a DA method, one can freely
!! specify the integer parameters, param_int(i) with i>2, and real
!! parameters, param_real(i) with i>1. In this routine one hands over
!! the parameters from the input arrays param_int and param_real to
!! the method-specific variables. One has to check the value of 
!! dim_pint and dim_preal for the maximum index that should be 
!! considered.
!! When implemeting for a new DA method, you should replace 'GLOBALTEMPLATE'
!! by the name of the method.
!!
!! __Revision history:__
!! * 2024-12 - Lars Nerger - Initial code for template based on ETKF
!! * Later revisions - see repository log
!!
SUBROUTINE PDAF_GLOBALTEMPLATE_init(subtype, param_int, dim_pint, param_real, dim_preal, &
     ensemblefilter, fixedbasis, verbose, outflag)

  USE PDAF_mod_filter, &
       ONLY: incremental, dim_ens, dim_lag, localfilter

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(inout) :: subtype             ! Sub-type of filter
  INTEGER, INTENT(in) :: dim_pint               ! Number of integer parameters
  INTEGER, INTENT(inout) :: param_int(dim_pint) ! Integer parameter array
  INTEGER, INTENT(in) :: dim_preal              ! Number of real parameters 
  REAL, INTENT(inout) :: param_real(dim_preal)  ! Real parameter array
  LOGICAL, INTENT(out) :: ensemblefilter        ! Is the chosen filter ensemble-based?
  LOGICAL, INTENT(out) :: fixedbasis            ! Does the filter run with fixed error-space basis?
  INTEGER, INTENT(in) :: verbose                ! Control screen output
  INTEGER, INTENT(inout):: outflag              ! Status flag


! **********************************************************
! *** Initialize variables from param_int and param_real ***
! **********************************************************

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++ TEMPLATE:                                          +++
! +++ Adapt the initialization of the internal parameter +++
! +++ variables from param_int and param_real            +++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Size of lag considered for smoother
  IF (dim_pint>=3) THEN
     IF (param_int(3) > 0) THEN
        dim_lag = param_int(3)
     ELSE
        dim_lag = 0
     END IF
  END IF

  ! Example initializing a real parameter
  IF (dim_preal>=2) THEN
     IF (param_int(2) > 0) THEN
!        myvar = param_real(2)
     ELSE
!        myvar = 0
     END IF
  END IF


! *************************************************************
! *** Initialize further variables specifying the DA method ***
! *************************************************************

  ! Define whether DA method is mode-based or ensemble-based
  ensemblefilter = .TRUE.

  ! Define whether filter is domain localized
  localfilter = 0

  ! Initialize flag for fixed-basis ensemble methods
  IF (subtype == 2 .OR. subtype == 3) THEN
     fixedbasis = .TRUE.
  ELSE
     fixedbasis = .FALSE.
  END IF


! *********************
! *** Screen output ***
! *********************

! +++++++++++++++++++++++++++++++++++++++++++++++++++
! +++ TEMPLATE:                                   +++
! +++ Adapt the screen output to give an overview +++
! +++ of the specific parameters of the DA method +++
! +++++++++++++++++++++++++++++++++++++++++++++++++++

  filter_pe2: IF (verbose > 0) THEN
  
     WRITE(*, '(/a, 4x, a)') 'PDAF', '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++      TEMPLATE FOR GLOBALTEMPLATE DA METHOD       +++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++ replace this string upon implementing the method +++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     ! *** General output ***
     WRITE (*, '(/a, 4x, a)') 'PDAF', 'GLOBALTEMPLATE configuration'
     WRITE (*, '(a, 9x, a, i1)') 'PDAF', 'sub-type = ', subtype
     IF (subtype == 0) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> GLOBALTEMPLATE subtype 0'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> GLOBALTEMPLATE subtype 1'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> offline mode'

        ! Reset subtype
        subtype = 0
     ELSE
        WRITE (*, '(/5x, a/)') 'PDAF-ERROR(2): No valid sub type!'
        outflag = 2
     END IF
     IF (incremental == 1) &
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> Perform incremental updating'
     IF (dim_lag > 0) &
          WRITE (*, '(a, 12x, a, i6)') 'PDAF', '--> Apply smoother up to lag:',dim_lag
     WRITE (*, '(a, 12x, a, i5)') 'PDAF', '--> ensemble size:', dim_ens

  END IF filter_pe2

END SUBROUTINE PDAF_GLOBALTEMPLATE_init
