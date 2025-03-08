! Copyright (c) 2004-2025 Lars Nerger
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
!> Module for GENOBS holding shared parameters and some helper routines
!!
!! This module declares the parameters that are used in GENOBS. 
!! Parameters that are specific for GENOBS are declared while some
!! other parameters are use-included from PDAF_mod_core. This allows
!! us to only include this module in the method-specific analysis routines.
!! In addition, subroutines are included that initialize these parameters.
!!
!!    ! This is a core routine of PDAF and !
!!    ! should not be changed by the user  !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code from restructuring
!! *  Other revisions - see repository log
!!
MODULE PDAF_GENOBS

  USE PDAF_mod_core, &
       ONLY: debug, dim_lag

  IMPLICIT NONE

! *** Integer parameters ***
  INTEGER :: seedset=1     ! Choice of seed set for random numbers
                           ! Valid are choices between 1 and 20


!-------------------------------------------------------------------------------
  
CONTAINS

!> PDAF-internal initialization of observation generation
!!
!! Initialization of NETF within PDAF. Performed are:\\
!!   - initialize filter-specific parameters\\
!!   - print screen information on filter configuration.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2014-05 - Paul Kirchgessner - Initial code based on code for ETKF
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_genobs_init(subtype, param_int, dim_pint, param_real, dim_preal, &
       ensemblefilter, fixedbasis, verbose, outflag)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: subtype                !< Sub-type of filter
    INTEGER, INTENT(in) :: dim_pint               !< Number of integer parameters
    INTEGER, INTENT(inout) :: param_int(dim_pint) !< Integer parameter array
    INTEGER, INTENT(in) :: dim_preal              !< Number of real parameters 
    REAL, INTENT(inout) :: param_real(dim_preal)  !< Real parameter array
    LOGICAL, INTENT(out) :: ensemblefilter        !< Is the chosen filter ensemble-based?
    LOGICAL, INTENT(out) :: fixedbasis            !< Does the filter run with fixed error-space basis?
    INTEGER, INTENT(in) :: verbose                !< Control screen output
    INTEGER, INTENT(inout):: outflag              !< Status flag

! *** Local variables ***
    INTEGER :: i                         ! Counter
    INTEGER :: subtype_dummy             ! Dummy variable to prevent compiler warning
    REAL :: param_real_dummy(dim_preal)  ! Dummy variable to prevent compiler warning


! *********************
! *** Screen output ***
! *********************

    writeout: IF (verbose == 1) THEN
       WRITE(*, '(/a)') 'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       WRITE(*, '(a)')  'PDAF    +++       PDAF Generator for synthetic observations       +++'
       WRITE(*, '(a)')  'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    END IF writeout


! ****************************
! *** INITIALIZE VARIABLES ***
! ****************************

    ! Initialize dummy variables to prevent compiler warnings
    subtype_dummy = subtype
    param_real_dummy = param_real

    DO i=3, dim_pint
       CALL PDAF_genobs_set_iparam(i, param_int(i), outflag)
    END DO

    ! Define whether filter is mode-based or ensemble-based
    ensemblefilter = .TRUE.

    ! Initialize flag for fixed-basis filters
    fixedbasis = .FALSE.

    ! Set status flag
    outflag = 0

END SUBROUTINE PDAF_genobs_init

!-------------------------------------------------------------------------------
!> Perform allocation of arrays for GENOBS.
!!
!! __Revision history:__
!! * 2019-01 - Lars Nerger - Initial code
!! * 2025-02 - Lars Nerger - Restructuring introducing generic PDAF_alloc
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_genobs_alloc(outflag)

    USE PDAF_mod_core, &
         ONLY: dim_ens, dim_p
    USE PDAF_mod_parallel, &
         ONLY: dim_ens_l
    USE PDAF_utils, &
         ONLY: PDAF_alloc

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout):: outflag      !< Status flag


! ******************************
! *** Allocate filter fields ***
! ******************************

    CALL PDAF_alloc(dim_p, dim_ens, dim_ens_l, 1, 0, &
         0, 0, outflag)

  END SUBROUTINE PDAF_genobs_alloc


!-------------------------------------------------------------------------------
!>  Print information on configuration of GENOBS
!!
!!  !  This is a core routine of PDAF and   !
!!  !   should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code by splitting from PDAF_genobs_init
!! *  Other revisions - see repository log
!!
  SUBROUTINE PDAF_genobs_config(subtype, verbose)

    USE PDAF_mod_core, &
         ONLY: dim_ens

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: subtype               !< Sub-type of filter
    INTEGER, INTENT(in)    :: verbose               !< Control screen output


! *********************
! *** Screen output ***
! *********************

    writeout: IF (verbose > 0) THEN

       WRITE (*, '(/a, 4x, a)') 'PDAF', 'GENOBS configuration'
       WRITE (*, '(a, 10x, a, i5)') 'PDAF', 'ensemble size:', dim_ens
       WRITE (*, '(a, 10x, a, i1)') 'PDAF', 'filter sub-type= ', subtype
       IF (subtype == 0) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> standard observation generation'
       END IF
       WRITE(*, '(a, 10x, a, i3)') &
            'PDAF', 'param_int(3) seedset=', seedset

    END IF writeout

  END SUBROUTINE PDAF_genobs_config


!-------------------------------------------------------------------------------
!> Set integer parameter specific for GENOBS
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Other revisions - see repository log
!!
  SUBROUTINE PDAF_genobs_set_iparam(id, value, flag)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)  :: id      !< Index of parameter
    INTEGER, INTENT(in)  :: value   !< Parameter value
    INTEGER, INTENT(out) :: flag    !< Status flag: 0 for no error


! ****************************
! *** INITIALIZE VARIABLES ***
! ****************************

    ! Initialize status flag
    flag = 0

    SELECT CASE(id) 
    CASE(1)
       CALL PDAF_reset_dim_p(value, flag)
    CASE(2)
       CALL PDAF_reset_dim_ens(value, flag)
    CASE(3)
       seedset = value
       IF (seedset<1 .OR. seedset>20) THEN
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(8): Invalid setting for random number seed set - param_int(3)!'
          flag = 8
       END IF
    CASE DEFAULT
       WRITE (*,'(/5x, a, i3/)') &
            'PDAF-WARNING: Invalid integer parameter index', id
    END SELECT

  END SUBROUTINE PDAF_genobs_set_iparam

!-------------------------------------------------------------------------------
!> Information output on options for GENOBS
!!
!! Subroutine to perform information output on options
!! available for the NETF.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! !__Revision history:__
!! * 2019-01 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_genobs_options()

    IMPLICIT NONE

! *********************
! *** Screen output ***
! *********************

    WRITE(*, '(/a)') 'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    WRITE(*, '(a)')  'PDAF    +++       PDAF Generator for synthetic observations       +++'
    WRITE(*, '(a)')  'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

    WRITE(*, '(/a, 5x, a)') 'PDAF', 'Available options for GENOBS:'

    WRITE(*, '(a, 5x, a)') 'PDAF', '--- Sub-types (Parameter subtype) ---'
    WRITE(*, '(a, 7x, a)') 'PDAF', '0: Standard implementation with ensemble integration'

    WRITE(*, '(a, 5x, a)') 'PDAF', '--- Integer parameters (Array param_int) ---'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(1): Dimension of state vector (>0), required'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(2): Ensemble size (>0), required'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(3): seedset'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'seed set index for random number generator, optional'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'valid are values between 1 and 20; default=1'

    WRITE(*, '(a, 5x, a)') 'PDAF', '--- Floating point parameters (Array param_real) ---'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_real(1): Forgetting factor (usually >0 and <=1), required, but not used'

    WRITE(*, '(a, 5x, a)') 'PDAF', '--- Further parameters ---'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'n_modeltasks: Number of parallel model integration tasks'
    WRITE(*, '(a, 11x, a)') &
         'PDAF', '=1 for GENOBS; not larger than total number of processors'
    WRITE(*, '(a, 11x, a)') 'PDAF', '=1 required for subtypes 2 and 3'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'screen: Control verbosity of PDAF'
    WRITE(*, '(a, 11x, a)') 'PDAF', '0: no outputs'
    WRITE(*, '(a, 11x, a)') 'PDAF', '1: basic output (default)'
    WRITE(*, '(a, 11x, a)') 'PDAF', '2: 1 plus timing output'
    WRITE(*, '(a, 11x, a)') 'PDAF', '3: 2 plus debug output'

    WRITE(*, '(a, 5x, a)') &
         'PDAF', '+++++++++ End of option overview for GENOBS  ++++++++++'

  END SUBROUTINE PDAF_genobs_options

END MODULE PDAF_GENOBS
