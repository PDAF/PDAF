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
!> Module for ESTKF holding shared parameters and some helper routines
!!
!! This module declares the parameters that are used in ESTKF. 
!! Parameters that are specific for ESTKF are declared while some
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
MODULE PDAF_ESTKF

  USE PDAF_mod_core, &
       ONLY: localfilter, debug, dim_lag

  IMPLICIT NONE

! *** Integer parameters ***
  INTEGER :: type_forget=0 !< Type of forgetting factor
                           !< (0): fixed; (1) global adaptive
  INTEGER :: type_trans=0  !< Type of ensemble transformation
                           !< For ESTKF/LESTKF:
                           !< (0) use deterministic Omega
                           !< (1) use random orthonormal Omega orthogonal to (1,...,1)^T
                           !< (2) use product of (0) with random orthonomal matrix with
                           !<     eigenvector (1,...,1)^T
  INTEGER :: type_sqrt=0   !< Type of sqrt of A in ESTKF/LESTKF
                           !< (0): symmetric sqrt; (1): Cholesky decomposition

! *** Real parameters ***
  REAL    :: forget=1.0    !< Forgetting factor


!-------------------------------------------------------------------------------
  
CONTAINS

!>  PDAF-internal initialization of ESTKF
!!
!! Initialization of ESTKF within PDAF. Performed are:
!! * initialize filter-specific parameters
!! * print screen information on filter configuration.
!!
!!  !  This is a core routine of PDAF and   !
!!  !   should not be changed by the user   !
!!
!! __Revision history:__
!! * 2011-09 - Lars Nerger - Initial code
!! *  Other revisions - see repository log
!!
  SUBROUTINE PDAF_estkf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
       ensemblefilter, fixedbasis, verbose, outflag)

    USE PDAF_mod_core, &
         ONLY: dim_lag
    USE PDAFobs, &
         ONLY: observe_ens

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: subtype               !< Sub-type of filter
    INTEGER, INTENT(in)    :: dim_pint              !< Number of integer parameters
    INTEGER, INTENT(inout) :: param_int(dim_pint)   !< Integer parameter array
    INTEGER, INTENT(in)    :: dim_preal             !< Number of real parameters 
    REAL, INTENT(inout)    :: param_real(dim_preal) !< Real parameter array
    LOGICAL, INTENT(out)   :: ensemblefilter        !< Is the chosen filter ensemble-based?
    LOGICAL, INTENT(out)   :: fixedbasis            !< Does the filter run with fixed error-space basis?
    INTEGER, INTENT(in)    :: verbose               !< Control screen output
    INTEGER, INTENT(inout) :: outflag               !< Status flag

! *** local variables ***
    INTEGER :: i                ! Counter


! *********************
! *** Screen output ***
! *********************

    writeout: IF (verbose > 0) THEN
       WRITE(*, '(/a, 4x, a)') 'PDAF' ,'++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       WRITE(*, '(a, 4x, a)')  'PDAF' ,'+++ Error Subspace Transform Kalman Filter (ESTKF) +++'
       WRITE(*, '(a, 4x, a)')  'PDAF' ,'+++                                                +++'
       WRITE(*, '(a, 4x, a)')  'PDAF' ,'+++  Nerger et al., Mon. Wea. Rev. 140 (2012) 2335 +++'
       WRITE(*, '(a, 4x, a)')  'PDAF' ,'+++           doi:10.1175/MWR-D-11-00102.1         +++'
       WRITE(*, '(a, 4x, a)')  'PDAF' ,'++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    END IF writeout


! ****************************
! *** INITIALIZE VARIABLES ***
! ****************************

    ! Set parameter default values
    ! (Other defaults are set in the module)
    observe_ens = .false.
    dim_lag = 0

    ! Parse provided parameters
    DO i=3, dim_pint
       CALL PDAF_estkf_set_iparam(i, param_int(i), outflag)
    END DO
    DO i=1, dim_preal
       CALL PDAF_estkf_set_rparam(i, param_real(i), outflag)
    END DO


    ! Define whether filter is mode-based or ensemble-based
    ensemblefilter = .TRUE.

    ! Define whether filter is a domain-local filter
    localfilter = 0

    ! Initialize flag for fixed-basis filters
    IF (subtype == 10 .OR. subtype == 11) THEN
       fixedbasis = .TRUE.
    ELSE
       fixedbasis = .FALSE.
    END IF


! *********************
! *** Check subtype ***
! *********************

    IF (.NOT.(subtype==0 .OR. subtype==10 .OR. subtype==11)) THEN
       WRITE (*, '(/5x, a/)') 'PDAF-ERROR(3): No valid subtype!'
       outflag = 3
    END IF

  END SUBROUTINE PDAF_estkf_init


!-------------------------------------------------------------------------------
!> Perform allocation of arrays for ESTKF.
!!
!! __Revision history:__
!! * 2011-09 - Lars Nerger - Initial code adapted from SEIK
!! * 2025-02 - Lars Nerger - Restructuring introducing generic PDAF_alloc
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_estkf_alloc(outflag)

    USE PDAF_mod_core, &
         ONLY: dim_ens, dim_p, dim_bias_p
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

    CALL PDAF_alloc(dim_p, dim_ens, dim_ens_l, dim_ens-1, dim_bias_p, &
         dim_lag, 0, outflag)

  END SUBROUTINE PDAF_estkf_alloc


!-------------------------------------------------------------------------------
!>  Print information on configuration of ESTKF
!!
!!  !  This is a core routine of PDAF and   !
!!  !   should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code by splitting from PDAF_estkf_init
!! *  Other revisions - see repository log
!!
  SUBROUTINE PDAF_estkf_config(subtype, verbose)

    USE PDAF_mod_core, &
         ONLY: dim_ens, dim_lag
    USE PDAFobs, &
         ONLY: observe_ens, type_obs_init

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: subtype               !< Sub-type of filter
    INTEGER, INTENT(in)    :: verbose               !< Control screen output


! *********************
! *** Screen output ***
! *********************

    writeout: IF (verbose > 0) THEN

       ! *** General output ***
       WRITE (*, '(/a, 4x, a)') 'PDAF', 'ESTKF configuration'
       WRITE (*, '(a, 10x, a, i5)') 'PDAF', 'ensemble size:', dim_ens
       WRITE (*, '(a, 10x, a, i1)') 'PDAF', 'filter sub-type= ', subtype
       IF (subtype == 0) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> Standard ESTKF'
       ELSE IF (subtype == 10) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> ESTKF with fixed error-space basis'
       ELSE IF (subtype == 11) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> ESTKF with fixed state covariance matrix'
       END IF
       WRITE(*, '(a, 10x, a, i3)') &
            'PDAF', 'param_int(5) type_forget=', type_forget
       IF (type_forget == 0) THEN
          WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> Use fixed forgetting factor:', forget
       ELSEIF (type_forget == 1) THEN
          WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> Use adaptive forgetting factor, default value=', forget
       ENDIF
       WRITE(*, '(a, 10x, a, i3)') &
            'PDAF', 'param_int(6) type_trans=', type_trans
       IF (type_trans == 0) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> Deterministic ensemble transformation (default)'
       ELSE IF (type_trans == 1) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> Transform ensemble with random orthonormal Omega'
       ELSE IF (type_trans == 2) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> Transform ensemble including product with random matrix'
       END IF
       WRITE(*, '(a, 10x, a, i3)') &
            'PDAF', 'param_int(7) type_sqrt=', type_sqrt
       IF (type_sqrt == 0) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> symmetric square root (default)'
       ELSE IF (type_sqrt == 1) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> Cholesky decomposition'
       END IF
       WRITE(*, '(a, 10x, a, l)') &
            'PDAF', 'param_int(8) observe_ens'
       IF (observe_ens) THEN
          WRITE(*, '(a, 12x, a)') 'PDAF', '--> 1: Apply H to ensemble states and compute innovation as mean (default)'
       ELSE
          WRITE(*, '(a, 12x, a)') 'PDAF', '--> 0: Apply H to ensemble mean to compute innovation'
       END IF
       WRITE(*, '(a, 10x, a, i3)') &
            'PDAF', 'param_int(9) type_obs_init=', type_obs_init
       IF (type_obs_init==0) THEN
          WRITE(*, '(a, 12x, a)') 'PDAF', '--> Initialize observations before PDAF prestep'
       ELSE IF (type_obs_init==1) THEN
          WRITE(*, '(a, 12x, a)') 'PDAF', '--> Initialize observations after PDAF prestep'
       END IF
       IF (dim_lag > 0) &
            WRITE (*, '(a, 12x, a, i6)') 'PDAF', '--> Apply smoother up to lag:',dim_lag
       WRITE(*, '(a, 10x, a, f10.3)') &
            'PDAF', 'param_real(1) forget=', forget

    END IF writeout

  END SUBROUTINE PDAF_estkf_config

!-------------------------------------------------------------------------------
!> Set integer parameter specific for ESTKF
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Other revisions - see repository log
!!
  SUBROUTINE PDAF_estkf_set_iparam(id, value, flag)

    USE PDAFobs, &
         ONLY: type_obs_init, observe_ens

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
       dim_lag = value
       IF (dim_lag<0) THEN
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(8): Invalid value for smoother lag - param_int(3)!'
          flag = 8
       END IF
    CASE(4)
       ! Not used
    CASE(5)
       type_forget = value
       IF (type_forget<0 .OR. type_forget>1) THEN
          WRITE (*, '(/5x, a/)') 'PDAF-ERROR(8): Invalid type of forgetting factor - param_int(5)!'
          flag = 8
       END IF
    CASE(6)
       type_trans = value
       IF (type_trans<0 .OR. type_trans>2) THEN
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(8): Invalid setting for ensemble transformation - param_int(6)!'
          flag = 8
       END IF
    CASE(7)
       type_sqrt = value
       IF (type_sqrt<0 .OR. type_sqrt>1) THEN
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(8): Invalid setting for square root type - param_int(7)!'
          flag = 8
       END IF
    CASE(8)
       if (value==0) THEN
          observe_ens = .false. ! Apply H to ensemble mean to compute residual
       ELSE
          observe_ens = .true.  ! Apply H to X, compute mean of HX and then residual
       END IF
    CASE(9)
       type_obs_init = value    ! Initialize obs (0) before or (1) after prepoststep
       IF (type_obs_init<0 .OR. type_obs_init>1) THEN
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(8): Invalid setting type_obs_init - param_int(9)!'
          flag = 8
       END IF
    CASE DEFAULT
       WRITE (*,'(/5x, a, i3/)') &
            'PDAF-WARNING: Invalid integer parameter index', id
    END SELECT

  END SUBROUTINE PDAF_estkf_set_iparam


!-------------------------------------------------------------------------------
!> Set real parameter specific for ESTKF
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Other revisions - see repository log
!!
  SUBROUTINE PDAF_estkf_set_rparam(id, value, flag)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)  :: id       !< Index of parameter
    REAL, INTENT(in)     :: value    !< Parameter value
    INTEGER, INTENT(out) :: flag     !< Status flag: 0 for no error


! ****************************
! *** INITIALIZE VARIABLES ***
! ****************************

    ! Initialize status flag
    flag = 0

    SELECT CASE(id) 
    CASE(1)
       forget = value
       IF (forget <= 0.0) THEN
          WRITE (*,'(/5x,a/)') &
               'PDAF-ERROR(7): Invalid value of forgetting factor - param_real(1)!'
          flag = 7
       END IF
    CASE DEFAULT
       WRITE (*,'(/5x, a, i3/)') &
            'PDAF-WARNING: Invalid real parameter index', id
    END SELECT

  END SUBROUTINE PDAF_estkf_set_rparam

!-------------------------------------------------------------------------------
!> Information output on options for ESTKF
!!
!! Subroutine to perform information output on options
!! available for the ESTKF filter.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __REVISION HISTORY:__
!! * 2011-08 - Lars Nerger - Initial code
!! *  Other revisions - see repository log
!!
  SUBROUTINE PDAF_estkf_options()

    IMPLICIT NONE

! *********************
! *** Screen output ***
! *********************

    WRITE(*, '(/a, 5x, a)') 'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    WRITE(*, '(a, 5x, a)')  'PDAF', '+++  Error Subspace Transform Kalman Filter (ESTKF) +++'
    WRITE(*, '(a, 5x, a)')  'PDAF', '+++                                                 +++'
    WRITE(*, '(a, 5x, a)')  'PDAF', '+++  Nerger et al., Mon. Wea. Rev. 140 (2012) 2335  +++'
    WRITE(*, '(a, 5x, a)')  'PDAF', '+++           doi:10.1175/MWR-D-11-00102.1          +++'
    WRITE(*, '(a, 5x, a)')  'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'

    WRITE(*, '(/a, 5x, a)') 'PDAF', 'Available options for ESTKF:'

    WRITE(*, '(a, 5x, a)') 'PDAF', '--- Sub-types (Parameter subtype) ---'
    WRITE(*, '(a, 7x, a)') 'PDAF', '0: Standard implementation with ensemble integration'
    WRITE(*, '(a, 7x, a)') 'PDAF', '10: Fixed error space basis'
    WRITE(*, '(a, 7x, a)') 'PDAF', '11: Fixed state covariance matrix'

    WRITE(*, '(a, 5x, a)') 'PDAF', '--- Integer parameters (Array param_int) ---'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(1): Dimension of state vector (>0), required'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(2): Ensemble size (>0), required'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(3): dim_lag'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'Size of smoothing lag (>=0), optional'
    WRITE(*, '(a, 12x, a)') 'PDAF', '0: no smoothing (default)'
    WRITE(*, '(a, 12x, a)') 'PDAF', '>0: apply smoother up to specified lag'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(4): not used'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(5) type_forget'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'Type of forgetting factor; optional'
    WRITE(*, '(a, 12x, a)') 'PDAF', '0: fixed forgetting factor (default)'
    WRITE(*, '(a, 12x, a)') 'PDAF', '1: adaptive forgetting factor (experimental)'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(6) type_trans'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'Type of ensemble transformation matrix; optional'
    WRITE(*, '(a, 12x, a)') 'PDAF', '0: deterministic Omega (default)'
    WRITE(*, '(a, 12x, a)') 'PDAF', '1: random orthonormal Omega orthogonal to (1,...,1)^T'
    WRITE(*, '(a, 12x, a)') &
         'PDAF', '2: use product of 0 with random orthonomal matrix with eigenvector (1,...,1)^T'
    WRITE(*, '(a, 14x, a)') &
         'PDAF', '(experimental; for random transformations, 0 or 1 are recommended)'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(7) type_sqrt'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'Type of transformation matrix square root; optional'
    WRITE(*, '(a, 12x, a)') 'PDAF', '0: symmetric square root (default)'
    WRITE(*, '(a, 12x, a)') 'PDAF', '1: Cholesky decomposition'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(8): observe_ens'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'Application of observation operator H, optional'
    WRITE(*, '(a, 12x, a)') 'PDAF', '0: Apply H to ensemble mean to compute innovation'
    WRITE(*, '(a, 12x, a)') 'PDAF', '1: Apply H to ensemble states; then compute innovation from their mean (default)'
    WRITE(*, '(a, 12x, a)') 'PDAF', '   param_int(8)=1 is the recomended choice for nonlinear H'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(9): type_obs_init'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'Initialize observations before or after call to prepoststep_pdaf'
    WRITE(*, '(a, 11x, a)') 'PDAF', '0: Initialize observations before call to prepoststep_pdaf'
    WRITE(*, '(a, 11x, a)') 'PDAF', '1: Initialize observations after call to prepoststep_pdaf (default)'


    WRITE(*, '(a, 5x, a)') 'PDAF', '--- Floating point parameters (Array param_real) ---'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_real(1): forget'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'Forgetting factor (usually >0 and <=1), required'

    WRITE(*, '(a, 5x, a)') 'PDAF', '--- Further parameters ---'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'n_modeltasks: Number of parallel model integration tasks'
    WRITE(*, '(a, 11x, a)') &
         'PDAF', '>=1 for subtype 0; not larger than total number of processors'
    WRITE(*, '(a, 11x, a)') 'PDAF', '=1 required for subtypes 10 and 11'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'screen: Control verbosity of PDAF'
    WRITE(*, '(a, 11x, a)') 'PDAF', '0: no outputs'
    WRITE(*, '(a, 11x, a)') 'PDAF', '1: basic output (default)'
    WRITE(*, '(a, 11x, a)') 'PDAF', '2: 1 plus timing output'
    WRITE(*, '(a, 11x, a)') 'PDAF', '3: 2 plus debug output'

    WRITE(*, '(a, 5x, a)') &
         'PDAF', '+++++++++ End of option overview for the ESTKF  ++++++++++'

  END SUBROUTINE PDAF_estkf_options


!-------------------------------------------------------------------------------
!> Display timing and memory information for ESTKF
!!
!! This routine displays the PDAF-internal timing and
!! memory information for the ESTKF filter.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2011-09 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_estkf_memtime(printtype)

    USE PDAF_timer, &
         ONLY: PDAF_time_tot
    USE PDAF_memcounting, &
         ONLY: PDAF_memcount_get, PDAF_memcount_get_global
    USE PDAF_mod_core, &
         ONLY: subtype_filter, offline_mode
    USE PDAF_mod_parallel, &
         ONLY: filterpe, mype_world, COMM_pdaf
    USE PDAFomi_obs_f, &
         ONLY: omi_was_used

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: printtype    !< Type of screen output:  
                                      !< (1) timings, (2) memory

! *** Local variables ***
    INTEGER :: i                        ! Counter
    REAL :: memcount_global(4)          ! Globally counted memory
    REAL :: time_omi                    ! Sum of timers for OMI-internal call-back routines


! ********************************
! *** Print screen information ***
! ********************************

    ptype: IF (printtype == 1) THEN

! **************************************
! *** Print basic timing information ***
! **************************************

       ! Generic part
       WRITE (*, '(//a, 21x, a)') 'PDAF', 'PDAF Timing information'
       WRITE (*, '(a, 10x, 45a)') 'PDAF', ('-', i=1, 45)
       WRITE (*, '(a, 18x, a, F11.3, 1x, a)') &
            'PDAF', 'Initialize PDAF:', pdaf_time_tot(1), 's'
       IF (.not.offline_mode) THEN
          IF (subtype_filter<10) THEN
             WRITE (*, '(a, 16x, a, F11.3, 1x, a)') 'PDAF', 'Ensemble forecast:', pdaf_time_tot(2), 's'
          ELSE
             WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'State forecast:', pdaf_time_tot(2), 's'
          END IF
       END IF

       IF (filterpe) THEN
          ! Filter-specific part
          WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'ESTKF analysis:', pdaf_time_tot(3), 's'

          ! Generic part B
          WRITE (*, '(a, 22x, a, F11.3, 1x, a)') 'PDAF', 'Prepoststep:', pdaf_time_tot(5), 's'
       END IF

    ELSE IF (printtype == 2) THEN ptype

! *****************************************
! *** Formerly: Print allocated memory  ***
! *****************************************

       WRITE (*, '(/a, 23x, a)') 'PDAF', 'PDAF Memory overview'
       WRITE (*, '(/a, 23x, a)') 'PDAF', 'Note: The memory overview is moved to printtype=10 and printtype=11'

    ELSE IF (printtype == 3) THEN ptype

! *******************************************************
! *** Print timing information for call-back routines ***
! *******************************************************

       ! Generic part
       WRITE (*, '(//a, 12x, a)') 'PDAF', 'PDAF Timing information - call-back routines'
       WRITE (*, '(a, 8x, 52a)') 'PDAF', ('-', i=1, 52)
       WRITE (*, '(a, 10x, a, 15x, F11.3, 1x, a)') 'PDAF', 'Initialize PDAF:', pdaf_time_tot(1), 's'
       WRITE (*, '(a, 12x, a, 17x, F11.3, 1x, a)') 'PDAF', 'init_ens_pdaf:', pdaf_time_tot(39), 's'
       IF (.not.offline_mode) THEN
          IF (subtype_filter<10) THEN
             WRITE (*, '(a, 10x, a, 13x, F11.3, 1x, a)') 'PDAF', 'Ensemble forecast:', pdaf_time_tot(2), 's'
          ELSE
             WRITE (*, '(a, 10x, a, 17x, F11.3, 1x, a)') 'PDAF', 'State forecast:', pdaf_time_tot(2), 's'
          END IF
          WRITE (*, '(a, 12x, a, 5x, F11.3, 1x, a)') 'PDAF', 'MPI communication in PDAF:', pdaf_time_tot(4), 's'
          WRITE (*, '(a, 12x, a, 9x, F11.3, 1x, a)') 'PDAF', 'distribute_state_pdaf:', pdaf_time_tot(40), 's'
          WRITE (*, '(a, 12x, a, 12x, F11.3, 1x, a)') 'PDAF', 'collect_state_pdaf:', pdaf_time_tot(41), 's'
          IF (.not.filterpe) WRITE (*, '(a, 7x, a)') 'PDAF', &
               'Note: for filterpe=F, the time (2) includes the wait time for the analysis step'
       END IF

       IF (filterpe) THEN
          ! Filter-specific part
          WRITE (*, '(a, 10x, a, 16x, F11.3, 1x, a)') 'PDAF', 'ESTKF analysis:', pdaf_time_tot(3), 's'
          WRITE (*, '(a, 12x, a, 6x, F11.3, 1x, a)') 'PDAF', 'PDAF-internal operations:', pdaf_time_tot(51), 's'

          IF (omi_was_used) THEN
             ! Output when using OMI

             time_omi = pdaf_time_tot(50) + pdaf_time_tot(48)
             IF (type_forget==1) &
                  time_omi = time_omi + pdaf_time_tot(49) 
             WRITE (*, '(a, 12x, a, 9x, F11.3, 1x, a)') 'PDAF', 'OMI-internal routines:', &
                  time_omi, 's'
             WRITE (*, '(a, 12x, a)') 'PDAF', 'Time in OMI observation module routines '
             WRITE (*, '(a, 14x, a, 8x, F11.3, 1x, a)') 'PDAF', 'init_dim_obs_pdafomi:', pdaf_time_tot(43), 's'
             WRITE (*, '(a, 14x, a, 14x, F11.3, 1x, a)') 'PDAF', 'obs_op_pdafomi:', pdaf_time_tot(44), 's'

!            WRITE (*, '(a, 12x, a, 11x, F11.3, 1x, a)') 'PDAF', 'Time in OMI-internal routines'
!            WRITE (*, '(a, 14x, a, 12x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_init_obs:', pdaf_time_tot(50), 's'
!            IF (type_forget==1) THEN
!               WRITE (*, '(a, 14x, a, 9x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_init_obsvar:', pdaf_time_tot(49), 's'
!            END IF
!            WRITE (*, '(a, 14x, a, 11x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_prodRinvA:', pdaf_time_tot(48), 's'
          ELSE
             ! Output when NOT using OMI
             WRITE (*, '(a, 12x, a, 13x, F11.3, 1x, a)') 'PDAF', 'init_dim_obs_pdaf:', pdaf_time_tot(43), 's'
             WRITE (*, '(a, 12x, a, 19x, F11.3, 1x, a)') 'PDAF', 'obs_op_pdaf:', pdaf_time_tot(44), 's'
             WRITE (*, '(a, 12x, a, 17x, F11.3, 1x, a)') 'PDAF', 'init_obs_pdaf:', pdaf_time_tot(50), 's'
             IF (type_forget==1) THEN
                WRITE (*, '(a, 12x, a, 14x, F11.3, 1x, a)') 'PDAF', 'init_obsvar_pdaf:', pdaf_time_tot(49), 's'
             END IF
             WRITE (*, '(a, 12x, a, 16x, F11.3, 1x, a)') 'PDAF', 'prodRinvA_pdaf:', pdaf_time_tot(48), 's'
          END IF

          ! Generic part B
          WRITE (*, '(a, 10x, a, 14x, F11.3, 1x, a)') 'PDAF', 'prepoststep_pdaf:', pdaf_time_tot(5), 's'
       END IF

    ELSE IF (printtype == 4) THEN ptype

! *********************************************
! *** Print second-level timing information ***
! *********************************************

       ! Generic part
       WRITE (*, '(//a, 21x, a)') 'PDAF', 'PDAF Timing information'
       WRITE (*, '(a, 10x, 45a)') 'PDAF', ('-', i=1, 45)
       WRITE (*, '(a, 21x, a, F11.3, 1x, a)') 'PDAF', 'Initialize PDAF (1):', pdaf_time_tot(1), 's'
       IF (.not.offline_mode) THEN
          IF (subtype_filter<10) THEN
             WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'Ensemble forecast (2):', pdaf_time_tot(2), 's'
          ELSE
             WRITE (*, '(a, 22x, a, F11.3, 1x, a)') 'PDAF', 'State forecast (2):', pdaf_time_tot(2), 's'
          END IF
          WRITE (*, '(a, 12x, a, F11.3, 1x, a)') 'PDAF', 'MPI communication in PDAF (4):', pdaf_time_tot(4), 's'
          IF (.not.filterpe) WRITE (*, '(a, 7x, a)') 'PDAF', &
               'Note: for filterpe=F, the time (2) includes the wait time for the analysis step'
       END IF

       IF (filterpe) THEN
          ! Filter-specific part
          WRITE (*, '(a, 22x, a, F11.3, 1x, a)') 'PDAF', 'ESTKF analysis (3):', pdaf_time_tot(3), 's'
          WRITE (*, '(a, 24x, a, F11.3, 1x, a)') 'PDAF', 'get mean state (9):', pdaf_time_tot(9), 's'
          WRITE (*, '(a, 18x, a, F11.3, 1x, a)') 'PDAF', 'prepare observations (6):', pdaf_time_tot(6), 's'
          WRITE (*, '(a, 24x, a, F11.3, 1x, a)') 'PDAF', 'init residual (10):', pdaf_time_tot(10), 's'
          WRITE (*, '(a, 25x, a, F11.3, 1x, a)') 'PDAF', 'compute Ainv (11):', pdaf_time_tot(11), 's'
          WRITE (*, '(a, 14x, a, F11.3, 1x, a)') 'PDAF', 'get state weight vector (12):', pdaf_time_tot(12), 's'
          IF (subtype_filter /= 11) THEN
             WRITE (*, '(a, 13x, a, F11.3, 1x, a)') 'PDAF', 'prepare ensemble weights (20):', pdaf_time_tot(20), 's'
             WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'transform ensemble (21):', pdaf_time_tot(21), 's'
             IF (dim_lag >0) &
                  WRITE (*, '(a, 20x, a, F11.3, 1x, a)') 'PDAF', 'perform smoothing (15):', pdaf_time_tot(15), 's'
          ELSE
             WRITE (*, '(a, 12x, a, F11.3, 1x, a)') 'PDAF', 'update state and ensemble (21):', pdaf_time_tot(21), 's'
          END IF

          ! Generic part B
          WRITE (*, '(a, 25x, a, F11.3, 1x, a)') 'PDAF', 'Prepoststep (5):', pdaf_time_tot(5), 's'
       END IF

    ELSE IF (printtype == 5) THEN ptype

! *****************************************
! *** Print detailed timing information ***
! *****************************************

       ! Generic part
       WRITE (*, '(//a, 21x, a)') 'PDAF', 'PDAF Timing information'
       WRITE (*, '(a, 10x, 45a)') 'PDAF', ('-', i=1, 45)
       WRITE (*, '(a, 21x, a, F11.3, 1x, a)') 'PDAF', 'Initialize PDAF (1):', pdaf_time_tot(1), 's'
       IF (.not.offline_mode) THEN
          IF (subtype_filter<10) THEN
             WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'Ensemble forecast (2):', pdaf_time_tot(2), 's'
          ELSE
             WRITE (*, '(a, 22x, a, F11.3, 1x, a)') 'PDAF', 'State forecast (2):', pdaf_time_tot(2), 's'
          END IF
          WRITE (*, '(a, 12x, a, F11.3, 1x, a)') 'PDAF', 'MPI communication in PDAF (19):', pdaf_time_tot(19), 's'
          IF (.not.filterpe) WRITE (*, '(a, 7x, a)') 'PDAF', &
               'Note: for filterpe=F, the time (2) includes the wait time for the analysis step'
       END IF

       IF (filterpe) THEN
          ! Filter-specific part
          WRITE (*, '(a, 22x, a, F11.3, 1x, a)') 'PDAF', 'ESTKF analysis (3):', pdaf_time_tot(3), 's'
          WRITE (*, '(a, 24x, a, F11.3, 1x, a)') 'PDAF', 'get mean state (9):', pdaf_time_tot(9), 's'
          WRITE (*, '(a, 18x, a, F11.3, 1x, a)') 'PDAF', 'prepare observations (6):', pdaf_time_tot(6), 's'
          WRITE (*, '(a, 24x, a, F11.3, 1x, a)') 'PDAF', 'init residual (10):', pdaf_time_tot(10), 's'
          WRITE (*, '(a, 25x, a, F11.3, 1x, a)') 'PDAF', 'compute Ainv (11):', pdaf_time_tot(11), 's'
          WRITE (*, '(a, 14x, a, F11.3, 1x, a)') 'PDAF', 'get state weight vector (12):', pdaf_time_tot(12), 's'
          IF (subtype_filter /= 11) THEN
             WRITE (*, '(a, 13x, a, F11.3, 1x, a)') 'PDAF', 'prepare ensemble weights (20):', pdaf_time_tot(20), 's'
             WRITE (*, '(a, 29x, a, F11.3, 1x, a)') 'PDAF', 'SQRT(Ainv) (32):', pdaf_time_tot(32), 's'
             WRITE (*, '(a, 29x, a, F11.3, 1x, a)') 'PDAF', 'init Omega (33):', pdaf_time_tot(33), 's'
             WRITE (*, '(a, 21x, a, F11.3, 1x, a)') 'PDAF', 'compute C^T OmegaT (34):', pdaf_time_tot(34), 's'
             WRITE (*, '(a, 23x, a, F11.3, 1x, a)') 'PDAF', 'complete weights (35):', pdaf_time_tot(35), 's'
             WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'transform ensemble (21):', pdaf_time_tot(21), 's'
             IF (dim_lag >0) &
                  WRITE (*, '(a, 20x, a, F11.3, 1x, a)') 'PDAF', 'perform smoothing (15):', pdaf_time_tot(15), 's'
          ELSE
             WRITE (*, '(a, 12x, a, F11.3, 1x, a)') 'PDAF', 'update state and ensemble (21):', pdaf_time_tot(21), 's'
          END IF

          ! Generic part B
          WRITE (*, '(a, 25x, a, F11.3, 1x, a)') 'PDAF', 'Prepoststep (5):', pdaf_time_tot(5), 's'
       END IF

    ELSE IF (printtype == 10) THEN ptype

! *******************************
! *** Print allocated memory  ***
! *******************************

       WRITE (*, '(/a, 23x, a)') 'PDAF', 'PDAF Memory overview'
       WRITE (*, '(a, 10x, 45a)') 'PDAF', ('-', i=1, 45)
       WRITE (*, '(a, 21x, a)') 'PDAF', 'Allocated memory  (MiB)'
       WRITE (*, '(a, 14x, a, 1x, f10.3, a)') &
            'PDAF', 'state and A:', pdaf_memcount_get(1, 'M'), ' MiB (persistent)'
       WRITE (*, '(a, 11x, a, 1x, f10.3, a)') &
            'PDAF', 'ensemble array:', pdaf_memcount_get(2, 'M'), ' MiB (persistent)'
       WRITE (*, '(a, 12x, a, 1x, f10.3, a)') &
            'PDAF', 'analysis step:', pdaf_memcount_get(3, 'M'), ' MiB (temporary)'
       IF (omi_was_used) THEN
          WRITE (*, '(a, 17x, a, 1x, f10.3, a)') &
            'PDAF', 'PDAF-OMI:', pdaf_memcount_get(6, 'M'), ' MiB (temporary)'
       END IF

    ELSE IF (printtype == 11) THEN ptype

! ****************************************
! *** Print globally allocated memory  ***
! ****************************************

       memcount_global(1) = pdaf_memcount_get_global(1, 'M', COMM_pdaf)
       memcount_global(2) = pdaf_memcount_get_global(2, 'M', COMM_pdaf)
       memcount_global(3) = pdaf_memcount_get_global(3, 'M', COMM_pdaf)
       memcount_global(4) = pdaf_memcount_get_global(6, 'M', COMM_pdaf)

       IF (mype_world==0) THEN
          WRITE (*, '(/a, 23x, a)') 'PDAF', 'PDAF Memory overview'
          WRITE (*, '(a, 10x, 45a)') 'PDAF', ('-', i=1, 45)
          WRITE (*, '(a, 17x, a)') 'PDAF', 'Globally allocated memory  (MiB)'
          WRITE (*, '(a, 14x, a, 1x, f12.3, a)') &
               'PDAF', 'state and A:', memcount_global(1), ' MiB (persistent)'
          WRITE (*, '(a, 11x, a, 1x, f12.3, a)') &
               'PDAF', 'ensemble array:', memcount_global(2), ' MiB (persistent)'
          WRITE (*, '(a, 12x, a, 1x, f12.3, a)') &
               'PDAF', 'analysis step:', memcount_global(3), ' MiB (temporary)'
          IF (omi_was_used) THEN
             WRITE (*, '(a, 17x, a, 1x, f12.3, a)') &
                  'PDAF', 'PDAF-OMI:', memcount_global(4), ' MiB (temporary)'
          END IF
       END IF

    END IF ptype

  END SUBROUTINE PDAF_estkf_memtime

END MODULE PDAF_ESTKF
