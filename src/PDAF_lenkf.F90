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
!> Module for LEnKF holding shared parameters and some helper routines
!!
!! This module declares the parameters that are used in LEnKF. 
!! Parameters that are specific for LEnKF are declared while some
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
MODULE PDAF_LEnKF

  USE PDAF_mod_core, &
       ONLY: debug, dim_lag

  IMPLICIT NONE

! *** Integer parameters ***
  INTEGER :: rank_ana_enkf=0 !< Rank to be considered for inversion of HPH in analysis of LEnKF

! *** Real parameters ***
  REAL    :: forget=1.0      !< Forgetting factor


!-------------------------------------------------------------------------------
  
CONTAINS

!>  PDAF-internal initialization of ETKF
!!
!! Initialization of ETKF within PDAF. Performed are:
!! * initialize filter-specific parameters
!! * print screen information on filter configuration.
!!
!!  !  This is a core routine of PDAF and   !
!!  !   should not be changed by the user   !
!!
!! __Revision history:__
!! *  2015-12 - Lars Nerger - Initial code by copying PDAF_enkf_init
!! *  Other revisions - see repository log
!!
  SUBROUTINE PDAF_lenkf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
       ensemblefilter, fixedbasis, verbose, outflag)

    USE PDAF_mod_core, &
         ONLY: localfilter, covarloc, dim_lag
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

    writeout: IF (verbose == 1) THEN
       WRITE(*, '(/a, 5x, a)') 'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       WRITE(*, '(a, 5x, a)') 'PDAF',  '+++    Localized Ensemble Kalman Filter (LEnKF)     +++'
       WRITE(*, '(a, 5x, a)') 'PDAF',  '+++                                                 +++'     
       WRITE(*, '(a, 5x, a)') 'PDAF',  '+++   Evensen, J. Geophys. Res. 99C (1994) 10143    +++'     
       WRITE(*, '(a, 5x, a)') 'PDAF',  '+++ using an ensemble of observations according to  +++'     
       WRITE(*, '(a, 5x, a)') 'PDAF',  '+++ Burgers et al., Mon. Wea. Rev. 126 (1998) 1719  +++'     
       WRITE(*, '(a, 5x, a)') 'PDAF',  '+++          This implementation follows            +++'
       WRITE(*, '(a, 5x, a)') 'PDAF',  '+++      Nerger et al., Tellus 57A (2005) 715       +++'
       WRITE(*, '(a, 5x, a)') 'PDAF',  '+++   The localization is covariance lozalization   +++'
       WRITE(*, '(a, 5x, a)') 'PDAF',  '+++        of PH^T and HPH^T as described in        +++'
       WRITE(*, '(a, 5x, a)') 'PDAF',  '+++   Houtekamer & Mitchell, MWR, 129 (2001) 123    +++'
       WRITE(*, '(a, 5x, a)') 'PDAF',  '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    END IF writeout


! ****************************
! *** INITIALIZE VARIABLES ***
! ****************************

    ! Set parameter default values - other defaults are set directly in the module
    observe_ens = .false.
    dim_lag = 0

    ! Parse provided parameters
    DO i=3, dim_pint
       CALL PDAF_lenkf_set_iparam(i, param_int(i), outflag)
    END DO
    DO i=1, dim_preal
       CALL PDAF_lenkf_set_rparam(i, param_real(i), outflag)
    END DO


    ! Define whether filter is mode-based or ensemble-based
    ensemblefilter = .TRUE.

    ! Define whether filter is a domain-local filter
    localfilter = 0

    ! Define whether the filter uses covariance localization
    covarloc = 1

    ! Initialize flag for fixed-basis filters
    fixedbasis = .FALSE.


! *********************
! *** Check subtype ***
! *********************

    IF (subtype /= 0) THEN
       WRITE (*, '(/5x, a/)') 'PDAF-ERROR(3): No valid subtype!'
       outflag = 3
    END IF

  END SUBROUTINE PDAF_lenkf_init


!-------------------------------------------------------------------------------
!> Perform allocation of arrays for LEnKF.
!!
!! __Revision history:__
!! * 2015-12 - Lars Nerger - Initial code by copying PDAF_enkf_alloc
!! * 2025-02 - Lars Nerger - Restructuring introducing generic PDAF_alloc
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_lenkf_alloc(outflag)

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

    CALL PDAF_alloc(dim_p, dim_ens, dim_ens_l, 1, dim_bias_p, &
         dim_lag, 0, outflag)

  END SUBROUTINE PDAF_lenkf_alloc


!-------------------------------------------------------------------------------
!>  Print information on configuration of LEnKF
!!
!!  !  This is a core routine of PDAF and   !
!!  !   should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code by splitting from PDAF_seik_init
!! *  Other revisions - see repository log
!!
  SUBROUTINE PDAF_lenkf_config(subtype, verbose)

    USE PDAF_mod_core, &
         ONLY: dim_ens
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

       WRITE (*, '(/a, 4x, a)') 'PDAF', 'LEnKF configuration'
       WRITE (*, '(a, 10x, a, i5)') 'PDAF', 'ensemble size:', dim_ens
       WRITE (*, '(a, 10x, a, i1)') 'PDAF', 'filter sub-type= ', subtype
       IF (subtype == 0) THEN
          WRITE (*, '(a, 14x, a)') 'PDAF', '--> local EnKF (analysis for small observation dimension)'
       END IF
       WRITE(*, '(a, 10x, a, i3)') &
            'PDAF', 'param_int(4) rank_ana_enkf=', rank_ana_enkf
       IF (rank_ana_enkf == 0) THEN
          WRITE (*, '(a, 12x, a)') &
               'PDAF', '---> analysis with direct inversion'
       ELSE
          WRITE (*, '(a, 12x, a, i5)') &
               'PDAF', '--->analysis with pseudo-inverse of HPH, rank=', rank_ana_enkf
       END IF
       WRITE (*, '(a, 10x, a, f5.2)') 'PDAF' ,'Use fixed forgetting factor:', forget
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
       WRITE(*, '(a, 10x, a, f10.3)') &
            'PDAF', 'param_real(1) forget=', forget

    END IF writeout

  END SUBROUTINE PDAF_lenkf_config


!-------------------------------------------------------------------------------
!> Set integer parameter specific for LEnKF
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Other revisions - see repository log
!!
  SUBROUTINE PDAF_lenkf_set_iparam(id, value, flag)

    USE PDAF_mod_core, &
         ONLY: dim_ens
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
       ! Not used
    CASE(4)
       IF (value==0 .AND. value < dim_ens) THEN
          rank_ana_enkf = value
       ELSE
          WRITE (*,'(/5x, a/)') &
             'PDAF-ERROR(8): Invalid setting of param_int(3)/rank_ana_enkf!'
          flag = 8
          rank_ana_enkf = 0 ! Just for safety: Fall back to default
       END IF
    CASE(5)
       ! Not used
    CASE(6)
       ! Not used
    CASE(7)
       ! Not used
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
               'PDAF-ERROR(10): Invalid setting type_obs_init - param_int(9)!'
          flag = 8
       END IF
    CASE DEFAULT
       WRITE (*,'(/5x, a, i3/)') &
            'PDAF-WARNING: Invalid integer parameter index', id
    END SELECT

  END SUBROUTINE PDAF_lenkf_set_iparam


!-------------------------------------------------------------------------------
!> Set real parameter specific for LEnKF
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Other revisions - see repository log
!!
  SUBROUTINE PDAF_lenkf_set_rparam(id, value, flag)

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

  END SUBROUTINE PDAF_lenkf_set_rparam

!-------------------------------------------------------------------------------
!> Information output on options for EnKF
!!
!! Subroutine to perform information output on options
!! available for the EnKF filter.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __REVISION HISTORY:__
!! * 2015-12 - Lars Nerger - Initial code by copying and adapting PDAF_enkf_options
!! *  Other revisions - see repository log
!!
  SUBROUTINE PDAF_lenkf_options()

    IMPLICIT NONE

! *********************
! *** Screen output ***
! *********************

    WRITE(*, '(/a, 5x, a)') 'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    WRITE(*, '(a, 5x, a)') 'PDAF',  '+++    Localized Ensemble Kalman Filter (LEnKF)     +++'
    WRITE(*, '(a, 5x, a)') 'PDAF',  '+++                                                 +++'     
    WRITE(*, '(a, 5x, a)') 'PDAF',  '+++   Evensen, J. Geophys. Res. 99C (1994) 10143    +++'     
    WRITE(*, '(a, 5x, a)') 'PDAF',  '+++ using an ensemble of observations according to  +++'     
    WRITE(*, '(a, 5x, a)') 'PDAF',  '+++ Burgers et al., Mon. Wea. Rev. 126 (1998) 1719  +++'     
    WRITE(*, '(a, 5x, a)') 'PDAF',  '+++          This implementation follows            +++'
    WRITE(*, '(a, 5x, a)') 'PDAF',  '+++      Nerger et al., Tellus 57A (2005) 715       +++'
    WRITE(*, '(a, 5x, a)') 'PDAF',  '+++   The localization is covariance lozalization   +++'
    WRITE(*, '(a, 5x, a)') 'PDAF',  '+++        of PH^T and HPH^T as described in        +++'
    WRITE(*, '(a, 5x, a)') 'PDAF',  '+++   Houtekamer & Mitchell, MWR, 129 (2001) 123    +++'
    WRITE(*, '(a, 5x, a)') 'PDAF',  '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'

    WRITE(*, '(/a, 5x, a)') 'PDAF', 'Available options for LEnKF:'

    WRITE(*, '(a, 5x, a)') 'PDAF', '--- Sub-types (Parameter subtype) ---'
    WRITE(*, '(a, 7x, a)') 'PDAF', '0: Standard EnKF analysis with covariance localization'

    WRITE(*, '(a, 5x, a)') 'PDAF', '--- Integer parameters (Array param_int) ---'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(1): Dimension of state vector (>0), required'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(2): Ensemble size (>0), required'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(3): not used'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(4): rank_ana_enkf'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'maximum rank for inversion of HPH^T, optional, default=0'
    WRITE(*, '(a, 12x, a)') 'PDAF', 'for =0, HPH is inverted by solving the representer equation'
    WRITE(*, '(a, 12x, a)') 'PDAF', 'allowed range is 0 to ensemble size - 1'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(5): not used'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(6): not used'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(7): not used'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(8): observe_ens'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'Application of observation operator H, optional'
    WRITE(*, '(a, 12x, a)') 'PDAF', '0: Apply H to ensemble mean to compute innovation'
    WRITE(*, '(a, 12x, a)') 'PDAF', '1: Apply H to ensemble states; then compute innovation from their mean (default)'
    WRITE(*, '(a, 12x, a)') 'PDAF', '   param_int(8)=1 is the recomended choice for nonlinear H'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(9): type_obs_init'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'Initialize observations before or after call to prepoststep_pdaf'
    WRITE(*, '(a, 12x, a)') 'PDAF', '0: Initialize observations before call to prepoststep_pdaf'
    WRITE(*, '(a, 12x, a)') 'PDAF', '1: Initialize observations after call to prepoststep_pdaf (default)'

    WRITE(*, '(a, 5x, a)') 'PDAF', '--- Floating point parameters (Array param_real) ---'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_real(1): forget'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'Forgetting factor (usually >0 and <=1), required'

    WRITE(*, '(a, 5x, a)') 'PDAF', '--- Further parameters ---'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'n_modeltasks: Number of parallel model integration tasks'
    WRITE(*, '(a, 11x, a)') &
         'PDAF', '(>=1; not larger than total number of processors)'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'screen: Control verbosity of PDAF'
    WRITE(*, '(a, 11x, a)') 'PDAF', '0: no outputs'
    WRITE(*, '(a, 11x, a)') 'PDAF', '1: basic output (default)'
    WRITE(*, '(a, 11x, a)') 'PDAF', '2: 1 plus timing output'
    WRITE(*, '(a, 11x, a)') 'PDAF', '3: 2 plus debug output'


    WRITE(*, '(a, 5x, a)') &
         'PDAF', '+++++++++ End of option overview for the LEnKF ++++++++++'

  END SUBROUTINE PDAF_lenkf_options


!-------------------------------------------------------------------------------
!> Display timing and memory information for LEnKF
!!
!! This routine displays the PDAF-internal timing and
!! memory information for the LEnKF.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2015-12 - Lars Nerger - Initial code by copying PDAF_enkf_memtime
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_lenkf_memtime(printtype)

    USE PDAF_timer, &
         ONLY: PDAF_time_tot
    USE PDAF_memcounting, &
         ONLY: PDAF_memcount_get, PDAF_memcount_get_global
    USE PDAF_mod_core, &
         ONLY: subtype_filter, offline_mode
    USE PDAF_enkf, &
         ONLY: rank_ana_enkf
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
    REAL :: memcount_global(3)          ! Globally counted memory
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
       WRITE (*, '(a, 18x, a, F11.3, 1x, a)') 'PDAF', 'Initialize PDAF:', pdaf_time_tot(1), 's'
       IF (.not.offline_mode) THEN
          IF (subtype_filter<10) THEN
             WRITE (*, '(a, 16x, a, F11.3, 1x, a)') 'PDAF', 'Ensemble forecast:', pdaf_time_tot(2), 's'
          ELSE
             WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'State forecast:', pdaf_time_tot(2), 's'
          END IF
       END IF

       IF (filterpe) THEN
          ! Filter-specific part
          WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'LEnKF analysis:', pdaf_time_tot(3), 's'

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
          WRITE (*, '(a, 10x, a, 16x, F11.3, 1x, a)') 'PDAF', 'LEnKF analysis:', pdaf_time_tot(3), 's'
          WRITE (*, '(a, 12x, a, 6x, F11.3, 1x, a)') 'PDAF', 'PDAF-internal operations:', pdaf_time_tot(51), 's'

          IF(omi_was_used) THEN
             ! Output when using OMI

             time_omi = pdaf_time_tot(50) + pdaf_time_tot(49) + pdaf_time_tot(46)
             WRITE (*, '(a, 12x, a, 9x, F11.3, 1x, a)') 'PDAF', 'OMI-internal routines:', &
                  time_omi, 's'
             WRITE (*, '(a, 12x, a)') 'PDAF', 'Time in OMI observation module routines '
             WRITE (*, '(a, 14x, a, 8x, F11.3, 1x, a)') 'PDAF', 'init_dim_obs_pdafomi:', pdaf_time_tot(43), 's'
             WRITE (*, '(a, 14x, a, 14x, F11.3, 1x, a)') 'PDAF', 'obs_op_pdafomi:', pdaf_time_tot(44), 's'
             WRITE (*, '(a, 14x, a, 6x, F11.3, 1x, a)') 'PDAF', 'localize_covar_pdafomi:', pdaf_time_tot(45), 's'

!            WRITE (*, '(a, 12x, a, 11x, F11.3, 1x, a)') 'PDAF', 'Time in OMI-internal routines'
!            WRITE (*, '(a, 14x, a, 2x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_add_obs_error_pdaf:', pdaf_time_tot(46), 's'
!            WRITE (*, '(a, 14x, a, 7x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_init_obs_pdaf:', pdaf_time_tot(50), 's'
!            WRITE (*, '(a, 14x, a, 2x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_init_obscovar_pdaf:', pdaf_time_tot(49), 's'
          ELSE
             ! Output when NOT using OMI

             WRITE (*, '(a, 12x, a, 13x, F11.3, 1x, a)') 'PDAF', 'init_dim_obs_pdaf:', pdaf_time_tot(43), 's'
             WRITE (*, '(a, 12x, a, 19x, F11.3, 1x, a)') 'PDAF', 'obs_op_pdaf:', pdaf_time_tot(44), 's'
             WRITE (*, '(a, 12x, a, 11x, F11.3, 1x, a)') 'PDAF', 'localize_covar_pdaf:', pdaf_time_tot(45), 's'
             WRITE (*, '(a, 12x, a, 12x, F11.3, 1x, a)') 'PDAF', 'add_obs_error_pdaf:', pdaf_time_tot(46), 's'
             WRITE (*, '(a, 12x, a, 17x, F11.3, 1x, a)') 'PDAF', 'init_obs_pdaf:', pdaf_time_tot(50), 's'
             WRITE (*, '(a, 12x, a, 12x, F11.3, 1x, a)') 'PDAF', 'init_obscovar_pdaf:', pdaf_time_tot(49), 's'
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
          WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'Ensemble forecast (2):', pdaf_time_tot(2), 's'
          WRITE (*, '(a, 13x, a, F11.3, 1x, a)') 'PDAF', 'MPI communication in PDAF (4):', pdaf_time_tot(4), 's'
          IF (.not.filterpe) WRITE (*, '(a, 7x, a)') 'PDAF', &
               'Note: for filterpe=F, the time (2) includes the wait time for the analysis step'
       END IF

       IF (filterpe) THEN
          ! Filter-specific part
          WRITE (*, '(a, 22x, a, F11.3, 1x, a)') 'PDAF', 'LEnKF analysis (3):', pdaf_time_tot(3), 's'
          WRITE (*, '(a, 24x, a, F11.3, 1x, a)') 'PDAF', 'get mean state (9):', pdaf_time_tot(9), 's'
          WRITE (*, '(a, 18x, a, F11.3, 1x, a)') 'PDAF', 'prepare observations (6):', pdaf_time_tot(6), 's'
          WRITE (*, '(a, 17x, a, F11.3, 1x, a)') 'PDAF', 'compute HPH+R and HP (10):', pdaf_time_tot(10), 's'
          WRITE (*, '(a, 18x, a, F11.3, 1x, a)') 'PDAF', 'localize HP_p and HPH (45):', pdaf_time_tot(45), 's'
          WRITE (*, '(a, 17x, a, F11.3, 1x, a)') 'PDAF', 'sample obs. ensemble (11):', pdaf_time_tot(11), 's'
          WRITE (*, '(a, 13x, a, F11.3, 1x, a)') 'PDAF', 'init innovation ensemble (12):', pdaf_time_tot(12), 's'
          WRITE (*, '(a, 13x, a, F11.3, 1x, a)') 'PDAF', 'compute transform matrix (13):', pdaf_time_tot(13), 's'
          WRITE (*, '(a, 14x, a, F11.3, 1x, a)') 'PDAF', 'ensemble transformation (14):', pdaf_time_tot(14), 's'

          ! Generic part
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
          WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'Ensemble forecast (2):', pdaf_time_tot(2), 's'
          WRITE (*, '(a, 13x, a, F11.3, 1x, a)') 'PDAF', 'MPI communication in PDAF (4):', pdaf_time_tot(4), 's'
          IF (.not.filterpe) WRITE (*, '(a, 7x, a)') 'PDAF', &
               'Note: for filterpe=F, the time (2) includes the wait time for the analysis step'
       END IF

       IF (filterpe) THEN
          ! Filter-specific part
          WRITE (*, '(a, 22x, a, F11.3, 1x, a)') 'PDAF', 'LEnKF analysis (3):', pdaf_time_tot(3), 's'
          WRITE (*, '(a, 24x, a, F11.3, 1x, a)') 'PDAF', 'get mean state (9):', pdaf_time_tot(9), 's'
          WRITE (*, '(a, 18x, a, F11.3, 1x, a)') 'PDAF', 'prepare observations (6):', pdaf_time_tot(6), 's'
          WRITE (*, '(a, 17x, a, F11.3, 1x, a)') 'PDAF', 'compute HPH+R and HP (10):', pdaf_time_tot(10), 's'
          WRITE (*, '(a, 33x, a, F11.3, 1x, a)') 'PDAF', 'HXpert (30):', pdaf_time_tot(30), 's'
          WRITE (*, '(a, 35x, a, F11.3, 1x, a)') 'PDAF', 'HP_p (31):', pdaf_time_tot(31), 's'
          WRITE (*, '(a, 36x, a, F11.3, 1x, a)') 'PDAF', 'HPH (32):', pdaf_time_tot(32), 's'
          WRITE (*, '(a, 18x, a, F11.3, 1x, a)') 'PDAF', 'localize HP_p and HPH (45):', pdaf_time_tot(45), 's'
          WRITE (*, '(a, 27x, a, F11.3, 1x, a)') 'PDAF', 'add matrix R (46):', pdaf_time_tot(46), 's'
          WRITE (*, '(a, 17x, a, F11.3, 1x, a)') 'PDAF', 'sample obs. ensemble (11):', pdaf_time_tot(11), 's'
          WRITE (*, '(a, 13x, a, F11.3, 1x, a)') 'PDAF', 'init innovation ensemble (12):', pdaf_time_tot(12), 's'
          WRITE (*, '(a, 13x, a, F11.3, 1x, a)') 'PDAF', 'compute transform matrix (13):', pdaf_time_tot(13), 's'
          IF (rank_ana_enkf == 0) THEN
             WRITE (*, '(a, 17x, a, F11.3, 1x, a)') 'PDAF', 'solve for representers (33):', pdaf_time_tot(33), 's'
             WRITE (*, '(a, 25x, a, F11.3, 1x, a)') 'PDAF', 'HX*representer (34):', pdaf_time_tot(34), 's'
          ELSE
             WRITE (*, '(a, 16x, a, F11.3, 1x, a)') 'PDAF', 'pseudo inverse of HPH^T (35):', pdaf_time_tot(35), 's'
             WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'compute representers (36):', pdaf_time_tot(36), 's'
             WRITE (*, '(a, 25x, a, F11.3, 1x, a)') 'PDAF', 'HX*representer (37):', pdaf_time_tot(37), 's'
          END IF
          WRITE (*, '(a, 14x, a, F11.3, 1x, a)') 'PDAF', 'ensemble transformation (14):', pdaf_time_tot(14), 's'

          ! Generic part
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


    ELSE IF (printtype == 11) THEN ptype

! ****************************************
! *** Print globally allocated memory  ***
! ****************************************

       memcount_global(1) = pdaf_memcount_get_global(1, 'M', COMM_pdaf)
       memcount_global(2) = pdaf_memcount_get_global(2, 'M', COMM_pdaf)
       memcount_global(3) = pdaf_memcount_get_global(3, 'M', COMM_pdaf)

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
       END IF

    END IF ptype

  END SUBROUTINE PDAF_lenkf_memtime

END MODULE PDAF_LEnKF
