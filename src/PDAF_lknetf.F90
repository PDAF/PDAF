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
!> Module for LKNETF holding shared parameters and some helper routines
!!
!! This module declares the parameters that are used in LKNETF. 
!! Parameters that are specific for LKNETF are declared while some
!! other parameters are use-included from PDAF_mod_filter. This allows
!! us to only include this module in the method-specific analysis routines.
!! In addition, subroutines are included that initialize these parameters.
!!
!!    ! This is a core routine of PDAF and !
!!    ! should not be changed by the user  !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code from restructuring
!! *  Later revisions - see repository log
!!
MODULE PDAF_LKNETF

  USE PDAF_mod_filter, &
       ONLY: localfilter, incremental, debug, dim_lag, &
       member_save, skewness, kurtosis

  IMPLICIT NONE

! *** Integer parameters ***
  INTEGER :: type_hyb=1    !< Type of hybrid weight: 
                           !< (0) fixed
                           !< (2) adaptive
  INTEGER :: type_forget=0 !< Type of forgetting factor
                           !< (0) inflate forecast ensemble
                           !< (1) inflate forecast ensemble only observed domains
                           !< (2) inflate analysis ensemble
                           !< (3) inflate analysis ensemble only observed domains
                           !< (5) globally addaptive forgetting factor
                           !< (6) locally adaptive forgetting factor
  INTEGER :: type_trans=0  !< Type of ensemble transformation
                           !< For LKNETF:
                           !< (0) use product with random orthonomal matrix with
                           !<     eigenvector (1,...,1)^T
                           !< (1) use deterministic transformation 

! *** Real parameters ***
  REAL :: forget=1.0       !< Forgetting factor
  REAL :: hyb_g=0.95       !< Hybrid weight for state in LKNEF (1.0 for LETKF; 0.0 for LKNETF)
  REAL :: hyb_k            !< Hybrid weight norm for using skewness and kurtosis

! *** Internal variables ***
  LOGICAL :: store_rndmat = .false.  ! Whether to recompute or store the random matrix
  LOGICAL :: inloop=.false. ! Whether the program is in the local analysis loop

!-------------------------------------------------------------------------------
  
CONTAINS

!>  PDAF-internal initialization of LKNETF
!!
!! Initialization of LKNETF within PDAF. Performed are:
!! * initialize filter-specific parameters
!! * print screen information on filter configuration.
!!
!!  !  This is a core routine of PDAF and   !
!!  !   should not be changed by the user   !
!!
!! __Revision history:__
!! * 2017-08 - Lars Nerger - Initial code based on code for LNETF
!! *  Later revisions - see repository log
!!
  SUBROUTINE PDAF_lknetf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
       ensemblefilter, fixedbasis, verbose, outflag)

    USE PDAF_mod_filter, &
         ONLY: incremental, dim_lag, dim_ens
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
    INTEGER :: flagsum          ! Sum of status flags


! ****************************
! *** INITIALIZE VARIABLES ***
! ****************************

    ! Set parameter default values
    ! (Other defaults are set in the module)
    incremental = 0
    observe_ens = .false.
    dim_lag = 0

    ! hyb_k needs to be set at runtime
    hyb_k = dim_ens

    ! Parse provided parameters
    flagsum = 0
    DO i=3, dim_pint
       CALL PDAF_lknetf_set_iparam(i, param_int(i), outflag)
       flagsum = flagsum+outflag
    END DO
    DO i=1, dim_preal
       CALL PDAF_lknetf_set_rparam(i, param_real(i), outflag)
       flagsum = flagsum+outflag
    END DO


    ! Define whether filter is mode-based or ensemble-based
    ensemblefilter = .TRUE.

    ! Define whether filter is domain localized
    localfilter = 1

    ! Initialize flag for fixed-basis filters
    fixedbasis = .FALSE.


! *********************
! *** Screen output ***
! *********************

    writeout: IF (verbose == 1) THEN

       WRITE(*, '(/a, 4x, a)') 'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '+++  Local Hybrid Kalman-Nonlinear Ensemble Transform Filter  +++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '+++                                                           +++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '+++                Domain-localized LKNETF by                 +++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '+++ L. Nerger, QJRMS, 148 (2022) 620-640, doi:10.1002/qj.4221 +++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

       IF (flagsum == 0) THEN

          ! *** General output ***
          WRITE (*, '(/a, 4x, a)') 'PDAF', 'LKNETF configuration'
          WRITE (*, '(a, 10x, a, i1)') 'PDAF', 'filter sub-type = ', subtype
          IF (subtype == 0) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> (HNK) 2-step LKNETF: NETF before LETKF'
          ELSE IF (subtype == 1) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> (HKN) 2-step LKNETF: LETKF before NETF'
          ELSE IF (subtype == 4) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> (HSync) LKNETF synchronous'
          ELSE IF (subtype == 5) THEN
          ELSE
             WRITE (*, '(/5x, a/)') 'PDAF-ERROR(2): No valid sub type!'
             outflag = 3
          END IF
          IF (type_trans == 0) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> Transform ensemble including product with random matrix'
          ELSE IF (type_trans == 1) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> Deterministic symmetric ensemble transformation'
          END IF
          IF (incremental == 1) &
               WRITE (*, '(a, 12x, a)') 'PDAF', '--> Perform incremental updating'
          IF (type_forget == 0) THEN
             WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> prior inflation, forgetting factor:', forget
          ELSEIF (type_forget == 1) THEN
             WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> prior inflation on observed domains, forgetting factor: ', forget
          ELSEIF (type_forget == 2) THEN
             WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> posterior inflation, forgetting factor:', forget
          ELSEIF (type_forget == 3) THEN
             WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> posterior inflation on observed domains, forgetting factor: ', forget
          ELSEIF (type_forget == 5) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> Use global adaptive forgetting factor in LETKF'
          ELSEIF (type_forget == 6) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> Use local adaptive forgetting factors in LETKF'
          ENDIF
          WRITE (*, '(a, 10x, a, i1)') 'PDAF', 'hybridization type = ', type_hyb
          WRITE (*, '(a, 12x, a, es10.2)') 'PDAF','--> hybrid weight input (gamma):', hyb_g
          IF (type_hyb==3 .OR. type_hyb==4) &
               WRITE (*, '(a, 12x, a, es10.2)') 'PDAF','--> hybrid norm (kappa):', hyb_k
          IF (type_hyb == 0) THEN
             WRITE(*, '(a, 12x, a, f8.3)') 'PDAF', '--> use gamma_fix: fixed hybrid weight', hyb_g
          ELSEIF (type_hyb == 1) THEN
             WRITE(*, '(a, 12x, a, f8.3)') 'PDAF', '--> use gamma_lin: (1 - N_eff/N_e)*', hyb_g
          ELSEIF (type_hyb == 2) THEN
             WRITE(*, '(a, 12x, a, f8.3)') 'PDAF', '--> use gamma_alpha: hybrid weight from N_eff/N>=', hyb_g
          ELSEIF (type_hyb == 3) THEN
             WRITE(*, '(a, 12x, a, f8.3, a, f8.3)') &
                  'PDAF', '--> use gamma_ska: 1 - min(s,k)/sqrt(', hyb_k, ') with N_eff/N>=', hyb_g
          ELSEIF (type_hyb == 4) THEN
             WRITE(*, '(a, 12x, a, f8.3, a, f8.3)') &
                  'PDAF', '--> use gamma_sklin: 1 - min(s,k)/sqrt(', hyb_k, ') >= 1-N_eff/N>=', hyb_g
          END IF
          WRITE (*, '(a, 12x, a, i5)') 'PDAF', '--> ensemble size N:', dim_ens
          IF (observe_ens) &
               WRITE (*, '(a, 12x, a, 1x, l)') 'PDAF', '--> observe_ens:', observe_ens
       ELSE
          WRITE (*, '(/5x, a/)') 'PDAF-ERROR: Invalid parameter setting - check prior output!'
       END IF

    END IF writeout

  END SUBROUTINE PDAF_lknetf_init


!-------------------------------------------------------------------------------
!> Perform allocation of arrays for LKNETF.
!!
!! __Revision history:__
!! * 2017-08 - Lars Nerger - Initial code based in code for LETKF
!! * 2025-02 - Lars Nerger - Restructuring introducing generic PDAF_alloc
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAF_lknetf_alloc(outflag)

    USE PDAF_mod_filter, &
         ONLY: dim_ens, dim_p, dim_bias_p
    USE PDAF_mod_filtermpi, &
         ONLY: dim_ens_l

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout):: outflag      !< Status flag


! ******************************
! *** Allocate filter fields ***
! ******************************

    CALL PDAF_alloc(dim_p, dim_ens, dim_ens_l, dim_ens, dim_bias_p, &
         dim_lag, 0, 1, outflag)

  END SUBROUTINE PDAF_lknetf_alloc


!-------------------------------------------------------------------------------
!> Set integer parameter specific for LKNETF
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
  SUBROUTINE PDAF_lknetf_set_iparam(id, value, flag)

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
       ! Not used
    CASE(5)
       type_forget = value
       IF (type_forget<0 .OR. type_forget>6 .OR. type_forget==4) THEN
          WRITE (*, '(/5x, a/)') 'PDAF-ERROR(8): Invalid type of forgetting factor - param_int(5)!'
          flag = 8
       END IF
    CASE(6)
       type_trans = value
       IF (type_trans<0 .OR. type_trans>1) THEN
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(8): Invalid setting for ensemble transformation - param_int(6)!'
          flag = 8
       END IF
    CASE(7)
       type_hyb = value
       IF (type_hyb<0 .OR. type_hyb>4) THEN
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(8): Invalid setting for type of hybrid weight - param_int(7)!'
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
            'PDAF-ERROR(8): Invalid integer parameter index'
       flag = 8
    END SELECT

  END SUBROUTINE PDAF_lknetf_set_iparam


!-------------------------------------------------------------------------------
!> Set real parameter specific for LKNETF
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
  SUBROUTINE PDAF_lknetf_set_rparam(id, value, flag)

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
    CASE(2)
       hyb_g = value
       IF (hyb_g < 0.0) THEN
          WRITE (*,'(/5x,a/)') &
               'PDAF-ERROR(8): Hybrid weight gamma should be >0 - param_real(2)!'
          flag = 8
       END IF
    CASE(3)
       hyb_k = value
       IF(hyb_k<0.0) THEN
          WRITE (*,'(/5x,a/)') &
               'PDAF-ERROR(8): Hybrid norm kappa should be >0 - param_real(3)!'
          flag = 8
       END IF
    CASE DEFAULT
       WRITE (*,'(/5x, a, i3/)') &
            'PDAF-WARNING: Invalid real parameter index', id
    END SELECT

  END SUBROUTINE PDAF_lknetf_set_rparam

!-------------------------------------------------------------------------------
!> Information output on options for LNETF
!!
!! Subroutine to perform information output on options
!! available for the LNETF.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __REVISION HISTORY:__
!! * 2018-07 - Lars Nerger - Initial code based on code for LNETF
!! *  Later revisions - see repository log
!!
  SUBROUTINE PDAF_lknetf_options()

    IMPLICIT NONE

! *********************
! *** Screen output ***
! *********************

    WRITE(*, '(/a, 5x, a)') 'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    WRITE(*, '(a, 5x, a)')  'PDAF', '+++  Local Hybrid Kalman-Nonlinear Ensemble Transform Filter  +++'
    WRITE(*, '(a, 5x, a)')  'PDAF', '+++                                                           +++'
    WRITE(*, '(a, 5x, a)')  'PDAF', '+++                Domain-localized LKNETF by                 +++'
    WRITE(*, '(a, 5x, a)')  'PDAF', '+++ L. Nerger, QJRMS, 148 (2022) 620-640, doi:10.1002/qj.4221 +++'
    WRITE(*, '(a, 5x, a)')  'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

    WRITE(*, '(/a, 5x, a)') 'PDAF', 'Available options for LKNETF:'

    WRITE(*, '(a, 5x, a)') 'PDAF', '--- Sub-types (Parameter subtype) ---'
    WRITE(*, '(a, 7x, a)') 'PDAF', '0: HNK: 2-step LKNETF with NETF before LETKF'
    WRITE(*, '(a, 7x, a)') 'PDAF', '1: HKN: 2-step LKNETF with LETKF before NETF'
    WRITE(*, '(a, 7x, a)') 'PDAF', '4: HSync: LKNETF synchronous'

    WRITE(*, '(a, 5x, a)') 'PDAF', '--- Integer parameters (Array param_int) ---'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(1): Dimension of state vector (>0), required'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(2): Ensemble size (>0), required'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(3): not used'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(4): not used'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(5): type_forget'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'Type of forgetting factor; optional'
    WRITE(*, '(a, 12x, a)') 'PDAF', '0: inflate forecast ensemble by 1/forget (default)'
    WRITE(*, '(a, 12x, a)') 'PDAF', '1: inflate forecast ensemble by 1/forget only observed domains'
    WRITE(*, '(a, 12x, a)') 'PDAF', '2: inflate analysis ensemble by 1/forget'
    WRITE(*, '(a, 12x, a)') 'PDAF', '3: inflate analysis ensemble by 1/forget only observed domains'
    WRITE(*, '(a, 12x, a)') 'PDAF', '5: adaptive forgetting factor for full domain in LETKF part'
    WRITE(*, '(a, 12x, a)') 'PDAF', '6: locally adaptive forgetting factor in LETKF part'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(6): type_trans'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'Type of ensemble transformation matrix; optional'
    WRITE(*, '(a, 12x, a)') 'PDAF', '0: random orthonormal matrix orthogonal to (1,...,1)^T (default)'
    WRITE(*, '(a, 12x, a)') 'PDAF', '1: deterministic transformation'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(7): type_hyb'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'Type of hybrid weight; optional'
    WRITE(*, '(a, 12x, a)') 'PDAF', '0: fixed value'
    WRITE(*, '(a, 12x, a)') 'PDAF', '1: gamma_lin: (1 - N_eff/N_e)*param_real(2) (default)'
    WRITE(*, '(a, 12x, a)') 'PDAF', '2: gamma_alpha: hybrid weight from N_eff/N>=param_real(2)'
    WRITE(*, '(a, 12x, a)') 'PDAF', '3: gamma_ska: 1 - min(s,k)/sqrt(param_real(3)) with N_eff/N>=param_real(2)'
    WRITE(*, '(a, 12x, a)') 'PDAF', '4: gamma_sklin: 1 - min(s,k)/sqrt(param_real(3)) >= 1-N_eff/N>=param_real(2)'


    WRITE(*, '(a, 5x, a)') 'PDAF', '--- Floating point parameters (Array param_real) ---'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_real(1): forget'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'Forgetting factor (usually >0 and <=1), required'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_real(2): hyb_g'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'prescribed hybrid weight gamma (usually >0 and <=1), optional, default=0.95'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_real(3): hyb_k'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'hybrid norm kappa (>0), optional, default=dim_ens'

    WRITE(*, '(a, 5x, a)') 'PDAF', '--- Further parameters ---'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'n_modeltasks: Number of parallel model integration tasks'
    WRITE(*, '(a, 11x, a)') &
         'PDAF', '>=1 for subtypes 0 and 1; not larger than total number of processors'
    WRITE(*, '(a, 11x, a)') 'PDAF', '=1 required for subtypes 2 and 3'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'screen: Control verbosity of PDAF'
    WRITE(*, '(a, 11x, a)') 'PDAF', '0: no outputs'
    WRITE(*, '(a, 11x, a)') 'PDAF', '1: basic output (default)'
    WRITE(*, '(a, 11x, a)') 'PDAF', '2: 1 plus timing output'
    WRITE(*, '(a, 11x, a)') 'PDAF', '3: 2 plus debug output'

    WRITE(*, '(a, 5x, a)') &
         'PDAF', '+++++++++ End of option overview for the LKNETF  ++++++++++'

  END SUBROUTINE PDAF_lknetf_options


!-------------------------------------------------------------------------------
!> Display timing and memory information for LKNETF
!!
!! This routine displays the PDAF-internal timing and
!! memory information for the LKNETF.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2022-03 - Lars Nerger - Initial code
!! * Later revisions - see svn log
!!
  SUBROUTINE PDAF_lknetf_memtime(printtype)

    USE PDAF_timer, &
         ONLY: PDAF_time_tot
    USE PDAF_memcounting, &
         ONLY: PDAF_memcount_get, PDAF_memcount_get_global
    USE PDAF_mod_filter, &
         ONLY: subtype_filter, offline_mode
    USE PDAF_mod_filtermpi, &
         ONLY: filterpe, mype_world, COMM_pdaf
    USE PDAFomi, &
         ONLY: omi_was_used
    USE PDAFlocal, &
         ONLY: pdaflocal_was_used

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
       WRITE (*, '(a, 18x, a, F11.3, 1x, a)') &
            'PDAF', 'Initialize PDAF:', pdaf_time_tot(1), 's'
       IF (.not.offline_mode) THEN
          IF (subtype_filter<2) THEN
             WRITE (*, '(a, 16x, a, F11.3, 1x, a)') 'PDAF', 'Ensemble forecast:', pdaf_time_tot(2), 's'
          ELSE
             WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'State forecast:', pdaf_time_tot(2), 's'
          END IF
       END IF

       IF (filterpe) THEN
          ! Filter-specific part
          WRITE (*, '(a, 18x, a, F11.3, 1x, a)') 'PDAF', 'LKNETF analysis:', pdaf_time_tot(3), 's'

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
       WRITE (*, '(a, 8x, 54a)') 'PDAF', ('-', i=1, 54)
       WRITE (*, '(a, 10x, a, 15x, F11.3, 1x, a)') 'PDAF', 'Initialize PDAF:', pdaf_time_tot(1), 's'
       WRITE (*, '(a, 12x, a, 17x, F11.3, 1x, a)') 'PDAF', 'init_ens_pdaf:', pdaf_time_tot(39), 's'
       IF (.not.offline_mode) THEN
          IF (subtype_filter<2) THEN
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
          WRITE (*, '(a, 10x, a, 15x, F11.3, 1x, a)') 'PDAF', 'LKNETF analysis:', pdaf_time_tot(3), 's'
          WRITE (*, '(a, 12x, a, 6x, F11.3, 1x, a)') 'PDAF', 'PDAF-internal operations:', pdaf_time_tot(51), 's'

          IF(omi_was_used) THEN
             ! Output when using OMI

             time_omi = pdaf_time_tot(46) + pdaf_time_tot(47) + pdaf_time_tot(48)
             IF (type_forget==1) &
                  time_omi = time_omi + pdaf_time_tot(50) + pdaf_time_tot(49) + pdaf_time_tot(52)
             WRITE (*, '(a, 12x, a, 9x, F11.3, 1x, a)') 'PDAF', 'OMI-internal routines:', &
                  time_omi, 's'
             WRITE (*, '(a, 12x, a, 11x, F11.3, 1x, a)') 'PDAF', 'init_n_domains_pdaf:', pdaf_time_tot(42), 's'
             WRITE (*, '(a, 12x, a, 15x, F11.3, 1x, a)') 'PDAF', 'init_dim_l_pdaf:', pdaf_time_tot(45), 's'
             IF (.NOT.pdaflocal_was_used) THEN
                WRITE (*, '(a, 12x, a, 16x, F11.3, 1x, a)') 'PDAF', 'g2l_state_pdaf:', pdaf_time_tot(10), 's'
                WRITE (*, '(a, 12x, a, 16x, F11.3, 1x, a)') 'PDAF', 'l2g_state_pdaf:', pdaf_time_tot(14), 's'
             END IF
             WRITE (*, '(a, 12x, a)') 'PDAF', 'Time in OMI observation module routines '
             WRITE (*, '(a, 14x, a, 8x, F11.3, 1x, a)') 'PDAF', 'init_dim_obs_pdafomi:', pdaf_time_tot(43), 's'
             WRITE (*, '(a, 14x, a, 14x, F11.3, 1x, a)') 'PDAF', 'obs_op_pdafomi:', pdaf_time_tot(44), 's'
             WRITE (*, '(a, 14x, a, 6x, F11.3, 1x, a)') 'PDAF', 'init_dim_obs_l_pdafomi:', pdaf_time_tot(38), 's'

!            WRITE (*, '(a, 12x, a, 11x, F11.3, 1x, a)') 'PDAF', 'Time in OMI-internal routines'
!            IF (type_forget==1) THEN
!               WRITE (*, '(a, 14x, a, 10x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_init_obs_f:', pdaf_time_tot(50), 's'
!               WRITE (*, '(a, 14x, a, 9x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_init_obsvar:', pdaf_time_tot(49), 's'
!            END IF
!            WRITE (*, '(a, 14x, a, 13x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_g2l_obs:', pdaf_time_tot(46), 's'
!            IF (type_forget==1) THEN
!               WRITE (*, '(a, 14x, a, 7x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_init_obsvar_l:', pdaf_time_tot(52), 's'
!            END IF
!            WRITE (*, '(a, 14x, a, 10x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_init_obs_l:', pdaf_time_tot(47), 's'
!            WRITE (*, '(a, 14x, a, 3x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_prodRinvA_l_(hyb):', pdaf_time_tot(48), 's'
!            WRITE (*, '(a, 14x, a, 2x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_likelihood_l_(hyb):', pdaf_time_tot(49), 's'
          ELSE
             ! Output when NOT using OMI

             WRITE (*, '(a, 12x, a, 11x, F11.3, 1x, a)') 'PDAF', 'init_n_domains_pdaf:', pdaf_time_tot(42), 's'
             WRITE (*, '(a, 12x, a, 11x, F11.3, 1x, a)') 'PDAF', 'init_dim_obs_f_pdaf:', pdaf_time_tot(43), 's'
             WRITE (*, '(a, 12x, a, 17x, F11.3, 1x, a)') 'PDAF', 'obs_op_f_pdaf:', pdaf_time_tot(44), 's'
             IF (type_forget==1) THEN
                WRITE (*, '(a, 12x, a, 15x, F11.3, 1x, a)') 'PDAF', 'init_obs_f_pdaf:', pdaf_time_tot(50), 's'
                WRITE (*, '(a, 12x, a, 14x, F11.3, 1x, a)') 'PDAF', 'init_obsvar_pdaf:', pdaf_time_tot(49), 's'
             END IF
             WRITE (*, '(a, 12x, a, 15x, F11.3, 1x, a)') 'PDAF', 'init_dim_l_pdaf:', pdaf_time_tot(45), 's'
             WRITE (*, '(a, 12x, a, 11x, F11.3, 1x, a)') 'PDAF', 'init_dim_obs_l_pdaf:', pdaf_time_tot(38), 's'
             WRITE (*, '(a, 12x, a, 16x, F11.3, 1x, a)') 'PDAF', 'g2l_state_pdaf:', pdaf_time_tot(10), 's'
             WRITE (*, '(a, 12x, a, 18x, F11.3, 1x, a)') 'PDAF', 'g2l_obs_pdaf:', pdaf_time_tot(46), 's'
             IF (type_forget==1) THEN
                WRITE (*, '(a, 12x, a, 12x, F11.3, 1x, a)') 'PDAF', 'init_obsvar_l_pdaf:', pdaf_time_tot(52), 's'
             END IF
             WRITE (*, '(a, 12x, a, 15x, F11.3, 1x, a)') 'PDAF', 'init_obs_l_pdaf:', pdaf_time_tot(47), 's'
             WRITE (*, '(a, 12x, a, 8x, F11.3, 1x, a)') 'PDAF', 'prodRinvA_l_(hyb_)pdaf:', pdaf_time_tot(48), 's'
             WRITE (*, '(a, 12x, a, 7x, F11.3, 1x, a)') 'PDAF', 'likelihood_l_(hyb_)pdaf:', pdaf_time_tot(49), 's'
             WRITE (*, '(a, 12x, a, 16x, F11.3, 1x, a)') 'PDAF', 'l2g_state_pdaf:', pdaf_time_tot(14), 's'
          END IF

          ! Generic part B
          WRITE (*, '(a, 10x, a, 14x, F11.3, 1x, a)') 'PDAF', 'prepoststep_pdaf:', pdaf_time_tot(5), 's'
       END IF

    ELSE IF (printtype == 4 .OR. printtype == 5) THEN ptype

! *********************************************
! *** Print second-level timing information ***
! *********************************************

       ! Generic part
       WRITE (*, '(//a, 23x, a)') 'PDAF', 'PDAF Timing information'
       WRITE (*, '(a, 8x, 52a)') 'PDAF', ('-', i=1, 52)
       WRITE (*, '(a, 10x, a, 11x, F11.3, 1x, a)') 'PDAF', 'Initialize PDAF (1):', pdaf_time_tot(1), 's'
       IF (.not.offline_mode) THEN
          IF (subtype_filter<2) THEN
             WRITE (*, '(a, 10x, a, 9x, F11.3, 1x, a)') 'PDAF', 'Ensemble forecast (2):', pdaf_time_tot(2), 's'
          ELSE
             WRITE (*, '(a, 10x, a, 12x, F11.3, 1x, a)') 'PDAF', 'State forecast (2):', pdaf_time_tot(2), 's'
          END IF
          WRITE (*, '(a, 12x, a, 1x, F11.3, 1x, a)') 'PDAF', 'MPI communication in PDAF (4):', pdaf_time_tot(4), 's'
          IF (.not.filterpe) WRITE (*, '(a, 7x, a)') 'PDAF', &
               'Note: for filterpe=F, the time (2) includes the wait time for the analysis step'
       END IF

       IF (filterpe) THEN
          ! Filter-specific part
          WRITE (*, '(a, 10x, a, 11x, F11.3, 1x, a)') 'PDAF', 'LKNETF analysis (3):', pdaf_time_tot(3), 's'
          WRITE (*, '(a, 12x, a, 6x, F11.3, 1x, a)') 'PDAF', 'prepare observations (6):', pdaf_time_tot(6), 's'
          WRITE (*, '(a, 12x, a, 8x, F11.3, 1x, a)') 'PDAF', 'compute mean state (9):', pdaf_time_tot(9), 's'
          WRITE (*, '(a, 12x, a, 7x, F11.3, 1x, a)') 'PDAF', 'global preparations (7):', pdaf_time_tot(7), 's'
          WRITE (*, '(a, 12x, a, 7x, F11.3, 1x, a)') 'PDAF', 'local analysis loop (8):', pdaf_time_tot(8), 's'
          WRITE (*, '(a, 14x, a, 10x, F11.3, 1x, a)') 'PDAF', 'global to local (10):', pdaf_time_tot(10), 's'
          WRITE (*, '(a, 14x, a, 4x, F11.3, 1x, a)') 'PDAF', 'localize observations (11):', pdaf_time_tot(11), 's'
          WRITE (*, '(a, 14x, a, 11x, F11.3, 1x, a)') 'PDAF', 'local analysis (12):', pdaf_time_tot(12), 's'
          WRITE (*, '(a, 14x, a, 10x, F11.3, 1x, a)') 'PDAF', 'local to global (14):', pdaf_time_tot(14), 's'
          IF (dim_lag >0) &
               WRITE (*, '(a, 14x, a, 8x, F11.3, 1x, a)') 'PDAF', 'perform smoothing (15):', pdaf_time_tot(15), 's'


          ! Generic part B
          WRITE (*, '(a, 12x, a, 13x, F11.3, 1x, a)') 'PDAF', 'Prepoststep (5):', pdaf_time_tot(5), 's'
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

  END SUBROUTINE PDAF_lknetf_memtime

END MODULE PDAF_LKNETF
