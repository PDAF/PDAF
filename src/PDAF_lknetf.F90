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
!! *  Other revisions - see repository log
!!
MODULE PDAF_LKNETF

  USE PDAF_mod_filter, &
       ONLY: localfilter, debug, dim_lag, &
       member_save, skewness, kurtosis

  IMPLICIT NONE

! *** Integer parameters ***
  INTEGER :: type_hyb=1    !< Type of hybrid weight: 
                           !< (0) fixed
                           !< (1) gamma_lin (default)
                           !< (2) gamma_alpha
                           !< (3) gamma_ska
                           !< (4) gamma_sklin
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
!! *  Other revisions - see repository log
!!
  SUBROUTINE PDAF_lknetf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
       ensemblefilter, fixedbasis, verbose, outflag)

    USE PDAF_mod_filter, &
         ONLY: dim_lag, dim_ens
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
       WRITE(*, '(/a, 4x, a)') 'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '+++  Local Hybrid Kalman-Nonlinear Ensemble Transform Filter  +++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '+++                                                           +++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '+++                Domain-localized LKNETF by                 +++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '+++ L. Nerger, QJRMS, 148 (2022) 620-640, doi:10.1002/qj.4221 +++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    END IF writeout


! ****************************
! *** INITIALIZE VARIABLES ***
! ****************************

    ! Set parameter default values
    ! (Other defaults are set in the module)
    observe_ens = .false.
    dim_lag = 0

    ! hyb_k needs to be set at runtime
    hyb_k = dim_ens

    ! Parse provided parameters
    DO i=3, dim_pint
       CALL PDAF_lknetf_set_iparam(i, param_int(i), outflag)
    END DO
    DO i=1, dim_preal
       CALL PDAF_lknetf_set_rparam(i, param_real(i), outflag)
    END DO


    ! Define whether filter is mode-based or ensemble-based
    ensemblefilter = .TRUE.

    ! Define whether filter is domain localized
    localfilter = 1

    ! Initialize flag for fixed-basis filters
    fixedbasis = .FALSE.


! *********************
! *** Check subtype ***
! *********************

    IF (subtype<0 .OR. subtype>4) THEN
       WRITE (*, '(/5x, a/)') 'PDAF-ERROR(3): No valid subtype!'
       outflag = 3
    END IF

  END SUBROUTINE PDAF_lknetf_init


!-------------------------------------------------------------------------------
!> Perform allocation of arrays for LKNETF.
!!
!! __Revision history:__
!! * 2017-08 - Lars Nerger - Initial code based in code for LETKF
!! * 2025-02 - Lars Nerger - Restructuring introducing generic PDAF_alloc
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_lknetf_alloc(outflag)

    USE PDAF_mod_filter, &
         ONLY: dim_ens, dim_p, dim_bias_p
    USE PDAF_mod_filtermpi, &
         ONLY: dim_ens_l
    USE PDAF_utils, &
         ONLY: PDAF_alloc

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout):: outflag      !< Status flag


! ******************************
! *** Allocate filter fields ***
! ******************************

    CALL PDAF_alloc(dim_p, dim_ens, dim_ens_l, dim_ens, dim_bias_p, &
         dim_lag, 0, outflag)

  END SUBROUTINE PDAF_lknetf_alloc


!-------------------------------------------------------------------------------
!>  Print information on configuration of LKNETF
!!
!!  !  This is a core routine of PDAF and   !
!!  !   should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code by splitting from PDAF_lknetf_init
!! *  Other revisions - see repository log
!!
  SUBROUTINE PDAF_lknetf_config(subtype, verbose)

    USE PDAF_mod_filter, &
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

       WRITE (*, '(/a, 4x, a)') 'PDAF', 'LKNETF configuration'
       WRITE (*, '(a, 10x, a, i5)') 'PDAF', 'ensemble size:', dim_ens
       WRITE (*, '(a, 10x, a, i1)') 'PDAF', 'filter sub-type= ', subtype
       IF (subtype == 0) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> (HNK) 2-step LKNETF: NETF before LETKF'
       ELSE IF (subtype == 1) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> (HKN) 2-step LKNETF: LETKF before NETF'
       ELSE IF (subtype == 4) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> (HSync) LKNETF synchronous'
       END IF
       IF (dim_lag > 0) &
            WRITE (*, '(a, 12x, a, i6)') 'PDAF', '--> Apply smoother up to lag:',dim_lag
       WRITE(*, '(a, 10x, a, i3)') &
            'PDAF', 'param_int(5) type_forget=', type_forget
       IF (type_forget == 0) THEN
          WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> prior inflation (default), forgetting factor:', forget
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
       WRITE(*, '(a, 10x, a, i3)') &
            'PDAF', 'param_int(6) type_trans=', type_trans
       IF (type_trans == 0) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> Transform ensemble including product with random matrix (default)'
       ELSE IF (type_trans == 1) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> Deterministic symmetric ensemble transformation'
       END IF
       WRITE(*, '(a, 10x, a, i3)') &
            'PDAF', 'param_int(7) type_hyb=', type_hyb
       IF (type_hyb == 0) THEN
          WRITE(*, '(a, 12x, a, f8.3)') 'PDAF', '--> use gamma_fix: fixed hybrid weight', hyb_g
       ELSEIF (type_hyb == 1) THEN
          WRITE(*, '(a, 12x, a, f8.3)') 'PDAF', '--> use gamma_lin (default): (1 - N_eff/N_e)*', hyb_g
       ELSEIF (type_hyb == 2) THEN
          WRITE(*, '(a, 12x, a, f8.3)') 'PDAF', '--> use gamma_alpha: hybrid weight from N_eff/N>=', hyb_g
       ELSEIF (type_hyb == 3) THEN
          WRITE(*, '(a, 12x, a, f8.3, a, f8.3)') &
               'PDAF', '--> use gamma_ska: 1 - min(s,k)/sqrt(', hyb_k, ') with N_eff/N>=', hyb_g
       ELSEIF (type_hyb == 4) THEN
          WRITE(*, '(a, 12x, a, f8.3, a, f8.3)') &
               'PDAF', '--> use gamma_sklin: 1 - min(s,k)/sqrt(', hyb_k, ') >= 1-N_eff/N>=', hyb_g
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
       WRITE(*, '(a, 10x, a, f10.3)') &
            'PDAF', 'param_real(1) forget=', forget
       WRITE(*, '(a, 10x, a, f10.3)') &
            'PDAF', 'param_real(2) hybrid weight input (gamma): hyb_g=', hyb_g
       WRITE(*, '(a, 10x, a, f10.3)') &
            'PDAF', 'param_real(3) hybrid norm (kappa): hyb_k=', hyb_k

    END IF writeout

  END SUBROUTINE PDAF_lknetf_config


!-------------------------------------------------------------------------------
!> Set integer parameter specific for LKNETF
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Other revisions - see repository log
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
!! *  Other revisions - see repository log
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
!! *  Other revisions - see repository log
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
!! * Other revisions - see repository log
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

!-------------------------------------------------------------------------------
!> Get tempering weight according to N_eff
!!
!! Routine to compute an adaptive tempering factor alpha
!! according to the effective sample size.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2019-08 - Lars Nerger
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_lknetf_alpha_neff(dim_ens, weights, hlimit, alpha)

    USE PDAF_diag, ONLY: PDAF_diag_effsample

    IMPLICIT NONE

! *** Arguments ***
!  Variable naming scheme:
!    suffix _p: Denotes a full variable on the PE-local domain
!    suffix _l: Denotes a local variable on the current analysis domain
    INTEGER, INTENT(in) :: dim_ens        !< Size of ensemble 
    REAL, INTENT(in) :: weights(dim_ens)  !< Weights
    REAL, INTENT(in) :: hlimit            !< Minimum of n_eff / N
    REAL, INTENT(inout) :: alpha          !< hybrid weight  

! *** Local variables ***
    INTEGER :: i
    REAL, ALLOCATABLE :: locw(:)
    REAL, ALLOCATABLE :: hweights(:)
    REAL :: a_step, tot_weight
    REAL :: nhlim, n_eff


! *****************************************
! *** Iteratively compute alpha so that ***
! ***      N_eff/dim_ens > hlimit       ***
! *****************************************

    ALLOCATE(locw(dim_ens))
    ALLOCATE(hweights(dim_ens))

    ! Get logarithm of weights
    DO i=1, dim_ens
       locw(i) = LOG(weights(i))
    END DO

    ! Initialize iterations
    alpha = 0.0
    a_step = 0.05
    nhlim = ABS(hlimit * REAL(dim_ens))

    aloop: DO

       ! scale 
       DO i = 1, dim_ens
          hweights(i) = EXP(locw(i) * (1-alpha))
       END DO
  
       ! Normalize weights
       tot_weight = 0.0
       DO i = 1, dim_ens
          tot_weight = tot_weight + hweights(i)
       END DO
       IF (tot_weight /= 0.0) THEN
          hweights = hweights / tot_weight

          ! Compute effective ensemble size
          CALL PDAF_diag_effsample(dim_ens, hweights, n_eff)

          IF (REAL(n_eff) >= nhlim) EXIT aloop
       END IF

       IF (alpha>=1.0) THEN
          alpha = 1.0
          EXIT aloop
       END IF

       alpha = alpha + a_step

    END DO aloop


! ********************
! *** Finishing up ***
! ********************

    DEALLOCATE(hweights, locw)

  END SUBROUTINE PDAF_lknetf_alpha_neff


!-------------------------------------------------------------------------------
!> Compute hybrid weight gamma for hybrid LKNETF
!!
!! This routine computes the hybrid weight gamma for the
!! two-step variants HNK and HKN. 
!! For this, it first computes the particle weights and
!! then it calls PDAF_lknetf_set_gamma to get gamma
!! according to the chosen hybridication weight type. 
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2018-01 - Lars Nerger - 
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_lknetf_compute_gamma(domain_p, step, dim_obs_l, dim_ens, &
       HX_l, HXbar_l, obs_l, type_hyb, hyb_g, hyb_k, &
       gamma, n_eff_out, skew_mabs, kurt_mabs, &
       U_likelihood_l, screen, flag)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    USE PDAF_timer, &
         ONLY: PDAF_timeit
    USE PDAF_memcounting, &
         ONLY: PDAF_memcount
    USE PDAF_mod_filter, &
         ONLY: obs_member, debug

    IMPLICIT NONE

! *** Arguments ***
!  Variable naming scheme:
!    suffix _p: Denotes a full variable on the PE-local domain
!    suffix _l: Denotes a local variable on the current analysis domain
    INTEGER, INTENT(in) :: domain_p              !< Current local analysis domain
    INTEGER, INTENT(in) :: step                  !< Current time step
    INTEGER, INTENT(in) :: dim_obs_l             !< Size of obs. vector on local ana. domain
    INTEGER, INTENT(in) :: dim_ens               !< Size of ensemble 
    REAL, INTENT(in) :: HX_l(dim_obs_l, dim_ens) !< local observed state ens.
    REAL, INTENT(in) :: HXbar_l(dim_obs_l)       !< local mean observed ensemble
    REAL, INTENT(in) :: obs_l(dim_obs_l)         !< Local observation vector
    INTEGER, INTENT(in) :: type_hyb              !< Type of hybrid weight
    REAL, INTENT(in) :: hyb_g                    !< Prescribed hybrid weight for state transformation
    REAL, INTENT(in) :: hyb_k                    !< Hybrid weight for covariance transformation
    REAL, INTENT(inout) :: gamma(1)              !< Hybrid weight for state transformation
    REAL, INTENT(inout) :: n_eff_out(1)          !< Effective ensemble size
    REAL, INTENT(inout) :: skew_mabs(1)          !< Mean absolute skewness
    REAL, INTENT(inout) :: kurt_mabs(1)          !< Mean absolute kurtosis
    INTEGER, INTENT(in) :: screen                !< Verbosity flag
    INTEGER, INTENT(inout) :: flag               !< Status flag

! *** External subroutines ***
!  (PDAF-internal names, real names are defined in the call to PDAF)
    EXTERNAL :: U_likelihood_l         !< Compute observation likelihood for an ensemble member

! *** local variables ***
    INTEGER :: i, member               ! Counters
    INTEGER, SAVE :: allocflag=0       ! Flag whether first time allocation is done
    REAL :: weight                     ! Ensemble weight (likelihood)
    REAL, ALLOCATABLE :: resid_i(:)    ! Residual
    REAL, ALLOCATABLE :: weights(:)    ! weight vector
    REAL :: total_weight               ! Sum of weight

!$OMP THREADPRIVATE(allocflag)


! **********************************************
! *** Compute particle weights as likelihood ***
! **********************************************

    ! Allocate weight vector
    ALLOCATE(weights(dim_ens))
    IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens)

    ! Allocate temporary array for obs-ens_i
    ALLOCATE(resid_i(dim_obs_l))
    IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l)

    CALL PDAF_timeit(22, 'new')
    ! Get residual as difference of observation and observed state for 
    ! each ensemble member only on domains where observations are availible

    IF (debug>0) &
         WRITE (*,*) '++ PDAF-debug: ', debug, &
         'PDAF_lknetf_compute_gamma -- call g2l_obs and likelihood_l', dim_ens, 'times'

    CALC_w: DO member = 1,dim_ens

       ! Store member index to make it accessible with PDAF_get_obsmemberid
       obs_member = member
     
       ! Calculate local residual  
       resid_i = obs_l - HX_l(:,member)

       IF (debug>0) THEN
          WRITE (*,*) '++ PDAF-debug: ', debug, &
               'PDAF_lknetf_compute_gamma -- member', member
          WRITE (*,*) '++ PDAF-debug PDAF_lknetf_compute_gamma:', debug, '  innovation d_l', resid_i
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lknetf_compute_gamma -- call likelihood_l'
       end IF

       ! Compute likelihood
       CALL PDAF_timeit(49, 'new')
       CALL U_likelihood_l(domain_p, step, dim_obs_l, obs_l, resid_i, weight)
       CALL PDAF_timeit(49, 'old')
       weights(member) = weight

    END DO CALC_w

    IF (debug>0) &
         WRITE (*,*) '++ PDAF-debug PDAF_lknetf_compute_gamma:', debug, '  raw non-hybrid weights', weights

    ! Normalize weights
    total_weight = 0.0
    DO i = 1, dim_ens
       total_weight = total_weight + weights(i)
    END DO

    IF (total_weight /= 0.0) THEN
       weights = weights / total_weight

       IF (debug>0) &
            WRITE (*,*) '++ PDAF-debug PDAF_lknetf_compute_gamma:', debug, '  normalized non-hybrid weights', weights
    ELSE
       ! ERROR: weights are zero
       WRITE(*,'(/5x,a/)') 'WARNING: Zero weights - reset to 1/dim_ens'
       weights = 1.0 / REAL(dim_ens)
    END IF

    CALL PDAF_timeit(22, 'old')


! *******************************
! *** Set hybrid weight gamma ***
! *******************************

    CALL PDAF_lknetf_set_gamma(domain_p, dim_obs_l, dim_ens, &
         HX_l, HXbar_l, weights, type_hyb, hyb_g, hyb_k, &
         gamma, n_eff_out, skew_mabs, kurt_mabs, &
         screen, flag)


! ********************
! *** Finishing up ***
! ********************

    DEALLOCATE(weights, resid_i)

    IF (allocflag == 0) allocflag = 1

  END SUBROUTINE PDAF_lknetf_compute_gamma

!-------------------------------------------------------------------------------
!> Compute hybrid weight gamma for hybrid LKNETF
!!
!! Compute hybrid weight for LKNETF
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2018-01 - Lars Nerger - 
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_lknetf_set_gamma(domain_p, dim_obs_l, dim_ens, &
       HX_l, HXbar_l, weights, type_hyb, hyb_g, hyb_k, &
       gamma, n_eff_out, maSkew, maKurt, &
       screen, flag)

! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

    USE PDAF_timer, &
         ONLY: PDAF_timeit
    USE PDAF_memcounting, &
         ONLY: PDAF_memcount
    USE PDAF_mod_filtermpi, &
         ONLY: mype
    USE PDAF_mod_filter, &
         ONLY: debug
    USE PDAF_diag, &
         ONLY: PDAF_diag_effsample, PDAF_diag_ensstats
#if defined (_OPENMP)
    USE omp_lib, &
         ONLY: omp_get_num_threads, omp_get_thread_num
#endif

    IMPLICIT NONE

! *** Arguments ***
!  Variable naming scheme:
!    suffix _p: Denotes a full variable on the PE-local domain
!    suffix _l: Denotes a local variable on the current analysis domain
    INTEGER, INTENT(in) :: domain_p               !< Current local analysis domain
    INTEGER, INTENT(in) :: dim_obs_l              !< Size of obs. vector on local ana. domain
    INTEGER, INTENT(in) :: dim_ens                !< Size of ensemble 
    REAL, INTENT(in) :: HX_l(dim_obs_l, dim_ens)  !< local observed state ens.
    REAL, INTENT(in) :: HXbar_l(dim_obs_l)        !< local mean observed ensemble
    REAL, INTENT(in) :: weights(dim_ens)          !< Weight vector
    INTEGER, INTENT(in) :: type_hyb               !< Type of hybrid weight
    REAL, INTENT(in) :: hyb_g                     !< Prescribed hybrid weight for state transformation
    REAL, INTENT(in) :: hyb_k                     !< Scale factor kappa (for type_hyb 3 and 4)
    REAL, INTENT(inout) :: gamma(1)               !< Hybrid weight for state transformation
    REAL, INTENT(inout) :: n_eff_out(1)           !< Effective ensemble size
    REAL, INTENT(inout) :: maSkew(1)              !< Mean absolute skewness
    REAL, INTENT(inout) :: maKurt(1)              !< Mean absolute kurtosis
    INTEGER, INTENT(in) :: screen                 !< Verbosity flag
    INTEGER, INTENT(inout) :: flag                !< Status flag

! *** local variables ***
    INTEGER :: i                              ! Counter
    INTEGER, SAVE :: allocflag=0              ! Flag whether first time allocation is done
    LOGICAL :: screenout = .TRUE.             ! Whether to print information to stdout
    REAL :: n_eff                             ! Effective sample size
    INTEGER , SAVE :: lastdomain=-1           !save domain index
    INTEGER, SAVE :: mythread, nthreads       ! Thread variables for OpenMP
    REAL, PARAMETER :: pi=3.14159265358979    ! Pi
    REAL :: skew, kurt                        ! Skewness and kurtosis of observed ensemble
    REAL :: kurt_limit                        ! Asymptotic value of kurtosis
    REAL :: gamma_kurt                        ! Gamma from kurtosis
    REAL :: gamma_skew                        ! Gamma from Skewness
    REAL :: gamma_Neff                        ! Gamma from effective sample size
    REAL :: gamma_stat                        ! Gamma from combining kurtosis and skewness

!$OMP THREADPRIVATE(mythread, nthreads, lastdomain, allocflag, screenout)


! *******************
! *** Preparation ***
! *******************

#if defined (_OPENMP)
    nthreads = omp_get_num_threads()
    mythread = omp_get_thread_num()
#else
    nthreads = 1
    mythread = 0
#endif

    ! Control screen output
    IF (lastdomain<domain_p .AND. lastdomain>-1) THEN
       screenout = .FALSE.
    ELSE
       screenout = .TRUE.

       ! In case of OpenMP, let only thread 0 write output to the screen
       IF (mythread>0) screenout = .FALSE.

       ! Output, only in case of OpenMP parallelization
    END IF

    ! Dummy initialization to prevent compiler warning
    gamma_Neff = hyb_g


    ! *********************
    ! *** Screen output ***
    ! *********************

    IF (mype == 0 .AND. screen > 0 .AND. screenout) THEN
       IF (type_hyb==0) THEN
          WRITE (*, '(a, 5x, a, f8.3)') 'PDAF', 'Set hybrid weight in LETKF to ', hyb_g
       ELSE
          WRITE (*, '(a, 5x, a)') 'PDAF', 'Compute adaptive hybrid weight'
       END IF

       ! First four methods are discussed in Nerger (2022)
       IF (type_hyb==1) THEN
          WRITE (*, '(a, 5x, a, f8.3)') 'PDAF', '--- gamma_lin: (1 - N_eff/N)*scale, scale=', hyb_g
       ELSEIF (type_hyb==2) THEN
          WRITE (*, '(a, 5x, a, f8.3)') 'PDAF', '--- gamma_alpha: hybrid weight from N_eff/N>=', hyb_g
       ELSEIF (type_hyb==3) THEN
          WRITE (*, '(a, 5x, a, f8.3, a, f8.2)') 'PDAF', &
               '--- gamma_ska: 1 - min(s,k)/sqrt(kappa) with N_eff/N>=', hyb_g, ', kappa=', hyb_k
       ELSEIF (type_hyb==4) THEN
          WRITE (*, '(a, 5x, a, a, f8.3)') 'PDAF', &
               '--- gamma_sklin: 1 - min(s,k)/sqrt(kappa) >= 1-N_eff/N', ', kappa', hyb_k

       ! Additional methods
       ELSEIF (type_hyb==23) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- quadratic dependence on 1 - N_eff/N'
       ELSEIF (type_hyb==24) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- square-root dependence on 1 - N_eff/N; 0 for N_eff=N'
       ELSEIF (type_hyb==25) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- square-root dependence on 1 - N_eff/N; 1 for N_eff=N'
       ELSEIF (type_hyb==26) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- sine dependence on N_eff/N with minimum constraint'
       ELSEIF (type_hyb==27) THEN
          WRITE (*, '(a, 5x, a,f8.3)') 'PDAF', '--- square-root dependence on N_eff/N, minimum limit', hyb_g
       ELSEIF (type_hyb==28) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- linear dependence 1 - 0.5 N_eff/N'
       ELSEIF (type_hyb==29) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- quadratic dependence 1 - 0.5 (N_eff/N)^2'
       ELSEIF (type_hyb==10) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - (1-limit)* skewness/sqrt(N)'
       ELSEIF (type_hyb==11) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - skewness/sqrt(N) with minimum limit'
       ELSEIF (type_hyb==12) THEN
          WRITE (*, '(a, 5x, a, f12.3)') 'PDAF', '--- dependence 1 - skewness/sqrt(N) with min. from N_eff/N>=', hyb_g
       ELSEIF (type_hyb==13) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - skewness/sqrt(N) with min. from 1-N_eff/N'
       ELSEIF (type_hyb==14) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- linear dependence on (N_eff-1)/(N-1)'
       ELSEIF (type_hyb==15) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- linear dependence on (N_eff-1)/(N-1), limit 0.95'
       ELSEIF (type_hyb==110) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - (1-limit)* kurtosis/sqrt(N)'
       ELSEIF (type_hyb==111) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - kurtosis/sqrt(N) with minimum limit'
       ELSEIF (type_hyb==112) THEN
          WRITE (*, '(a, 5x, a, f12.3)') 'PDAF', '--- dependence 1 - kurtosis/sqrt(N) with min. from N_eff/N>=', hyb_g
       ELSEIF (type_hyb==113) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - kurtosis/sqrt(N) with min. from 1-N_eff/N'
       ELSEIF (type_hyb==212) THEN
          WRITE (*, '(a, 5x, a, f12.3)') 'PDAF', '--- dependence 1 - kurtnorm/sqrt(N) with min. from N_eff/N>=', hyb_g
       ELSEIF (type_hyb==213) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - kurtnorm/sqrt(N) with min. from 1-N_eff/N'
       ELSEIF (type_hyb==312) THEN
          WRITE (*, '(a, 5x, a, f12.3)') 'PDAF', '--- dependence 1 - ensstat/sqrt(N) with min. from N_eff/N>=', hyb_g
       ELSEIF (type_hyb==313) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - ensstat/sqrt(N) with min. from 1-N_eff/N'
       ELSEIF (type_hyb==413) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- min of 1 - ensstat/sqrt(N) and 1-N_eff/N'
       ELSEIF (type_hyb==513) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - ensstat/sqrt(N) with min. from 1-(N_eff-1)/(N-1)'
       ELSEIF (type_hyb==613) THEN
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- dependence 1 - ensstat/sqrt(N) with min. from 1-(N_eff-1)/(N-1) B'
       END IF
    END IF


    ! **********************************************
    ! *** Compute particle weights as likelihood ***
    ! **********************************************

    CALL PDAF_timeit(54, 'new')
    ! Get residual as difference of observation and observed state for 
    ! each ensemble member only on domains where observations are availible

    ! Compute adaptive hybrid weights
    IF (type_hyb==2 .OR. type_hyb==3 .OR. type_hyb==12 .OR. type_hyb==112 &
         .OR. type_hyb==212 .OR. type_hyb==312) THEN
       IF (debug>0) &
            WRITE (*,*) '++ PDAF-debug: ', debug, &
            'PDAF_lknetf_compute_gamma -- Determine gamma for N_eff/N>=limit iteratively'
       CALL PDAF_lknetf_alpha_neff(dim_ens, weights, hyb_g, gamma(1))
       gamma_Neff = gamma(1)

       IF (debug>0 .AND. type_hyb/=2) &
            WRITE (*,*) '++ PDAF-debug PDAF_lknetf_compute_gamma:', debug, '  gamma for N_eff/N>=limit', gamma_Neff
    END IF

    CALL PDAF_timeit(54, 'old')


    ! Compute effective ensemble size
    CALL PDAF_diag_effsample(dim_ens, weights, n_eff)
    n_eff_out(1) = n_eff

    ! Compute ensemble statistics for observed ensemble
    DO i = 1, dim_obs_l
       CALL PDAF_diag_ensstats(dim_obs_l, dim_ens, i, &
            HXbar_l, HX_l, skew, kurt, flag)
       maSkew = maSkew + ABS(skew)
       maKurt = maKurt + ABS(kurt)
    END DO
    maSkew = maSkew / dim_obs_l
    maKurt = maKurt / dim_obs_l

    IF (debug>0) &
         WRITE (*,*) '++ PDAF-debug PDAF_lknetf_compute_gamma:', debug, '  MAS, MAK', maSkew, maKurt


! *********************************************
! *** Compute adaptive hybridization weight ***
! *********************************************

    IF (type_hyb>0) THEN
       ! Compute adaptive hybrid weight for state update

       ! First four methods are discussed in Nerger (2022)
       IF (type_hyb==1) THEN
          gamma(1) = (1.0 - n_eff/REAL(dim_ens))*(hyb_g)
       ELSEIF (type_hyb==2) THEN
          gamma(1) = gamma_Neff
       ELSEIF (type_hyb==3) THEN
          gamma_skew = 1.0 - maSkew(1)/SQRT(hyb_k)
          IF (gamma_skew<0.0) gamma_skew = 0.0
          gamma_kurt = 1.0 - maKurt(1)/hyb_k
          IF (gamma_kurt<0.0) gamma_kurt = 0.0
          gamma_stat = MIN(gamma_skew, gamma_kurt)
          gamma(1) = MAX(gamma_stat, gamma_Neff)
       ELSEIF (type_hyb==4) THEN
          gamma_skew = 1.0 - maSkew(1)/SQRT(hyb_k)
          IF (gamma_skew<0.0) gamma_skew = 0.0
          gamma_kurt = 1.0 - maKurt(1)/hyb_k
          IF (gamma_kurt<0.0) gamma_kurt = 0.0
          gamma_stat = MIN(gamma_skew, gamma_kurt)
          gamma_Neff = 1.0 - n_eff/REAL(dim_ens)
          gamma(1) = MAX(gamma_stat, gamma_Neff)

       ! Additional methods - experimental
       ELSEIF (type_hyb==23) THEN
          gamma(1) = ((1.0 - n_eff/REAL(dim_ens))**2)*(hyb_g)
       ELSEIF (type_hyb==24 .OR. type_hyb==25) THEN
          IF (1.0 - n_eff/REAL(dim_ens) > 0.0) THEN
             gamma(1) = SQRT(1.0 - n_eff/REAL(dim_ens))*(hyb_g)
          ELSE
             IF (type_hyb==24) THEN
                gamma(1) = 0.0
             ELSEIF (type_hyb==25) THEN
                gamma(1) = 1.0
             ENDIF
          ENDIF
       ELSEIF (type_hyb==26) THEN
          gamma(1) = n_eff/REAL(dim_ens)
          gamma(1) = 1.0 - SIN(gamma(1)*pi)*(1.0-hyb_g)
          IF (gamma(1)<hyb_g) gamma(1) = hyb_g
          IF (gamma(1)>1.0) gamma(1) = 1.0
       ELSEIF (type_hyb==27) THEN
          IF (1.0 - n_eff/REAL(dim_ens) > 0.0) THEN
             gamma(1) = SQRT(1.0 - n_eff/REAL(dim_ens))
          ELSE
             gamma(1) = 0.0
          ENDIF
          IF (gamma(1)<hyb_g) gamma(1) = hyb_g
       ELSEIF (type_hyb==28) THEN
          gamma(1) = (1.0 - 0.5*n_eff/REAL(dim_ens))*(hyb_g)
       ELSEIF (type_hyb==29) THEN
          gamma(1) = (1.0 - 0.5*(n_eff/REAL(dim_ens))**2)*(hyb_g)
       ELSEIF (type_hyb==10) THEN
          gamma(1) = 1.0 - maSkew(1)/SQRT(REAL(dim_ens))*(1.0-hyb_g)
          IF (gamma(1) < (hyb_g)) gamma(1) = hyb_g
       ELSEIF (type_hyb==11) THEN
          gamma(1) = 1.0 - maSkew(1)/SQRT(REAL(dim_ens))
          IF (gamma(1) < (hyb_g)) gamma(1) = hyb_g
       ELSEIF (type_hyb==12) THEN
          gamma_skew = 1.0 - maSkew(1)/SQRT(REAL(dim_ens))
          IF (gamma_skew<0.0) gamma_skew = 0.0
          gamma(1) = MAX(gamma_skew, gamma_Neff)
       ELSEIF (type_hyb==13) THEN
          gamma_skew = 1.0 - maSkew(1)/SQRT(REAL(dim_ens))
          IF (gamma_skew<0.0) gamma_skew = 0.0
          gamma_Neff = 1.0 - n_eff/REAL(dim_ens)
          gamma(1) = MAX(gamma_skew, gamma_Neff)
       ELSEIF (type_hyb==14) THEN
          gamma(1) = (1.0 - (n_eff-1.0)/REAL(dim_ens-1))*(hyb_g)
          IF (gamma(1)<0.0) gamma(1) = 0.0
       ELSEIF (type_hyb==15) THEN
          gamma(1) = (1.0 - (n_eff-1.0)/REAL(dim_ens-1))*(hyb_g)
          IF (gamma(1)<0.0) gamma(1) = 0.0
          IF (gamma(1)>0.95) gamma(1)=0.95
       ELSEIF (type_hyb==110) THEN
          kurt_limit = 0.5*(REAL(dim_ens-3))/(REAL(dim_ens-2))*maSkew(1)*maSkew(1) - REAL(dim_ens)/0.5-3.0
          gamma(1) = 1.0 - maKurt(1)/kurt_limit*(1.0-hyb_g)
          IF (gamma(1) < (hyb_g)) gamma(1) = hyb_g
       ELSEIF (type_hyb==111) THEN
          kurt_limit = 0.5*(REAL(dim_ens-3))/(REAL(dim_ens-2))*maSkew(1)*maSkew(1) - REAL(dim_ens)/0.5-3.0
          gamma(1) = 1.0 - maKurt(1)/kurt_limit
          IF (gamma(1) < (hyb_g)) gamma(1) = hyb_g
       ELSEIF (type_hyb==112) THEN
          kurt_limit = 0.5*(REAL(dim_ens-3))/(REAL(dim_ens-2))*maSkew(1)*maSkew(1) - REAL(dim_ens)/0.5-3.0
          gamma_kurt = 1.0 - maKurt(1)/kurt_limit
          IF (gamma_kurt<0.0) gamma_kurt = 0.0
          gamma(1) = MAX(gamma_kurt, gamma_Neff)
       ELSEIF (type_hyb==113) THEN
          kurt_limit = 0.5*(REAL(dim_ens-3))/(REAL(dim_ens-2))*maSkew(1)*maSkew(1) - REAL(dim_ens)/0.5-3.0
          gamma_kurt = 1.0 - maKurt(1)/kurt_limit
          IF (gamma_kurt<0.0) gamma_kurt = 0.0
          gamma_Neff = 1.0 - n_eff/REAL(dim_ens)
          gamma(1) = MAX(gamma_kurt, gamma_Neff)
       ELSEIF (type_hyb==212) THEN
          kurt_limit = REAL(dim_ens)
          gamma_kurt = 1.0 - maKurt(1)/kurt_limit
          IF (gamma_kurt<0.0) gamma_kurt = 0.0
          gamma(1) = MAX(gamma_kurt, gamma_Neff)
       ELSEIF (type_hyb==213) THEN
          kurt_limit = REAL(dim_ens)
          gamma_kurt = 1.0 - maKurt(1)/kurt_limit
          IF (gamma_kurt<0.0) gamma_kurt = 0.0
          gamma_Neff = 1.0 - n_eff/REAL(dim_ens)
          gamma(1) = MAX(gamma_kurt, gamma_Neff)
       ELSEIF (type_hyb==312) THEN
          gamma_skew = 1.0 - maSkew(1)/SQRT(REAL(dim_ens))
          IF (gamma_skew<0.0) gamma_skew = 0.0
          gamma_kurt = 1.0 - maKurt(1)/REAL(dim_ens)
          IF (gamma_kurt<0.0) gamma_kurt = 0.0
          gamma_stat = MIN(gamma_skew, gamma_kurt)
          gamma(1) = MAX(gamma_stat, gamma_Neff)
       ELSEIF (type_hyb==313) THEN
          gamma_skew = 1.0 - maSkew(1)/SQRT(REAL(dim_ens))
          IF (gamma_skew<0.0) gamma_skew = 0.0
          gamma_kurt = 1.0 - maKurt(1)/REAL(dim_ens)
          IF (gamma_kurt<0.0) gamma_kurt = 0.0
          gamma_stat = MIN(gamma_skew, gamma_kurt)
          gamma_Neff = 1.0 - n_eff/REAL(dim_ens)
          gamma(1) = MAX(gamma_stat, gamma_Neff)
       ELSEIF (type_hyb==413) THEN
          gamma_skew = 1.0 - maSkew(1)/SQRT(REAL(dim_ens))
          IF (gamma_skew<0.0) gamma_skew = 0.0
          gamma_kurt = 1.0 - maKurt(1)/REAL(dim_ens)
          IF (gamma_kurt<0.0) gamma_kurt = 0.0
          gamma_stat = MIN(gamma_skew, gamma_kurt)
          gamma_Neff = 1.0 - n_eff/REAL(dim_ens)
          IF (gamma_Neff<0.0) gamma_Neff=0.0
          gamma(1) = MIN(gamma_stat, gamma_Neff)
       ELSEIF (type_hyb==513) THEN
          gamma_skew = 1.0 - maSkew(1)/SQRT(REAL(dim_ens))
          IF (gamma_skew<0.0) gamma_skew = 0.0
          gamma_kurt = 1.0 - maKurt(1)/REAL(dim_ens)
          IF (gamma_kurt<0.0) gamma_kurt = 0.0
          gamma_stat = MIN(gamma_skew, gamma_kurt)
          gamma_Neff = 1.0 - (n_eff-1.0)/REAL(dim_ens-1)
          gamma(1) = MAX(gamma_stat, gamma_Neff)
       ELSEIF (type_hyb==613) THEN
          gamma_skew = 1.0 - maSkew(1)/SQRT(REAL(dim_ens))
          IF (gamma_skew<0.0) gamma_skew = 0.0
          gamma_kurt = 1.0 - maKurt(1)/REAL(dim_ens)
          IF (gamma_kurt<0.0) gamma_kurt = 0.0
          gamma_stat = MIN(gamma_skew, gamma_kurt)
          gamma_Neff = 1.0 - (n_eff-1.0)/REAL(dim_ens-1)
          gamma(1) = MAX(gamma_stat, gamma_Neff)
          IF (gamma(1)>0.95) gamma(1)=0.95
       END IF

    ELSE

       ! fixed hybrid weight
       gamma(1) = hyb_g

    END IF


! ********************
! *** Finishing up ***
! ********************

    IF (allocflag == 0) allocflag = 1

    lastdomain = domain_p

  END SUBROUTINE PDAF_lknetf_set_gamma

!-------------------------------------------------------------------------------
!> reset gamma value of LKNETF
!!
!! This routine resets the hybrid weight value
!! gamma in the LKNETF.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Resision history
!! * 2019-09 - Lars Nerger - initial code
!! * Other revisions - see repository log
!!
SUBROUTINE PDAF_lknetf_reset_gamma(gamma_in)

  IMPLICIT NONE

! *** Argument ***
  REAL, INTENT(in) :: gamma_in     !< Prescribed hybrid weight


! *** Set hybrid weights ***
  hyb_g = gamma_in

END SUBROUTINE PDAF_lknetf_reset_gamma

END MODULE PDAF_LKNETF
