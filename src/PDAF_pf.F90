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
!> Module for PF holding shared parameters and some helper routines
!!
!! This module declares the parameters that are used in PF. 
!! Parameters that are specific for PF are declared while some
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
MODULE PDAF_PF

  USE PDAF_mod_filter, &
       ONLY: incremental, debug, dim_lag

  IMPLICIT NONE

! *** Integer parameters ***

  INTEGER :: type_resample=1 !< Resampling type 
                             !< (1) probabilistic resamping
                             !< (2) stochastic universal resampling
                             !< (3) residual resampling
  INTEGER :: type_forget=0   !< Type of forgetting factor
                             !< (0) inflate forecast ensemble
                             !< (2) inflate analysis ensemble
  INTEGER :: type_noise=0    !< Type of perturbing noise in PF
                             !< (0) no noise added
                             !< (1) constant variance
                             !< (2) amplitude relative to ensemble std.
  INTEGER :: type_winf=0     !< Type of weights inflation for PF
                             !< (0): none; (1) inflate for N_eff/N > limit_winf

! *** Real parameters ***
  REAL :: forget=1.0         !< Forgetting factor
  REAL :: noise_amp=0.0      !< Amplitude of noise perturbing particles
  REAL :: limit_winf=0.0     !< Limit to weights inflation

!-------------------------------------------------------------------------------
  
CONTAINS

!>  PDAF-internal initialization of PF
!!
!! Initialization of PF within PDAF. Performed are:
!! * initialize filter-specific parameters
!! * print screen information on filter configuration.
!!
!!  !  This is a core routine of PDAF and   !
!!  !   should not be changed by the user   !
!!
!! __Revision history:__
!! * 2019-05 - Lars Nerger - Initial code based on code for NETF
!! *  Later revisions - see repository log
!!
  SUBROUTINE PDAF_PF_init(subtype, param_int, dim_pint, param_real, dim_preal, &
       ensemblefilter, fixedbasis, verbose, outflag)

    USE PDAF_mod_filter, &
         ONLY: localfilter, dim_lag, globalobs
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
    observe_ens = .false.
    dim_lag = 0

! Parse provided parameters
    flagsum = 0
    DO i=3, dim_pint
       CALL PDAF_pf_set_iparam(i, param_int(i), outflag)
       flagsum = flagsum+outflag
    END DO
    DO i=1, dim_preal
       CALL PDAF_pf_set_rparam(i, param_real(i), outflag)
       flagsum = flagsum+outflag
    END DO

    ! Define whether filter is mode-based or ensemble-based
    ensemblefilter = .TRUE.

    ! Define whether filter is a domain-local filter
    localfilter = 0

    ! Define that filter needs global observations (used for OMI)
    globalobs = 1

    ! Initialize flag for fixed-basis filters
    fixedbasis = .FALSE.


! *********************
! *** Screen output ***
! *********************

    writeout: IF (verbose == 1) THEN

       WRITE(*, '(/a)') 'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       WRITE(*, '(a)')  'PDAF    +++           Particle Filter with resampling             +++'
       WRITE(*, '(a)')  'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

       IF (flagsum == 0) THEN

          ! *** General output ***
          WRITE (*, '(/a, 4x, a)') 'PDAF', 'PF configuration'
          WRITE (*, '(a, 11x, a, i1)') 'PDAF', 'filter sub-type = ', subtype
          IF (subtype == 0) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> PF '
          ELSE
             WRITE (*, '(/5x, a/)') 'PDAF-ERROR(2): No valid sub type!'
             outflag = 3
          END IF
          IF (type_resample == 1) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> Resample using probabilistic resampling'
          ELSE IF (type_resample == 2) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> Resample using stochastic universal resampling'
          ELSE IF (type_resample == 3) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> Resample using residual resampling'
          END IF
          IF (type_forget == 0) THEN
             WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> prior inflation, forgetting factor:', forget
          ELSEIF (type_forget == 2) THEN
             WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> posterior inflation, forgetting factor:', forget
          ENDIF
          IF (type_noise == 0) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> no noise added to particles'
          ELSEIF (type_noise == 1) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> use noise of constant variance'
          ELSEIF (type_noise == 2) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> use noise with amplitude relative to ensemble standard deviation'
          END IF
          WRITE (*, '(a, 12x, a, f8.3)') 'PDAF', '--> noise amplitude/factor', noise_amp
          IF (type_winf == 1) THEN
             WRITE (*, '(a, 12x, a, f8.3)') 'PDAF', '--> inflate particle weights so that N_eff/N> ', limit_winf
          END IF
          IF (observe_ens) &
               WRITE (*, '(a, 12x, a, 1x, l)') 'PDAF', '--> observe_ens:', observe_ens
       ELSE
          WRITE (*, '(/5x, a/)') 'PDAF-ERROR: Invalid parameter setting - check prior output!'
       END IF

    END IF writeout

  END SUBROUTINE PDAF_PF_init


!-------------------------------------------------------------------------------
!> Perform allocation of arrays for PF.
!!
!! __Revision history:__
!! * 2019-05 - Lars Nerger - Initial code based on NETF
!! * 2025-02 - Lars Nerger - Restructuring introducing generic PDAF_alloc
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAF_pf_alloc(outflag)

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

    CALL PDAF_alloc(dim_p, dim_ens, dim_ens_l, 1, dim_bias_p, &
         dim_lag, 0, 0, outflag)

  END SUBROUTINE PDAF_pf_alloc


!-------------------------------------------------------------------------------
!> Set integer parameter specific for PF
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
  SUBROUTINE PDAF_pf_set_iparam(id, value, flag)

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
       type_resample = value
       IF (type_resample < 1 .OR. type_resample > 3) THEN
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(8): Invalid value for resampling type - param_int(3)!'
          flag = 8
       END IF
    CASE(4)
       type_noise = value
       IF (type_noise < 0 .OR. type_noise > 2) THEN
          WRITE (*, '(/5x, a/)') 'PDAF-ERROR(8): Invalid choise of noise type - param_int(4)!'
          flag = 8
       END IF
    CASE(5)
       type_forget = value
       IF (.NOT.(type_forget == 0 .OR. type_forget == 2)) THEN
          WRITE (*, '(/5x, a/)') 'PDAF-ERROR(8): Invalid type of forgetting factor - param_int(5)!'
          flag = 8
       END IF
    CASE(6)
       type_winf = value
       IF (.NOT.(type_winf == 0 .OR. type_winf == 1)) THEN
          WRITE (*, '(/5x, a/)') 'PDAF-ERROR(8): Invalid type weights inflation - param_int(6)!'
          flag = 8
       END IF
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
               'PDAF-ERROR(8): Invalid setting type_obs_init - param_int(9)!'
          flag = 8
       END IF
    CASE DEFAULT
       WRITE (*,'(/5x, a, i3/)') &
            'PDAF-WARNING: Invalid integer parameter index', id
    END SELECT

  END SUBROUTINE PDAF_pf_set_iparam


!-------------------------------------------------------------------------------
!> Set real parameter specific for PF
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
  SUBROUTINE PDAF_pf_set_rparam(id, value, flag)

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
       noise_amp = value
       IF(noise_amp<0.0) THEN
          WRITE (*,'(/5x,a/)') &
               'PDAF-ERROR(8): Noise amplitude cannot be negative - param_real(1)!'
          flag = 8
       END IF
    CASE(2)
       forget = value
       IF (forget <= 0.0) THEN
          WRITE (*,'(/5x,a/)') &
               'PDAF-ERROR(7): Invalid value of forgetting factor - param_real(2)!'
          flag = 7
       END IF
    CASE(3)
       limit_winf = value
       IF (limit_winf < 0.0) THEN
          WRITE (*,'(/5x,a/)') &
               'PDAF-ERROR(8): Invalid limit for weight inflation - param_real(3)!'
          flag = 8
       END IF
    CASE DEFAULT
       WRITE (*,'(/5x, a, i3/)') &
            'PDAF-WARNING: Invalid real parameter index', id
    END SELECT

  END SUBROUTINE PDAF_pf_set_rparam

END MODULE PDAF_PF
