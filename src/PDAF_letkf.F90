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
!> Module for LETKF holding shared parameters and some helper routines
!!
!! This module declares the parameters that are used in LETKF. 
!! Parameters that are specific for LETKF are declared while some
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
MODULE PDAF_letkf

  USE PDAF_mod_filter, &
       ONLY: localfilter, incremental, debug, dim_lag, member_save

  IMPLICIT NONE

! *** Integer parameters ***
  INTEGER :: type_forget=0 !< Type of forgetting factor
                           !< (0): fixed; (1) global adaptive; (2) local adaptive
  INTEGER :: type_trans=0  !< Type of ensemble transformation
                           !< For ESTKF/LETKF:
                           !< (0) use deterministic Omega
                           !< (1) use random orthonormal Omega orthogonal to (1,...,1)^T
                           !< (2) use product of (0) with random orthonomal matrix with
                           !<     eigenvector (1,...,1)^T

! *** Real parameters ***
  REAL    :: forget=1.0    !< Forgetting factor
  REAL    :: forget_l      !< Forgetting factor in local analysis loop

! *** Internal variable ***
  LOGICAL :: inloop=.false. ! Whether the program is in the local analysis loop


!$OMP THREADPRIVATE(forget_l)

!-------------------------------------------------------------------------------
  
CONTAINS

!>  PDAF-internal initialization of LETKF
!!
!! Initialization of LETKF within PDAF. Performed are:
!! * initialize filter-specific parameters
!! * print screen information on filter configuration.
!!
!!  !  This is a core routine of PDAF and   !
!!  !   should not be changed by the user   !
!!
!! __Revision history:__
!! * 2009-07 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
  SUBROUTINE PDAF_letkf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
       ensemblefilter, fixedbasis, verbose, outflag)

    USE PDAF_mod_filter, &
         ONLY: dim_ens, dim_bias_p
    USE PDAFobs, &
         ONLY: observe_ens

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout) :: subtype             !< Sub-type of filter
    INTEGER, INTENT(in) :: dim_pint               !< Number of integer parameters
    INTEGER, INTENT(inout) :: param_int(dim_pint) !< Integer parameter array
    INTEGER, INTENT(in) :: dim_preal              !< Number of real parameters 
    REAL, INTENT(inout) :: param_real(dim_preal)  !< Real parameter array
    LOGICAL, INTENT(out) :: ensemblefilter        !< Is the chosen filter ensemble-based?
    LOGICAL, INTENT(out) :: fixedbasis            !< Does the filter run with fixed error-space basis?
    INTEGER, INTENT(in) :: verbose                !< Control screen output
    INTEGER, INTENT(inout):: outflag              !< Status flag

! *** local variables ***
    INTEGER :: i                ! Counter
    INTEGER :: flagsum          ! Sum of status flags


! ****************************
! *** INITIALIZE VARIABLES ***
! ****************************

    ! Set parameter default values - other defaults are set directly in the module
    incremental = 0
    observe_ens = .true.
    dim_lag = 0

    ! Parse provided parameters
    flagsum = 0
    DO i=3, dim_pint
       CALL PDAF_letkf_set_iparam(i, param_int(i), outflag)
       flagsum = flagsum+outflag
    END DO
    DO i=1, dim_preal
       CALL PDAF_letkf_set_rparam(i, param_real(i), outflag)
       flagsum = flagsum+outflag
    END DO


    ! Define whether filter is mode-based or ensemble-based
    ensemblefilter = .TRUE.

    ! Define whether filter is domain localized
    localfilter = 1

    ! Initialize flag for fixed-basis filters
    IF (subtype == 2 .OR. subtype == 3) THEN
       fixedbasis = .TRUE.
    ELSE
       fixedbasis = .FALSE.
    END IF


! *********************
! *** Screen output ***
! *********************

    writeout: IF (verbose == 1) THEN

       WRITE(*, '(/a, 4x, a)') 'PDAF', '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '+++  Local Ensemble Transform Kalman Filter (LETKF)  +++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '+++                                                  +++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '+++            Domain-localized LETKF by             +++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '+++      Hunt et al., Physica D 230 (2007) 112       +++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

       IF (flagsum== 0 ) THEN

          ! *** General output ***
          WRITE (*, '(/a, 4x, a)') 'PDAF', 'LETKF configuration'
          WRITE (*, '(a, 10x, a, i1)') 'PDAF', 'filter sub-type = ', subtype
          IF (subtype == 0) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> LETKF using T-matrix'
          ELSE IF (subtype == 1) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> LETKF following Hunt et al. (2007)'
          ELSE IF (subtype == 2) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> LETKF with fixed error-space basis'
          ELSE IF (subtype == 3) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> LETKF with fixed state covariance matrix'
          ELSE
             WRITE (*, '(/5x, a/)') 'PDAF-ERROR(3): No valid subtype!'
             outflag = 3
          END IF
          IF (type_trans == 0) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> Deterministic symmetric ensemble transformation'
          ELSE IF (type_trans == 2) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> Transform ensemble including product with random matrix'
          END IF
          IF (incremental == 1) &
               WRITE (*, '(a, 12x, a)') 'PDAF', '--> Perform incremental updating'
          IF (type_forget == 0) THEN
             WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> Use fixed forgetting factor:', forget
          ELSEIF (type_forget == 1) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> Use global adaptive forgetting factor'
          ELSEIF (type_forget == 2) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> Use local adaptive forgetting factors'
          ENDIF
          IF (dim_lag > 0) &
               WRITE (*, '(a, 12x, a, i6)') 'PDAF', '--> Apply smoother up to lag:',dim_lag
          IF (dim_bias_p > 0) WRITE (*, '(a, 12x, a)') 'PDAF', '--> perform bias correction'
          WRITE (*, '(a, 12x, a, i5)') 'PDAF', '--> ensemble size:', dim_ens
          IF (observe_ens) &
               WRITE (*, '(a, 12x, a, 1x, l)') 'PDAF', '--> observe_ens:', observe_ens
       ELSE
          WRITE (*, '(/5x, a/)') 'PDAF-ERROR: Invalid parameter setting - check prior output!'
       END IF

    END IF writeout

  END SUBROUTINE PDAF_letkf_init


!-------------------------------------------------------------------------------
!> Perform allocation of arrays for LETKF.
!!
!! __Revision history:__
!! * 2010-08 - Lars Nerger - Initial code from splitting PDAF_letkf_init
!! * 2025-02 - Lars Nerger - Restructuring introducing generic PDAF_alloc
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAF_letkf_alloc(outflag)

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

  END SUBROUTINE PDAF_letkf_alloc


!-------------------------------------------------------------------------------
!> Set integer parameter specific for LETKF
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
  SUBROUTINE PDAF_letkf_set_iparam(id, value, flag)

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
       incremental = value
       IF (incremental /= 0) THEN
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(10): LETKF does not yet support incremental updating!'
          flag = 10
       END IF
    CASE(5)
       type_forget = value
       IF (type_forget<0 .OR. type_forget>2) THEN
          WRITE (*, '(/5x, a/)') 'PDAF-ERROR(8): Invalid type of forgetting factor - param_int(5)!'
          flag = 8
       END IF
    CASE(6)
       type_trans = value
       IF (type_trans<0 .OR. type_trans>2) THEN
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(10): Invalid setting for ensemble transformation - param_int(6)!'
          flag = 8
       END IF
    CASE(7)
       ! Not used
    CASE(8)
       if (value==0) THEN
          observe_ens = .false. ! Apply H to ensemble mean to compute residual
       ELSE
          observe_ens = .true.  ! Apply H to X, compute mean of HX and then residual (the default for LETKF in PDAF2.3)
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

  END SUBROUTINE PDAF_letkf_set_iparam


!-------------------------------------------------------------------------------
!> Set real parameter specific for LETKF filter
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
  SUBROUTINE PDAF_letkf_set_rparam(id, value, flag)

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
       IF (localfilter == 0) THEN
          forget = value
       ELSE
          IF (inloop) THEN
             forget_l = value
          ELSE
             forget = value
          END IF
       END IF
       IF (forget <= 0.0) THEN
          WRITE (*,'(/5x,a/)') &
               'PDAF-ERROR(7): Invalid value of forgetting factor - param_real(1)!'
          flag = 7
       END IF
    CASE DEFAULT
       WRITE (*,'(/5x, a, i3/)') &
            'PDAF-WARNING: Invalid real parameter index', id
    END SELECT

  END SUBROUTINE PDAF_letkf_set_rparam

END MODULE PDAF_LETKF
