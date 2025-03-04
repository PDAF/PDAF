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
!> Module for 3DVAR holding shared parameters and some helper routines
!!
!! This module declares the parameters that are used in 3DVAR. 
!! Parameters that are specific for 3DVAR are declared while some
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
MODULE PDAF_3DVAR

  USE PDAF_mod_filter, &
       ONLY: localfilter, debug, dim_lag
  USE PDAF_estkf, &
       ONLY: type_forget, type_trans, type_sqrt, forget
  USE PDAF_lestkf, &
       ONLY: type_forget_l => type_forget, type_trans_l => type_trans, &
       type_sqrt_l => type_sqrt, forget_l => forget

  IMPLICIT NONE

! *** Integer parameters ***

  INTEGER :: type_opt=0            !< Type of minimizer for 3DVar
                                   !< (0) LBFGS, (1) CG+, (-1) steepest descent
  INTEGER :: dim_cvec=0            !< Size of control vector (fixed part)
  INTEGER :: dim_cvec_ens=0        !< Size of control vector (ensemble part)
  INTEGER :: m_lbfgs_var=5         !< Parameter 'm' of LBFGS
  INTEGER :: method_cgplus_var=2   !< Parameter 'method' of CG+
  INTEGER :: irest_cgplus_var=1    !< Parameter 'irest' of CG+
  INTEGER :: maxiter_cg_var=200    !< Parameter 'maxiter' of CG

! *** Real parameters ***
  REAL :: beta_3dvar=0.5           !< Hybrid weight for hybrid 3D-Var
  REAL :: eps_cg_var = 1.0e-6      !< Parameter 'EPS' of  CG
  REAL :: eps_cgplus_var = 1.0e-5  !< Parameter 'EPS' of CG+
  REAL :: pgtol_lbfgs_var=1.0e-5   !< Parameter 'pgtol' of LBFGS
  REAL :: factr_lbfgs_var=1.0e7    !< Parameter 'factr' of LBFGS


!-------------------------------------------------------------------------------
  
CONTAINS

!>  PDAF-internal initialization of 3DVAR
!!
!! Initialization of 3DVAR within PDAF. Performed are:
!! * initialize filter-specific parameters
!! * print screen information on filter configuration.
!!
!!  !  This is a core routine of PDAF and   !
!!  !   should not be changed by the user   !
!!
!! __Revision history:__
!! * 2021-03 - Lars Nerger - Initial code
!! *  Other revisions - see repository log
!!
  SUBROUTINE PDAF_3dvar_init(subtype, param_int, dim_pint, param_real, dim_preal, &
       ensemblefilter, fixedbasis, verbose, outflag)

    USE PDAF_mod_filter, &
         ONLY: dim_ens, localfilter, dim_lag
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
       WRITE(*, '(/a, 4x, a)') 'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '+++                      3D-Var                     +++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '+++                                                 +++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '+++      3D-Var variants implemented following      +++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '+++      Bannister, Q. J. Royal Meteorol. Soc.,     +++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '+++     143 (2017) 607-633, doi:10.1002/qj.2982     +++'
       WRITE(*, '(a, 4x, a)')  'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    END IF writeout


! ****************************
! *** INITIALIZE VARIABLES ***
! ****************************

    ! Set parameter default values - other defaults are set directly in the module
    observe_ens = .false.
    dim_lag = 0

    ! Parse provided parameters
    DO i=3, dim_pint
       CALL PDAF_3dvar_set_iparam(i, param_int(i), outflag)
    END DO
    DO i=1, dim_preal
       CALL PDAF_3dvar_set_rparam(i, param_real(i), outflag)
    END DO

    IF (subtype==0 .AND. dim_ens > 1) THEN
       WRITE (*, '(/5x, a/)') 'PDAF-ERROR(6): 3D-Var must be run with ensemble size = 1!'
       outflag = 6
    END IF


    ! Some special conditions
    IF (dim_pint<4) THEN
       IF (subtype==0 .OR. subtype==4 .OR. subtype==6 .OR. subtype==7) THEN
          WRITE (*, '(/5x, a/)') 'PDAF-ERROR(3): Missing specification of control vector dimension - param_int(4)!'
          outflag = 3
       END IF
    END IF

    IF (dim_pint<5) THEN
       IF (subtype==1 .OR. subtype==4) THEN
          dim_cvec_ens = dim_ens
       END IF
    END IF


    ! Define whether filter is domain localized
    localfilter = 0

    ! Define whether filter is mode-based or ensemble-based
    ensemblefilter = .TRUE.
 
    ! Initialize flag for fixed-basis filters
    IF (subtype == 2 .OR. subtype == 3) THEN
       fixedbasis = .TRUE.
    ELSE
       fixedbasis = .FALSE.
    END IF


! *********************
! *** Check subtype ***
! *********************

    IF (subtype<0 .OR. subtype>8 .OR. subtype==2 .OR. subtype==3 .OR. subtype==5) THEN
       WRITE (*, '(/5x, a/)') 'PDAF-ERROR(3): No valid subtype!'
       outflag = 3
    END IF

  END SUBROUTINE PDAF_3dvar_init


!-------------------------------------------------------------------------------
!> Perform allocation of arrays for 3DVAR.
!!
!! __Revision history:__
!! * 2019-05 - Lars Nerger - Initial code based on NETF
!! * 2025-02 - Lars Nerger - Restructuring introducing generic PDAF_alloc
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_3dvar_alloc(subtype, outflag)

    USE PDAF_mod_filter, &
         ONLY: dim_ens, dim_p, dim_bias_p
    USE PDAF_mod_filtermpi, &
         ONLY: dim_ens_l
    USE PDAF_utils, &
         ONLY: PDAF_alloc

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(inout):: outflag      !< Status flag
    INTEGER, INTENT(in) :: subtype        !< Sub-type of filter

! *** Local variable ***
    INTEGER :: dim_es                     ! Dimension of error space (size of Ainv)


! ******************************
! *** Allocate filter fields ***
! ******************************

    IF (subtype == 0) THEN
       dim_es = 1
    ELSE
       dim_es = dim_ens-1
    END IF

    CALL PDAF_alloc(dim_p, dim_ens, dim_ens_l, dim_es, dim_bias_p, &
         dim_lag, 0, outflag)

  END SUBROUTINE PDAF_3dvar_alloc


!-------------------------------------------------------------------------------
!>  Print information on configuration of 3DVAR
!!
!!  !  This is a core routine of PDAF and   !
!!  !   should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code by splitting from PDAF_3dvar_init
!! *  Other revisions - see repository log
!!
  SUBROUTINE PDAF_3dvar_config(subtype, verbose)

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

       WRITE (*, '(/a, 4x, a)') 'PDAF', '3DVAR configuration'
       WRITE (*, '(a, 10x, a, i1)') 'PDAF', 'filter sub-type= ', subtype
       IF (subtype == 0) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> 3DVAR incremental with control variable transform'
          WRITE (*, '(a, 12x, a, i7)') 'PDAF', '--> size of control vector', dim_cvec
       ELSEIF (subtype == 1) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> ensemble 3DVAR using LESTKF for ensemble transformation'
          WRITE (*, '(a, 12x, a, i7)') 'PDAF', '--> size of control vector', dim_cvec_ens
       ELSEIF (subtype == 4) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> ensemble 3DVAR using ESTKF for ensemble transformation'
          WRITE (*, '(a, 12x, a, i7)') 'PDAF', '--> size of control vector', dim_cvec_ens
       ELSEIF (subtype == 6) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> hybrid 3DVAR using LESTKF for ensemble transformation'
          WRITE (*, '(a, 12x, a, f10.3)') 'PDAF', '--> hybrid weight', beta_3dvar
          WRITE (*, '(a, 12x, a, i7)') 'PDAF', '--> total size of control vector', dim_cvec_ens + dim_cvec
          WRITE (*, '(a, 12x, a, 2i7)') 'PDAF', '--> size of ensemble and parameterized parts', dim_cvec_ens, dim_cvec
       ELSEIF (subtype == 7) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> hybrid 3DVAR using ESTKF for ensemble transformation'
          WRITE (*, '(a, 12x, a, f10.3)') 'PDAF', '--> hybrid weight', beta_3dvar
          WRITE (*, '(a, 12x, a, i7)') 'PDAF', '--> total size of control vector', dim_cvec_ens + dim_cvec
          WRITE (*, '(a, 12x, a, 2i7)') 'PDAF', '--> size of ensemble and parameterized parts', dim_cvec_ens, dim_cvec
       END IF
       IF (dim_lag > 0) &
            WRITE (*, '(a, 12x, a, i6)') 'PDAF', '--> Apply smoother up to lag:',dim_lag
       WRITE(*, '(a, 10x, a, i3)') &
            'PDAF', 'param_int(3) Solver: type_opt=', type_opt
       IF (type_opt == 1) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> LBFGS (default)'
       ELSE IF (type_opt == 2) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> CG+'
       ELSE IF (type_opt == 3) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> direct implementation of CG'
       ELSE IF (type_opt == 4) THEN
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> direct implementation of CG with decomposed control vector'
       END IF
       IF (subtype==0 .OR. subtype>=6) &
            WRITE(*, '(a, 10x, a, i3)') &
            'PDAF', 'param_int(4) dim_cvec=', dim_cvec
       IF (subtype>0) &
            WRITE(*, '(a, 10x, a, i3)') &
            'PDAF', 'param_int(5) dim_cvec_ens=', dim_cvec_ens
       IF (type_opt == 1) THEN
          WRITE(*, '(a, 10x, a, i3)') &
               'PDAF', 'param_int(6) solver-specific parameter: m=', m_lbfgs_var
       ELSEIF (type_opt == 2) THEN
          WRITE(*, '(a, 10x, a, i3)') &
               'PDAF', 'param_int(6) solver-specific parameter: CG method=', method_cgplus_var
          IF (method_cgplus_var==1) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> Fletcher-Reeves'
          ELSE IF (method_cgplus_var==2) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> Polak-Ribiere'
          ELSE IF (method_cgplus_var==3) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> positive Polak-Ribiere'
          END IF
          WRITE(*, '(a, 10x, a, i7)') &
               'PDAF', 'param_int(7) solver-specific parameter: number of restarts=', irest_cgplus_var
       ELSEIF (type_opt == 3 .OR. type_opt==4) THEN
          WRITE(*, '(a, 10x, a, i7)') &
               'PDAF', 'param_int(6) solver-specific parameter: maximum number of iterations=', maxiter_cg_var
       END IF
       IF (subtype>0) THEN
          WRITE(*, '(a, 10x, a)') &
               'PDAF', 'param_int(8) observe_ens'
          IF (observe_ens) THEN
             WRITE(*, '(a, 12x, a)') 'PDAF', '--> 1: Apply H to ensemble states and compute innovation as mean (default)'
          ELSE
             WRITE(*, '(a, 12x, a)') 'PDAF', '--> 0: Apply H to ensemble mean to compute innovation'
          END IF
       END IF
       WRITE(*, '(a, 10x, a, i3)') &
            'PDAF', 'param_int(9) type_obs_init=', type_obs_init
       IF (type_obs_init==0) THEN
          WRITE(*, '(a, 12x, a)') 'PDAF', '--> Initialize observations before PDAF prestep'
       ELSE IF (type_obs_init==1) THEN
          WRITE(*, '(a, 12x, a)') 'PDAF', '--> Initialize observations after PDAF prestep'
       END IF

       IF (subtype>0) THEN
          IF (subtype==1 .OR. subtype==4) THEN
             WRITE(*, '(a, 10x, a)') &
                  'PDAF', '___ Parameters for LESTKF ___'
          ELSE
             WRITE(*, '(a, 10x, a)') &
                  'PDAF', '___ Parameters for ESTKF ___'
          END IF
          WRITE (*, '(a, 12x, a, i5)') 'PDAF', '---> ensemble size:', dim_ens
          WRITE(*, '(a, 10x, a, i3)') &
               'PDAF', 'param_int(11) type_forget=', type_forget
          IF (subtype==4 .OR. subtype==7) THEN
             IF (type_forget == 0) THEN
                WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> Use fixed forgetting factor:', forget
             ELSEIF (type_forget == 1) THEN
                WRITE (*, '(a, 12x, a)') 'PDAF', '--> Use adaptive forgetting factor'
             ENDIF
          ELSEIF (subtype==1 .OR. subtype==6) THEN
             IF (type_forget == 0) THEN
                WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> Use fixed forgetting factor:', forget
             ELSEIF (type_forget == 1) THEN
                WRITE (*, '(a, 12x, a)') 'PDAF', '--> Use global adaptive forgetting factor'
             ELSEIF (type_forget == 2) THEN
                WRITE (*, '(a, 12x, a)') 'PDAF', '--> Use local adaptive forgetting factors'
             ENDIF
          END IF
          WRITE(*, '(a, 10x, a, i3)') &
               'PDAF', 'param_int(12) type_trans=', type_trans
          IF (type_trans == 0) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> Deterministic ensemble transformation (default)'
          ELSE IF (type_trans == 1) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> Transform ensemble with random orthonormal Omega'
          ELSE IF (type_trans == 2) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> Transform ensemble including product with random matrix'
          END IF
          WRITE(*, '(a, 10x, a, i3)') &
               'PDAF', 'param_int(13) type_sqrt=', type_sqrt
          IF (type_sqrt == 0) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> symmetric square root (default)'
          ELSE IF (type_sqrt == 1) THEN
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> Cholesky decomposition'
          END IF
       END IF

       IF (subtype>0) &
            WRITE(*, '(a, 10x, a, f10.3)') &
            'PDAF', 'param_real(1) forget=', forget
       IF (subtype>=6) &
            WRITE(*, '(a, 10x, a, f10.3)') &
            'PDAF', 'param_real(2) hybrid weight in hyb3DVar: beta_3dvar=', beta_3dvar
       IF (type_opt == 1) THEN
          WRITE(*, '(a, 10x, a, es10.3)') &
               'PDAF', 'param_real(3) solver-specific parameter: pgtol=', pgtol_lbfgs_var
          WRITE(*, '(a, 10x, a, es10.3)') &
               'PDAF', 'param_real(4) solver-specific parameter: factr=', factr_lbfgs_var
       ELSEIF (type_opt == 2) THEN
          WRITE(*, '(a, 10x, a, es10.3)') &
               'PDAF', 'param_real(3) solver-specific parameter: eps=', eps_cgplus_var
       ELSEIF (type_opt == 3 .OR. type_opt == 4) THEN
          WRITE(*, '(a, 10x, a, es10.3)') &
               'PDAF', 'param_real(3) solver-specific parameter: eps=', eps_cg_var
       END IF

    END IF writeout

  END SUBROUTINE PDAF_3dvar_config


!-------------------------------------------------------------------------------
!> Set integer parameter specific for 3DVAR
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Other revisions - see repository log
!!
  SUBROUTINE PDAF_3dvar_set_iparam(id, value, flag)

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
       type_opt = value
       IF (.NOT.(type_opt==1 .OR. type_opt==2 .OR. type_opt==3 &
            .OR. type_opt==12 .OR. type_opt==13)) THEN
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(8): Invalid type of solver - param_int(3)!'
          flag = 8
       END IF
    CASE(4)
       dim_cvec = value
       IF (dim_cvec < 0) THEN
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(8): dim_cvec cannot be below 1 - param_int(3)!'
          flag = 10
       END IF
    CASE(5)
       dim_cvec_ens = value
       IF (dim_cvec < 0) THEN
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(8): dim_cvec_ens cannot be below 1 - param_int(3)!'
          flag = 10
       END IF
    CASE(6)
       m_lbfgs_var = value
       method_cgplus_var = value
       maxiter_cg_var = value
    CASE(7)
       irest_cgplus_var = value
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
    CASE(10)
       ! Not used
    CASE(11)
       type_forget = value
       type_forget_l = value
       IF (type_forget<0 .OR. type_forget>2) THEN
          WRITE (*, '(/5x, a/)') 'PDAF-ERROR(8): Invalid type of forgetting factor - param_int(11)!'
          flag = 8
       END IF
    CASE(12)
       type_trans = value
       type_trans_l = value
       IF (type_trans<0 .OR. type_trans>2) THEN
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(8): Invalid setting for ensemble transformation - param_int(12)!'
          flag = 8
       END IF
    CASE(13)
       type_sqrt = value
       type_sqrt_l = value
       IF (type_sqrt<0 .OR. type_sqrt>1) THEN
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(8): Invalid setting for square root type - param_int(13)!'
          flag = 8
       END IF
    CASE DEFAULT
       WRITE (*,'(/5x, a, i3/)') &
            'PDAF-WARNING: Invalid integer parameter index', id
    END SELECT

  END SUBROUTINE PDAF_3dvar_set_iparam


!-------------------------------------------------------------------------------
!> Set real parameter specific for 3DVAR
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Other revisions - see repository log
!!
  SUBROUTINE PDAF_3dvar_set_rparam(id, value, flag)

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
       forget_l = value
       IF (forget <= 0.0) THEN
          WRITE (*,'(/5x,a/)') &
               'PDAF-ERROR(7): Invalid value of forgetting factor - param_real(1)!'
          flag = 7
       END IF
    CASE(2)
       beta_3dvar = value
    CASE(3)
       eps_cg_var = value
       eps_cgplus_var = value
       pgtol_lbfgs_var = value
    CASE(4)
       factr_lbfgs_var = value
    CASE DEFAULT
       WRITE (*,'(/5x, a, i3/)') &
            'PDAF-WARNING: Invalid real parameter index', id
    END SELECT

  END SUBROUTINE PDAF_3dvar_set_rparam

!-------------------------------------------------------------------------------
!> Information output on options for 3DVAR
!!
!! Subroutine to perform information output on options
!! available for the Ensemble Transform Kalman filter.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2011-08 - Lars Nerger - Initial code
!! *Other revisions - see repository log
!!
  SUBROUTINE PDAF_3dvar_options()

    IMPLICIT NONE

    WRITE(*, '(/a, 5x, a)') 'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    WRITE(*, '(a, 5x, a)')  'PDAF', '+++                      3D-Var                     +++'
    WRITE(*, '(a, 5x, a)')  'PDAF', '+++                                                 +++'
    WRITE(*, '(a, 5x, a)')  'PDAF', '+++      3D-Var variants implemented following      +++'
    WRITE(*, '(a, 5x, a)')  'PDAF', '+++      Bannister, Q. J. Royal Meteorol. Soc.,     +++'
    WRITE(*, '(a, 5x, a)')  'PDAF', '+++     143 (2017) 607-633, doi:10.1002/qj.2982     +++'
    WRITE(*, '(a, 5x, a)')  'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'

    WRITE(*, '(/a, 5x, a)') 'PDAF', 'Available options for 3D-Var:'

    WRITE(*, '(a, 5x, a)') 'PDAF', '--- Sub-types (Parameter subtype) ---'
    WRITE(*, '(a, 7x, a)') 'PDAF', '0: incremental 3D-Var with parameterized covariance matrix'
    WRITE(*, '(a, 7x, a)') 'PDAF', '1: 3D ensemble Var using LESTKF for ensemble transformation'
    WRITE(*, '(a, 7x, a)') 'PDAF', '4: 3D ensemble Var using ESTKF for ensemble transformation'
    WRITE(*, '(a, 7x, a)') 'PDAF', '6: hybrid 3D-Var using LESTKF for ensemble transformation'
    WRITE(*, '(a, 7x, a)') 'PDAF', '7: hybrid 3D-Var using ESTKF for ensemble transformation'

    WRITE(*, '(a, 5x, a)') 'PDAF', '--- Integer parameters (Array param_int) ---'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(1): Dimension of state vector (>0), required'
    WRITE(*, '(a, 7x, a)') 'PDAF', 'param_int(2): Ensemble size (>0), required'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(3): type_opt'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'Select optimization method (solver), required'
    WRITE(*, '(a, 12x, a)') 'PDAF', '0: LBFGS (default)'
    WRITE(*, '(a, 12x, a)') 'PDAF', '1: CG+'
    WRITE(*, '(a, 12x, a)') 'PDAF', '2: direct implementation of CG'
    WRITE(*, '(a, 12x, a)') 'PDAF', '3: direct implementation of CG with decomposed control vector'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(4): size of parameterized control vector (for 3D-Var and hybrid 3D-Var), required'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(5): size of ensemble control vector (required for ensemble and hybrid 3D-Var), '
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(6): Solver-specific parameter, optional'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'LBFGS: parameter m (default=5)'
    WRITE(*, '(a, 16x, a)') 'PDAF', 'Number of corrections used in limited memory matrix; 3<=m<=20'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'CG+: parameter method (default=2)'
    WRITE(*, '(a, 16x, a)') 'PDAF', '(1) Fletcher-Reeves, (2) Polak-Ribiere, (3) positive Polak-Ribiere'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'CG: maximum number of iterations (default=200)'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(7): Solver-specific parameter, optional'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'CG+: parameter irest (default=1)'
    WRITE(*, '(a, 16x, a)') 'PDAF', '(0) no restarts; (n>0) restart every n steps'
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
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(10): not used'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', '___Options for ESTKF/LESTKF for En3DVar/hyb3DVar___'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(11) type_forget'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'Type of forgetting factor; optional'
    WRITE(*, '(a, 12x, a)') 'PDAF', '0: fixed forgetting factor (default)'
    WRITE(*, '(a, 12x, a)') 'PDAF', '1: adaptive forgetting factor (experimental)'
    WRITE(*, '(a, 12x, a)') 'PDAF', '2: locally adaptive forgetting factor (experimental)'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(12) type_trans'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'Type of ensemble transformation matrix; optional'
    WRITE(*, '(a, 12x, a)') 'PDAF', '0: deterministic Omega (default)'
    WRITE(*, '(a, 12x, a)') 'PDAF', '1: random orthonormal Omega orthogonal to (1,...,1)^T'
    WRITE(*, '(a, 12x, a)') &
         'PDAF', '2: use product of 0 with random orthonomal matrix with eigenvector (1,...,1)^T'
    WRITE(*, '(a, 14x, a)') &
         'PDAF', '(experimental; for random transformations, 0 or 1 are recommended)'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_int(13) type_sqrt'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'Type of transformation matrix square root; optional'
    WRITE(*, '(a, 12x, a)') 'PDAF', '0: symmetric square root (default)'
    WRITE(*, '(a, 12x, a)') 'PDAF', '1: Cholesky decomposition'

    WRITE(*, '(a, 5x, a)') 'PDAF', '--- Floating point parameters (Array param_real) ---'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_real(1): Forgetting factor (usually >0 and <=1), required;'
    WRITE(*, '(a, 11x, a)') 'PDAF', '(only used for ensemble and hybrid 3D-Var)'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_real(2): hybrid weight beta, optional (only for hybrid 3D-Var)'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'range >=0.0 and <=1.0, =1.0 for pure ensemble 3D-var  (default=0.5)'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_real(3): Solver-specific parameter, optional'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'LBFGS: Limit for stopping iterations (pgtol, default=1.0e-5)'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'CG+: convergence parameter eps (default=1.0e-5)'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'CG: convergence parameter eps (default=1.0e-6)'
    WRITE(*, '(a, 7x, a)') &
         'PDAF', 'param_real(4): Solver-specific parameter, optional'
    WRITE(*, '(a, 11x, a)') 'PDAF', 'LBFGS: Tolerance in termination test (factr, default=1.0e+7)'

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
         'PDAF', '+++++++++ End of option overview for 3DVAR ++++++++++'

  END SUBROUTINE PDAF_3dvar_options


!-------------------------------------------------------------------------------
!> Display timing and memory information for 3DVAR
!!
!! This routine displays the PDAF-internal timing and
!! memory information for the 3DVAR.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2008-09 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_3dvar_memtime(printtype)

    USE PDAF_timer, &
         ONLY: PDAF_time_tot
    USE PDAF_memcounting, &
         ONLY: PDAF_memcount_get, PDAF_memcount_get_global
    USE PDAF_mod_filter, &
         ONLY: subtype_filter, offline_mode, dim_lag
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
          IF (subtype_filter==0) THEN
             WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'State forecast:', pdaf_time_tot(2), 's'
          ELSE
             WRITE (*, '(a, 16x, a, F11.3, 1x, a)') 'PDAF', 'Ensemble forecast:', pdaf_time_tot(2), 's'
          END IF
       END IF

       IF (filterpe) THEN
          ! Filter-specific part
          WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', '3DVAR analysis:', pdaf_time_tot(3), 's'

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
          IF (subtype_filter==0) THEN
             WRITE (*, '(a, 10x, a, 16x, F11.3, 1x, a)') 'PDAF', 'State forecast:', pdaf_time_tot(2), 's'
          ELSE
             WRITE (*, '(a, 10x, a, 13x, F11.3, 1x, a)') 'PDAF', 'Ensemble forecast:', pdaf_time_tot(2), 's'
          END IF
          WRITE (*, '(a, 12x, a, 5x, F11.3, 1x, a)') 'PDAF', 'MPI communication in PDAF:', pdaf_time_tot(4), 's'
          WRITE (*, '(a, 12x, a, 9x, F11.3, 1x, a)') 'PDAF', 'distribute_state_pdaf:', pdaf_time_tot(40), 's'
          WRITE (*, '(a, 12x, a, 12x, F11.3, 1x, a)') 'PDAF', 'collect_state_pdaf:', pdaf_time_tot(41), 's'
          IF (.not.filterpe) WRITE (*, '(a, 7x, a)') 'PDAF', &
               'Note: for filterpe=F, the time (2) includes the wait time for the analysis step'
       END IF

       IF (filterpe) THEN
! Filter-specific part

          IF (subtype_filter==0) THEN
             WRITE (*, '(a, 10x, a, 16x, F11.3, 1x, a)') 'PDAF', '3DVAR analysis:', pdaf_time_tot(3), 's'
          ELSEIF (subtype_filter>0 .AND. subtype_filter<4) THEN
             WRITE (*, '(a, 10x, a, 14x, F11.3, 1x, a)') 'PDAF', 'En3DVAR analysis:', pdaf_time_tot(3), 's'
          ELSE
             WRITE (*, '(a, 10x, a, 13x, F11.3, 1x, a)') 'PDAF', 'Hyb3DVAR analysis:', pdaf_time_tot(3), 's'
          END IF
          WRITE (*, '(a, 12x, a, 6x, F11.3, 1x, a)') 'PDAF', 'PDAF-internal operations:', pdaf_time_tot(51), 's'

          IF(omi_was_used) THEN
             ! Output when using OMI

             time_omi = pdaf_time_tot(50) + pdaf_time_tot(48)
             IF(subtype_filter==1 .OR. subtype_filter==6) THEN
                time_omi = time_omi + pdaf_time_tot(46) + pdaf_time_tot(47)
                IF (type_forget==1) &
                     time_omi = time_omi + pdaf_time_tot(49) + pdaf_time_tot(52)
             END IF
             WRITE (*, '(a, 12x, a, 9x, F11.3, 1x, a)') 'PDAF', 'OMI-internal routines:', &
                  time_omi, 's'
             WRITE (*, '(a, 12x, a, 24x, F11.3, 1x, a)') 'PDAF', 'Solver:', pdaf_time_tot(54), 's'
             IF (subtype_filter==0) THEN
                WRITE (*, '(a, 12x, a, 22x, F11.3, 1x, a)') 'PDAF', 'cvt_pdaf:', pdaf_time_tot(60), 's'
                WRITE (*, '(a, 12x, a, 18x, F11.3, 1x, a)') 'PDAF', 'cvt_adj_pdaf:', pdaf_time_tot(62), 's'
             ELSE
                WRITE (*, '(a, 12x, a, 18x, F11.3, 1x, a)') 'PDAF', 'cvt_ens_pdaf:', pdaf_time_tot(61), 's'
                WRITE (*, '(a, 12x, a, 14x, F11.3, 1x, a)') 'PDAF', 'cvt_ens_adj_pdaf:', pdaf_time_tot(63), 's'
                IF(subtype_filter==1 .OR. subtype_filter==6) THEN
                   WRITE (*, '(a, 12x, a)') 'PDAF', 'Timers in LESTKF only'
                   WRITE (*, '(a, 14x, a, 9x, F11.3, 1x, a)') 'PDAF', 'init_n_domains_pdaf:', pdaf_time_tot(42), 's'
                   WRITE (*, '(a, 14x, a, 13x, F11.3, 1x, a)') 'PDAF', 'init_dim_l_pdaf:', pdaf_time_tot(45), 's'
                   IF (.NOT.pdaflocal_was_used) THEN
                      WRITE (*, '(a, 14x, a, 14x, F11.3, 1x, a)') 'PDAF', 'g2l_state_pdaf:', pdaf_time_tot(10), 's'
                      WRITE (*, '(a, 14x, a, 14x, F11.3, 1x, a)') 'PDAF', 'l2g_state_pdaf:', pdaf_time_tot(14), 's'
                   END IF
                END IF
             END IF

             WRITE (*, '(a, 12x, a)') 'PDAF', 'Time in OMI observation module routines '
             WRITE (*, '(a, 14x, a, 8x, F11.3, 1x, a)') 'PDAF', 'init_dim_obs_pdafomi:', pdaf_time_tot(43), 's'
             WRITE (*, '(a, 14x, a, 14x, F11.3, 1x, a)') 'PDAF', 'obs_op_pdafomi:', pdaf_time_tot(44), 's'
             WRITE (*, '(a, 14x, a, 10x, F11.3, 1x, a)') 'PDAF', 'obs_op_lin_pdafomi:', pdaf_time_tot(64), 's'
             WRITE (*, '(a, 14x, a, 10x, F11.3, 1x, a)') 'PDAF', 'obs_op_adj_pdafomi:', pdaf_time_tot(65), 's'
             IF(subtype_filter==1 .OR. subtype_filter==6) &
                  WRITE (*, '(a, 14x, a, 6x, F11.3, 1x, a)') 'PDAF', 'init_dim_obs_l_pdafomi:', pdaf_time_tot(38), 's'

!            WRITE (*, '(a, 12x, a)') 'PDAF', 'Time in OMI-internal routines'
!            IF(subtype_filter==1 .OR. subtype_filter==6) THEN
!               IF (type_forget==1) THEN
!                  WRITE (*, '(a, 14x, a, 9x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_init_obsvar:', pdaf_time_tot(49), 's'
!               END IF
!               WRITE (*, '(a, 14x, a, 13x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_g2l_obs:', pdaf_time_tot(46), 's'
!               IF (type_forget==1) THEN
!                  WRITE (*, '(a, 14x, a, 7x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_init_obsvar_l:', pdaf_time_tot(52), 's'
!               END IF
!               WRITE (*, '(a, 14x, a, 10x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_init_obs_l:', pdaf_time_tot(47), 's'
!            END IF
! 
!            WRITE (*, '(a, 14x, a, 12x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_init_obs:', pdaf_time_tot(50), 's'
!            WRITE (*, '(a, 14x, a, 11x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_prodRinvA:', pdaf_time_tot(48), 's'
          ELSE
             ! Output when NOT using OMI
             WRITE (*, '(a, 12x, a, 24x, F11.3, 1x, a)') 'PDAF', 'Solver:', pdaf_time_tot(54), 's'
             WRITE (*, '(a, 12x, a, 13x, F11.3, 1x, a)') 'PDAF', 'init_dim_obs_pdaf:', pdaf_time_tot(43), 's'
             WRITE (*, '(a, 12x, a, 19x, F11.3, 1x, a)') 'PDAF', 'obs_op_pdaf:', pdaf_time_tot(44), 's'
             WRITE (*, '(a, 12x, a, 17x, F11.3, 1x, a)') 'PDAF', 'init_obs_pdaf:', pdaf_time_tot(50), 's'
             WRITE (*, '(a, 12x, a, 16x, F11.3, 1x, a)') 'PDAF', 'prodRinvA_pdaf:', pdaf_time_tot(48), 's'
             IF (subtype_filter==0) THEN
                WRITE (*, '(a, 12x, a, 22x, F11.3, 1x, a)') 'PDAF', 'cvt_pdaf:', pdaf_time_tot(60), 's'
                WRITE (*, '(a, 12x, a, 15x, F11.3, 1x, a)') 'PDAF', 'obs_op_lin_pdaf:', pdaf_time_tot(64), 's'
                WRITE (*, '(a, 12x, a, 18x, F11.3, 1x, a)') 'PDAF', 'cvt_adj_pdaf:', pdaf_time_tot(62), 's'
                WRITE (*, '(a, 12x, a, 15x, F11.3, 1x, a)') 'PDAF', 'obs_op_adj_pdaf:', pdaf_time_tot(65), 's'
             ELSE
                WRITE (*, '(a, 12x, a, 18x, F11.3, 1x, a)') 'PDAF', 'cvt_ens_pdaf:', pdaf_time_tot(61), 's'
                WRITE (*, '(a, 12x, a, 11x, F11.3, 1x, a)') 'PDAF', 'obs_ens_op_lin_pdaf:', pdaf_time_tot(64), 's'
                WRITE (*, '(a, 12x, a, 14x, F11.3, 1x, a)') 'PDAF', 'cvt_ens_adj_pdaf:', pdaf_time_tot(63), 's'
                WRITE (*, '(a, 12x, a, 11x, F11.3, 1x, a)') 'PDAF', 'obs_ens_op_adj_pdaf:', pdaf_time_tot(65), 's'
                IF(subtype_filter==1 .OR. subtype_filter==6) THEN
                   WRITE (*, '(a, 10x, a)') 'PDAF', 'Timers in LESTKF only'
                   WRITE (*, '(a, 12x, a, 11x, F11.3, 1x, a)') 'PDAF', 'init_n_domains_pdaf:', pdaf_time_tot(42), 's'
                   WRITE (*, '(a, 12x, a, 11x, F11.3, 1x, a)') 'PDAF', 'init_dim_obs_f_pdaf:', pdaf_time_tot(43), 's'
                   IF (type_forget==1) THEN
                      WRITE (*, '(a, 12x, a, 15x, F11.3, 1x, a)') 'PDAF', 'init_obs_f_pdaf:', pdaf_time_tot(50), 's'
                      WRITE (*, '(a, 12x, a, 14x, F11.3, 1x, a)') 'PDAF', 'init_obsvar_pdaf:', pdaf_time_tot(49), 's'
                   END IF
                   WRITE (*, '(a, 12x, a, 15x, F11.3, 1x, a)') 'PDAF', 'init_dim_l_pdaf:', pdaf_time_tot(45), 's'
                   WRITE (*, '(a, 12x, a, 11x, F11.3, 1x, a)') 'PDAF', 'init_dim_obs_l_pdaf:', pdaf_time_tot(38), 's'
                   WRITE (*, '(a, 12x, a, 16x, F11.3, 1x, a)') 'PDAF', 'g2l_state_pdaf:', pdaf_time_tot(10), 's'
                   IF (type_forget==1) THEN
                      WRITE (*, '(a, 12x, a, 12x, F11.3, 1x, a)') 'PDAF', 'init_obsvar_l_pdaf:', pdaf_time_tot(52), 's'
                   END IF
                   WRITE (*, '(a, 12x, a, 15x, F11.3, 1x, a)') 'PDAF', 'init_obs_l_pdaf:', pdaf_time_tot(47), 's'
                   WRITE (*, '(a, 12x, a, 16x, F11.3, 1x, a)') 'PDAF', 'l2g_state_pdaf:', pdaf_time_tot(14), 's'
                END IF
             END IF
          END IF

          ! Generic part B
          WRITE (*, '(a, 10x, a, 15x, F11.3, 1x, a)') 'PDAF', 'Prepoststep (5):', pdaf_time_tot(5), 's'
       END IF
    ELSE IF (printtype == 4) THEN ptype

! *********************************************
! *** Print second-level timing information ***
! *********************************************

       ! Generic part
       WRITE (*, '(//a, 23x, a)') 'PDAF', 'PDAF Timing information'
       WRITE (*, '(a, 8x, 52a)') 'PDAF', ('-', i=1, 52)
       WRITE (*, '(a, 10x, a, 11x, F11.3, 1x, a)') &
            'PDAF', 'Initialize PDAF (1):', pdaf_time_tot(1), 's'
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
          IF (subtype_filter==0) THEN
             WRITE (*, '(a, 10x, a, 12x, F11.3, 1x, a)') 'PDAF', '3DVAR analysis (3):', pdaf_time_tot(3), 's'
          ELSEIF(subtype_filter==1 .OR. subtype_filter==6) THEN
             WRITE (*, '(a, 10x, a, 10x, F11.3, 1x, a)') 'PDAF', 'En3DVAR analysis (3):', pdaf_time_tot(3), 's'
          ELSEIF(subtype_filter==4 .OR. subtype_filter==7) THEN
             WRITE (*, '(a, 10x, a, 9x, F11.3, 1x, a)') 'PDAF', 'hyb3DVAR analysis (3):', pdaf_time_tot(3), 's'
          END IF
          IF (subtype_filter>0) &
               WRITE (*, '(a, 12x, a, 5x, F11.3, 1x, a)') 'PDAF', 'compute ensemble mean (9):', pdaf_time_tot(9), 's'
          WRITE (*, '(a, 12x, a, 6x, F11.3, 1x, a)') 'PDAF', 'prepare observations (6):', pdaf_time_tot(6), 's'
          WRITE (*, '(a, 12x, a, 10x, F11.3, 1x, a)') 'PDAF', 'init innovation (13):', pdaf_time_tot(13), 's'
          WRITE (*, '(a, 12x, a, 9x, F11.3, 1x, a)') 'PDAF', 'run optimization (52):', pdaf_time_tot(52), 's'
          WRITE (*, '(a, 14x, a, 5x, F11.3, 1x, a)') 'PDAF', 'compute J and grad J (53):', pdaf_time_tot(53), 's'
          WRITE (*, '(a, 14x, a, 11x, F11.3, 1x, a)') 'PDAF', 'execute solver (54):', pdaf_time_tot(54), 's'
          WRITE (*, '(a, 12x, a, 6x, F11.3, 1x, a)') 'PDAF', 'update state vector (19):', pdaf_time_tot(19), 's'
          IF(subtype_filter==1 .OR. subtype_filter==6) THEN
             WRITE (*, '(a, 10x, a)') 'PDAF', 'Timers in ESTKF only'
             WRITE (*, '(a, 12x, a, 10x, F11.3, 1x, a)') 'PDAF', 'init innovation (10):', pdaf_time_tot(10), 's'
             WRITE (*, '(a, 12x, a, 5x, F11.3, 1x, a)') 'PDAF', 'compute inverse of A (11):', pdaf_time_tot(11), 's'
             WRITE (*, '(a, 12x, a, 2x, F11.3, 1x, a)') 'PDAF', 'get state weight vector (12):', pdaf_time_tot(12), 's'
             WRITE (*, '(a, 12x, a, 1x, F11.3, 1x, a)') 'PDAF', 'prepare ensemble weights (20):', pdaf_time_tot(20), 's'
             WRITE (*, '(a, 12x, a, 7x, F11.3, 1x, a)') 'PDAF', 'transform ensemble (21):', pdaf_time_tot(21), 's'
             IF (dim_lag >0) &
                  WRITE (*, '(a, 20x, a, F11.3, 1x, a)') 'PDAF', 'perform smoothing (15):', pdaf_time_tot(15), 's'
          ELSEIF(subtype_filter==4 .OR. subtype_filter==7) THEN
             WRITE (*, '(a, 10x, a)') 'PDAF', 'Timers in LESTKF only'
             WRITE (*, '(a, 12x, a, 7x, F11.3, 1x, a)') 'PDAF', 'global preparations (7):', pdaf_time_tot(7), 's'
             WRITE (*, '(a, 12x, a, 7x, F11.3, 1x, a)') 'PDAF', 'local analysis loop (8):', pdaf_time_tot(8), 's'
             WRITE (*, '(a, 14x, a, 10x, F11.3, 1x, a)') 'PDAF', 'global to local (10):', pdaf_time_tot(10), 's'
             WRITE (*, '(a, 14x, a, 4x, F11.3, 1x, a)') 'PDAF', 'localize observations (11):', pdaf_time_tot(11), 's'
             WRITE (*, '(a, 14x, a, 11x, F11.3, 1x, a)') 'PDAF', 'local analysis (12):', pdaf_time_tot(12), 's'
             WRITE (*, '(a, 14x, a, 10x, F11.3, 1x, a)') 'PDAF', 'local to global (14):', pdaf_time_tot(14), 's'
             IF (dim_lag >0) &
                  WRITE (*, '(a, 14x, a, 8x, F11.3, 1x, a)') 'PDAF', 'perform smoothing (15):', pdaf_time_tot(15), 's'
          END IF

          ! Generic part B
          WRITE (*, '(a, 10x, a, 15x, F11.3, 1x, a)') 'PDAF', 'Prepoststep (5):', pdaf_time_tot(5), 's'
       END IF

    ELSE IF (printtype == 5) THEN ptype

! *****************************************
! *** Print detailed timing information ***
! *****************************************

       ! Generic part
       WRITE (*, '(//a, 21x, a)') 'PDAF', 'PDAF Timing information'
       WRITE (*, '(a, 8x, 52a)') 'PDAF', ('-', i=1, 52)
       WRITE (*, '(a, 10x, a, 11x, F11.3, 1x, a)') &
            'PDAF', 'Initialize PDAF (1):', pdaf_time_tot(1), 's'
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
          IF (subtype_filter==0) THEN
             WRITE (*, '(a, 10x, a, 12x, F11.3, 1x, a)') 'PDAF', '3DVAR analysis (3):', pdaf_time_tot(3), 's'
          ELSEIF(subtype_filter==1 .OR. subtype_filter==6) THEN
             WRITE (*, '(a, 10x, a, 10x, F11.3, 1x, a)') 'PDAF', 'En3DVAR analysis (3):', pdaf_time_tot(3), 's'
          ELSEIF(subtype_filter==4 .OR. subtype_filter==7) THEN
             WRITE (*, '(a, 10x, a, 9x, F11.3, 1x, a)') 'PDAF', 'hyb3DVAR analysis (3):', pdaf_time_tot(3), 's'
          END IF
          IF (subtype_filter>0) &
               WRITE (*, '(a, 12x, a, 5x, F11.3, 1x, a)') 'PDAF', 'compute ensemble mean (9):', pdaf_time_tot(9), 's'
          WRITE (*, '(a, 12x, a, 6x, F11.3, 1x, a)') 'PDAF', 'prepare observations (6):', pdaf_time_tot(6), 's'
          WRITE (*, '(a, 12x, a, 10x, F11.3, 1x, a)') 'PDAF', 'init innovation (13):', pdaf_time_tot(13), 's'
          WRITE (*, '(a, 12x, a, 9x, F11.3, 1x, a)') 'PDAF', 'run optimization (52):', pdaf_time_tot(52), 's'
          WRITE (*, '(a, 14x, a, 5x, F11.3, 1x, a)') 'PDAF', 'compute J and grad J (53):', pdaf_time_tot(53), 's'
          WRITE (*, '(a, 16x, a, 4x, F11.3, 1x, a)') 'PDAF', 'compute cost function (55):', pdaf_time_tot(55), 's'
          WRITE (*, '(a, 16x, a, 7x, F11.3, 1x, a)') 'PDAF', 'J observation part (56):', pdaf_time_tot(56), 's'
          WRITE (*, '(a, 16x, a, 8x, F11.3, 1x, a)') 'PDAF', 'J background part (57):', pdaf_time_tot(57), 's'
          WRITE (*, '(a, 16x, a, 11x, F11.3, 1x, a)') 'PDAF', 'compute grad J (58):', pdaf_time_tot(58), 's'
          IF (type_opt==3 .OR. type_opt==13) &
               WRITE (*, '(a, 14x, a, 2x, F11.3, 1x, a)') 'PDAF', 'compute Hessian times d (59):', pdaf_time_tot(59), 's'
          WRITE (*, '(a, 14x, a, 11x, F11.3, 1x, a)') 'PDAF', 'execute solver (54):', pdaf_time_tot(54), 's'
          WRITE (*, '(a, 12x, a, 6x, F11.3, 1x, a)') 'PDAF', 'update state vector (19):', pdaf_time_tot(19), 's'
          IF(subtype_filter==1 .OR. subtype_filter==6) THEN
             WRITE (*, '(a, 10x, a)') 'PDAF', 'Timers in ESTKF only'
             WRITE (*, '(a, 12x, a, 10x, F11.3, 1x, a)') 'PDAF', 'init innovation (10):', pdaf_time_tot(10), 's'
             WRITE (*, '(a, 12x, a, 5x, F11.3, 1x, a)') 'PDAF', 'compute inverse of A (11):', pdaf_time_tot(11), 's'
             WRITE (*, '(a, 14x, a, 12x, F11.3, 1x, a)') 'PDAF', 'complete Ainv (31):', pdaf_time_tot(31), 's'
             WRITE (*, '(a, 12x, a, 2x, F11.3, 1x, a)') 'PDAF', 'get state weight vector (12):', pdaf_time_tot(12), 's'
             WRITE (*, '(a, 12x, a, 1x, F11.3, 1x, a)') 'PDAF', 'prepare ensemble weights (20):', pdaf_time_tot(20), 's'
             WRITE (*, '(a, 14x, a, 15x, F11.3, 1x, a)') 'PDAF', 'SQRT(Ainv) (32):', pdaf_time_tot(32), 's'
             WRITE (*, '(a, 14x, a, 15x, F11.3, 1x, a)') 'PDAF', 'init Omega (33):', pdaf_time_tot(33), 's'
             WRITE (*, '(a, 14x, a, 7x, F11.3, 1x, a)') 'PDAF', 'compute C^T OmegaT (34):', pdaf_time_tot(34), 's'
             WRITE (*, '(a, 14x, a, 9x, F11.3, 1x, a)') 'PDAF', 'complete weights (35):', pdaf_time_tot(35), 's'
             WRITE (*, '(a, 12x, a, 7x, F11.3, 1x, a)') 'PDAF', 'transform ensemble (21):', pdaf_time_tot(21), 's'
             IF (dim_lag >0) &
                  WRITE (*, '(a, 20x, a, F11.3, 1x, a)') 'PDAF', 'perform smoothing (15):', pdaf_time_tot(15), 's'

          ELSEIF(subtype_filter==4 .OR. subtype_filter==7) THEN
             WRITE (*, '(a, 10x, a)') 'PDAF', 'Timers in LESTKF only'
             WRITE (*, '(a, 12x, a, 7x, F11.3, 1x, a)') 'PDAF', 'global preparations (7):', pdaf_time_tot(7), 's'
             WRITE (*, '(a, 14x, a, 15x, F11.3, 1x, a)') 'PDAF', 'init Omega (33):', pdaf_time_tot(33), 's'
             WRITE (*, '(a, 12x, a, 7x, F11.3, 1x, a)') 'PDAF', 'local analysis loop (8):', pdaf_time_tot(8), 's'
             WRITE (*, '(a, 14x, a, 10x, F11.3, 1x, a)') 'PDAF', 'global to local (10):', pdaf_time_tot(10), 's'
             WRITE (*, '(a, 14x, a, 4x, F11.3, 1x, a)') 'PDAF', 'localize observations (11):', pdaf_time_tot(11), 's'
             WRITE (*, '(a, 16x, a, 6x, F11.3, 1x, a)') 'PDAF', 'init_dim_obs_l_pdaf (38):', pdaf_time_tot(38), 's'
             WRITE (*, '(a, 16x, a, 13x, F11.3, 1x, a)') 'PDAF', 'g2l_obs_pdaf (46):', pdaf_time_tot(46), 's'
             WRITE (*, '(a, 16x, a, 15x, F11.3, 1x, a)') 'PDAF', 'init_obs_l (47):', pdaf_time_tot(47), 's'
             WRITE (*, '(a, 14x, a, 11x, F11.3, 1x, a)') 'PDAF', 'local analysis (12):', pdaf_time_tot(12), 's'
             WRITE (*, '(a, 16x, a, 5x, F11.3, 1x, a)') 'PDAF', 'compute inverse of A (16):', pdaf_time_tot(16), 's'
             WRITE (*, '(a, 18x, a, 10x, F11.3, 1x, a)') 'PDAF', 'init innovation (20):', pdaf_time_tot(20), 's'
             WRITE (*, '(a, 18x, a, 12x, F11.3, 1x, a)') 'PDAF', 'complete Ainv (21):', pdaf_time_tot(21), 's'
             WRITE (*, '(a, 16x, a, 1x, F11.3, 1x, a)') 'PDAF', 'compute ensemble weights (17):', pdaf_time_tot(17), 's'
             WRITE (*, '(a, 18x, a, 2x, F11.3, 1x, a)') 'PDAF', 'get state weight vector (22):', pdaf_time_tot(22), 's'
             WRITE (*, '(a, 18x, a, F11.3, 1x, a)') 'PDAF', 'complete transform matrix (23):', pdaf_time_tot(23), 's'
             WRITE (*, '(a, 16x, a, 7x, F11.3, 1x, a)') 'PDAF', 'transform ensemble (18):', pdaf_time_tot(18), 's'
             WRITE (*, '(a, 14x, a, 10x, F11.3, 1x, a)') 'PDAF', 'local to global (14):', pdaf_time_tot(14), 's'
             IF (dim_lag >0) &
                  WRITE (*, '(a, 14x, a, 8x, F11.3, 1x, a)') 'PDAF', 'perform smoothing (15):', pdaf_time_tot(15), 's'

          END IF
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

  END SUBROUTINE PDAF_3dvar_memtime

END MODULE PDAF_3DVAR
