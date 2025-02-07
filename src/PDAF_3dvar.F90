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
!! *  Later revisions - see repository log
!!
MODULE PDAF_3DVAR

  USE PDAF_mod_filter, &
       ONLY: filterstr, incremental, debug, localfilter, dim_lag
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
!! *  Later revisions - see repository log
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
    INTEGER :: flagsum          ! Sum of status flags


! ****************************
! *** INITIALIZE VARIABLES ***
! ****************************

    incremental = 0
    observe_ens = .false.

    ! Settings for variational part
    type_opt = 0
    dim_cvec = 0
    dim_cvec_ens = 0
    m_lbfgs_var = 5
    method_cgplus_var = 2
    irest_cgplus_var = 1
    maxiter_cg_var = 200
    beta_3dvar = 0.5
    eps_cg_var = 1.0e-6
    eps_cgplus_var = 1.0e-5
    pgtol_lbfgs_var = 1.0e-5
    factr_lbfgs_var  =1.0e7

    ! Settings for ensemble filter
    type_forget = 0
    type_trans = 0
    dim_lag = 0
    forget = 1.0
  

    ! Parse provided parameters
    flagsum = 0
    DO i=3, dim_pint
       CALL PDAF_3dvar_set_iparam(i, param_int(i), outflag)
       flagsum = flagsum+outflag
    END DO
    DO i=1, dim_preal
       CALL PDAF_3dvar_set_rparam(i, param_real(i), outflag)
       flagsum = flagsum+outflag
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

       IF (flagsum== 0 ) THEN

          ! *** General output ***
          WRITE (*, '(/a, 4x, a)') 'PDAF', '3DVAR configuration'
          WRITE (*, '(a, 9x, a, i1)') 'PDAF', 'filter sub-type = ', subtype
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
          ELSE
             WRITE (*, '(/5x, a/)') 'PDAF-ERROR(2): No valid sub type!'
             outflag = 2
          END IF
          IF (incremental == 1) &
               WRITE (*, '(a, 12x, a)') 'PDAF', '--> Perform incremental updating'

          IF (subtype>0) THEN
             IF (type_trans == 0) THEN
                WRITE (*, '(a, 12x, a)') 'PDAF', '--> Deterministic ensemble transformation'
             ELSE IF (type_trans == 1) THEN
                WRITE (*, '(a, 12x, a)') 'PDAF', '--> Transform ensemble with random orthonormal Omega'
             ELSE IF (type_trans == 2) THEN
                WRITE (*, '(a, 12x, a)') 'PDAF', '--> Transform ensemble including product with random matrix'
             END IF
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
             WRITE (*, '(a, 12x, a, i5)') 'PDAF', '--> ensemble size:', dim_ens
             IF (observe_ens) &
                  WRITE (*, '(a, 12x, a, 1x, l)') 'PDAF', '--> observe_ens:', observe_ens
          END IF
       ELSE
          WRITE (*, '(/5x, a/)') 'PDAF-ERROR: Invalid parameter setting - check prior output!'
       END IF

    END IF writeout

  END SUBROUTINE PDAF_3dvar_init


!-------------------------------------------------------------------------------
!> Perform allocation of arrays for 3DVAR.
!!
!! __Revision history:__
!! * 2019-05 - Lars Nerger - Initial code based on NETF
!! * 2025-02 - Lars Nerger - Restructuring introducing generic PDAF_alloc
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAF_3dvar_alloc(subtype, outflag)

    USE PDAF_mod_filter, &
         ONLY: dim_ens, dim_p, dim_bias_p
    USE PDAF_mod_filtermpi, &
         ONLY: dim_ens_l

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
       dim_ens = dim_ens-1
    END IF

    CALL PDAF_alloc(dim_p, dim_ens, dim_ens_l, dim_ens, dim_bias_p, &
         dim_lag, 0, incremental, outflag)

  END SUBROUTINE PDAF_3dvar_alloc


!-------------------------------------------------------------------------------
!> Set integer parameter specific for 3DVAR
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
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
       IF (type_opt<1 .OR. type_opt>3) THEN
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
!! *  Later revisions - see repository log
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

END MODULE PDAF_3DVAR
