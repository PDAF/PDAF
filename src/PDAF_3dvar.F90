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
!! * 2025-02 - Lars Nerger - Initial code
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

  INTEGER :: type_opt=0           !< Type of minimizer for 3DVar
                                  !< (0) LBFGS, (1) CG+, (-1) steepest descent
  INTEGER :: dim_cvec=0           !< Size of control vector (fixed part)
  INTEGER :: dim_cvec_ens=0       !< Size of control vector (ensemble part)
  INTEGER :: m_lbfgs_var=5        !< Parameter 'm' of LBFGS
  INTEGER :: method_cgplus_var=2  !< Parameter 'method' of CG+
  INTEGER :: irest_cgplus_var=1   !< Parameter 'irest' of CG+
  INTEGER :: maxiter_cg_var=200   !< Parameter 'maxiter' of CG

! *** Real parameters ***
  REAL :: beta_3dvar=0.5           !< Hybrid weight for hybrid 3D-Var
  REAL :: eps_cg_var = 1.0e-6      !< Parameter 'EPS' of  CG
  REAL :: eps_cgplus_var = 1.0e-5  !< Parameter 'EPS' of CG+
  REAL :: pgtol_lbfgs_var=1.0e-5   !< Parameter 'pgtol' of LBFGS
  REAL :: factr_lbfgs_var=1.0e7    !< Parameter 'factr' of LBFGS


!-------------------------------------------------------------------------------
  
CONTAINS

!> PDAF_3dvar_set_iparam --- Set integer parameter specific for 3DVAR
!!
!! 3DVAR-specific initialization of integer parameter
!!
!!    ! This is a core routine of PDAF and !
!!    ! should not be changed by the user  !
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

write (*,*) 'set_iparam: id', id,' value', value

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
       WRITE (*,'(/5x, a/)') &
            'PDAF-ERROR(8): Invalid integer parameter index'
       flag = 8
    END SELECT

  END SUBROUTINE PDAF_3dvar_set_iparam


!-------------------------------------------------------------------------------
!> PDAF_3dvar_set_rparam --- Set real parameter specific for 3DVAR
!!
!! 3DVAR-specific initialization of real parameter
!!
!!    ! This is a core routine of PDAF and !
!!    ! should not be changed by the user  !
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

    write (*,*) 'set_rparam: id', id,' value', value

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
       WRITE (*,'(/5x, a/)') &
            'PDAF-ERROR(10): Invalid real parameter index'
       flag = 10
    END SELECT

  END SUBROUTINE PDAF_3dvar_set_rparam

END MODULE PDAF_3DVAR
