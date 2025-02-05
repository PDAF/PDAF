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
!> Module for LNETF holding shared parameters and some helper routines
!!
!! This module declares the parameters that are used in LNETF. 
!! Parameters that are specific for LNETF are declared while some
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
MODULE PDAF_LNETF

  USE PDAF_mod_filter, &
       ONLY: filterstr, incremental, debug, localfilter, dim_lag, &
       member_save

  IMPLICIT NONE

! *** Integer parameters ***
  INTEGER :: type_forget=0 !< Type of forgetting factor
                           !< (0) inflate forecast ensemble
                           !< (1) inflate forecast ensemble only observed domains
                           !< (2) inflate analysis ensemble
                           !< (3) inflate analysis ensemble only observed domains
  INTEGER :: type_trans=0  !< Type of ensemble transformation
                           !< For LNETF:
                           !< (0) use product with random orthonomal matrix with
                           !<     eigenvector (1,...,1)^T
                           !< (1) use deterministic transformation 
  INTEGER :: type_noise=0  !< Type of perturbing noise in PF
                           !< (0) no noise added
                           !< (1) constant variance
                           !< (2) amplitude relative to ensemble std.
  INTEGER :: type_winf=0   !< Type of weights inflation for LNETF
                           !< (0) none
                           !< (1) inflate for N_eff/N > limit_winf

! *** Real parameters ***
  REAL :: forget=1.0       !< Forgetting factor
  REAL :: noise_amp=0.0    !< Amplitude of noise perturbing particles
  REAL :: limit_winf=0.0   !< Limit to weights inflation


! *** Internal variable ***
  LOGICAL :: inloop=.false. ! Whether the program is in the local analysis loop

!-------------------------------------------------------------------------------
  
CONTAINS

!> PDAF_lnetf_set_iparam --- Set integer parameter specific for LNETF
!!
!! LNETF-specific initialization of integer parameter
!!
!!    ! This is a core routine of PDAF and !
!!    ! should not be changed by the user  !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
  SUBROUTINE PDAF_lnetf_set_iparam(id, value, flag)

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
       dim_lag = value
       IF (dim_lag<0) THEN
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(8): Invalid value for smoother lag - param_int(3)!'
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
       IF (type_forget<0 .OR. type_forget>3) THEN
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
       type_winf = value
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
       WRITE (*,'(/5x, a/)') &
            'PDAF-ERROR(8): Invalid integer parameter index'
       flag = 8
    END SELECT

  END SUBROUTINE PDAF_lnetf_set_iparam


!-------------------------------------------------------------------------------
!> PDAF_lnetf_set_rparam --- Set real parameter specific for LNETF
!!
!! LNETF-specific initialization of real parameter
!!
!!    ! This is a core routine of PDAF and !
!!    ! should not be changed by the user  !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
  SUBROUTINE PDAF_lnetf_set_rparam(id, value, flag)

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
       IF (forget <= 0.0) THEN
          WRITE (*,'(/5x,a/)') &
               'PDAF-ERROR(7): Invalid value of forgetting factor - param_real(1)!'
          flag = 7
       END IF
    CASE(2)
       limit_winf = value
       IF (limit_winf < 0.0) THEN
          WRITE (*,'(/5x,a/)') &
               'PDAF-ERROR(8): Invalid limit for weight inflation - param_real(2)!'
          flag = 8
       END IF
    CASE(3)
       noise_amp = value
       IF(noise_amp<0.0) THEN
          WRITE (*,'(/5x,a/)') &
               'PDAF-ERROR(8): Noise amplitude cannot be negative - param_real(3)!'
          flag = 8
       END IF
    CASE DEFAULT
       WRITE (*,'(/5x, a/)') &
            'PDAF-ERROR(10): Invalid real parameter index'
       flag = 10
    END SELECT

  END SUBROUTINE PDAF_lnetf_set_rparam

END MODULE PDAF_LNETF
