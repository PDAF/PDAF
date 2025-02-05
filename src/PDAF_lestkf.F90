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
!> Module for LESTKF holding shared parameters and some helper routines
!!
!! This module declares the parameters that are used in LESTKF. 
!! Parameters that are specific for LESTKF are declared while some
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
MODULE PDAF_lestkf

  USE PDAF_mod_filter, &
       ONLY: filterstr, incremental, debug, localfilter, dim_lag, &
       member_save

  IMPLICIT NONE

! *** Integer parameters ***
  INTEGER :: type_forget=0 !< Type of forgetting factor
                           !< (0): fixed; (1) global adaptive; (2) local adaptive
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
  REAL    :: forget_l      !< Forgetting factor in local analysis loop

! *** Internal variable ***
  LOGICAL :: inloop=.false. ! Whether the program is in the local analysis loop


!$OMP THREADPRIVATE(forget_l)

!-------------------------------------------------------------------------------
  
CONTAINS

!> PDAF_lestkf_set_iparam --- Set integer parameter specific for LESTKF filter
!!
!! LESTKF-specific initialization of integer parameters
!!
!!    ! This is a core routine of PDAF and !
!!    ! should not be changed by the user  !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
  SUBROUTINE PDAF_lestkf_set_iparam(id, value, flag)

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
       incremental = value
       IF (incremental /= 0) THEN
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(10): LESTKF does not yet support incremental updating!'
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
       type_sqrt = value
       IF (type_sqrt<0 .OR. type_sqrt>1) THEN
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(10): Invalid setting for square root type - param_int(7)!'
          flag = 8
       END IF
    CASE(8)
       if (value==0) THEN
          observe_ens = .false. ! Apply H to ensemble mean to compute residual
       ELSE
          observe_ens = .true.  ! Apply H to X, compute mean of HX and then residual (the default for LESTKF in PDAF2.3)
       END IF
    CASE(9)
       type_obs_init = value    ! Initialize obs (0) before or (1) after prepoststep
       IF (type_obs_init<0 .OR. type_obs_init>1) THEN
          WRITE (*,'(/5x, a/)') &
               'PDAF-ERROR(10): Invalid setting type_obs_init - param_int(9)!'
          flag = 8
       END IF
    CASE DEFAULT
       WRITE (*,'(/5x, a/)') &
            'PDAF-ERROR(10): Invalid integer parameter index'
       flag = 10
    END SELECT

  END SUBROUTINE PDAF_lestkf_set_iparam


!-------------------------------------------------------------------------------
!> PDAF_lestkf_set_rparam --- Set real parameter specific for LESTKF filter
!!
!! LESTKF-specific initialization of real parameters
!!
!!    ! This is a core routine of PDAF and
!!      should not be changed by the user  !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
  SUBROUTINE PDAF_lestkf_set_rparam(id, value, flag)

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
       IF (localfilter == 0) THEN
          forget = value
       ELSE
          IF (inloop) THEN
             forget_l = value
          ELSE
             forget = value
write (*,*) 'set forget=', forget
          END IF
       END IF
       IF (forget <= 0.0) THEN
          WRITE (*,'(/5x,a/)') &
               'PDAF-ERROR(7): Invalid value of forgetting factor - param_real(1)!'
          flag = 7
       END IF
    CASE DEFAULT
       WRITE (*,'(/5x, a/)') &
            'PDAF-ERROR(10): Invalid real parameter index'
       flag = 10
    END SELECT

  END SUBROUTINE PDAF_lestkf_set_rparam

END MODULE PDAF_LESTKF
