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
       ONLY: dim_ens, localfilter, rank, dim_lag
  USE PDAF_3dvar, &
       ONLY: incremental, type_opt, dim_cvec, dim_cvec_ens, &
       beta_3dvar, forget, type_forget, type_trans, eps_cg_var, &
       maxiter_cg_var, eps_cgplus_var, method_cgplus_var, irest_cgplus_var, &
       factr_lbfgs_var, pgtol_lbfgs_var, m_lbfgs_var, &
       PDAF_3dvar_set_iparam, PDAF_3dvar_set_rparam
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

!   ! Initialize variable to prevent compiler warning
!   forget = param_real(1)
! 
!   ! choice of optimizer
!   IF (dim_pint>=3) THEN
!      type_opt = param_int(3)
!      IF (type_opt==0) THEN
!         WRITE (*, '(/5x, a/)') 'PDAF-ERROR(4): Incorrect choice of solver!'
!         outflag = 4
!      END IF
!   END IF

  ! Some special conditions
  IF (dim_pint<4) THEN
     IF (subtype==0 .OR. subtype==4 .OR. subtype==6 .OR. subtype==7) THEN
        WRITE (*, '(/5x, a/)') 'PDAF-ERROR(3): Missing specification of control vector dimension!'
        outflag = 3
     END IF
  END IF

  IF (dim_pint<5) THEN
     IF (subtype==1 .OR. subtype==4) THEN
        dim_cvec_ens = dim_ens
     END IF
  END IF

!   IF (dim_pint>=6) THEN
!      m_lbfgs_var = param_int(6)
!      method_cgplus_var = param_int(6)
!      maxiter_cg_var = param_int(6)
!   END IF
! 
!   IF (dim_pint>=7) THEN
!      irest_cgplus_var = param_int(7)
!   END IF
! 
!   IF (dim_preal>=2) THEN
!      beta_3dvar = param_real(2)
!   END IF

!   IF (dim_preal>=3) THEN
!      eps_cg_var = param_real(3)
!      eps_cgplus_var = param_real(3)
!      pgtol_lbfgs_var = param_real(3)
!   END IF
! 
!   IF (dim_preal>=4) THEN
!      factr_lbfgs_var = param_real(4)
!   END IF

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
