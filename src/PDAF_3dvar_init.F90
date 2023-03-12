! Copyright (c) 2004-2021 Lars Nerger
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
!$Id$
!BOP
!
! !ROUTINE: PDAF_3dvar_init --- PDAF-internal initialization of 3DVAR
!
! !INTERFACE:
SUBROUTINE PDAF_3dvar_init(subtype, param_int, dim_pint, param_real, dim_preal, &
     ensemblefilter, fixedbasis, verbose, outflag)

! !DESCRIPTION:
! Initialization of 3DVAR within PDAF. Performed are:\\
!   - initialize filter-specific parameters\\
!   - print screen information on filter configuration.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2021-03 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filter, &
       ONLY: incremental, dim_ens, type_opt, dim_cvec, dim_cvec_ens, &
       beta_3dvar, localfilter, forget, &
       type_forget, dim_bias_p, type_trans, dim_lag

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: subtype                ! Sub-type of filter
  INTEGER, INTENT(in) :: dim_pint               ! Number of integer parameters
  INTEGER, INTENT(inout) :: param_int(dim_pint) ! Integer parameter array
  INTEGER, INTENT(in) :: dim_preal              ! Number of real parameters 
  REAL, INTENT(inout) :: param_real(dim_preal)  ! Real parameter array
  LOGICAL, INTENT(out) :: ensemblefilter        ! Is the chosen filter ensemble-based?
  LOGICAL, INTENT(out) :: fixedbasis            ! Does the filter run with fixed error-space basis?
  INTEGER, INTENT(in) :: verbose                ! Control screen output
  INTEGER, INTENT(inout):: outflag              ! Status flag

! !CALLING SEQUENCE:
! Called by: PDAF_init_filters
!EOP

! *** local variables ***


! ****************************
! *** INITIALIZE VARIABLES ***
! ****************************

  dim_lag = 0
  type_forget = 0
  incremental = 0
  type_trans = 0
  dim_bias_p = 0

  ! Define whether filter is domain localized
  localfilter = 0

  IF (subtype==0 .AND. dim_ens > 1) THEN
     WRITE (*, '(/5x, a/)') 'PDAF-ERROR(6): 3D-Var must be run with ensemble size = 1!'
     outflag = 6
  END IF

  ! Initialize variable to prevent compiler warning
  forget = param_real(1)

  ! choice of optimizer
  IF (dim_pint>=3) THEN
     type_opt = param_int(3)
  END IF

  IF (dim_pint>=4) THEN
     dim_cvec = param_int(4)
  ELSE
     IF (subtype==0 .OR. subtype==4 .OR. subtype==6 .OR. subtype==7) THEN
        WRITE (*, '(/5x, a/)') 'PDAF-ERROR(3): Missing specification of control vector dimension!'
        outflag = 3
     END IF
  END IF

  IF (dim_pint>=5) THEN
     dim_cvec_ens = param_int(5)
  ELSE
     IF (subtype==1 .OR. subtype==4) THEN
        dim_cvec_ens = dim_ens
     END IF
  END IF

  IF (dim_preal>=2) THEN
     beta_3dvar = param_real(2)
  END IF

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

  filter_pe2: IF (verbose > 0) THEN
  
     WRITE(*, '(/a, 4x, a)') 'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++                      3D-Var                     +++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++                                                 +++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++      3D-Var variants implemented following      +++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++      Bannister, Q. J. Royal Meteorol. Soc.,     +++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++     143 (2017) 607-633, doi:10.1002/qj.2982     +++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'

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
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> offline mode - choose method by the PDAF_assimilate/_put_state routine'
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
!      IF (type_forget == 0) THEN
!         WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> Use fixed forgetting factor:', forget
     IF (subtype == 1 .OR. subtype == 2) &
          WRITE (*, '(a, 12x, a, i5)') 'PDAF', '--> ensemble size:', dim_ens

  END IF filter_pe2

END SUBROUTINE PDAF_3dvar_init
