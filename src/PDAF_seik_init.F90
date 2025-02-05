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
!>  PDAF-internal initialization of SEIK filter
!!
!! Initialization of SEIK within PDAF. Performed are:
!! * initialize filter-specific parameters
!! * print screen information on filter configuration.
!!
!!  !  This is a core routine of PDAF and   !
!!  !   should not be changed by the user   !
!!
!! __Revision history:__
!! * 2003-08 - Lars Nerger - Initial code
!! *  Later revisions - see repository log
!!
SUBROUTINE PDAF_seik_init(subtype, param_int, dim_pint, param_real, dim_preal, &
     ensemblefilter, fixedbasis, verbose, outflag)

  USE PDAF_mod_filter, &
       ONLY: dim_ens, localfilter, rank
  USE PDAF_seik, &
       ONLY: incremental, forget, type_forget, type_trans, type_sqrt, &
       PDAF_seik_set_iparam, PDAF_seik_set_rparam
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
  incremental = 0
  observe_ens = .false.
  forget = 1.0

  ! Parse provided parameters
  flagsum = 0
  DO i=3, dim_pint
     CALL PDAF_seik_set_iparam(i, param_int(i), outflag)
     flagsum = flagsum+outflag
  END DO
  DO i=1, dim_preal
     CALL PDAF_seik_set_rparam(i, param_real(i), outflag)
     flagsum = flagsum+outflag
  END DO

  ! *** Special setting
  IF (subtype==3) type_sqrt = 1 ! For fixed covariance we always use Cholesky decomposition


  ! Rank of initial covariance matrix
  rank = dim_ens - 1

  ! Define whether filter is mode-based or ensemble-based
  ensemblefilter = .TRUE.

  ! Define whether filter is a domain-local filter
  localfilter = 0

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
  
     WRITE(*, '(/a, 4x, a)') 'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++                  SEIK Filter                    +++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++                                                 +++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++ Pham et al., C. R. Acad. Sci. II, 326(1998) 255 +++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++    and Pham, Mon. Wea. Rev. 129 (2001) 1194     +++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++          This implementation follows            +++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++      Nerger et al., Tellus 57A (2005) 715       +++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     IF (flagsum== 0 ) THEN

        ! *** General output ***
        WRITE (*, '(/a, 4x, a)') 'PDAF', 'SEIK configuration'
        WRITE (*, '(a, 10x, a, i1)') 'PDAF', 'filter sub-type = ', subtype
        IF (subtype == 0) THEN
           WRITE (*, '(a, 12x, a)') 'PDAF', '--> Standard SEIK'
        ELSE IF (subtype == 1) THEN
           WRITE (*, '(a, 12x, a)') 'PDAF', '--> Standard SEIK - old formulation'
        ELSE IF (subtype == 2) THEN
           WRITE (*, '(a, 12x, a)') 'PDAF', '--> SEIK with fixed error-space basis'
        ELSE IF (subtype == 3) THEN
           WRITE (*, '(a, 12x, a)') 'PDAF', '--> SEIK with fixed state covariance matrix'
        ELSE IF (subtype == 4) THEN
           WRITE (*, '(a, 12x, a)') 'PDAF', '--> SEIK with ensemble transformation'
        ELSE
           WRITE (*, '(/5x, a/)') 'PDAF-ERROR(3): No valid subtype!'
           outflag = 3
        END IF
        IF (type_trans == 0) THEN
           WRITE (*, '(a, 12x, a)') 'PDAF', '--> Transform ensemble with deterministic Omega'
        ELSE IF (type_trans == 1) THEN
           WRITE (*, '(a, 12x, a)') 'PDAF', '--> Transform ensemble with random orthonormal Omega'
        ELSE IF (type_trans == 2) THEN
           WRITE (*, '(a, 12x, a)') 'PDAF', '--> Transform ensemble with product Omega'
        END IF
        IF (incremental == 1) &
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> Perform incremental updating'
        IF (type_forget == 0) THEN
           WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> Use fixed forgetting factor:', forget
        ELSEIF (type_forget == 1) THEN
           WRITE (*, '( a, 12x, a)') 'PDAF', '--> Use adaptive forgetting factor'
        ENDIF
        WRITE (*, '(a, 12x, a, i5)') 'PDAF', '--> ensemble size:', dim_ens
        IF (observe_ens) &
             WRITE (*, '(a, 12x, a, 1x, l)') 'PDAF', '--> observe_ens:', observe_ens
     ELSE
        WRITE (*, '(/5x, a/)') 'PDAF-ERROR: Invalid parameter setting - check prior output!'
     END IF

  END IF writeout

END SUBROUTINE PDAF_seik_init
