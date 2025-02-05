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
!>  PDAF-internal initialization of SEEK filter
!!
!! Initialization of SEEK within PDAF. Performed are:
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
SUBROUTINE PDAF_seek_init(subtype, param_int, dim_pint, param_real, dim_preal, &
     ensemblefilter, fixedbasis, verbose, outflag)

  USE PDAF_mod_filter, &
       ONLY: dim_ens, localfilter, dim_eof, offline_mode
  USE PDAF_seek, &
       ONLY: incremental, forget, int_rediag, epsilon, &
       PDAF_seek_set_iparam, PDAF_seek_set_rparam
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
  int_rediag = 1
  epsilon = 0.1

  ! Parse provided parameters
  flagsum = 0
  DO i=3, dim_pint
     CALL PDAF_seek_set_iparam(i, param_int(i), outflag)
     flagsum = flagsum+outflag
  END DO
  DO i=1, dim_preal
     CALL PDAF_seek_set_rparam(i, param_real(i), outflag)
     flagsum = flagsum+outflag
  END DO


  ! For fixed basis SEEK do not perform rediagonalization
  IF (subtype == 2 .OR. subtype == 3) THEN
     int_rediag = 0
  END IF

  ! Special for SEEK: Initialize number of modes
  dim_eof = dim_ens

  ! Define whether filter is a domain-local filter
  localfilter = 0

  ! Define whether filter is mode-based or ensemble-based
  ensemblefilter = .FALSE.

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
  
     WRITE(*, '(/a, 5x, a)') 'PDAF', '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     WRITE(*, '(a, 5x, a)')  'PDAF', '+++                  SEEK Filter                   +++'
     WRITE(*, '(a, 5x, a)')  'PDAF', '+++                                                +++'
     WRITE(*, '(a, 5x, a)')  'PDAF', '+++    Pham et al., J. Mar. Syst. 16 (1998) 323    +++'
     WRITE(*, '(a, 5x, a)')  'PDAF', '+++          This implementation follows           +++'
     WRITE(*, '(a, 5x, a)')  'PDAF', '+++      Nerger et al., Tellus 57A (2005) 715      +++'
     WRITE(*, '(a, 5x, a)')  'PDAF', '+++                                                +++'     
     WRITE(*, '(a, 5x, a)')  'PDAF', '+++ NOTE: The SEEK filter in PDAF is deprecated    +++'     
     WRITE(*, '(a, 5x, a)')  'PDAF', '+++       as of Version 1.14. It will be removed   +++'
     WRITE(*, '(a, 5x, a)')  'PDAF', '+++       in the future.                           +++'     
     WRITE(*, '(a, 5x, a)')  'PDAF', '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     IF (flagsum== 0 ) THEN

        ! *** General output ***
        WRITE (*, '(/a, 4x, a)') 'PDAF', 'SEEK configuration'
        WRITE (*, '(a, 10x, a, i1)') 'PDAF', 'filter sub-type = ', subtype
        IF (subtype == 0) THEN
           WRITE (*, '(a, 12x, a)') 'PDAF', '--> Standard SEEK with unit modes'
        ELSE IF (subtype == 1) THEN
           WRITE (*, '(a, 12x, a)') 'PDAF', '--> SEEK with non-unit modes'
        ELSE IF (subtype == 2) THEN
           WRITE (*, '(a, 12x, a)') 'PDAF', '--> fixed basis filter with update of matrix U'
           WRITE (*, '(a, 12x, a)') 'PDAF', '--> no re-diagonalization of VUV^T'
        ELSE IF (subtype == 3) THEN
           WRITE (*, '(a, 12x, a)') 'PDAF', '--> fixed basis filter & no update of matrix U'
           WRITE (*, '(a, 12x, a)') 'PDAF', '--> no re-diagonalization of VUV^T'
        ELSE
           WRITE (*, '(/5x, a/)') 'PDAF-ERROR(3): No valid subtype!'
           outflag = 3
        END IF
        IF (incremental == 1) &
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> Perform incremental updating'
        IF (.not.offline_mode) THEN
           IF ((int_rediag > 0) .AND. ((subtype /= 2) .OR. (subtype /= 3))) THEN
              IF (int_rediag == 1) THEN
                 WRITE (*, '(a, 10x, a, i4, a)') 'PDAF', 'Re-diag at each analysis step'
              ELSE
                 WRITE (*, '(a, 10x, a, i4, a)') 'PDAF', 'Re-diag at each ', int_rediag, &
                      '-th analysis step'
              END IF
           END IF
        ELSE
           IF (int_rediag == 1) THEN
              WRITE (*, '(a, 5x, a)') 'PDAF', 'Perform re-diagonalization'
           ELSE
              WRITE (*, '(a, 5x, a)') 'PDAF', 'No re-diagonalization'
           END IF
        END IF
        WRITE (*, '(a, 12x, a, i5)') 'PDAF', '--> number of EOFs:', dim_eof
     ELSE
        WRITE (*, '(/5x, a/)') 'PDAF-ERROR: Invalid parameter setting - check prior output!'
     END IF

  END IF writeout

END SUBROUTINE PDAF_seek_init
