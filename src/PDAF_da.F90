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
!> Module providing shared variables for DA-methods
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
MODULE PDAF_DA
  
  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: PDAF_DA_SEIK = 1
  INTEGER, PARAMETER :: PDAF_DA_LSEIK = 3
  INTEGER, PARAMETER :: PDAF_DA_ENKF = 2
  INTEGER, PARAMETER :: PDAF_DA_LENKF = 8
  INTEGER, PARAMETER :: PDAF_DA_ETKF = 4
  INTEGER, PARAMETER :: PDAF_DA_LETKF = 5
  INTEGER, PARAMETER :: PDAF_DA_ESTKF = 6
  INTEGER, PARAMETER :: PDAF_DA_LESTKF = 7
  INTEGER, PARAMETER :: PDAF_DA_NETF = 9
  INTEGER, PARAMETER :: PDAF_DA_LNETF = 10
  INTEGER, PARAMETER :: PDAF_DA_LKNETF = 11
  INTEGER, PARAMETER :: PDAF_DA_PF = 12
  INTEGER, PARAMETER :: PDAF_DA_ENSRF = 13
  INTEGER, PARAMETER :: PDAF_DA_GENOBS = 100
  INTEGER, PARAMETER :: PDAF_DA_3DVAR = 200

CONTAINS


!-------------------------------------------------------------------------------
!> Print information on available filter types
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!

  SUBROUTINE PDAF_print_filter_types(verbose)

    IMPLICIT NONE

! *** Argument ***
    INTEGER, INTENT(in) :: verbose

    ! *** Print list of available filter types
    IF (verbose > 0) THEN

       WRITE (*,'(a,5x,a)') 'PDAF','Available named filtertypes'
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_LESTKF', PDAF_DA_LESTKF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_ESTKF ', PDAF_DA_ESTKF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_LETKF ', PDAF_DA_LETKF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_ETKF  ', PDAF_DA_ETKF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_LENKF ', PDAF_DA_LENKF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_ENKF  ', PDAF_DA_ENKF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_LSEIK ', PDAF_DA_LSEIK
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_SEIK  ', PDAF_DA_SEIK
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_ENSRF ', PDAF_DA_ENSRF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_LNETF ', PDAF_DA_LNETF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_NETF  ', PDAF_DA_NETF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_PF    ', PDAF_DA_PF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_LKNETF', PDAF_DA_LKNETF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_GENOBS', PDAF_DA_GENOBS
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_3DVAR ', PDAF_DA_3DVAR
       WRITE (*,'(a,5x,a,i4)') 'PDAF','Note: DA methods can be specified by their value or the named parameter'
       WRITE (*,'(a,5x,a)') 'PDAF','+++++++++ End of list of filtertypes +++++++++'

    END IF

  END SUBROUTINE PDAF_print_filter_types


!-------------------------------------------------------------------------------
!> Print information on available types of DA methods
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2025-03 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!

  SUBROUTINE PDAF_print_DA_types(verbose)

    IMPLICIT NONE

! *** Argument ***
    INTEGER, INTENT(in) :: verbose

    ! *** Print list of available filter types
    IF (verbose > 0) THEN

       WRITE (*,'(a,5x,a)') 'PDAF','Available named DA methods'
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_LESTKF', PDAF_DA_LESTKF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_ESTKF ', PDAF_DA_ESTKF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_LETKF ', PDAF_DA_LETKF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_ETKF  ', PDAF_DA_ETKF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_LENKF ', PDAF_DA_LENKF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_ENKF  ', PDAF_DA_ENKF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_LSEIK ', PDAF_DA_LSEIK
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_SEIK  ', PDAF_DA_SEIK
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_ENSRF ', PDAF_DA_ENSRF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_LNETF ', PDAF_DA_LNETF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_NETF  ', PDAF_DA_NETF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_PF    ', PDAF_DA_PF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_LKNETF', PDAF_DA_LKNETF
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_GENOBS', PDAF_DA_GENOBS
       WRITE (*,'(a,8x,a,i4)') 'PDAF','PDAF_DA_3DVAR ', PDAF_DA_3DVAR
       WRITE (*,'(a,5x,a,i4)') 'PDAF','Note: DA methods can be specified by their value or the named parameter'
       WRITE (*,'(a,5x,a)') 'PDAF','+++++++++ End of list of DA method types +++++++++'

    END IF

  END SUBROUTINE PDAF_print_DA_types

END MODULE PDAF_DA
