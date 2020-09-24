! Copyright (c) 2004-2020 Lars Nerger
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
!$Id: PDAF-D_set_smootherens.F90 374 2020-02-26 12:49:56Z lnerger $
!BOP
!
! !ROUTINE: PDAF_set_timeavg --- Activate time-averaging for ensemble
SUBROUTINE PDAF_set_timeavg(averagetime, status)

! !DESCRIPTION:
! This routine is use to set or unset whether the analysis step
! of uses a time-averaged ensembe to compute the assimilation
! increment (thus using the cross-covariance matrix mean(X')*HX'.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2020-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_mod_filter, &
       ONLY: timeAvg
  USE PDAF_mod_filtermpi, &
       ONLY: mype_world

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)        :: averagetime  ! >0 to activate averaging; =0 to deactivate 
  INTEGER, INTENT(out)       :: status  ! Status flag, 
                                        ! 0: no error, 1: maxlag too large

! !CALLING SEQUENCE:
! Called by: PDAF ensemble initialization routine
!EOP

  
! *******************
! *** Set pointer ***
! *******************

  IF (averagetime>0) THEN
     write (*,*) 'PDAF  timeavg: activate, mype_world', mype_world
     timeAvg = .true.
  ELSE
     write (*,*) 'PDAF  timeavg: deactivate, mype_world', mype_world
     timeAvg = .false.
  END IF

  status = 0

END SUBROUTINE PDAF_set_timeavg
