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
! !ROUTINE: PDAF_options_filters --- Interface routine for information output
!
! !INTERFACE:
SUBROUTINE PDAF_options_filters(type_filter)

! !DESCRIPTION:
! This subroutine builds the interface for calling
! the screen output routine for the overview of
! options for the selected filter.

! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2011-08 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mpi
  USE PDAF_mod_filtermpi, &
       ONLY: mype_world, MPIerr, COMM_pdaf

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: type_filter     ! Type of filter

! !CALLING SEQUENCE:
! Called by: PDAF_init
!EOP

  ! Determine parallel rank of process
  CALL MPI_Comm_rank(COMM_pdaf, mype_world, MPIerr)


! *** Call output routine for specified filter type
  IF (mype_world==0) THEN
     ! Output on process 0 only

     IF (type_filter == 0) THEN
        CALL PDAF_seek_options()
     ELSE IF (type_filter == 1) THEN
        CALL PDAF_seik_options()
     ELSE IF (type_filter == 2) THEN
        CALL PDAF_enkf_options()
     ELSE IF (type_filter == 3) THEN
        CALL PDAF_lseik_options()
     ELSE IF (type_filter == 4) THEN
        CALL PDAF_etkf_options()
     ELSE IF (type_filter == 5) THEN
        CALL PDAF_letkf_options()
     ELSE IF (type_filter == 6) THEN
        CALL PDAF_estkf_options()
     ELSE IF (type_filter == 7) THEN
        CALL PDAF_lestkf_options()
     ELSE IF (type_filter == 8) THEN
        CALL PDAF_lenkf_options()
     ELSE IF (type_filter == 9) THEN
        CALL PDAF_netf_options()
     ELSE IF (type_filter == 10) THEN
        CALL PDAF_lnetf_options()
     ELSE IF (type_filter == 12) THEN
        CALL PDAF_pf_options()
     ELSE IF (type_filter == 100) THEN
        CALL PDAF_genobs_options()
     ELSE IF (type_filter == 200) THEN
        CALL PDAF_3dvar_options()
     ELSE
        WRITE (*,'(a, 5x, a)') 'PDAF', 'No options overview available for the selected filter!'
     END IF

  END IF

END SUBROUTINE PDAF_options_filters
