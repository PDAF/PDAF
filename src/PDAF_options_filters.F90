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
!
!> Interface routine for information output
!!
!! This subroutine builds the interface for calling
!! the screen output routine for the overview of
!! options for the selected filter.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2011-08 - Lars Nerger - Initial code
!! * Later revisions - see svn log
!!
SUBROUTINE PDAF_options_filters(type_filter)

  USE mpi
  USE PDAF_mod_filtermpi, &
       ONLY: mype_world, MPIerr, COMM_pdaf
  USE PDAF_seek, ONLY: PDAF_seek_options
  USE PDAF_seik, ONLY: PDAF_seik_options
  USE PDAF_lseik, ONLY: PDAF_lseik_options
  USE PDAF_enkf, ONLY: PDAF_enkf_options
  USE PDAF_lenkf, ONLY: PDAF_lenkf_options
  USE PDAF_estkf, ONLY: PDAF_estkf_options
  USE PDAF_lestkf, ONLY: PDAF_lestkf_options
  USE PDAF_etkf, ONLY: PDAF_etkf_options
  USE PDAF_letkf, ONLY: PDAF_letkf_options
  USE PDAF_netf, ONLY: PDAF_netf_options
  USE PDAF_lnetf, ONLY: PDAF_lnetf_options
  USE PDAF_lknetf, ONLY: PDAF_lknetf_options
  USE PDAF_pf, ONLY: PDAF_pf_options
  USE PDAF_genobs, ONLY: PDAF_genobs_options
  USE PDAF_3dvar, ONLY: PDAF_3dvar_options

  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: type_filter     !< Type of filter


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
     ELSE IF (type_filter == 11) THEN
        CALL PDAF_lknetf_options()
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
