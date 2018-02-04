! Copyright (c) 2004-2016 Lars Nerger
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
!$Id: PDAF-D_seek_alloc.F90 1681 2016-12-11 12:43:58Z lnerger $
!BOP
!
! !ROUTINE: PDAF_seek_alloc --- PDAF-internal initialization of SEEK filter
!
! !INTERFACE:
SUBROUTINE PDAF_seek_alloc(subtype, outflag)

! !DESCRIPTION:
! Perform allocation of arrays for SEEK.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2010-08 - Lars Nerger - Initial code from splitting PDAF_seek_init
! Later revisions - see svn log
!
! !USES:
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filter, &
       ONLY: screen, incremental, dim_eof, dim_p, state, &
       state_inc, eofU, eofV
  USE PDAF_mod_filtermpi, &
       ONLY: mype, mype_model, mype_couple, statetask, filterpe, &
       dim_eof_l, task_id, COMM_couple, MPI_COMM_NULL

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: subtype        ! Sub-type of filter
  INTEGER, INTENT(out):: outflag        ! Status flag

! !CALLING SEQUENCE:
! Called by: PDAF_alloc_filters
! Calls: PDAF_memcount
!EOP

! *** local variables ***
  INTEGER :: allocstat    ! Status flag for allocate


! ******************************
! *** Allocate filter fields ***
! ******************************

  on_filterpe: IF (filterpe) THEN
     ! Allocate all arrays and full mode matrix on Filter-PEs

     ALLOCATE(state(dim_p), stat = allocstat)
     IF (allocstat /= 0) THEN
        WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of STATE'
        outflag = 20
     END IF
     ! count allocated memory
     CALL PDAF_memcount(1, 'r', dim_p)

     IF (incremental == 1) THEN
        ALLOCATE(state_inc(dim_p), stat = allocstat)
        IF (allocstat /= 0) THEN
           WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of STATE_INC'
           outflag = 20
        END IF
        ! count allocated memory
        CALL PDAF_memcount(1, 'r', dim_p)

        state_inc = 0.0
     END IF

     ALLOCATE(eofU(dim_eof, dim_eof), stat = allocstat)
     IF (allocstat /= 0) THEN
        WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of eofU'
        outflag = 20
     END IF
     ! count allocated memory
     CALL PDAF_memcount(1, 'r', dim_eof**2)

     ! Allocate full mode matrix on filter-PEs
     ALLOCATE(eofV(dim_p, dim_eof), stat = allocstat)
     IF (allocstat /= 0) THEN
        WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of eofV'
        outflag = 20
     END IF
     ! count allocated memory
     CALL PDAF_memcount(2, 'r', dim_p * dim_eof)

     IF (screen > 2) &
          WRITE (*,*) 'PDAF: seek_alloc - allocate eofV of size ', &
          dim_eof, ' on pe(f) ', mype

  ELSE on_filterpe
     ! Model-PEs that are not Filter-PEs only need an array for the local modes
     ! if they participate in the coupling communication

     ! Allocate partial ensemble on model-only PEs
     IF (COMM_couple /= MPI_COMM_NULL) THEN
        ALLOCATE(eofV(dim_p, dim_eof_l), stat = allocstat)
        IF (allocstat /= 0) THEN
           WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of eofV on model-pe'
           outflag = 20
        END IF
        ! count allocated memory
        CALL PDAF_memcount(2, 'r', dim_p * dim_eof_l)

        IF (screen > 2) WRITE (*,*) 'PDAF: seek_alloc - allocate eofV of size ', &
             dim_eof_l, ' on pe(m) ', mype_model, ' of model task ',task_id

        ! Some of the model-PEs may integrate the central state
        ifSEEK2: IF (mype_couple+1 == statetask) THEN
           ALLOCATE(state(dim_p), stat = allocstat)
           IF (allocstat /= 0) THEN
              WRITE (*,'(5x, a)') 'PDAF-ERROR(20): error in allocation of STATE'
              outflag = 20
           END IF

           IF (screen > 2) WRITE (*,*) 'PDAF: seek_alloc - allocate state of size ', &
                dim_p, ' on pe(m) ', mype_model, ' of model task ',task_id
        ENDIF ifSEEK2
        ! count allocated memory
        CALL PDAF_memcount(2, 'r', dim_p)
     END IF

  END IF on_filterpe

END SUBROUTINE PDAF_seek_alloc
