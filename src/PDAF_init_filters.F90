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
! !ROUTINE: PDAF_init_filters --- internal interface to filter initializations
!
! !INTERFACE:
SUBROUTINE PDAF_init_filters(type_filter, subtype, param_int, dim_pint, param_real, &
     dim_preal, filterstr, ensemblefilter, fixedbasis, screen, flag)

! !DESCRIPTION:
! Interface routine to the filter-specific initialization
! routines.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2010-08 - Lars Nerger - Initial code for restructuring PDAF
! Later revisions - see svn log
!
! !USES:
  USE mpi
  USE PDAF_mod_filtermpi, &
       ONLY: MPIerr, COMM_pdaf, mype_world

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: type_filter     ! Type of filter
  INTEGER, INTENT(in) :: subtype         ! Sub-type of filter
  INTEGER, INTENT(in) :: dim_pint        ! Number of integer parameters
  INTEGER, INTENT(inout) :: param_int(dim_pint) ! Integer parameter array
  INTEGER, INTENT(in) :: dim_preal       ! Number of real parameters 
  REAL, INTENT(inout) :: param_real(dim_preal)  ! Real parameter array
  CHARACTER(len=10), INTENT(out) :: filterstr   ! Name of filter algorithm
  LOGICAL, INTENT(out) :: ensemblefilter ! Is the chosen filter ensemble-based?
  LOGICAL, INTENT(out) :: fixedbasis     ! Does the filter run with fixed error-space basis?
  INTEGER, INTENT(in)  ::  screen        ! Control screen output
  INTEGER, INTENT(inout):: flag          ! Status flag

! !CALLING SEQUENCE:
! Called by: PDAF_init
! Calls: PDAF_seek_init
! Calls: PDAF_seik_init
! Calls: PDAF_enkf_init
! Calls: PDAF_lseik_init
! Calls: PDAF_etkf_init
! Calls: PDAF_letkf_init
! Calls: PDAF_estkf_init
! Calls: PDAF_lestkf_init
!EOP

! *** local variables ***
  INTEGER :: verbose      ! Control verbosity of info routine


! ***************************
! *** Set writing process ***  
! ***************************

  IF (screen > 0) THEN
        ! Define a single process that writes the information
     CALL MPI_Comm_rank(COMM_pdaf, mype_world, MPIerr)

     IF (mype_world == 0) THEN
        verbose = 1
     ELSE
        verbose = 0
     END IF
  ELSE
     verbose = 0
  END IF


! *******************************************************
! *** 1. Initialize the string identifying the filter ***
! *** 2. Call filter-specific initialization routine  ***
! *******************************************************



  checkflag: IF (flag == 0) THEN
  
     IF (verbose == 1) WRITE (*, '(/a)') 'PDAF: Initialize filter'

     IF (type_filter == 0) THEN

        filterstr = 'SEEK'

        CALL PDAF_seek_init(subtype, param_int, dim_pint, param_real, dim_preal, &
             ensemblefilter, fixedbasis, verbose, flag)

     ELSE IF (type_filter == 1) THEN

        filterstr = 'SEIK'

        CALL PDAF_seik_init(subtype, param_int, dim_pint, param_real, dim_preal, &
             ensemblefilter, fixedbasis, verbose, flag)

     ELSE IF (type_filter == 2) THEN

        filterstr = 'ENKF'

        CALL PDAF_enkf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
             ensemblefilter, fixedbasis, verbose, flag)

     ELSE IF (type_filter == 3) THEN

        filterstr = 'LSEIK'

        CALL PDAF_lseik_init(subtype, param_int, dim_pint, param_real, dim_preal, &
             ensemblefilter, fixedbasis, verbose, flag)

     ELSE IF (type_filter == 4) THEN

        filterstr = 'ETKF'

        CALL PDAF_etkf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
             ensemblefilter, fixedbasis, verbose, flag)

     ELSE IF (type_filter == 5) THEN

        filterstr = 'LETKF'

        CALL PDAF_letkf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
             ensemblefilter, fixedbasis, verbose, flag)

     ELSE IF (type_filter == 6) THEN

        filterstr = 'ESTKF'

        CALL PDAF_estkf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
             ensemblefilter, fixedbasis, verbose, flag)

     ELSE IF (type_filter == 7) THEN

        filterstr = 'LESTKF'

        CALL PDAF_lestkf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
             ensemblefilter, fixedbasis, verbose, flag)

     ELSE IF (type_filter == 8) THEN

        filterstr = 'LENKF'

        CALL PDAF_lenkf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
             ensemblefilter, fixedbasis, verbose, flag)
     ELSE IF (type_filter == 9) THEN

        filterstr = 'NETF'

        CALL PDAF_netf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
             ensemblefilter, fixedbasis, verbose, flag)
     ELSE IF (type_filter == 10) THEN

        filterstr = 'LNETF'

        CALL PDAF_lnetf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
             ensemblefilter, fixedbasis, verbose, flag)

     ELSE IF (type_filter == 11) THEN

        WRITE (*,'(/5x,a/)') 'PDAF-ERROR(1): No valid filter type specified! GENOBS is filtertype=100 from PDAF V2'
        flag = 1        

     ELSE IF (type_filter == 12) THEN

        filterstr = 'PF'

        CALL PDAF_pf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
             ensemblefilter, fixedbasis, verbose, flag)

     ELSE IF (type_filter == 100) THEN

        filterstr = 'GENOBS'

        CALL PDAF_genobs_init(subtype, param_int, dim_pint, param_real, dim_preal, &
             ensemblefilter, fixedbasis, verbose, flag)

     ELSE IF (type_filter == 200) THEN

        filterstr = '3DVAR'

        CALL PDAF_3dvar_init(subtype, param_int, dim_pint, param_real, dim_preal, &
             ensemblefilter, fixedbasis, verbose, flag)
     ELSE

        WRITE (*,'(/5x,a/)') 'PDAF-ERROR(1): No valid filter type specified!'
        flag = 1

     ENDIF
  END IF checkflag

END SUBROUTINE PDAF_init_filters
