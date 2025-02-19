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
!> Collection of routines interfacing to method-specifc routines
!! 
!! These routines are interfaces to the DA-method specific
!! routines. The are called by PDAF_init, PDAF_print_info,
!! PDAF_set_iparam or PDAF_set_rparam.
!!
!! !  These are core routines of PDAF and
!!    should not be changed by the user   !
!! 
MODULE PDAF_utils_filters

CONTAINS

!> Interface routine to the method-specific initialization routines.
!!
!! __Revision history:__
!! * 2010-08 - Lars Nerger - Initial code for restructuring PDAF
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_init_filters(type_filter, subtype, param_int, dim_pint, param_real, &
       dim_preal, filterstr, ensemblefilter, fixedbasis, screen, flag)

    USE mpi
    USE PDAF_mod_filtermpi, &
         ONLY: MPIerr, COMM_pdaf, mype_world
    USE PDAF_DA
    USE PDAF_seik, ONLY: PDAF_seik_init
    USE PDAF_lseik, ONLY: PDAF_lseik_init
    USE PDAF_enkf, ONLY: PDAF_enkf_init
    USE PDAF_lenkf, ONLY: PDAF_lenkf_init
    USE PDAF_etkf, ONLY: PDAF_etkf_init
    USE PDAF_letkf, ONLY: PDAF_letkf_init
    USE PDAF_estkf, ONLY: PDAF_estkf_init
    USE PDAF_lestkf, ONLY: PDAF_lestkf_init
    USE PDAF_netf, ONLY: PDAF_netf_init
    USE PDAF_lnetf, ONLY: PDAF_lnetf_init
    USE PDAF_lknetf, ONLY: PDAF_lknetf_init
    USE PDAF_pf, ONLY: PDAF_pf_init
    USE PDAF_genobs, ONLY: PDAF_genobs_init
    USE PDAF_3dvar, ONLY: PDAF_3dvar_init

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: type_filter            !< Type of filter
    INTEGER, INTENT(inout) :: subtype             !< Sub-type of filter
    INTEGER, INTENT(in) :: dim_pint               !< Number of integer parameters
    INTEGER, INTENT(inout) :: param_int(dim_pint) !< Integer parameter array
    INTEGER, INTENT(in) :: dim_preal              !< Number of real parameters 
    REAL, INTENT(inout) :: param_real(dim_preal)  !< Real parameter array
    CHARACTER(len=10), INTENT(out) :: filterstr   !< Name of filter algorithm
    LOGICAL, INTENT(out) :: ensemblefilter        !< Is the chosen filter ensemble-based?
    LOGICAL, INTENT(out) :: fixedbasis            !< Does the filter run with fixed error-space basis?
    INTEGER, INTENT(in)  ::  screen               !< Control screen output
    INTEGER, INTENT(inout):: flag                 !< Status flag

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

       IF (type_filter == PDAF_DA_SEIK) THEN

          filterstr = 'SEIK'

          CALL PDAF_seik_init(subtype, param_int, dim_pint, param_real, dim_preal, &
               ensemblefilter, fixedbasis, verbose, flag)

       ELSE IF (type_filter == PDAF_DA_ENKF) THEN

          filterstr = 'ENKF'

          CALL PDAF_enkf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
               ensemblefilter, fixedbasis, verbose, flag)

       ELSE IF (type_filter == PDAF_DA_LSEIK) THEN

          filterstr = 'LSEIK'

          CALL PDAF_lseik_init(subtype, param_int, dim_pint, param_real, dim_preal, &
               ensemblefilter, fixedbasis, verbose, flag)

       ELSE IF (type_filter == PDAF_DA_ETKF) THEN

          filterstr = 'ETKF'

          CALL PDAF_etkf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
               ensemblefilter, fixedbasis, verbose, flag)

       ELSE IF (type_filter == PDAF_DA_LETKF) THEN

          filterstr = 'LETKF'

          CALL PDAF_letkf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
               ensemblefilter, fixedbasis, verbose, flag)

       ELSE IF (type_filter == PDAF_DA_ESTKF) THEN

          filterstr = 'ESTKF'

          CALL PDAF_estkf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
               ensemblefilter, fixedbasis, verbose, flag)

       ELSE IF (type_filter == PDAF_DA_LESTKF) THEN

          filterstr = 'LESTKF'

          CALL PDAF_lestkf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
               ensemblefilter, fixedbasis, verbose, flag)

       ELSE IF (type_filter == PDAF_DA_LENKF) THEN

          filterstr = 'LENKF'

          CALL PDAF_lenkf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
               ensemblefilter, fixedbasis, verbose, flag)
       ELSE IF (type_filter == PDAF_DA_NETF) THEN

          filterstr = 'NETF'

          CALL PDAF_netf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
               ensemblefilter, fixedbasis, verbose, flag)
       ELSE IF (type_filter == PDAF_DA_LNETF) THEN

          filterstr = 'LNETF'

          CALL PDAF_lnetf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
               ensemblefilter, fixedbasis, verbose, flag)

       ELSE IF (type_filter == PDAF_DA_LKNETF) THEN

          filterstr = 'LKNETF'

          CALL PDAF_lknetf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
               ensemblefilter, fixedbasis, verbose, flag)

       ELSE IF (type_filter == PDAF_DA_PF) THEN

          filterstr = 'PF'

          CALL PDAF_pf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
               ensemblefilter, fixedbasis, verbose, flag)

       ELSE IF (type_filter == PDAF_DA_GENOBS) THEN

          filterstr = 'GENOBS'

          CALL PDAF_genobs_init(subtype, param_int, dim_pint, param_real, dim_preal, &
               ensemblefilter, fixedbasis, verbose, flag)

       ELSE IF (type_filter == PDAF_DA_3DVAR) THEN

          filterstr = '3DVAR'

          CALL PDAF_3dvar_init(subtype, param_int, dim_pint, param_real, dim_preal, &
               ensemblefilter, fixedbasis, verbose, flag)
       ELSE

          WRITE (*,'(/5x,a/)') 'PDAF-ERROR(1): No valid filter type specified!'
          flag = 1

       ENDIF
    END IF checkflag

  END SUBROUTINE PDAF_init_filters

!-------------------------------------------------------------------------------
!> Interface routine to the method-specific allocation routines.
!!
!! __Revision history:__
!! 2010-08 - Lars Nerger - Initial code from restructuring PDAF
!! Other revisions - see repository log
!!
  SUBROUTINE PDAF_alloc_filters(filterstr, subtype, flag)

    USE PDAF_seik, ONLY: PDAF_seik_alloc
    USE PDAF_lseik, ONLY: PDAF_lseik_alloc
    USE PDAF_enkf, ONLY: PDAF_enkf_alloc
    USE PDAF_etkf, ONLY: PDAF_etkf_alloc
    USE PDAF_letkf, ONLY: PDAF_letkf_alloc
    USE PDAF_estkf, ONLY: PDAF_estkf_alloc
    USE PDAF_lestkf, ONLY: PDAF_lestkf_alloc
    USE PDAF_lenkf, ONLY: PDAF_lenkf_alloc
    USE PDAF_netf, ONLY: PDAF_netf_alloc
    USE PDAF_lnetf, ONLY: PDAF_lnetf_alloc
    USE PDAF_lknetf, ONLY: PDAF_lknetf_alloc
    USE PDAF_pf, ONLY: PDAF_pf_alloc
    USE PDAF_genobs, ONLY: PDAF_genobs_alloc
    USE PDAF_3dvar, ONLY: PDAF_3dvar_alloc

    IMPLICIT NONE

! *** Arguments ***
    CHARACTER(len=10), INTENT(in) :: filterstr ! Name of filter algorithm
    INTEGER, INTENT(in) :: subtype             ! Sub-type of filter
    INTEGER, INTENT(inout)::flag               ! Status flag


! ***********************************************
! *** Call filter-specific allocation routine ***
! ***********************************************

    checkflag: IF (flag == 0) THEN
       IF (TRIM(filterstr) == 'SEIK') THEN
          CALL PDAF_seik_alloc(flag)
       ELSE IF (TRIM(filterstr) == 'ENKF') THEN
          CALL PDAF_enkf_alloc(flag)
       ELSE IF (TRIM(filterstr) == 'LSEIK') THEN
          CALL PDAF_lseik_alloc(flag)
       ELSE IF (TRIM(filterstr) == 'ETKF') THEN
          CALL PDAF_etkf_alloc(flag)
       ELSE IF (TRIM(filterstr) == 'LETKF') THEN
          CALL PDAF_letkf_alloc(flag)
       ELSE IF (TRIM(filterstr) == 'ESTKF') THEN
          CALL PDAF_estkf_alloc(flag)
       ELSE IF (TRIM(filterstr) == 'LESTKF') THEN
          CALL PDAF_lestkf_alloc(flag)
       ELSE IF (TRIM(filterstr) == 'LENKF') THEN
          CALL PDAF_lenkf_alloc(flag)
       ELSE IF (TRIM(filterstr) == 'NETF') THEN
          CALL PDAF_netf_alloc(flag)
       ELSE IF (TRIM(filterstr) == 'LNETF') THEN
          CALL PDAF_lnetf_alloc(flag)
       ELSE IF (TRIM(filterstr) == 'LKNETF') THEN
          CALL PDAF_lknetf_alloc(flag)
       ELSE IF (TRIM(filterstr) == 'PF') THEN
          CALL PDAF_pf_alloc(flag)
       ELSE IF (TRIM(filterstr) == 'GENOBS') THEN
          CALL PDAF_genobs_alloc(flag)
       ELSE IF (TRIM(filterstr) == '3DVAR') THEN
          CALL PDAF_3dvar_alloc(subtype, flag)
       ENDIF
    END IF checkflag

  END SUBROUTINE PDAF_alloc_filters


!-------------------------------------------------------------------------------
!> Print configuration info of the active filter
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_configinfo_filters(subtype, verbose)

    USE PDAF_mod_filter, ONLY: filterstr
    USE PDAF_seik, ONLY: PDAF_seik_config
    USE PDAF_lseik, ONLY: PDAF_lseik_config
    USE PDAF_enkf, ONLY: PDAF_enkf_config
    USE PDAF_lenkf, ONLY: PDAF_lenkf_config
    USE PDAF_estkf, ONLY: PDAF_estkf_config
    USE PDAF_lestkf, ONLY: PDAF_lestkf_config
    USE PDAF_etkf, ONLY: PDAF_etkf_config
    USE PDAF_letkf, ONLY: PDAF_letkf_config
    USE PDAF_netf, ONLY: PDAF_netf_config
    USE PDAF_lnetf, ONLY: PDAF_lnetf_config
    USE PDAF_lknetf, ONLY: PDAF_lknetf_config
    USE PDAF_pf, ONLY: PDAF_pf_config
    USE PDAF_genobs, ONLY: PDAF_genobs_config
    USE PDAF_3dvar, ONLY: PDAF_3dvar_config

    IMPLICIT NONE

! *** Arguments ***

    INTEGER, INTENT(inout) :: subtype               !< Sub-type of filter
    INTEGER, INTENT(in)    :: verbose               !< Control screen output


! ********************************
! *** Print screen information ***
! ********************************

    IF (TRIM(filterstr) == 'SEIK') THEN
       CALL PDAF_seik_config(subtype, verbose)
    ELSE IF (TRIM(filterstr) == 'LSEIK') THEN
       CALL PDAF_lseik_config(subtype, verbose)
    ELSE IF (TRIM(filterstr) == 'ENKF') THEN
       CALL PDAF_enkf_config(subtype, verbose)
    ELSE IF (TRIM(filterstr) == 'LENKF') THEN
       CALL PDAF_lenkf_config(subtype, verbose)
    ELSE IF (TRIM(filterstr) == 'ETKF') THEN
       CALL PDAF_etkf_config(subtype, verbose)
    ELSE IF (TRIM(filterstr) == 'LETKF') THEN
       CALL PDAF_letkf_config(subtype, verbose)
    ELSE IF (TRIM(filterstr) == 'ESTKF') THEN
       CALL PDAF_estkf_config(subtype, verbose)
    ELSE IF (TRIM(filterstr) == 'LESTKF') THEN
       CALL PDAF_lestkf_config(subtype, verbose)
    ELSE IF (TRIM(filterstr) == 'NETF') THEN
       CALL PDAF_netf_config(subtype, verbose)
    ELSE IF (TRIM(filterstr) == 'LNETF') THEN
       CALL PDAF_lnetf_config(subtype, verbose)
    ELSE IF (TRIM(filterstr) == 'LKNETF') THEN
       CALL PDAF_lknetf_config(subtype, verbose)
    ELSE IF (TRIM(filterstr) == 'PF') THEN
       CALL PDAF_pf_config(subtype, verbose)
    ELSE IF (TRIM(filterstr) == 'GENOBS') THEN
       CALL PDAF_genobs_config(subtype, verbose)
    ELSE IF (TRIM(filterstr) == '3DVAR') THEN
       CALL PDAF_3dvar_config(subtype, verbose)
     END IF

   END SUBROUTINE PDAF_configinfo_filters

!-------------------------------------------------------------------------------
!> Interface routine for information output
!!
!! This subroutine builds the interface for calling
!! the screen output routine for the overview of
!! options for the selected DA method.
!!
!! __Revision history:__
!! * 2011-08 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_options_filters(type_filter)

    USE mpi
    USE PDAF_mod_filtermpi, &
         ONLY: mype_world, MPIerr, COMM_pdaf
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

       IF (type_filter == 1) THEN
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


!-------------------------------------------------------------------------------
!> Print information for PDAF (timing and memory) to screen
!!
!! This routine displays the information from PDAF.
!! Possible are to display the timing information and
!! allocated memory.
!!
!! __Revision history:__
!! * 2008-09 - Lars Nerger - Initial code
!! * Other revisions - see repository log
!!
  SUBROUTINE PDAF_print_info_filters(printtype)

    USE PDAF_mod_filter, ONLY: filterstr
    USE PDAF_seik, ONLY: PDAF_seik_memtime
    USE PDAF_lseik, ONLY: PDAF_lseik_memtime
    USE PDAF_enkf, ONLY: PDAF_enkf_memtime
    USE PDAF_lenkf, ONLY: PDAF_lenkf_memtime
    USE PDAF_estkf, ONLY: PDAF_estkf_memtime
    USE PDAF_lestkf, ONLY: PDAF_lestkf_memtime
    USE PDAF_etkf, ONLY: PDAF_etkf_memtime
    USE PDAF_letkf, ONLY: PDAF_letkf_memtime
    USE PDAF_netf, ONLY: PDAF_netf_memtime
    USE PDAF_lnetf, ONLY: PDAF_lnetf_memtime
    USE PDAF_lknetf, ONLY: PDAF_lknetf_memtime
    USE PDAF_pf, ONLY: PDAF_pf_memtime
    USE PDAF_3dvar, ONLY: PDAF_3dvar_memtime

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: printtype     !< Type of screen output:  
                                         !< (1) general timings
                                         !< (3) timers focused on call-back routines (recommended)
                                         !< (4,5) detailed and very detailed timers (to analyze filters)
                                         !< (10) allocated memory of calling MPI task
                                         !< (11) globally used memory (needs to e called by all processes)


! ********************************
! *** Print screen information ***
! ********************************

    IF (TRIM(filterstr) == 'SEIK') THEN
       CALL PDAF_seik_memtime(printtype)
    ELSE IF (TRIM(filterstr) == 'ENKF') THEN
       CALL PDAF_enkf_memtime(printtype)
    ELSE IF (TRIM(filterstr) == 'LSEIK') THEN
       CALL PDAF_lseik_memtime(printtype)
    ELSE IF (TRIM(filterstr) == 'ETKF') THEN
       CALL PDAF_etkf_memtime(printtype)
    ELSE IF (TRIM(filterstr) == 'LETKF') THEN
       CALL PDAF_letkf_memtime(printtype)
    ELSE IF (TRIM(filterstr) == 'ESTKF') THEN
       CALL PDAF_estkf_memtime(printtype)
    ELSE IF (TRIM(filterstr) == 'LESTKF') THEN
       CALL PDAF_lestkf_memtime(printtype)
    ELSE IF (TRIM(filterstr) == 'LENKF') THEN
       CALL PDAF_lenkf_memtime(printtype)
    ELSE IF (TRIM(filterstr) == 'NETF') THEN
       CALL PDAF_netf_memtime(printtype)
    ELSE IF (TRIM(filterstr) == 'LNETF') THEN
       CALL PDAF_lnetf_memtime(printtype)
    ELSE IF (TRIM(filterstr) == 'LKNETF') THEN
       CALL PDAF_lknetf_memtime(printtype)
    ELSE IF (TRIM(filterstr) == 'PF') THEN
       CALL PDAF_pf_memtime(printtype)
    ELSE IF (TRIM(filterstr) == '3DVAR') THEN
       CALL PDAF_3dvar_memtime(printtype)
    END IF

  END SUBROUTINE PDAF_print_info_filters

!-------------------------------------------------------------------------------
!> PDAF_set_iparam_filters --- Set integer parameter for PDAF
!!
!! Routine allowing the user to set the value of a
!! method-specific integer parameter. This routine
!! builds the interface to the specific routine
!! provided by each DA method.
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Other revisions - see repository log
!!
  SUBROUTINE PDAF_set_iparam_filters(id, value, flag)

    USE PDAF_mod_filter, &
         ONLY: filterstr
    USE PDAF_seik, &
         ONLY: PDAF_seik_set_iparam
    USE PDAF_enkf, &
         ONLY: PDAF_enkf_set_iparam
    USE PDAF_lseik, &
         ONLY: PDAF_lseik_set_iparam
    USE PDAF_etkf, &
         ONLY: PDAF_etkf_set_iparam
    USE PDAF_letkf, &
         ONLY: PDAF_letkf_set_iparam
    USE PDAF_estkf, &
         ONLY: PDAF_estkf_set_iparam
    USE PDAF_lestkf, &
         ONLY: PDAF_lestkf_set_iparam
    USE PDAF_lenkf, &
         ONLY: PDAF_lenkf_set_iparam
    USE PDAF_netf, &
         ONLY: PDAF_netf_set_iparam
    USE PDAF_lnetf, &
         ONLY: PDAF_lnetf_set_iparam
    USE PDAF_lknetf, &
         ONLY: PDAF_lknetf_set_iparam
    USE PDAF_pf, &
         ONLY: PDAF_pf_set_iparam
    USE PDAF_3dvar, &
         ONLY: PDAF_3dvar_set_iparam
    USE PDAF_genobs, &
         ONLY: PDAF_genobs_set_iparam

    IMPLICIT NONE

! Arguments
    INTEGER, INTENT(in) :: id       !< Index of parameter
    INTEGER, INTENT(in) :: value    !< Parameter value
    INTEGER, INTENT(out) :: flag    !< Status flag: 0 for no error


! ********************************
! *** Print screen information ***
! ********************************

    IF (TRIM(filterstr) == 'SEIK') THEN
       CALL PDAF_seik_set_iparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'ENKF') THEN
       CALL PDAF_enkf_set_iparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'LSEIK') THEN
       CALL PDAF_lseik_set_iparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'ETKF') THEN
       CALL PDAF_etkf_set_iparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'LETKF') THEN
       CALL PDAF_letkf_set_iparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'ESTKF') THEN
       CALL PDAF_estkf_set_iparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'LESTKF') THEN
       CALL PDAF_lestkf_set_iparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'LENKF') THEN
       CALL PDAF_lenkf_set_iparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'NETF') THEN
       CALL PDAF_netf_set_iparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'LNETF') THEN
       CALL PDAF_lnetf_set_iparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'LKNETF') THEN
       CALL PDAF_lknetf_set_iparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'PF') THEN
       CALL PDAF_pf_set_iparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == '3DVAR') THEN
       CALL PDAF_3dvar_set_iparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'GENOBS') THEN
       CALL PDAF_genobs_set_iparam(id, value, flag)
    ELSE
       WRITE (*,*) 'PDAF-ERROR: invalid DA method - likely PDAF is not yet initialized' 
    END IF

  END SUBROUTINE PDAF_set_iparam_filters

!-------------------------------------------------------------------------------
!> PDAF_set_rparam --- Set real parameter for PDAF
!!
!! This routine lets the user set the value of a
!! method-specific real parameter. This routine
!! builds the interface to the specific routine
!! provided by each DA method.
!!
!! __Revision history:__
!! * 2025-02 - Lars Nerger - Initial code
!! *  Other revisions - see repository log
!!
  SUBROUTINE PDAF_set_rparam_filters(id, value, flag)

    USE PDAF_mod_filter, &
         ONLY: filterstr
    USE PDAF_seik, &
         ONLY: PDAF_seik_set_rparam
    USE PDAF_enkf, &
         ONLY: PDAF_enkf_set_rparam
    USE PDAF_lseik, &
         ONLY: PDAF_lseik_set_rparam
    USE PDAF_etkf, &
         ONLY: PDAF_etkf_set_rparam
    USE PDAF_letkf, &
         ONLY: PDAF_letkf_set_rparam
    USE PDAF_estkf, &
         ONLY: PDAF_estkf_set_rparam
    USE PDAF_lestkf, &
         ONLY: PDAF_lestkf_set_rparam
    USE PDAF_lenkf, &
         ONLY: PDAF_lenkf_set_rparam
    USE PDAF_netf, &
         ONLY: PDAF_netf_set_rparam
    USE PDAF_lnetf, &
         ONLY: PDAF_lnetf_set_rparam
    USE PDAF_lknetf, &
         ONLY: PDAF_lknetf_set_rparam
    USE PDAF_pf, &
         ONLY: PDAF_pf_set_rparam
    USE PDAF_3dvar, &
         ONLY: PDAF_3dvar_set_rparam

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)  :: id       !< Index of parameter
    REAL, INTENT(in)     :: value    !< Parameter value
    INTEGER, INTENT(out) :: flag     !< Status flag: 0 for no error


! ********************************
! *** Print screen information ***
! ********************************

    IF (TRIM(filterstr) == 'SEIK') THEN
       CALL PDAF_seik_set_rparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'ENKF') THEN
       CALL PDAF_enkf_set_rparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'LSEIK') THEN
       CALL PDAF_lseik_set_rparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'ETKF') THEN
       CALL PDAF_etkf_set_rparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'LETKF') THEN
       CALL PDAF_letkf_set_rparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'ESTKF') THEN
       CALL PDAF_estkf_set_rparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'LESTKF') THEN
       CALL PDAF_lestkf_set_rparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'LENKF') THEN
       CALL PDAF_lenkf_set_rparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'NETF') THEN
       CALL PDAF_netf_set_rparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'LNETF') THEN
       CALL PDAF_lnetf_set_rparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'LKNETF') THEN
       CALL PDAF_lknetf_set_rparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'PF') THEN
       CALL PDAF_pf_set_rparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == '3DVAR') THEN
       CALL PDAF_3dvar_set_rparam(id, value, flag)
    ELSE IF (TRIM(filterstr) == 'GENOBS') THEN
       ! There are no real parameters in GENOBS
    ELSE
       WRITE (*,*) 'PDAF-ERROR: invalid DA method - likely PDAF is not yet initialized' 
    END IF

  END SUBROUTINE PDAF_set_rparam_filters

END MODULE PDAF_utils_filters
