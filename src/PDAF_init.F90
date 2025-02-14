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
!> Initialize PDAF
!!
!! Initialization of PDAF. Performed are:\\
!!   * Initialization of filter independent parameters\\
!!   * Call to filter-specific routine for parameter initialization\\
!!   * Initialization of PDAF-internal parallelization\\
!!   * Call to filter-specific routine for allocation of arrays\\
!!   * Call to user-routine for ensemble/mode initialization.
!!
!! !  This is a core routine of PDAF and
!!    should not be changed by the user   !
!!
!! __Revision history:__
!! * 2003-08 - Lars Nerger - Initial code
!! *  Other revisions - see repository log
!!
SUBROUTINE PDAF_init(filtertype, subtype, stepnull, param_int, dim_pint, &
     param_real, dim_preal, COMM_model, COMM_filter, COMM_couple, &
     task_id, n_modeltasks, in_filterpe, U_init_ens, in_screen, &
     outflag)

  USE mpi
  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount_ini
  USE PDAF_mod_filter, &
       ONLY: dim_ens, dim_eof, dim_p, flag, &
       screen, step, step_obs, type_filter, filterstr, &
       subtype_filter, ensemblefilter, state, Ainv, ens, &
       debug, offline_mode
  USE PDAF_mod_filtermpi, &
       ONLY: mype, filterpe, PDAF_init_parallel, COMM_pdaf, &
       isset_comm_pdaf
  USE PDAF_utils_filters, &
       ONLY: PDAF_init_filters, PDAF_alloc_filters, PDAF_options_filters

  IMPLICIT NONE

! *** Arguments ***
  ! For valid and default values see PDAF_mod_filter.F90
  INTEGER, INTENT(in) :: filtertype     !< Type of filter
  INTEGER, INTENT(in) :: subtype        !< Sub-type of filter
  INTEGER, INTENT(in) :: stepnull       !< Initial time step of assimilation
  INTEGER, INTENT(in) :: dim_pint       !< Number of integer parameters
  INTEGER, INTENT(inout) :: param_int(dim_pint) !< Integer parameter array
  INTEGER, INTENT(in) :: dim_preal      !< Number of real parameter 
  REAL, INTENT(inout) :: param_real(dim_preal) !< Real parameter array
  INTEGER, INTENT(in) :: COMM_model     !< Model communicator
  INTEGER, INTENT(in) :: COMM_couple    !< Coupling communicator
  INTEGER, INTENT(in) :: COMM_filter    !< Filter communicator
  INTEGER, INTENT(in) :: task_id        !< Id of my ensemble task
  INTEGER, INTENT(in) :: n_modeltasks   !< Number of parallel model tasks
  LOGICAL, INTENT(in) :: in_filterpe    !< Is my PE a filter-PE?
  INTEGER, INTENT(in) :: in_screen      !< Control screen output:
                                        !< (0) none, (1) some, default, (2) extensive
  INTEGER, INTENT(out):: outflag        !< Status flag, 0: no error, error codes:
                                   !< -1: Call with subtype=-1 for info display
                                   !<  1: No valid filter type
                                   !<  2: No valid sub type
                                   !<  3: Invalid dim_pint
                                   !<  4: Invalid dim_preal
                                   !<  5: Invalid state dimension
                                   !<  6: Invalid ensemble size
                                   !<  7: Invalid value for forgetting factor
                                   !<  8: Invalid other integer parameter value
                                   !<  9: Invalid other real parameter value
                                   !< 10: Invalid setting for incremental updating
                                   !< 20: error in allocation of array at PDAF init

! *** External subroutines ***
! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_init_ens  !< User-supplied routine for ensemble initialization

! *** local variables ***
  INTEGER :: i                     ! Counter
  LOGICAL :: fixedbasis            ! Does the filter run with fixed error-space basis?


! ********************************************
! *** INITIALIZE VARIABLES FOR ALL FILTERS ***
! ********************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_init -- START'

  ! set number of timers
  CALL PDAF_timeit(66, 'ini')

  ! Initialize memory counters
  CALL PDAF_memcount_ini(4)

  ! Call timer
  CALL PDAF_timeit(1, 'new')

  ! Set PDAF communicator if not set externally
  IF (.NOT. isset_comm_pdaf) THEN
     COMM_pdaf = MPI_COMM_WORLD

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_init:', debug, 'Use MPI_COMM_WORLD for COMM_PDAF'
  ELSE
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_init:', debug, 'Use user-defined communicator for COMM_PDAF'
  END IF

  ! Print version information
  CALL PDAF_print_version()

  info_assim: IF (subtype < 0) THEN

     ! ***********************************************************************
     ! *** For negative subtype only display information on filter options ***
     ! ***********************************************************************

     CALL PDAF_options_filters(filtertype)

     ! Set status flag
     flag = -1

  ELSE info_assim

     ! *** Check size of parameter arrays
     IF (dim_pint < 2) THEN
        WRITE (*,'(/5x,a/)') &
             'PDAF-ERROR(3): Invalid size of array of integer parameters!'
        flag = 3
     END IF
     IF (dim_preal < 1) THEN
        WRITE (*,'(/5x,a/)') &
             'PDAF-ERROR(4): Invalid size of array of real parameters!'
        flag = 4
     END IF

     ! *** Initialize variables for all PEs
     type_filter    = filtertype   ! Set filter type
     subtype_filter = subtype      ! Set sub-type of filter
     step_obs       = stepnull     ! Set initial time step
     step           = step_obs + 1 ! stepping index
     filterpe       = in_filterpe  ! Whether my PE is a PE of the filter
     screen         = in_screen    ! Control verbosity
     if (subtype==5) offline_mode = .true.   ! Set offline mode

     dim_p          = param_int(1) ! PE-local state dimension
     IF (param_int(1) < 1) THEN
        WRITE (*,'(/5x,a/)') &
             'PDAF-ERROR(5): Invalid state dimension!'
        flag = 5
     END IF

     ! Ensemble size
     dim_ens = param_int(2)
     IF (param_int(2) < 1) THEN
        WRITE (*,'(/5x,a/)') 'PDAF-ERROR(6): Invalid ensemble size!'
        flag = 6
     END IF

     IF (debug>0 .AND. flag==0) THEN
       WRITE (*,*) '++ PDAF-debug PDAF_init:', debug, 'param_int of size', dim_pint, &
       'values:', param_int(1:dim_pint)
       WRITE (*,*) '++ PDAF-debug PDAF_init:', debug, 'param_real of size', dim_preal, &
       'values:', param_real(1:dim_preal)
       WRITE (*,*) '++ PDAF-debug PDAF_init:', debug, &
          '  Note: If REAL values appear incorrect, please check if you provide them with the correct precision'
    END IF


! ********************************************
! *** Initialize filter-specific variables ***
! ********************************************

     IF (flag == 0) THEN
        CALL PDAF_init_filters(type_filter, subtype_filter, param_int, dim_pint, param_real, &
             dim_preal, filterstr, ensemblefilter, fixedbasis, screen, flag)
     END IF

  
! **********************************
! *** Initialize parallelization ***
! **********************************

     IF (flag == 0) THEN
        CALL PDAF_init_parallel(dim_ens, ensemblefilter, fixedbasis, &
             COMM_model, COMM_filter, COMM_couple, &
             n_modeltasks, task_id, screen, flag)
     END IF


! ********************************************
! *** Filter-specific allocation of arrays ***
! *** and screen output                    ***
! ********************************************

     IF (flag == 0) THEN
        CALL PDAF_alloc_filters(filterstr, subtype_filter, flag)
     END IF


! **********************************
! *** Initialize ensemble matrix ***
! **********************************

     filter_pe3: IF (filterpe .AND. flag == 0) THEN

        IF (mype == 0 .AND. screen > 0) &
             WRITE (*, '(/a)') 'PDAF: Call ensemble initialization'

        CALL PDAF_timeit(39, 'new')

        typef: IF (ensemblefilter) THEN
           ! *** Initialize ensemble of ensemble-based filter      ***
           ! *** EnKF/SEIK/LSEIK/ETKF/LETKF                        ***
           CALL U_init_ens(type_filter, dim_p, dim_ens, state, Ainv, &
                ens, flag)

           IF (debug>0) THEN
              DO i = 1, dim_ens
                 WRITE (*,*) '++ PDAF-debug PDAF_init:', debug, 'ensemble member', i, &
                      ' values (1:min(dim_p,6)):', ens(1:min(dim_p,6),i)
              END DO
           END IF
        ELSE
           ! *** Mode-based filter (SEEK)                          ***
           ! *** Initialize rank reduced covariance matrix         ***
           ! *** factors Ainv and ens and estimated initial state ***
           CALL U_init_ens(type_filter, dim_p, dim_eof, state, Ainv, &
                ens, flag)

           IF (debug>0) THEN
              DO i = 1, dim_ens
                 WRITE (*,*) '++ PDAF-debug PDAF_init:', debug, 'covar mode', i, &
                      ' values (1:min(dim_p,6)):', ens(1:min(dim_p,6),i)
              END DO
              WRITE (*,*) '++ PDAF-debug PDAF_init:', debug, 'mode weights (1:min(dim_eof,10)):', &
                   Ainv(1:min(dim_eof, 10),1:min(dim_eof, 10))
           END IF
        END IF typef

        CALL PDAF_timeit(39, 'old')

     END IF filter_pe3

  END IF info_assim


! ********************
! *** FINISHING UP ***
! ********************

  ! Store internal status flag
  outflag = flag

  CALL PDAF_timeit(1, 'old')

  IF (mype == 0 .AND. filterpe .AND. screen > 0) &
       WRITE (*, '(/a)') 'PDAF: Initialization completed'
  IF (mype == 0 .AND. filterpe .AND. screen > 1) &
       WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
       'PDAF', '--- duration of PDAF initialization:', PDAF_time_temp(1), 's'

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_init -- END'

END SUBROUTINE PDAF_init
