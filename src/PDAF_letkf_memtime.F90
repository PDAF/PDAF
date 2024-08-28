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
!$Id$
!BOP
!
! !ROUTINE: PDAF_letkf_memtime --- Display timing and memory information for LETKF
!
! !INTERFACE:
SUBROUTINE PDAF_letkf_memtime(printtype)

! !DESCRIPTION:
! This routine displays the PDAF-internal timing and
! memory information for the LETKF.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2008-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_timer, &
       ONLY: PDAF_time_tot
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount_get, PDAF_memcount_get_global
  USE PDAF_mod_filter, &
       ONLY: subtype_filter, offline_mode, dim_lag, type_forget
  USE PDAF_mod_filtermpi, &
       ONLY: filterpe, mype_world, COMM_pdaf
  USE PDAFomi, &
       ONLY: omi_was_used
  USE PDAFlocal, &
       ONLY: pdaflocal_was_used

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: printtype    ! Type of screen output:  
                                      ! (1) timings, (2) memory
!EOP

! *** Local variables ***
  INTEGER :: i                        ! Counter
  REAL :: memcount_global(3)          ! Globally counted memory
  REAL :: time_omi                    ! Sum of timers for OMI-internal call-back routines


! ********************************
! *** Print screen information ***
! ********************************

  ptype: IF (printtype == 1) THEN

! **************************************
! *** Print basic timing information ***
! **************************************

     ! Generic part
     WRITE (*, '(//a, 21x, a)') 'PDAF', 'PDAF Timing information'
     WRITE (*, '(a, 10x, 45a)') 'PDAF', ('-', i=1, 45)
     WRITE (*, '(a, 18x, a, F11.3, 1x, a)') &
          'PDAF', 'Initialize PDAF:', pdaf_time_tot(1), 's'
     IF (.not.offline_mode) THEN
        IF (subtype_filter<2) THEN
           WRITE (*, '(a, 16x, a, F11.3, 1x, a)') 'PDAF', 'Ensemble forecast:', pdaf_time_tot(2), 's'
        ELSE
           WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'State forecast:', pdaf_time_tot(2), 's'
        END IF
     END IF

     IF (filterpe) THEN
        ! Filter-specific part
        WRITE (*, '(a, 19x, a, F11.3, 1x, a)') 'PDAF', 'LETKF analysis:', pdaf_time_tot(3), 's'

        ! Generic part B
        WRITE (*, '(a, 22x, a, F11.3, 1x, a)') 'PDAF', 'Prepoststep:', pdaf_time_tot(5), 's'
     END IF

  ELSE IF (printtype == 2) THEN ptype

! *****************************************
! *** Formerly: Print allocated memory  ***
! *****************************************

     WRITE (*, '(/a, 23x, a)') 'PDAF', 'PDAF Memory overview'
     WRITE (*, '(/a, 23x, a)') 'PDAF', 'Note: The memory overview is moved to printtype=10 and printtype=11'

  ELSE IF (printtype == 3) THEN ptype

! *******************************************************
! *** Print timing information for call-back routines ***
! *******************************************************

     ! Generic part
     WRITE (*, '(//a, 12x, a)') 'PDAF', 'PDAF Timing information - call-back routines'
     WRITE (*, '(a, 8x, 52a)') 'PDAF', ('-', i=1, 52)
     WRITE (*, '(a, 10x, a, 15x, F11.3, 1x, a)') 'PDAF', 'Initialize PDAF:', pdaf_time_tot(1), 's'
     WRITE (*, '(a, 12x, a, 17x, F11.3, 1x, a)') 'PDAF', 'init_ens_pdaf:', pdaf_time_tot(39), 's'
     IF (.not.offline_mode) THEN
        IF (subtype_filter<2) THEN
           WRITE (*, '(a, 10x, a, 13x, F11.3, 1x, a)') 'PDAF', 'Ensemble forecast:', pdaf_time_tot(2), 's'
        ELSE
           WRITE (*, '(a, 10x, a, 17x, F11.3, 1x, a)') 'PDAF', 'State forecast:', pdaf_time_tot(2), 's'
        END IF
        WRITE (*, '(a, 12x, a, 5x, F11.3, 1x, a)') 'PDAF', 'MPI communication in PDAF:', pdaf_time_tot(19), 's'
        WRITE (*, '(a, 12x, a, 9x, F11.3, 1x, a)') 'PDAF', 'distribute_state_pdaf:', pdaf_time_tot(40), 's'
        WRITE (*, '(a, 12x, a, 12x, F11.3, 1x, a)') 'PDAF', 'collect_state_pdaf:', pdaf_time_tot(41), 's'
        IF (.not.filterpe) WRITE (*, '(a, 7x, a)') 'PDAF', &
             'Note: for filterpe=F, the time (2) includes the wait time for the analysis step'
     END IF

     IF (filterpe) THEN
        ! Filter-specific part
        WRITE (*, '(a, 10x, a, 16x, F11.3, 1x, a)') 'PDAF', 'LETKF analysis:', pdaf_time_tot(3), 's'
        WRITE (*, '(a, 12x, a, 6x, F11.3, 1x, a)') 'PDAF', 'PDAF-internal operations:', pdaf_time_tot(51), 's'

        IF(omi_was_used) THEN
           ! Output when using OMI

           time_omi = pdaf_time_tot(46) + pdaf_time_tot(47) + pdaf_time_tot(48)
           IF (type_forget==1) &
                time_omi = time_omi + pdaf_time_tot(50) + pdaf_time_tot(49) + pdaf_time_tot(52)
           WRITE (*, '(a, 12x, a, 9x, F11.3, 1x, a)') 'PDAF', 'OMI-internal routines:', &
                time_omi, 's'
           WRITE (*, '(a, 12x, a, 11x, F11.3, 1x, a)') 'PDAF', 'init_n_domains_pdaf:', pdaf_time_tot(42), 's'
           WRITE (*, '(a, 12x, a, 15x, F11.3, 1x, a)') 'PDAF', 'init_dim_l_pdaf:', pdaf_time_tot(45), 's'
           IF (.NOT.pdaflocal_was_used) THEN
              WRITE (*, '(a, 12x, a, 16x, F11.3, 1x, a)') 'PDAF', 'g2l_state_pdaf:', pdaf_time_tot(15), 's'
              WRITE (*, '(a, 12x, a, 16x, F11.3, 1x, a)') 'PDAF', 'l2g_state_pdaf:', pdaf_time_tot(16), 's'
           END IF
           WRITE (*, '(a, 12x, a)') 'PDAF', 'Time in OMI observation module routines '
           WRITE (*, '(a, 14x, a, 8x, F11.3, 1x, a)') 'PDAF', 'init_dim_obs_pdafomi:', pdaf_time_tot(43), 's'
           WRITE (*, '(a, 14x, a, 14x, F11.3, 1x, a)') 'PDAF', 'obs_op_pdafomi:', pdaf_time_tot(44), 's'
           WRITE (*, '(a, 14x, a, 6x, F11.3, 1x, a)') 'PDAF', 'init_dim_obs_l_pdafomi:', pdaf_time_tot(9), 's'

!            WRITE (*, '(a, 12x, a, 11x, F11.3, 1x, a)') 'PDAF', 'Time in OMI-internal routines'
!            IF (type_forget==1) THEN
!               WRITE (*, '(a, 14x, a, 10x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_init_obs_f:', pdaf_time_tot(50), 's'
!               WRITE (*, '(a, 14x, a, 9x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_init_obsvar:', pdaf_time_tot(49), 's'
!            END IF
!            WRITE (*, '(a, 14x, a, 13x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_g2l_obs:', pdaf_time_tot(46), 's'
!            IF (type_forget==1) THEN
!               WRITE (*, '(a, 14x, a, 7x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_init_obsvar_l:', pdaf_time_tot(52), 's'
!            END IF
!            WRITE (*, '(a, 14x, a, 10x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_init_obs_l:', pdaf_time_tot(47), 's'
!            WRITE (*, '(a, 14x, a, 9x, F11.3, 1x, a)') 'PDAF', 'PDAFomi_prodRinvA_l:', pdaf_time_tot(48), 's'
        ELSE
           ! Output when NOT using OMI

           ! Filter-specific part
           WRITE (*, '(a, 10x, a, 16x, F11.3, 1x, a)') 'PDAF', 'LETKF analysis:', pdaf_time_tot(3), 's'
           WRITE (*, '(a, 12x, a, 6x, F11.3, 1x, a)') 'PDAF', 'PDAF-internal operations:', pdaf_time_tot(51), 's'
           WRITE (*, '(a, 12x, a, 11x, F11.3, 1x, a)') 'PDAF', 'init_n_domains_pdaf:', pdaf_time_tot(42), 's'
           WRITE (*, '(a, 12x, a, 11x, F11.3, 1x, a)') 'PDAF', 'init_dim_obs_f_pdaf:', pdaf_time_tot(43), 's'
           WRITE (*, '(a, 12x, a, 17x, F11.3, 1x, a)') 'PDAF', 'obs_op_f_pdaf:', pdaf_time_tot(44), 's'
           IF (type_forget==1) THEN
              WRITE (*, '(a, 12x, a, 15x, F11.3, 1x, a)') 'PDAF', 'init_obs_f_pdaf:', pdaf_time_tot(50), 's'
              WRITE (*, '(a, 12x, a, 14x, F11.3, 1x, a)') 'PDAF', 'init_obsvar_pdaf:', pdaf_time_tot(49), 's'
           END IF
           WRITE (*, '(a, 12x, a, 15x, F11.3, 1x, a)') 'PDAF', 'init_dim_l_pdaf:', pdaf_time_tot(45), 's'
           WRITE (*, '(a, 12x, a, 11x, F11.3, 1x, a)') 'PDAF', 'init_dim_obs_l_pdaf:', pdaf_time_tot(9), 's'
           WRITE (*, '(a, 12x, a, 16x, F11.3, 1x, a)') 'PDAF', 'g2l_state_pdaf:', pdaf_time_tot(15), 's'
           WRITE (*, '(a, 12x, a, 18x, F11.3, 1x, a)') 'PDAF', 'g2l_obs_pdaf:', pdaf_time_tot(46), 's'
           IF (type_forget==1) THEN
              WRITE (*, '(a, 12x, a, 12x, F11.3, 1x, a)') 'PDAF', 'init_obsvar_l_pdaf:', pdaf_time_tot(52), 's'
           END IF
           WRITE (*, '(a, 12x, a, 15x, F11.3, 1x, a)') 'PDAF', 'init_obs_l_pdaf:', pdaf_time_tot(47), 's'
           WRITE (*, '(a, 12x, a, 14x, F11.3, 1x, a)') 'PDAF', 'prodRinvA_l_pdaf:', pdaf_time_tot(48), 's'
           WRITE (*, '(a, 12x, a, 16x, F11.3, 1x, a)') 'PDAF', 'l2g_state_pdaf:', pdaf_time_tot(16), 's'
        END IF

        ! Generic part B
        WRITE (*, '(a, 10x, a, 14x, F11.3, 1x, a)') 'PDAF', 'prepoststep_pdaf:', pdaf_time_tot(5), 's'
     END IF

  ELSE IF (printtype == 4) THEN ptype

! *********************************************
! *** Print second-level timing information ***
! *********************************************

     ! Generic part
     WRITE (*, '(//a, 23x, a)') 'PDAF', 'PDAF Timing information'
     WRITE (*, '(a, 8x, 52a)') 'PDAF', ('-', i=1, 52)
     WRITE (*, '(a, 10x, a, 11x, F11.3, 1x, a)') 'PDAF', 'Initialize PDAF (1):', pdaf_time_tot(1), 's'
     IF (.not.offline_mode) THEN
        IF (subtype_filter<2) THEN
           WRITE (*, '(a, 10x, a, 9x, F11.3, 1x, a)') 'PDAF', 'Ensemble forecast (2):', pdaf_time_tot(2), 's'
        ELSE
           WRITE (*, '(a, 10x, a, 12x, F11.3, 1x, a)') 'PDAF', 'State forecast (2):', pdaf_time_tot(2), 's'
        END IF
        WRITE (*, '(a, 12x, a, F11.3, 1x, a)') 'PDAF', 'MPI communication in PDAF (19):', pdaf_time_tot(19), 's'
        IF (.not.filterpe) WRITE (*, '(a, 7x, a)') 'PDAF', &
             'Note: for filterpe=F, the time (2) includes the wait time for the analysis step'
     END IF

     IF (filterpe) THEN
        ! Filter-specific part
        WRITE (*, '(a, 10x, a, 12x, F11.3, 1x, a)') 'PDAF', 'LETKF analysis (3):', pdaf_time_tot(3), 's'
        WRITE (*, '(a, 12x, a, 7x, F11.3, 1x, a)') 'PDAF', 'global preparations (4):', pdaf_time_tot(4), 's'
        WRITE (*, '(a, 14x, a, 7x, F11.3, 1x, a)') 'PDAF', 'compute mean state (11):', pdaf_time_tot(11), 's'
        WRITE (*, '(a, 12x, a, 7x, F11.3, 1x, a)') 'PDAF', 'local analysis loop (6):', pdaf_time_tot(6), 's'
        WRITE (*, '(a, 14x, a, 2x, F11.3, 1x, a)') 'PDAF', 'search local obs. domain (9):', pdaf_time_tot(9), 's'
        WRITE (*, '(a, 14x, a, 10x, F11.3, 1x, a)') 'PDAF', 'global to local (15):', pdaf_time_tot(15), 's'
        WRITE (*, '(a, 14x, a, 12x, F11.3, 1x, a)') 'PDAF', 'local analysis (7):', pdaf_time_tot(7), 's'
        WRITE (*, '(a, 14x, a, 10x, F11.3, 1x, a)') 'PDAF', 'local to global (16):', pdaf_time_tot(16), 's'
        IF (dim_lag >0) &
             WRITE (*, '(a, 14x, a, 8x, F11.3, 1x, a)') 'PDAF', 'perform smoothing (17):', pdaf_time_tot(17), 's'


        ! Generic part B
        WRITE (*, '(a, 12x, a, 13x, F11.3, 1x, a)') 'PDAF', 'Prepoststep (5):', pdaf_time_tot(5), 's'
     END IF

  ELSE IF (printtype == 5) THEN ptype

! *****************************************
! *** Print detailed timing information ***
! *****************************************

     ! Generic part
     WRITE (*, '(//a, 23x, a)') 'PDAF', 'PDAF Timing information'
     WRITE (*, '(a, 8x, 52a)') 'PDAF', ('-', i=1, 52)
     WRITE (*, '(a, 10x, a, 11x, F11.3, 1x, a)') &
          'PDAF', 'Initialize PDAF (1):', pdaf_time_tot(1), 's'
     IF (.not.offline_mode) THEN
        IF (subtype_filter<2) THEN
           WRITE (*, '(a, 10x, a, 9x, F11.3, 1x, a)') 'PDAF', 'Ensemble forecast (2):', pdaf_time_tot(2), 's'
        ELSE
           WRITE (*, '(a, 10x, a, 12x, F11.3, 1x, a)') 'PDAF', 'State forecast (2):', pdaf_time_tot(2), 's'
        END IF
        WRITE (*, '(a, 12x, a, F11.3, 1x, a)') 'PDAF', 'MPI communication in PDAF (19):', pdaf_time_tot(19), 's'
        IF (.not.filterpe) WRITE (*, '(a, 7x, a)') 'PDAF', &
             'Note: for filterpe=F, the time (2) includes the wait time for the analysis step'
     END IF

     IF (filterpe) THEN
        ! Filter-specific part
        WRITE (*, '(a, 10x, a, 12x, F11.3, 1x, a)') 'PDAF', 'LETKF analysis (3):', pdaf_time_tot(3), 's'
        WRITE (*, '(a, 12x, a, 7x, F11.3, 1x, a)') 'PDAF', 'global preparations (4):', pdaf_time_tot(4), 's'
        WRITE (*, '(a, 14x, a, 7x, F11.3, 1x, a)') 'PDAF', 'compute mean state (11):', pdaf_time_tot(11), 's'
        WRITE (*, '(a, 12x, a, 7x, F11.3, 1x, a)') 'PDAF', 'local analysis loop (6):', pdaf_time_tot(6), 's'
        WRITE (*, '(a, 14x, a, 2x, F11.3, 1x, a)') 'PDAF', 'search local obs. domain (9):', pdaf_time_tot(9), 's'
        WRITE (*, '(a, 14x, a, 10x, F11.3, 1x, a)') 'PDAF', 'global to local (15):', pdaf_time_tot(15), 's'
        WRITE (*, '(a, 14x, a, 12x, F11.3, 1x, a)') 'PDAF', 'local analysis (7):', pdaf_time_tot(7), 's'
        WRITE (*, '(a, 16x, a, 4x, F11.3, 1x, a)') 'PDAF', 'compute new inverse A (10):', pdaf_time_tot(10), 's'
        WRITE (*, '(a, 18x, a, 21x, F11.3, 1x, a)') 'PDAF', 'HZ_l (30):', pdaf_time_tot(30), 's'
        WRITE (*, '(a, 18x, a, 12x, F11.3, 1x, a)') 'PDAF', 'complete Uinv (31):', pdaf_time_tot(31), 's'
        WRITE (*, '(a, 16x, a, 12x, F11.3, 1x, a)') 'PDAF', 'init residual (12):', pdaf_time_tot(12), 's'
        WRITE (*, '(a, 16x, a, 2x, F11.3, 1x, a)') 'PDAF', 'get state weight vector (13):', pdaf_time_tot(13), 's'
        IF (subtype_filter /= 3) THEN
           WRITE (*, '(a, 16x, a, 5x, F11.3, 1x, a)') 'PDAF', 'prepare ens. weights (20):', pdaf_time_tot(20), 's'
           WRITE (*, '(a, 16x, a, 4x, F11.3, 1x, a)') 'PDAF', 'store ensemble matrix (21):', pdaf_time_tot(21), 's'
           WRITE (*, '(a, 16x, a, 10x, F11.3, 1x, a)') 'PDAF', 'update ensemble (22):', pdaf_time_tot(22), 's'
           WRITE (*, '(a, 14x, a, 10x, F11.3, 1x, a)') 'PDAF', 'local to global (16):', pdaf_time_tot(16), 's'
           IF (dim_lag >0) &
                WRITE (*, '(a, 14x, a, 8x, F11.3, 1x, a)') 'PDAF', 'perform smoothing (17):', pdaf_time_tot(17), 's'
        ELSE
           WRITE (*, '(a, 16x, a, F11.3, 1x, a)') 'PDAF', 'update state and ensemble (14):', pdaf_time_tot(14), 's'
        END IF

        ! Generic part B
        WRITE (*, '(a, 12x, a, 13x, F11.3, 1x, a)') 'PDAF', 'Prepoststep (5):', pdaf_time_tot(5), 's'
     END IF


  ELSE IF (printtype == 10) THEN ptype

! *******************************
! *** Print allocated memory  ***
! *******************************

     WRITE (*, '(/a, 23x, a)') 'PDAF', 'PDAF Memory overview'
     WRITE (*, '(a, 10x, 45a)') 'PDAF', ('-', i=1, 45)
     WRITE (*, '(a, 21x, a)') 'PDAF', 'Allocated memory  (MiB)'
     WRITE (*, '(a, 14x, a, 1x, f10.3, a)') &
          'PDAF', 'state and A:', pdaf_memcount_get(1, 'M'), ' MiB (persistent)'
     WRITE (*, '(a, 11x, a, 1x, f10.3, a)') &
          'PDAF', 'ensemble array:', pdaf_memcount_get(2, 'M'), ' MiB (persistent)'
     WRITE (*, '(a, 12x, a, 1x, f10.3, a)') &
          'PDAF', 'analysis step:', pdaf_memcount_get(3, 'M'), ' MiB (temporary)'


  ELSE IF (printtype == 11) THEN ptype

! ****************************************
! *** Print globally allocated memory  ***
! ****************************************

     memcount_global(1) = pdaf_memcount_get_global(1, 'M', COMM_pdaf)
     memcount_global(2) = pdaf_memcount_get_global(2, 'M', COMM_pdaf)
     memcount_global(3) = pdaf_memcount_get_global(3, 'M', COMM_pdaf)

     IF (mype_world==0) THEN
        WRITE (*, '(/a, 23x, a)') 'PDAF', 'PDAF Memory overview'
        WRITE (*, '(a, 10x, 45a)') 'PDAF', ('-', i=1, 45)
        WRITE (*, '(a, 17x, a)') 'PDAF', 'Globally allocated memory  (MiB)'
        WRITE (*, '(a, 14x, a, 1x, f12.3, a)') &
             'PDAF', 'state and A:', memcount_global(1), ' MiB (persistent)'
        WRITE (*, '(a, 11x, a, 1x, f12.3, a)') &
             'PDAF', 'ensemble array:', memcount_global(2), ' MiB (persistent)'
        WRITE (*, '(a, 12x, a, 1x, f12.3, a)') &
             'PDAF', 'analysis step:', memcount_global(3), ' MiB (temporary)'
     END IF

  END IF ptype


END SUBROUTINE PDAF_letkf_memtime
