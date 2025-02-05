! Copyright (c) 2014-2024 Paul Kirchgessner
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
!>  PDAF-internal initialization of NETF
!!
!! Initialization of NETF within PDAF. Performed are:
!! * initialize filter-specific parameters
!! * print screen information on filter configuration.
!!
!!  !  This is a core routine of PDAF and   !
!!  !   should not be changed by the user   !
!!
!! __Revision history:__
!! *  2014-05 - Paul Kirchgessner - Initial code based on code for ETKF
!! *  Later revisions - see repository log
!!
SUBROUTINE PDAF_NETF_init(subtype, param_int, dim_pint, param_real, dim_preal, &
     ensemblefilter, fixedbasis, verbose, outflag)

  USE PDAF_mod_filter, &
       ONLY: localfilter, dim_lag, globalobs
  USE PDAF_netf, &
       ONLY: incremental, forget, type_forget, type_trans, &
       type_noise, noise_amp, type_winf, limit_winf, &
       PDAF_netf_set_iparam, PDAF_netf_set_rparam
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
  dim_lag = 0
  type_noise = 0
  type_winf = 0
  limit_winf = 0.0
  noise_amp = 0.0

  ! Parse provided parameters
  flagsum = 0
  DO i=3, dim_pint
     CALL PDAF_netf_set_iparam(i, param_int(i), outflag)
     flagsum = flagsum+outflag
  END DO
  DO i=1, dim_preal
     CALL PDAF_netf_set_rparam(i, param_real(i), outflag)
     flagsum = flagsum+outflag
  END DO


  ! Define whether filter is mode-based or ensemble-based
  ensemblefilter = .TRUE.

  ! Define whether filter is a domain-local filter
  localfilter = 0

  ! Define that filter needs global observations (used for OMI)
  globalobs = 1

  ! Initialize flag for fixed-basis filters
  fixedbasis = .FALSE.


! *********************
! *** Screen output ***
! *********************

  writeout: IF (verbose == 1) THEN
  
     WRITE(*, '(/a)') 'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     WRITE(*, '(a)')  'PDAF    +++      Nonlinear Ensemble Transform Filter (NETF)       +++'
     WRITE(*, '(a)')  'PDAF    +++                                                       +++'
     WRITE(*, '(a)')  'PDAF    +++                         by                            +++'
     WRITE(*, '(a)')  'PDAF    +++ J. Toedter, B. Ahrens, Mon. Wea. Rev. 143 (2015) 1347 +++'
     WRITE(*, '(a)')  'PDAF    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     IF (flagsum == 0) THEN

        ! *** General output ***
        WRITE (*, '(/a, 4x, a)') 'PDAF', 'NETF configuration'
        WRITE (*, '(a, 11x, a, i1)') 'PDAF', 'filter sub-type = ', subtype
        IF (subtype == 0) THEN
           WRITE (*, '(a, 12x, a)') 'PDAF', '--> NETF '
        ELSE
           WRITE (*, '(/5x, a/)') 'PDAF-ERROR(3): No valid subtype!'
           outflag = 3
        END IF
        IF (type_trans == 0) THEN
           WRITE (*, '(a, 12x, a)') 'PDAF', '--> Transform ensemble including product with random matrix'
        ELSE IF (type_trans == 1) THEN
           WRITE (*, '(a, 12x, a)') 'PDAF', '--> Deterministic symmetric ensemble transformation'
        END IF
        IF (incremental == 1) &
             WRITE (*, '(a, 12x, a)') 'PDAF', '--> Perform incremental updating'
        IF (type_forget == 0) THEN
           WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> prior inflation, forgetting factor:', forget
        ELSEIF (type_forget == 2) THEN
           WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> posterior inflation, forgetting factor:', forget
        ENDIF
        IF (type_noise == 0) THEN
           WRITE (*, '(a, 12x, a)') 'PDAF', '--> no noise added to particles'
        ELSEIF (type_noise == 1) THEN
           WRITE (*, '(a, 12x, a)') 'PDAF', '--> use noise of constant variance'
        ELSEIF (type_noise == 2) THEN
           WRITE (*, '(a, 12x, a)') 'PDAF', '--> use noise with amplitude relative to ensemble standard deviation'
        END IF
        WRITE (*, '(a, 12x, a, f8.3)') 'PDAF', '--> noise amplitude/factor', noise_amp
        IF (type_winf == 1) THEN
           WRITE (*, '(a, 12x, a, f8.3)') 'PDAF', '--> inflate particle weights so that N_eff/N> ', limit_winf
        END IF
        IF (observe_ens) &
             WRITE (*, '(a, 12x, a, 1x, l)') 'PDAF', '--> observe_ens:', observe_ens
     ELSE
        WRITE (*, '(/5x, a/)') 'PDAF-ERROR: Invalid parameter setting - check prior output!'
     END IF

  END IF writeout

END SUBROUTINE PDAF_NETF_init
