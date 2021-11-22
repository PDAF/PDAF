!$Id: PDAFomi_obs_op.F90 552 2020-11-18 14:05:35Z lnerger $

!> PDAF-OMI TEMPLATE for observation operators
!!
!! This module contains a template routine for implementing
!! an observation operator using variables prepared in
!! init_dim_obs.
!!
!! Observation operators provided by PDAF can be found in
!! src/PDAFomi_obs_op.F90
!!
MODULE obs_op_pdafomi_TEMPLATE

  USE PDAFomi_obs_f, ONLY: obs_f, PDAFomi_gather_obsstate, debug

CONTAINS

!-------------------------------------------------------------------------------
!> Observation operator TEMPLATE
!!
!! The routine has to fill the part of the full observation 
!! vector OBS_F_ALL that represents the current observation
!! type.
!! 
!! The routine first applies the observation operator
!! for the current observation type and then calls
!! PDAFomi_gather_obsstate to gather the observation over
!! all processes and fills OBS_F_ALL.
!!
!! Thus, the observation operator has as least two parts
!! 1. process-local computation of observations of
!!    the current type (array OSTATE_P)
!! 2. gathering of full observation vector and filling
!!    of OBS_F_ALL by calling PDAFomi_gather_obsstate
!!
!! The operations in step 1 can include interpolation, 
!! averaging, or computing the observation values as
!! a function of different elements of the state vector STATE_P
!!
!! Please help to advance PDAF-OMI by providing new 
!! observation operator for inclusion into PDAF-OMI
!! under the LGPL license.
!!
!! __Revision history:__
!! * 2020-11 - Lars Nerger - Initial code for template
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_obs_op_TEMPLATE(thisobs, nrows, state_p, obs_f_all)

    IMPLICIT NONE

! *** Arguments ***
    ! *** The interface can be adapted to the user's needs. 
    ! *** However, OBS_F_ALL is always required.
    ! *** THISOBS and STATE_P are usually required, but one could
    ! *** implement a multi-step observation operator that separates
    ! *** the functional computation from interpolations
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    INTEGER, INTENT(in) :: nrows           !< Number of values to be averaged
    REAL, INTENT(in)    :: state_p(:)      !< PE-local model state (dim_p)
    REAL, INTENT(inout) :: obs_f_all(:)    !< Full observed state for all observation types (nobs_f_all)

! *** Local variables ***
    INTEGER :: i, row                      ! Counter
    REAL, ALLOCATABLE :: ostate_p(:)       ! local observed part of state vector
    REAL :: rrows                          ! Real-value for nrows


! *********************************************
! *** Perform application of observation    ***
! *** operator on vector state_p            ***
! *********************************************

    doassim: IF (thisobs%doassim == 1) THEN

       ! Print debug information
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_TEMPLATE -- START'
          WRITE (*,*) '++ OMI-debug: ', debug, '  PDAFomi_obs_op_TEMPLATE -- DESCRIBE OBS OPERATOR'
          WRITE (*,*) '++ OMI-debug obs_op_TEMPLATE:', debug, 'thisobs%dim_obs_p', thisobs%dim_obs_p
          WRITE (*,*) '++ OMI-debug obs_op_TEMPLATE:', debug, 'number of points to average', nrows
          WRITE (*,*) '++ OMI-debug obs_op_TEMPLATE:', debug, 'thisobs%id_obs_p', thisobs%id_obs_p
       END IF

       ! Consistency check
       IF (.NOT.ALLOCATED(thisobs%id_obs_p)) THEN
          WRITE (*,*) 'ERROR: PDAFomi_obs_op_TEMPLATE - thisobs%id_obs_p is not allocated'
       END IF


       ! ************************************************************
       ! *** Process-local: Initialize observed part state vector ***
       ! ************************************************************

       IF (thisobs%dim_obs_p>0) THEN
          ALLOCATE(ostate_p(thisobs%dim_obs_p))
       ELSE
          ALLOCATE(ostate_p(1))
       END IF


       ! *** Initialize ostate_p e.g. by                      ***
       ! *** - selecting elements of the state vector state_p ***
       ! ***      according to thisobs%id_obs_p               ***
       ! *** - deriving an observed variable as a function    ***
       ! ***      of elements of the state vector             ***

       ! Example for computing ostate_p by averaging nrows grid point values
!        rrows = REAL(nrows)
!        DO i = 1, thisobs%dim_obs_p
!           ostate_p(i) = 0.0
!           DO row = 1, nrows
!              ostate_p(i) = ostate_p(i) + state_p(thisobs%id_obs_p(row,i))
!           END DO
!           ostate_p(i) = ostate_p(i) / rrows
!        ENDDO


       ! *************************************************
       ! *** Global: Gather full observed state vector ***
       ! ***            THIS IS MANDATORY!             ***
       ! *************************************************

       CALL PDAFomi_gather_obsstate(thisobs, ostate_p, obs_f_all)

       ! *** Clean up
       DEALLOCATE(ostate_p)

       ! Print debug information
       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_TEMPLATE -- END'

    END IF doassim

  END SUBROUTINE PDAFomi_obs_op_TEMPLATE

END MODULE obs_op_pdafomi_TEMPLATE
