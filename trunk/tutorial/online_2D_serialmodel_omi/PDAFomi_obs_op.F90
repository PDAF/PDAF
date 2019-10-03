! Copyright (c) 2004-2019 Lars Nerger
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
! !MODULE:
MODULE PDAFomi_obs_op
!
! !DESCRIPTION:
! This module contains generic routines for several observation
! operators to be used after preparation with INIT_DIM_OBS_F
! The operators are:
!
! obs_op_f_gridpoint
!        Observation operator for data at grid points. The routine
!        selects values of the state vector according to an index array
! obs_op_f_gridavg
!        Observation operator for the case the observations are the
!        average of grid point values. The routine computes these
!        averages according to an index array. 

  INTERFACE obs_op_f
     MODULE PROCEDURE obs_op_f_gridpoint
     MODULE PROCEDURE obs_op_f_gridavg
  END INTERFACE

CONTAINS

!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: obs_op_f_gridpoint --- observation operator for data at grid points
!
! !INTERFACE:
  SUBROUTINE obs_op_f_gridpoint(dim_p, nobs_f_all, nobs_p_one, nobs_f_one, &
       id_obs_p_one, state_p, obs_f_all, offset_obs)

! !DESCRIPTION:
! Application of observation operator for the case that 
! model variables are observerved at model grid points. 
! For this case INIT_DIM_OBS_F will prepare an index 
! array ID_OBS_P_OBS containing the information which 
! elements of the  PE-local state vector contain the
! observed values.
!
! The routine is called by all filter processes. It first
! selects the observed elements for a PE-local domain. 
! Afterwards, the values are gathered into the full vector.
!
! The routine has to fill the part of the full observation 
! vector OBS_F_ALL that represents the current observation
! type. Its offset in the full observation vector is specified
! by OFFSET_OBS. Upon exit from the routine OFFSET_OBS has to
! be incremented by the number of observations filled in.
!
! !REVISION HISTORY:
! 2019-06 - Lars Nerger - Initial code from restructuring observation routines
! Later revisions - see svn log
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: dim_p                    ! PE-local dimension of state
    INTEGER, INTENT(in) :: nobs_f_all               ! Length of obs. vector for all observations
    INTEGER, INTENT(in) :: nobs_p_one               ! PE-local number observations of current observation type
    INTEGER, INTENT(in) :: nobs_f_one               ! Full number observations of current observation type
    INTEGER, INTENT(in) :: id_obs_p_one(1, nobs_p_one) ! Index of current observations in PE-local state vector
    REAL, INTENT(in)    :: state_p(dim_p)           ! PE-local model state
    REAL, INTENT(inout) :: obs_f_all(nobs_f_all)    ! Full observed state for all observation types
    INTEGER, INTENT(inout) :: offset_obs            ! Offset of current observation in overall observation vector
!EOP

! *** Local variables ***
    INTEGER :: i                       ! Counter
    REAL, ALLOCATABLE :: m_state_p(:)  ! local observed part of state vector
    INTEGER :: status                  ! status flag


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

    if (nobs_p_one>0) then
       ALLOCATE(m_state_p(nobs_p_one))
    else
       ALLOCATE(m_state_p(1))
    end if

    ! *** PE-local: Initialize observed part state vector
    DO i = 1, nobs_p_one
       m_state_p(i) = state_p(id_obs_p_one(1, i)) 
    ENDDO

    ! *** Gather observation vector - SST part from cnt_obs+1 in m_state_f ***
    CALL PDAF_gather_obs_f_flex(nobs_p_one, nobs_f_one, m_state_p, &
         obs_f_all(offset_obs+1), status)

    ! Increment offset in observaton vector
    offset_obs = offset_obs + nobs_f_one

    DEALLOCATE(m_state_p)

  END SUBROUTINE obs_op_f_gridpoint


!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: obs_op_f_gridavg --- observation operator for averaging grid point values
!
! !INTERFACE:
  SUBROUTINE obs_op_f_gridavg(dim_p, nobs_f_all, nobs_p_one, nobs_f_one, nrows, &
       id_obs_p_one, state_p, obs_f_all, offset_obs)

! !DESCRIPTION:
! Application of observation operator for the case that 
! the observation value is given as the average of model
! grid point values.
! For this case INIT_DIM_OBS_F will prepare an index array 
! that contains several rows holding the indices of the state
! vector elements which are to be averaged to represent an 
! observation.
!
! The routine is called by all filter processes, 
! and the operation has to be performed by each 
! these processes for its PE-local domain before the 
! information from all PEs is gathered.
!
! The routine has to fill the part of the full observation 
! vector OBS_F_ALL that represents the current observation
! type. Its offset in the full observation vector is specified
! by OFFSET_OBS. Upon exit from the routine OFFSET_OBS has to
! be incremented by the number of observations filled in.
!
! !REVISION HISTORY:
! 2019-06 - Lars Nerger - Initial code from restructuring observation routines
! Later revisions - see svn log
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: dim_p                    ! PE-local dimension of state
    INTEGER, INTENT(in) :: nobs_f_all               ! Length of obs. vector for all observations
    INTEGER, INTENT(in) :: nobs_p_one               ! PE-local number observations of current observation type
    INTEGER, INTENT(in) :: nobs_f_one               ! Full number observations of current observation type
    INTEGER, INTENT(in) :: nrows                    ! Number of values to be averaged
    INTEGER, INTENT(in) :: id_obs_p_one(nrows, nobs_p_one) ! Index of current observations in PE-local state vector
    REAL, INTENT(in)    :: state_p(dim_p)           ! PE-local model state
    REAL, INTENT(inout) :: obs_f_all(nobs_f_all)    ! Full observed state for all observation types
    INTEGER, INTENT(inout) :: offset_obs            ! Offset of current observation in overall observation vector
!EOP

! *** Local variables ***
    INTEGER :: i, row                  ! Counter
    REAL, ALLOCATABLE :: m_state_p(:)  ! local observed part of state vector
    REAL :: rrows                      ! Real-value for nrows
    INTEGER :: status                  ! status flag


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

    if (nobs_p_one>0) then
       ALLOCATE(m_state_p(nobs_p_one))
    else
       ALLOCATE(m_state_p(1))
    end if

    rrows = REAL(nrows)

    ! *** PE-local: Initialize observed part state vector by averaging
    DO i = 1, nobs_p_one
       m_state_p(i) = 0.0
       DO row = 1, nrows
          m_state_p(i) = m_state_p(i) + state_p(id_obs_p_one(row,i))
       END DO
       m_state_p(i) = m_state_p(i) / rrows
    ENDDO

    ! *** Gather observation vector - SST part from cnt_obs+1 in m_state_f ***
    CALL PDAF_gather_obs_f_flex(nobs_p_one, nobs_f_one, m_state_p, &
         obs_f_all(offset_obs+1), status)

    ! Increment offset in observaton vector
    offset_obs = offset_obs + nobs_f_one

    DEALLOCATE(m_state_p)

  END SUBROUTINE obs_op_f_gridavg

END MODULE PDAFomi_obs_op
