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
!        Observation operator for the case that the observations are the
!        average of grid point values. The routine computes these
!        averages according to an index array. 
! obs_op_f_interp_lin
!        Observation operator for the case that the observations are
!        linear interpolation from the grid points. The interpolation
!        coefficients are pre-computed.
!
! Helper routines for the operators:
! get_interp_coeff_tri
!        Routine to compute interpolation coefficients for triangular
!        interpolation from barycentric coordinates.


  INTERFACE obs_op_f
     MODULE PROCEDURE obs_op_f_gridpoint
     MODULE PROCEDURE obs_op_f_gridavg
     MODULE PROCEDURE obs_op_f_interp_lin
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

    ! *** Gather observation vector - part from cnt_obs+1 in obs_f_all ***
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

    ! *** Gather observation vector - part from cnt_obs+1 in obs_f_all ***
    CALL PDAF_gather_obs_f_flex(nobs_p_one, nobs_f_one, m_state_p, &
         obs_f_all(offset_obs+1), status)

    ! Increment offset in observaton vector
    offset_obs = offset_obs + nobs_f_one

    DEALLOCATE(m_state_p)

  END SUBROUTINE obs_op_f_gridavg


!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: obs_op_f_interp_lin --- observation operator for linear interpolation
!
! !INTERFACE:
  SUBROUTINE obs_op_f_interp_lin(dim_p, nobs_f_all, nobs_p_one, nobs_f_one, nrows, &
       id_obs_p_one, coeff_p_one, state_p, obs_f_all, offset_obs)

! !DESCRIPTION:
! Application of observation operator for the case that the
! observation value is given as the interpolation using
! pre-computed coefficients. 
! For this case, INIT_DIM_OBS_F will prepare the index array
! ID_OBS_P_ONE and the array COEFF_P_ONE of interpolation
! coefficients. ID_OBS_P_ONE contains several rows holding
! the indices of the state vector elements which are to be
! interpolated to represent an observation. COEFF_P_ONE
! contains the interpolation coefficients and can be prepared
! using a help routine like get_interp_coeff_tri.
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
! 2019-12 - Lars Nerger - Initial code
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
    REAL, INTENT(in)    :: coeff_p_one(nrows, nobs_p_one) ! interpolation coefficients for PE-local obs.
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

    ! *** PE-local: Initialize observed part state vector by weighted averaging
    DO i = 1, nobs_p_one
       m_state_p(i) = 0.0
       DO row = 1, nrows
          m_state_p(i) = m_state_p(i) + coeff_p_one(row,i)*state_p(id_obs_p_one(row,i))
       END DO
       m_state_p(i) = m_state_p(i) / rrows
    ENDDO

    ! *** Gather observation vector - part from cnt_obs+1 in obs_f_all ***
    CALL PDAF_gather_obs_f_flex(nobs_p_one, nobs_f_one, m_state_p, &
         obs_f_all(offset_obs+1), status)

    ! Increment offset in observaton vector
    offset_obs = offset_obs + nobs_f_one

    DEALLOCATE(m_state_p)

  END SUBROUTINE obs_op_f_interp_lin


!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: get_interp_coeff_tri --- Initialize interpolation coefficients in triangle
!
! !INTERFACE:
  SUBROUTINE get_interp_coeff_tri(c1, c2, c3, oc, icoeff)

! !DESCRIPTION:
! The routine computes the coefficients for triangular interpolation
! as barycentric coordinates.
! The computation is done for one observation given the 
! observation coordianates (OCOORDS) as well as the coordinates
! of 3 grid points (C1, C2, C3).
!
! !REVISION HISTORY:
! 2019-12 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    REAL, INTENT(in)    :: c1(2)       ! Coordinates of grid point 1
    REAL, INTENT(in)    :: c2(2)       ! Coordinates of grid point 2
    REAL, INTENT(in)    :: c3(2)       ! Coordinates of grid point 3
    REAL, INTENT(in)    :: oc(2)       ! Coordinates of observation
    REAL, INTENT(inout) :: icoeff(3)   ! Interpolation coefficients
!EOP

! *** Local variables ***
    REAL :: denum    ! denumerator


! ******************************************
! *** Compute interpolation coefficients ***
! *** as barycentric coordinates         ***
! ******************************************

    ! common denumerator for coefficients 1 and 2
    denum = (c2(2) - c3(2)) * (c1(1) - c3(1)) + (c3(1) - c2(1)) * (c1(2) - c3(2))

    ! compute coefficients
    icoeff(1) = (c2(2) - c3(2)) * (oc(1) - c3(1)) + (c3(1) - c2(1)) * (oc(2) - c3(2))
    icoeff(1) = icoeff(1) / denum

    icoeff(2) = (c3(2) - c1(2)) * (oc(1) - c3(1)) + (c1(1) - c3(1)) * (oc(2) - c3(2))
    icoeff(2) = icoeff(2) / denum

    icoeff(3) = 1.0 - icoeff(1) - icoeff(2)

  END SUBROUTINE get_interp_coeff_tri

END MODULE PDAFomi_obs_op
