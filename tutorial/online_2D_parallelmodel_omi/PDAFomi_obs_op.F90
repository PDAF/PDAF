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
! get_interp_coeff_lin1D
!        Routine to comput linear interpo;lation in 1D
! get_interp_coeff_lin
!        Routine to compute interpolation coefficients for linear
!        interpolations (linear, bi-linear, tri-linear)


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
  SUBROUTINE get_interp_coeff_tri(gpc, oc, icoeff)

! !DESCRIPTION:
! The routine computes the coefficients for triangular interpolation
! as barycentric coordinates.
! The computation is done for one observation given the 
! observation coordinates (OC) as well as the coordinates of the 
! grid points (GPC). In GPC each row contains the coordinates
! for one grid point.
!
! !REVISION HISTORY:
! 2019-12 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    REAL, INTENT(in)    :: gpc(3,2)    ! Coordinates of grid points
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
    denum = (gpc(2,2) - gpc(3,2)) * (gpc(1,1) - gpc(3,1)) + (gpc(3,1) - gpc(2,1)) * (gpc(1,2) - gpc(3,2))

    ! compute coefficients
    icoeff(1) = (gpc(2,2) - gpc(3,2)) * (oc(1) - gpc(3,1)) + (gpc(3,1) - gpc(2,1)) * (oc(2) - gpc(3,2))
    icoeff(1) = icoeff(1) / denum

    icoeff(2) = (gpc(3,2) - gpc(1,2)) * (oc(1) - gpc(3,1)) + (gpc(1,1) - gpc(3,1)) * (oc(2) - gpc(3,2))
    icoeff(2) = icoeff(2) / denum

    icoeff(3) = 1.0 - icoeff(1) - icoeff(2)

  END SUBROUTINE get_interp_coeff_tri


!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: get_interp_coeff_lin1D --- Initialize linear interpolation coefficients in 1D
!
! !INTERFACE:
  SUBROUTINE get_interp_coeff_lin1D(gpc, oc, icoeff)

! !DESCRIPTION:
! The routine computes the coefficients for linear interpolation
! in 1 dimensions.
! The computation is done for one observation given the 
! observation coordinates (OC) as well as the coordinates of the 
! grid points (GPC). 
!
! !REVISION HISTORY:
! 2019-12 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    REAL, INTENT(in)    :: gpc(2)      ! Coordinates of grid points
    REAL, INTENT(in)    :: oc          ! Coordinates of observation
    REAL, INTENT(inout) :: icoeff(2)   ! Interpolation coefficients
!EOP


! ****************************************************************
! *** Compute linear interpolation coefficients in 1 dimension ***
! ****************************************************************

    icoeff(1) = (gpc(2) - oc) / (gpc(2) - gpc(1))
    icoeff(2) = (oc - gpc(1)) / (gpc(2) - gpc(1))

  END SUBROUTINE get_interp_coeff_lin1D


!-------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: get_interp_coeff_lin --- Initialize linear interpolation coefficients
!
! !INTERFACE:
  SUBROUTINE get_interp_coeff_lin(num_gp, n_dim, gpc, oc, icoeff)

! !DESCRIPTION:
! The routine computes the coefficients for linear interpolation
! in 1, 2, or 3 dimensions.
! The computation is done for one observation given the 
! observation coordinates (OC) as well as the coordinates of the 
! grid points (GPC). In GPC each row contains the coordinates
! for one grid point.  
!
! Setup of GPC:
! The first index is specifies the grid point, while the second the
! coordinate
! - For n_dim=X only the first X coordinate values are used
! - The ordering of the coordinates and coefficient is the following:
!                       (7)------(8) 
!                       /|       /|    with
!                     (5)+-----(6)|       - column 1
!                      | |      | |       / column 2
!                      |(3)-----+(4)      | column 3
!                      |/       |/
!                     (1) ---- (2)
!   thus gpc(1,1)= gpc(2,1), gpc(1,2)=gpc(3,2), gpc(1,3)=gpc(5,3)
! - For bi-linear interpolation only the coordinates for grid
!   points 1, 2, and 3 are used to compute the coefficients
! - For tri-linear interpolation only the coordinates for grid
!   points 1, 2, 3, and 5 are used to compute the coefficients
!   (for bi-linear interpolation gpc only needs to have length 3
!   for tri-linear the length 5)
! - num_gp=2 for n_dim=1; num_gp=4 for n_dim=2; num_gp=8 for n_dim=3
!
! !REVISION HISTORY:
! 2019-12 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: num_gp         ! Length of ICOEFF
    INTEGER, INTENT(in) :: n_dim          ! Number of dimensions in interpolation
    REAL, INTENT(in)    :: gpc(:,:)       ! Coordinates of grid points
    REAL, INTENT(in)    :: oc(:)          ! Coordinates of observation
    REAL, INTENT(inout) :: icoeff(num_gp) ! Interpolation coefficients
!EOP

! *** Local variables ***
    REAL :: denum    ! denumerator


    IF (n_dim == 1) THEN

! ****************************************************************
! *** Compute linear interpolation coefficients in 1 dimension ***
! ****************************************************************
       
       ! Checks
       IF (num_gp /= 2) WRITE (*,'(a,3x,a)') &
            'PDAFomi', 'ERROR: get_interp_coeff_lin - NUM_GP=2 required!'
       IF (gpc(2,1) == gpc(1,1))  WRITE (*,'(a,3x,a)') &
            'PDAFomi', 'ERROR: get_interp_coeff_lin - wrong setting of coordinates!'

       ! Compute coefficients
       icoeff(1) = (gpc(2,1) - oc(1)) / (gpc(2,1) - gpc(1,1))
       icoeff(2) = (oc(1) - gpc(1,1)) / (gpc(2,1) - gpc(1,1))

    ELSE IF (n_dim == 2) THEN

! ********************************************************
! *** Compute coefficients for bi-linear interpolation ***
! *** Order of coefficients:  (3) ---- (4)             ***
! ***                          |        |              ***
! ***                         (1) ---- (2)             ***
! ********************************************************

       ! Checks
       IF (num_gp /= 4) WRITE (*,'(a,5x,a)') &
            'PDAF', 'ERROR: get_interp_coeff_lin - NUM_GP=4 required!'
       IF (gpc(2,1) == gpc(1,1) .OR. gpc(3,2) == gpc(1,2))  WRITE (*,'(a,3x,a)') &
            'PDAFomi', 'ERROR: get_interp_coeff_lin - wrong setting of coordinates!'

       ! Compute coefficients
       denum = (gpc(2,1) - gpc(1,1)) * (gpc(3,2) - gpc(1,2))

       icoeff(1) = (gpc(2,1) - oc(1)) * (gpc(3,2) - oc(2)) / denum
       icoeff(2) = (oc(1) - gpc(1,1)) * (gpc(3,2) - oc(2)) / denum
       icoeff(3) = (gpc(2,1) - oc(1)) * (oc(2) - gpc(1,2)) / denum
       icoeff(4) = (oc(1) - gpc(1,1)) * (oc(2) - gpc(1,2)) / denum

    ELSE IF (n_dim == 3) THEN

! *********************************************************
! *** Compute coefficients for tri-linear interpolation ***
! *** Order of coefficients:    (7)------(8)            ***
! ***                           /|       /|             ***
! ***                         (5)+-----(6)|             ***
! ***                          | |      | |             ***
! ***                          |(3)-----+(4)            ***
! ***                          |/       |/              ***
! ***                         (1) ---- (2)              ***
! *********************************************************

       ! Checks
       IF (num_gp /= 8) WRITE (*,'(a,5x,a)') &
            'PDAF', 'ERROR: get_interp_coeff_lin - NUM_GP=8 required!'
       IF (gpc(2,1) == gpc(1,1) .OR. gpc(3,2) == gpc(1,2) .OR. gpc(5,3) == gpc(1,3)) &
            WRITE (*,'(a,3x,a)') &
            'PDAFomi', 'ERROR: get_interp_coeff_lin - wrong setting of coordinates!'

       ! Compute coefficients
       denum = (gpc(2,1) - gpc(1,1)) * (gpc(3,2) - gpc(1,2)) * (gpc(5,3) - gpc(1,3))

       icoeff(1) = (gpc(2,1) - oc(1)) * (gpc(3,2) - oc(2)) * (gpc(5,3) - oc(3)) / denum
       icoeff(2) = (oc(1) - gpc(1,1)) * (gpc(3,2) - oc(2)) * (gpc(5,3) - oc(3)) / denum
       icoeff(3) = (gpc(2,1) - oc(1)) * (oc(2) - gpc(1,2)) * (gpc(5,3) - oc(3)) / denum
       icoeff(4) = (oc(1) - gpc(1,1)) * (oc(2) - gpc(1,2)) * (gpc(5,3) - oc(3)) / denum

       icoeff(5) = (gpc(2,1) - oc(1)) * (gpc(3,2) - oc(2)) * (oc(3) - gpc(1,3)) / denum
       icoeff(6) = (oc(1) - gpc(1,1)) * (gpc(3,2) - oc(2)) * (oc(3) - gpc(1,3)) / denum
       icoeff(7) = (gpc(2,1) - oc(1)) * (oc(2) - gpc(1,2)) * (oc(3) - gpc(1,3)) / denum
       icoeff(8) = (oc(1) - gpc(1,1)) * (oc(2) - gpc(1,2)) * (oc(3) - gpc(1,3)) / denum

    END IF

  END SUBROUTINE get_interp_coeff_lin

END MODULE PDAFomi_obs_op
