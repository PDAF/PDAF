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

!> PDAF-OMI observation operators
!!
!! This module contains generic routines for several observation
!! operators to be used after preparation with init_dim_obs_f
!!
!! The observation operators are:
!!
!! * PDAFomi_obs_op_gridpoint\n
!!        Observation operator for data at grid points. The routine
!!        selects values of the state vector according to an index array
!! * PDAFomi_obs_op_gridavg\n
!!        Observation operator for the case that the observations are the
!!        average of grid point values. The routine computes these
!!        averages according to an index array. 
!! * PDAFomi_obs_op_interp_lin\n
!!        Observation operator for the case that the observations are
!!        linear interpolated from the grid points. The interpolation
!!        coefficients are pre-computed.
!! * PDAFomi_obs_op_gatheronly\n
!!        Observation operator for the case of strongly coupled assimilation
!!        to gather an observation which only exists in other compartments.
!!
!! Adjoint observation operators:
!!
!! * PDAFomi_obs_op_adj_gridpoint\n
!!        Adjoint observation operator for data at grid points. The routine
!!        selects values of the state vector according to an index array
!! * PDAFomi_obs_op_adj_gridavg\n
!!        Adjoint observation operator for the case that the observations
!!        are the average of grid point values. The routine computes these
!!        averages according to an index array. 
!! * PDAFomi_obs_op_adj_interp_lin\n
!!        Adjoint observation operator for the case that the observations
!!        are linear interpolatied from the grid points. The interpolation
!!        coefficients are pre-computed.
!! * PDAFomi_obs_op_adj_gatheronly\n
!!        Adjoint observation operator for the case of strongly coupled assimilation
!!        to gather an observation which only exists in other compartments.
!!
!! Helper routines for the operators:
!! * PDAFomi_get_interp_coeff_tri \n
!!        Routine to compute interpolation coefficients for triangular
!!        interpolation from barycentric coordinates.
!! * PDAFomi_get_interp_coeff_lin1D \n
!!        Routine to comput linear interpolation in 1D
!! * PDAFomi_get_interp_coeff_lin \n
!!        Routine to compute interpolation coefficients for linear
!!        interpolations (linear, bi-linear, tri-linear)
!!
MODULE PDAFomi_obs_op

  USE PDAFomi_obs_f, ONLY: obs_f, PDAFomi_gather_obsstate, debug

CONTAINS

!-------------------------------------------------------------------------------
!> observation operator for data at grid points
!!
!! Application of observation operator for the case that 
!! model variables are observerved at model grid points. 
!!
!! For this case INIT_DIM_OBS_F will prepare the index 
!! array thisobs%id_obs_p containing the information which 
!! elements of the  PE-local state vector contain the
!! observed values.
!!
!! The routine is called by all filter processes. It first
!! selects the observed elements for a PE-local domain. 
!! Afterwards, the values are gathered into the full vector
!! using PDAFomi_gather_obsstate.
!!
!! The routine has to fill the part of the full observation 
!! vector OBS_F_ALL that represents the current observation
!! type. The routine first applied the observation operator
!! for the current observation type and the calls
!! PDAFomi_gather_obsstate to gather the observation over
!! all processes and fills OBS_F_ALL.
!!
!! The routine has to be called by all filter processes.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_obs_op_gridpoint(thisobs, state_p, obs_f_all)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    REAL, INTENT(in)    :: state_p(:)      !< PE-local model state (dim_p)
    REAL, INTENT(inout) :: obs_f_all(:)    !< Full observed state for all observation types (nobs_f_all)

! *** Local variables ***
    INTEGER :: i                           ! Counter
    REAL, ALLOCATABLE :: ostate_p(:)       ! local observed part of state vector


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

    doassim: IF (thisobs%doassim == 1) THEN

       ! Consistency check
       IF (.NOT.ALLOCATED(thisobs%id_obs_p)) THEN
          WRITE (*,*) 'ERROR: PDAFomi_obs_op_gridpoint - thisobs%id_obs_p is not allocated'
       END IF

       ! Print debug information
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_gridpoint -- START'
          WRITE (*,*) '++ OMI-debug: ', debug, '  PDAFomi_obs_op_gridpoint -- Process-local selection'
          WRITE (*,*) '++ OMI-debug obs_op_gridpoint:', debug, 'thisobs%dim_obs_p', thisobs%dim_obs_p
          WRITE (*,*) '++ OMI-debug obs_op_gridpoint:', debug, 'thisobs%id_obs_p', thisobs%id_obs_p
       END IF

       ! *** PE-local: Initialize observed part state vector

       IF (thisobs%dim_obs_p>0) THEN
          ALLOCATE(ostate_p(thisobs%dim_obs_p))
       ELSE
          ALLOCATE(ostate_p(1))
       END IF

       DO i = 1, thisobs%dim_obs_p
          ostate_p(i) = state_p(thisobs%id_obs_p(1, i)) 
       ENDDO

       ! *** Global: Gather full observed state vector
       CALL PDAFomi_gather_obsstate(thisobs, ostate_p, obs_f_all)

       ! *** Clean up
       DEALLOCATE(ostate_p)

       ! Print debug information
       IF (debug>0) &
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_gridpoint -- END'

    END IF doassim

  END SUBROUTINE PDAFomi_obs_op_gridpoint




!-------------------------------------------------------------------------------
!> Observation operator for averaging grid point values
!!
!! Application of observation operator for the case that 
!! the observation value is given as the average of model
!! grid point values.
!! For this case INIT_DIM_OBS_F will prepare the index
!! array thisobs%id_obs_p that contains several rows holding
!! the indices of the state vector elements which are to
!! be averaged to represent an observation.
!!
!! The routine is called by all filter processes, 
!! and the operation has to be performed by each 
!! these processes for its PE-local domain before the 
!! information from all PEs is gathered.
!!
!! The routine has to fill the part of the full observation 
!! vector OBS_F_ALL that represents the current observation
!! type. The routine first applied the observation operator
!! for the current observation type and the calls
!! PDAFomi_gather_obsstate to gather the observation over
!! all processes and fills OBS_F_ALL.
!!
!! The routine has to be called by all filter processes.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_obs_op_gridavg(thisobs, nrows, state_p, obs_f_all)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    INTEGER, INTENT(in) :: nrows           !< Number of values to be averaged
    REAL, INTENT(in)    :: state_p(:)      !< PE-local model state (dim_p)
    REAL, INTENT(inout) :: obs_f_all(:)    !< Full observed state for all observation types (nobs_f_all)

! *** Local variables ***
    INTEGER :: i, row                      ! Counter
    REAL, ALLOCATABLE :: ostate_p(:)       ! local observed part of state vector
    REAL :: rrows                          ! Real-value for nrows


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

    doassim: IF (thisobs%doassim == 1) THEN

       ! Print debug information
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_gridavg -- START'
          WRITE (*,*) '++ OMI-debug: ', debug, '  PDAFomi_obs_op_gridavg -- Process-local averaging'
          WRITE (*,*) '++ OMI-debug obs_op_gridavg:', debug, 'thisobs%dim_obs_p', thisobs%dim_obs_p
          WRITE (*,*) '++ OMI-debug obs_op_gridavg:', debug, 'number of points to average', nrows
          WRITE (*,*) '++ OMI-debug obs_op_gridavg:', debug, 'thisobs%id_obs_p', thisobs%id_obs_p
       END IF

       ! Consistency check
       IF (.NOT.ALLOCATED(thisobs%id_obs_p)) THEN
          WRITE (*,*) 'ERROR: PDAFomi_obs_op_gridavg - thisobs%id_obs_p is not allocated'
       END IF

       ! *** PE-local: Initialize observed part state vector by averaging

       IF (thisobs%dim_obs_p>0) THEN
          ALLOCATE(ostate_p(thisobs%dim_obs_p))
       ELSE
          ALLOCATE(ostate_p(1))
       END IF

       rrows = REAL(nrows)

       DO i = 1, thisobs%dim_obs_p
          ostate_p(i) = 0.0
          DO row = 1, nrows
             ostate_p(i) = ostate_p(i) + state_p(thisobs%id_obs_p(row,i))
          END DO
          ostate_p(i) = ostate_p(i) / rrows
       ENDDO

       ! *** Global: Gather full observed state vector
       CALL PDAFomi_gather_obsstate(thisobs, ostate_p, obs_f_all)

       ! *** Clean up
       DEALLOCATE(ostate_p)

       ! Print debug information
       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_gridavg -- END'

    END IF doassim

  END SUBROUTINE PDAFomi_obs_op_gridavg




!-------------------------------------------------------------------------------
!> Observation operator for linear interpolation
!!
!! Application of observation operator for the case that the
!! observation value is given as the interpolation using
!! pre-computed coefficients. 
!!
!! For this case INIT_DIM_OBS_F will prepare the index 
!! array thisobs%id_obs_p containing the information which 
!! elements of the  PE-local state vector contain the
!! observed values. Further the array thisobs%icoeff_p is
!! prepared which contains the interpolation coefficients. 
!! This can be prepared using a help routine like
!! get_interp_coeff_tri.
!!
!! The routine is called by all filter processes. It first
!! selects the observed elements for a PE-local domain. 
!! Afterwards, the values are gathered into the full vector
!! using PDAFomi_gather_obsstate.
!!
!! The routine has to fill the part of the full observation 
!! vector OBS_F_ALL that represents the current observation
!! type. The routine first applied the observation operator
!! for the current observation type and the calls
!! PDAFomi_gather_obsstate to gather the observation over
!! all processes and fills OBS_F_ALL.
!!
!! The routine has to be called by all filter processes.
!!
!! __Revision history:__
!! * 2019-12 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_obs_op_interp_lin(thisobs, nrows, state_p, obs_f_all)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    INTEGER, INTENT(in) :: nrows           !< Number of values to be averaged
    REAL, INTENT(in)    :: state_p(:)      !< PE-local model state (dim_p)
    REAL, INTENT(inout) :: obs_f_all(:)    !< Full observed state for all observation types (nobs_f_all)

! *** Local variables ***
    INTEGER :: i, row                      ! Counters
    REAL, ALLOCATABLE :: ostate_p(:)       ! local observed part of state vector


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

    doassim: IF (thisobs%doassim == 1) THEN

       ! Print debug information
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_interp_lin -- START'
          WRITE (*,*) '++ OMI-debug: ', debug, '  PDAFomi_obs_op_interp_lin -- Process-local interpolation'
          WRITE (*,*) '++ OMI-debug obs_op_interp_lin:', debug, 'thisobs%dim_obs_p', thisobs%dim_obs_p
          WRITE (*,*) '++ OMI-debug obs_op_interp_lin:', debug, 'number of points in interpolation', nrows
          WRITE (*,*) '++ OMI-debug obs_op_interp_lin:', debug, 'thisobs%id_obs_p', thisobs%id_obs_p
          WRITE (*,*) '++ OMI-debug obs_op_interp_lin:', debug, 'thisobs%icoeff_p', thisobs%icoeff_p
       END IF

       ! Check if required arrays are allocated (assuming that they are initialzed in this case)
       IF (.NOT.ALLOCATED(thisobs%id_obs_p)) THEN
          WRITE (*,*) 'ERROR: PDAFomi_obs_op_interp_lin - thisobs%id_obs_p is not allocated'
       END IF
       IF (.NOT.ALLOCATED(thisobs%icoeff_p)) THEN
          WRITE (*,*) 'ERROR: PDAFomi_obs_op_interp_lin - thisobs%icoeff_p is not allocated'
       END IF

       ! *** PE-local: Initialize observed part state vector by weighted averaging

       IF (thisobs%dim_obs_p>0) THEN
          ALLOCATE(ostate_p(thisobs%dim_obs_p))
       ELSE
          ALLOCATE(ostate_p(1))
       END IF

       DO i = 1, thisobs%dim_obs_p
          ostate_p(i) = 0.0
          DO row = 1, nrows
             ostate_p(i) = ostate_p(i) + thisobs%icoeff_p(row,i)*state_p(thisobs%id_obs_p(row,i))
          END DO
       ENDDO

       ! *** Global: Gather full observed state vector
       CALL PDAFomi_gather_obsstate(thisobs, ostate_p, obs_f_all)

       ! *** Clean up
       DEALLOCATE(ostate_p)

       ! Print debug information
       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_interp_lin -- END'

    END IF doassim

  END SUBROUTINE PDAFomi_obs_op_interp_lin




!-------------------------------------------------------------------------------
!> observation operator for the case that observations belong to other compartment
!!
!! Application of observation operator for the case that 
!! model variables are observerved in another compartment
!! only. Thus DIM_OBS_P of the current compartment is 0.
!! Accordingly, this observation operator only performs
!! the gather operation to obtain the full observations.
!!
!! The routine has to fill the part of the full observation 
!! vector OBS_F_ALL that represents the current observation
!! type. The routine first applied the observation operator
!! for the current observation type and the calls
!! PDAFomi_gather_obsstate to gather the observation over
!! all processes and fills OBS_F_ALL.
!!
!! The routine has to be called by all filter processes.
!!
!! __Revision history:__
!! * 2020-04 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_obs_op_gatheronly(thisobs, state_p, obs_f_all)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    REAL, INTENT(in)    :: state_p(:)      !< PE-local model state (dim_p)
    REAL, INTENT(inout) :: obs_f_all(:)    !< Full observed state for all observation types (nobs_f_all)

! *** Local variables ***
    REAL, ALLOCATABLE :: ostate_p(:)       ! local observed part of state vector
    REAL :: rdummy                         ! dummy variable to prevent compiler warning


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

    ! Initialize dummy to prevent compiler warning
    rdummy = state_p(1)

    doassim: IF (thisobs%doassim == 1) THEN

       ! Print debug information
       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_gatheronly -- START'

       ! *** PE-local: Nothing to be done!

       ALLOCATE(ostate_p(1))
       ostate_p = 0.0

       ! *** Global: Gather full observed state vector
       CALL PDAFomi_gather_obsstate(thisobs, ostate_p, obs_f_all)

       ! *** Clean up
       DEALLOCATE(ostate_p)

       ! Print debug information
       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_gatheronly -- END'

    END IF doassim

  END SUBROUTINE PDAFomi_obs_op_gatheronly




!-------------------------------------------------------------------------------
!> Helper routine: Initialize interpolation coefficients in triangle
!!
!! The routine computes the coefficients for triangular interpolation
!! as barycentric coordinates.
!! The computation is done for one observation given the 
!! observation coordinates (OC) as well as the coordinates of the 
!! grid points (GPC). In GPC each row contains the coordinates
!! for one grid point. Thus the first index determines the grid point,
!! while the second the coordinates of this grid point
!!
!! __Revision history:__
!! * 2019-12 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_get_interp_coeff_tri(gpc, oc, icoeff)

    IMPLICIT NONE

! *** Arguments ***
    REAL, INTENT(in)    :: gpc(:,:)    !< Coordinates of grid points; dim(3,2)
                                       ! 3 rows; each containing lon and lat coordinates
    REAL, INTENT(in)    :: oc(:)       !< Coordinates of observation; dim(2)
    REAL, INTENT(inout) :: icoeff(:)   !< Interpolation coefficients; dim(3)

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

  END SUBROUTINE PDAFomi_get_interp_coeff_tri



!-------------------------------------------------------------------------------
!> Helper routine: Initialize linear interpolation coefficients in 1D
!!
!! The routine computes the coefficients for linear interpolation
!! in 1 dimensions.
!! The computation is done for one observation given the 
!! observation coordinates (OC) as well as the coordinates of the 
!! grid points (GPC). 
!!
!! __Revision history:__
!! * 2019-12 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_get_interp_coeff_lin1D(gpc, oc, icoeff)

    IMPLICIT NONE

! *** Arguments ***
    REAL, INTENT(in)    :: gpc(:)      !< Coordinates of grid points (dim=2)
    REAL, INTENT(in)    :: oc          !< Coordinates of observation
    REAL, INTENT(inout) :: icoeff(:)   !< Interpolation coefficients (dim=2)


! ****************************************************************
! *** Compute linear interpolation coefficients in 1 dimension ***
! ****************************************************************

    icoeff(1) = (gpc(2) - oc) / (gpc(2) - gpc(1))
    icoeff(2) = (oc - gpc(1)) / (gpc(2) - gpc(1))

  END SUBROUTINE PDAFomi_get_interp_coeff_lin1D



!-------------------------------------------------------------------------------
!> Helper routine: Initialize linear interpolation coefficients
!!
!! The routine computes the coefficients for linear interpolation
!! in 1, 2, or 3 dimensions.
!! The computation is done for one observation given the 
!! observation coordinates (OC) as well as the coordinates of the 
!! grid points (GPC). In GPC each row contains the coordinates
!! for one grid point.  
!!
!! Setup of GPC:
!! The first index is specifies the grid point, while the second the
!! coordinate
!! * For n_dim=X only the first X coordinate values are used
!! * The ordering of the coordinates and coefficient is the following:
!!
!!                       (7)------(8) 
!!                       /|       /|    with
!!                     (5)+-----(6)|       - column 1
!!                      | |      | |       / column 2
!!                      |(3)-----+(4)      | column 3
!!                      |/       |/
!!                     (1) ---- (2)
!!
!!   thus gpc(1,1)/=gpc(2,1), gpc(1,2)/=gpc(3,2), gpc(1,3)/=gpc(5,3)
!!   but gpc(1,1)=gpc(3,1)=gpc(5,1), gpc(1,2)=gpc(2,2)=gpc(5,2), 
!!   gpc(1,3)=gpc(2,3)=gpc(3,3)
!! * For bi-linear interpolation only the coordinates for grid
!!   points 1, 2, and 3 are used to compute the coefficients
!! * For tri-linear interpolation only the coordinates for grid
!!   points 1, 2, 3, and 5 are used to compute the coefficients
!! * (for bi-linear interpolation gpc only needs to have length 3
!!   for tri-linear the length 5)
!! * num_gp=2 for n_dim=1; num_gp=4 for n_dim=2; num_gp=8 for n_dim=3
!!   is required
!!
!! __Revision history:__
!! * 2019-12 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_get_interp_coeff_lin(num_gp, n_dim, gpc, oc, icoeff)

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: num_gp         !< Length of icoeff
    INTEGER, INTENT(in) :: n_dim          !< Number of dimensions in interpolation
    REAL, INTENT(in)    :: gpc(:,:)       !< Coordinates of grid points
    REAL, INTENT(in)    :: oc(:)          !< Coordinates of observation
    REAL, INTENT(inout) :: icoeff(:)      !< Interpolation coefficients (num_gp)

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
            'PDAFomi', 'ERROR: get_interp_coeff_lin - NUM_GP=4 required!'
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
            'PDAFomi', 'ERROR: get_interp_coeff_lin - NUM_GP=8 required!'
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

  END SUBROUTINE PDAFomi_get_interp_coeff_lin

!-------------------------------------------------------------------------------
!> Adjoint observation operator for data at grid points
!!
!! Application of adjoint observation operator for the case 
!! that model variables are observerved at model grid points. 
!!
!! For this case INIT_DIM_OBS_F will prepare the index 
!! array thisobs%id_obs_p containing the information which 
!! elements of the  PE-local state vector contain the
!! observed values.
!!
!! The routine is called by all filter processes. It first
!! selects the observed elements for a PE-local domain. 
!! Afterwards, the values are gathered into the full vector
!! using PDAFomi_gather_obsstate.
!!
!! The routine has to fill the part of the full observation 
!! vector OBS_F_ALL that represents the current observation
!! type. The routine first applied the observation operator
!! for the current observation type and the calls
!! PDAFomi_gather_obsstate to gather the observation over
!! all processes and fills OBS_F_ALL.
!!
!! The routine has to be called by all filter processes.
!!
!! __Revision history:__
!! * 2021-04 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_obs_op_adj_gridpoint(thisobs, obs_f_all, state_p)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    REAL, INTENT(in)    :: obs_f_all(:)    !< Full observed state for all observation types (nobs_f_all)
    REAL, INTENT(inout) :: state_p(:)      !< PE-local model state (dim_p)

! *** Local variables ***
    INTEGER :: i                           ! Counter


! **************************************************
! *** Perform application of adjoint observation ***
! *** operator H^T on vector or matrix column    ***
! **************************************************

    doassim: IF (thisobs%doassim == 1) THEN

       ! Consistency check
       IF (.NOT.ALLOCATED(thisobs%id_obs_p)) THEN
          WRITE (*,*) 'ERROR: PDAFomi_obs_op_adj_gridpoint - thisobs%id_obs_p is not allocated'
       END IF

       ! Print debug information
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_adj_gridpoint -- START'
          WRITE (*,*) '++ OMI-debug: ', debug, '  PDAFomi_obs_op_adj_gridpoint -- Process-local selection'
          WRITE (*,*) '++ OMI-debug obs_op_adj_gridpoint:', debug, 'thisobs%dim_obs_p', thisobs%dim_obs_p
          WRITE (*,*) '++ OMI-debug obs_op_adj_gridpoint:', debug, 'thisobs%id_obs_p', thisobs%id_obs_p
          WRITE (*,*) '++ OMI-debug obs_op_adj_gridpoint:', debug, 'thisobs%off_obs_f', thisobs%off_obs_f
       END IF

       ! *** PE-local: Apply adjoint observation operator

       DO i = 1, thisobs%dim_obs_p
          state_p(thisobs%id_obs_p(1, i)) &
               = state_p(thisobs%id_obs_p(1, i)) + obs_f_all(thisobs%off_obs_f+i)
       ENDDO

       ! Print debug information
       IF (debug>0) &
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_adj_gridpoint -- END'

    END IF doassim

  END SUBROUTINE PDAFomi_obs_op_adj_gridpoint




!-------------------------------------------------------------------------------
!> Adjoint observation operator for averaging grid point values
!!
!! Application of adjoint observation operator for the case 
!! that the observation value is given as the average of
!! model grid point values.
!!
!! For this case INIT_DIM_OBS_F will prepare the index
!! array thisobs%id_obs_p that contains several rows holding
!! the indices of the state vector elements which are to
!! be averaged to represent an observation.
!!
!! The routine is called by all filter processes, 
!! and the operation has to be performed by each 
!! these processes for its PE-local domain before the 
!! information from all PEs is gathered.
!!
!! The routine has to fill the part of the full observation 
!! vector OBS_F_ALL that represents the current observation
!! type. The routine first applied the observation operator
!! for the current observation type and the calls
!! PDAFomi_gather_obsstate to gather the observation over
!! all processes and fills OBS_F_ALL.
!!
!! The routine has to be called by all filter processes.
!!
!! __Revision history:__
!! * 2021-12 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_obs_op_adj_gridavg(thisobs, nrows, obs_f_all, state_p)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    INTEGER, INTENT(in) :: nrows           !< Number of values to be averaged
    REAL, INTENT(in)    :: obs_f_all(:)    !< Full observed state for all observation types (nobs_f_all)
    REAL, INTENT(inout) :: state_p(:)      !< PE-local model state (dim_p)

! *** Local variables ***
    INTEGER :: i, row                      ! Counter
    REAL :: rrows                          ! Real-value for nrows


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

    doassim: IF (thisobs%doassim == 1) THEN

       ! Print debug information
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_adj_gridavg -- START'
          WRITE (*,*) '++ OMI-debug: ', debug, '  PDAFomi_obs_op_adj_gridavg -- Process-local averaging'
          WRITE (*,*) '++ OMI-debug obs_op_adj_gridpoint:', debug, 'thisobs%dim_obs_p', thisobs%dim_obs_p
          WRITE (*,*) '++ OMI-debug obs_op_adj_gridpoint:', debug, 'number of points to average', nrows
          WRITE (*,*) '++ OMI-debug obs_op_adj_gridpoint:', debug, 'thisobs%id_obs_p', thisobs%id_obs_p
       END IF

       ! Consistency check
       IF (.NOT.ALLOCATED(thisobs%id_obs_p)) THEN
          WRITE (*,*) 'ERROR: PDAFomi_obs_op_adj_gridavg - thisobs%id_obs_p is not allocated'
       END IF

       ! *** PE-local: Initialize observed part state vector by averaging

       rrows = REAL(nrows)

       DO i = 1, thisobs%dim_obs_p
          DO row = 1, nrows
             state_p(thisobs%id_obs_p(row, i)) &
                  = state_p(thisobs%id_obs_p(row, i)) + obs_f_all(thisobs%off_obs_f+i) / rrows
          end DO
       ENDDO

       ! Print debug information
       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_adj_gridavg -- END'

    END IF doassim

  END SUBROUTINE PDAFomi_obs_op_adj_gridavg




!-------------------------------------------------------------------------------
!> Observation operator for linear interpolation
!!
!! Application of adjoint observation operator for the case
!! that the observation value is given as the interpolation
!! using pre-computed coefficients. 
!!
!! For this case INIT_DIM_OBS_F will prepare the index 
!! array thisobs%id_obs_p containing the information which 
!! elements of the  PE-local state vector contain the
!! observed values. Further the array thisobs%icoeff_p is
!! prepared which contains the interpolation coefficients. 
!! This can be prepared using a help routine like
!! get_interp_coeff_tri.
!!
!! The routine is called by all filter processes. It first
!! selects the observed elements for a PE-local domain. 
!! Afterwards, the values are gathered into the full vector
!! using PDAFomi_gather_obsstate.
!!
!! The routine has to fill the part of the full observation 
!! vector OBS_F_ALL that represents the current observation
!! type. The routine first applied the observation operator
!! for the current observation type and the calls
!! PDAFomi_gather_obsstate to gather the observation over
!! all processes and fills OBS_F_ALL.
!!
!! The routine has to be called by all filter processes.
!!
!! __Revision history:__
!! * 2021-12 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_obs_op_adj_interp_lin(thisobs, nrows, obs_f_all, state_p)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    INTEGER, INTENT(in) :: nrows           !< Number of values to be averaged
    REAL, INTENT(in)    :: obs_f_all(:)    !< Full observed state for all observation types (nobs_f_all)
    REAL, INTENT(inout) :: state_p(:)      !< PE-local model state (dim_p)

! *** Local variables ***
    INTEGER :: i, row                      ! Counters


! **************************************************
! *** Perform application of adjoint observation ***
! *** operator H^T on vector or matrix column    ***
! **************************************************

    doassim: IF (thisobs%doassim == 1) THEN

       ! Print debug information
       IF (debug>0) THEN
          WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_adj_interp_lin -- START'
          WRITE (*,*) '++ OMI-debug: ', debug, '  PDAFomi_obs_op_adj_interp_lin -- Process-local interpolation'
          WRITE (*,*) '++ OMI-debug obs_op_adj_interp_lin:', debug, 'thisobs%dim_obs_p', thisobs%dim_obs_p
          WRITE (*,*) '++ OMI-debug obs_op_adj_interp_lin:', debug, 'number of points in interpolation', nrows
          WRITE (*,*) '++ OMI-debug obs_op_adj_interp_lin:', debug, 'thisobs%id_obs_p', thisobs%id_obs_p
          WRITE (*,*) '++ OMI-debug obs_op_adj_interp_lin:', debug, 'thisobs%icoeff_p', thisobs%icoeff_p
       END IF

       ! Check if required arrays are allocated (assuming that they are initialzed in this case)
       IF (.NOT.ALLOCATED(thisobs%id_obs_p)) THEN
          WRITE (*,*) 'ERROR: PDAFomi_obs_op_adj_interp_lin - thisobs%id_obs_p is not allocated'
       END IF
       IF (.NOT.ALLOCATED(thisobs%icoeff_p)) THEN
          WRITE (*,*) 'ERROR: PDAFomi_obs_op_adj_interp_lin - thisobs%icoeff_p is not allocated'
       END IF

       ! *** PE-local: Apply adjoint observation operator

       DO i = 1, thisobs%dim_obs_p
          DO row = 1, nrows
             state_p(thisobs%id_obs_p(row, i)) &
                  = state_p(thisobs%id_obs_p(row, i)) + thisobs%icoeff_p(row,i)*obs_f_all(thisobs%off_obs_f+i)
          end DO
       ENDDO

       ! Print debug information
       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_adj_interp_lin -- END'

    END IF doassim

  END SUBROUTINE PDAFomi_obs_op_adj_interp_lin




!-------------------------------------------------------------------------------
!> adjoint observation operator for the case that observations belong
!! to other compartment
!!
!! Application of adjoint observation operator for the case that 
!! model variables are observerved in another compartment
!! only. Thus DIM_OBS_P of the current compartment is 0.
!! Accordingly, this observation operator only performs
!! the gather operation to obtain the full observations.
!!
!! The routine does nothing, since there are no observations
!! on the compartment
!!
!! The routine has to be called by all filter processes.
!!
!! __Revision history:__
!! * 2021-12 - Lars Nerger - Initial code from restructuring observation routines
!! * Later revisions - see repository log
!!
  SUBROUTINE PDAFomi_obs_op_adj_gatheronly(thisobs, obs_f_all, state_p)

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    REAL, INTENT(in)    :: state_p(:)      !< PE-local model state (dim_p)
    REAL, INTENT(inout) :: obs_f_all(:)    !< Full observed state for all observation types (nobs_f_all)

! *** Local variables ***


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

    doassim: IF (thisobs%doassim == 1) THEN

       ! Print debug information
       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_gatheronly -- START'

       ! *** Nothing to be done!

       ! Print debug information
       IF (debug>0) &
            WRITE (*,*) '++ OMI-debug: ', debug, 'PDAFomi_obs_op_gatheronly -- END'

    END IF doassim

  END SUBROUTINE PDAFomi_obs_op_adj_gatheronly

END MODULE PDAFomi_obs_op
