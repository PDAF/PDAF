!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz and Guowei He (Forschungszentrum Juelich GmbH)
!
!This file is part of TerrSysMP-PDAF
!
!TerrSysMP-PDAF is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TerrSysMP-PDAF is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU LesserGeneral Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public License
!along with TerrSysMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------
!distribute_stateinc_pdaf.F90: TerrSysMP-PDAF implementation of routine
!                              'distribute_stateinc_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: distribute_stateinc_pdaf.F90 1383 2013-05-03 12:26:53Z lnerger $
!BOP
!
! !ROUTINE: distribute_stateinc_pdaf --- Add analysis increment to model fields
!
! !INTERFACE:
SUBROUTINE distribute_stateinc_pdaf(dim_p, state_inc_p, new_forecast, steps)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/EnKF/SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! This subroutine is called during the forecast 
! phase of the filter from PDAF\_incremental
! supplying the analysis state increment.
! The routine has to compute the fraction of 
! the increment to be added to the model state 
! at each time step. Further, it has to transform 
! the increment vector into increments of the 
! fields of the model (typically available 
! trough a module).
!
! The routine is executed by each process that 
! is participating in the model integrations.
!
! !REVISION HISTORY:
! 2006-08 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! Dimension of PE-local state
  REAL, INTENT(in)    :: state_inc_p(dim_p) ! PE-local state vector
  INTEGER, INTENT(in) :: new_forecast    ! Flag for first call of each forecast
  INTEGER, INTENT(in) :: steps           ! number of time steps in forecast

! !CALLING SEQUENCE:
! Called by: PDAF_incremental   (as U_dist_stateinc)
!EOP  


! *******************************
! *** Prepare state increment ***
! *******************************

  IF (new_forecast > 0) THEN
     ! At begin of each forecast phase distribute full increment to
     ! all processes and compute increment per update step.
     ! (E.g., at each time step)
  ENDIF


! *************************************
! *** Add increment to model fields ***
!**************************************


END SUBROUTINE distribute_stateinc_pdaf
