!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz, Guowei He and Mukund Pondkule (Forschungszentrum Juelich GmbH)
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
!g2l_obs_pdaf.F90: TerrSysMP-PDAF implementation of routine
!                  'g2l_obs_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: g2l_obs_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: g2l_obs_pdaf --- Restrict an obs. vector to local analysis domain
!
! !INTERFACE:
SUBROUTINE g2l_obs_pdaf(domain_p, step, dim_obs_f, dim_obs_l, mstate_f, &
     mstate_l)

  ! !DESCRIPTION:
  ! User-supplied routine for PDAF.
  ! Used in the filters: LSEIK/LETKF/LESTKF
  !
  ! The routine is called during the analysis step
  ! on each of the local analysis domains.
  ! It has to restrict the full vector of all 
  ! observations required for the loop of localized 
  ! analyses on the PE-local domain to the current 
  ! local analysis domain.
  !
  ! !REVISION HISTORY:
  ! 2013-02 - Lars Nerger - Initial code
  ! Later revisions - see svn log
  !
  ! !USES:
  USE mod_assimilation, &
       ONLY: distance, local_dims_obs, obs_index_p, &
       dim_obs_p, global_to_local, obs_index_l, m_id_f
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, npes_filter, comm_filter

  IMPLICIT NONE

  ! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain_p   ! Current local analysis domain
  INTEGER, INTENT(in) :: step       ! Current time step
  INTEGER, INTENT(in) :: dim_obs_f  ! Dimension of full PE-local obs. vector
  INTEGER, INTENT(in) :: dim_obs_l  ! Local dimension of observation vector
  REAL, INTENT(in)    :: mstate_f(dim_obs_f)   ! Full PE-local obs. vector
  REAL, INTENT(out)   :: mstate_l(dim_obs_l)   ! Obs. vector on local domain

  ! !CALLING SEQUENCE:
  ! Called by: PDAF_lseik_analysis   (as U_g2l_obs)
  ! Called by: PDAF_lestkf_analysis  (as U_g2l_obs)
  ! Called by: PDAF_letkf_analysis   (as U_g2l_obs)
  !EOP

  ! *** local variables ***
  INTEGER :: i      ! Counter

  ! *******************************************************
  ! *** Perform localization of some observation vector *** 
  ! *** to the current local analysis domain.           ***
  ! *******************************************************

  ! Intialize Obs. vector on local domain
  mstate_l(:) = 0.0

  ! kuw
  do i=1,dim_obs_l
    !mstate_l(i) = mstate_f(obs_index_l(i)) 
    mstate_l(i) = mstate_f(m_id_f(obs_index_l(i)))
  end do

END SUBROUTINE g2l_obs_pdaf
