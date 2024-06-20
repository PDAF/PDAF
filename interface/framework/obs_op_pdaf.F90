!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz and Guowei He (Forschungszentrum Juelich GmbH)
!
!This file is part of TSMP-PDAF
!
!TSMP-PDAF is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TSMP-PDAF is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU LesserGeneral Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public License
!along with TSMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------
!obs_op_pdaf.F90: TSMP-PDAF implementation of routine
!                 'obs_op_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: obs_op_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: obs_op_pdaf --- Implementation of observation operator
!
! !INTERFACE:
SUBROUTINE obs_op_pdaf(step, dim_p, dim_obs_p, state_p, m_state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/EnKF/ETKF/ESTKF
!
! The routine is called during the analysis step.
! It has to perform the operation of the
! observation operator acting on a state vector.
! For domain decomposition, the action is on the
! PE-local sub-domain of the state and has to
! provide the observed sub-state for the PE-local
! domain.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
   USE mod_assimilation, &
        ONLY: obs_index_p, &
#ifndef CLMSA
#ifndef OBS_ONLY_CLM
        sc_p, &
#endif
#endif
        obs_interp_indices_p, &
        obs_interp_weights_p
   use mod_tsmp, &
       only: obs_interp_switch, &
       soilay, &
       soilay_fortran, &
       nz_glob, &
       da_crns_depth_tol, &
       crns_flag
!      tcycle

   USE, INTRINSIC :: iso_c_binding

#if defined CLMSA
   USE enkf_clm_mod, & 
        ONLY : clm_varsize, clm_paramarr, clmupdate_swc, clmupdate_T
#endif
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step               ! Currrent time step
  INTEGER, INTENT(in) :: dim_p              ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs_p          ! Dimension of observed state
  REAL, INTENT(in)    :: state_p(dim_p)     ! PE-local model state
  REAL, INTENT(out) :: m_state_p(dim_obs_p) ! PE-local observed state
  integer :: i, j, k
  integer :: icorner
  logical :: lpointobs       !If true: no special observation; use point observation
! !CALLING SEQUENCE:
! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT   (as U_obs_op)
! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
!EOP

! *** local variables ***

! hcp test with hardcoding variable declaration
real(8), dimension(:), allocatable :: soide !soil depth
!real(8), dimension(0:12), parameter :: &
! soide=(/0.d0,  0.02d0,  0.05d0,  0.1d0,  0.17d0, 0.3d0,  0.5d0, &
!                0.8d0,   1.3d0,   2.d0,  3.d0, 5.d0,  12.d0/) !soil depth

real(8) :: tot, avesm, avesm_temp, Dp
integer :: nsc
! end of hcp 


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

! If no special observation operator is compiled, use point observations
lpointobs = .true.

#if defined CLMSA
if (clmupdate_T.EQ.1) then

  lpointobs = .false.

  DO i = 1, dim_obs_p
    ! Equation for LST computation from Kustas2009, Eq(7)
    ! http://dx.doi.org/10.1016/j.agrformet.2009.05.016
     m_state_p(i) &
    = (exp(-0.5*clm_paramarr(obs_index_p(i))) &
                     *state_p(obs_index_p(i))**4 & 
       +(1.-exp(-0.5*clm_paramarr(obs_index_p(i)))) &
                     *state_p(clm_varsize+obs_index_p(i))**4)**0.25
  END DO
!  write(*,*) 'now is cycle ', tcycle
!  write(*,*) 'obs_p =', obs_p(:)
!  write(*,*) 'model LST', m_state_p(:)
!  write(*,*) 'TG', state_p(obs_index_p(:))
!  write(*,*) 'TV', state_p(clm_varsize+obs_index_p(:))
endif
#endif


#ifndef CLMSA
#ifndef OBS_ONLY_CLM
 if (crns_flag.EQ.1) then
    !Schroen et al HESS 2017 modelled CRNS averaging
    lpointobs = .false.
     call C_F_POINTER(soilay,soilay_fortran,[nz_glob])
     Allocate(soide(0:nz_glob))
     soide(0)=0.d0
     do i=1,nz_glob
       soide(i)=soide(i-1)+soilay_fortran(nz_glob-i+1) 
     enddo
     do i = 1, dim_obs_p

       !Initial average soil moisture for 1st iteration
       avesm=0.d0
       do j=1,nz_glob
            avesm=avesm+(soide(j)-soide(j-1))*state_p(sc_p(j,i))/soide(nz_glob)
       enddo
       avesm_temp=0.d0

       !iteration
       do while (abs(avesm-avesm_temp)/avesm .GE. da_crns_depth_tol)
          !Averaging, conventional profile, Schroen et al HESS 2017 Eq. (3)
          avesm_temp=avesm
          Dp=0.058d0/(avesm+0.0829d0)

          !Sum weight*soil_moisture
          avesm=0.d0; nsc=nz_glob
          do j=1,nz_glob
             if ((soide(j-1).LT.Dp).AND.(Dp.LE.soide(j))) then
               nsc=j
             endif
          enddo
          do j=1, nsc-1
              avesm=avesm+(1.d0-0.5d0*(soide(j)+soide(j-1))/Dp)*(soide(j)-soide(j-1)) &
                    *state_p(sc_p(j,i))/Dp
          enddo
          avesm=avesm+(1.d0-0.5d0*(Dp+soide(nsc-1))/Dp)*(Dp-soide(nsc-1)) &
             *state_p(sc_p(nsc,i))/Dp

          !Sum weight
          tot=0.d0
          do j=1, nsc-1
              tot =   tot+(1.d0-0.5d0*(soide(j)+soide(j-1))/Dp)*(soide(j)-soide(j-1)) &
                                                    /Dp
          enddo
          tot  =  tot+(1.d0-0.5d0*(Dp+soide(nsc-1))/Dp)*(Dp-soide(nsc-1)) &
                                               /Dp

          avesm=avesm/tot
       enddo
       m_state_p(i)=avesm
     enddo
     deallocate(soide)
 end if
#endif
#endif

 if(obs_interp_switch == 1) then

      lpointobs = .false.

      do i = 1, dim_obs_p

          m_state_p(i) = 0
          do icorner = 1, 4
              m_state_p(i) = m_state_p(i) + state_p(obs_interp_indices_p(i,icorner)) * obs_interp_weights_p(i,icorner)
          enddo

      enddo

  end if

  if(lpointobs) then

  DO i = 1, dim_obs_p
     m_state_p(i) = state_p(obs_index_p(i))
  END DO
      
  end if

END SUBROUTINE obs_op_pdaf
