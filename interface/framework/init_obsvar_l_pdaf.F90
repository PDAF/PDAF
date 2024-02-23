!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz, Guowei He and Mukund Pondkule (Forschungszentrum Juelich GmbH)
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
!init_obsvar_l_pdaf.F90: TSMP-PDAF implementation of routine
!                       'init_obsvar_l_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: init_obsvar_l_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
! !ROUTINE: init_obsvar_l_pdaf --- Get local mean observation error variance
!
! !INTERFACE:
SUBROUTINE init_obsvar_l_pdaf(domain_p, step, dim_obs_l, obs_l, meanvar_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! This routine will only be called, if the 
! local adaptive forgetting factor feature 
! is used. Please note that this is an 
! experimental feature.
!
! The routine is called in the loop over all
! local analysis domains during each analysis
! by the routine PDAF\_set\_forget\_local that 
! estimates a local adaptive forgetting factor.
! The routine has to initialize the mean observation 
! error variance for the current local analysis 
! domain.  (See init_obsvar() for a global variant.)
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
   USE mod_assimilation, &
        ONLY:rms_obs

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: domain_p      ! Current local analysis domain
  INTEGER, INTENT(in) :: step          ! Current time step
  INTEGER, INTENT(in) :: dim_obs_l     ! Local dimension of observation vector
  REAL, INTENT(in) :: obs_l(dim_obs_l) ! Local observation vector
  REAL, INTENT(out)   :: meanvar_l     ! Mean local observation error variance

! !CALLING SEQUENCE:
! Called by: PDAF_set_forget_local    (as U_init_obsvar_l)
!EOP


! ***********************************
! *** Compute local mean variance ***
! ***********************************

   meanvar_l = rms_obs ** 2
!  meanvar_l = ?

  ! Due to domain decomposition in our case the mean variance is computed
  ! for the full domain using the function MPI_Allreduce
!!$  ! Due to domain decomposition in our case the mean variance is computed
!!$  ! for the full domain using the function MPI_Allreduce
!!$#ifndef CLMSA
!!$  if (model .eq. tag_model_parflow) then
!!$     meanvar_p = 0
!!$     sum_p = 0
!!$     count = 0
!!$     do i = 1, dim_obs_p
!!$        if(pressure_obserr_p(i) /= 0) then
!!$           sum_p = sum_p + pressure_obserr_p(i)
!!$           count = count + 1 
!!$        endif   
!!$     enddo
!!$     ! averaging the sum of observation errors with total no of non-zero observations
!!$     meanvar_p = sum_p/count
!!$     ! summing the average of observation errors and communicating it back to each rank
!!$     call MPI_Allreduce(meanvar_p, meanvar, 1, MPI_REAL8, MPI_SUM, COMM_filter, MPIerr)
!!$     ! to get the mean dividing the mean observation error by size of processors
!!$     meanvar = meanvar/npes_filter 
!!$  end if
!!$#endif
!!$
!!$#if defined CLMSA
!!$  if(model .eq. tag_model_clm) then
!!$     meanvar_p = 0
!!$     sum_p = 0
!!$     count = 0
!!$     do i = 1, dim_obs_p
!!$        if(clm_obserr_p(i) /= 0) then
!!$           sum_p = sum_p + clm_obserr_p(i)
!!$           count = count + 1 
!!$        endif   
!!$     enddo
!!$     ! averaging the sum of observation errors with total no of non-zero observations
!!$     meanvar_p = sum_p/count
!!$     ! summing the average of observation errors and communicating it back to each rank
!!$     call MPI_Allreduce(meanvar_p, meanvar, 1, MPI_REAL8, MPI_SUM, COMM_filter, MPIerr)
!!$     ! to get the mean dividing the mean observation error by size of processors
!!$     meanvar = meanvar/npes_filter    
!!$  end if
!!$#endif

END SUBROUTINE init_obsvar_l_pdaf
