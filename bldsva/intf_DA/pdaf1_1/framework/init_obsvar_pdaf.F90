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
!init_obsvar_pdaf.F90: TerrSysMP-PDAF implementation of routine
!                     'init_obsvar_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: init_obsvar_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: init_obsvar_pdaf --- Get mean observation error variance
!
! !INTERFACE:
SUBROUTINE init_obsvar_pdaf(step, dim_obs_p, obs_p, meanvar)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF

! This routine will only be called, if the adaptive
! forgetting factor feature is used. Please note that
! this is an experimental feature.
!
! The routine is called in global filters (like SEIK)
! during the analysis or in local filters (e.g. LSEIK)
! before the loop over local analysis domains 
! by the routine PDAF\_set\_forget that estimates an 
! adaptive forgetting factor.  The routine has to 
! initialize the mean observation error variance.  
! For global filters this should be the global mean,
! while for local filters it should be the mean for the
! PE-local  sub-domain.  (See init\_obsvar\_l_pdaf()
! for a localized variant for local filters.)
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
!USE mpi
USE mod_assimilation, &
    ONLY: rms_obs, pressure_obserr_p, clm_obserr_p
USE mod_parallel_pdaf, &
    ONLY: COMM_filter, MPIerr, MPI_REAL8, MPI_SUM, npes_filter 
USE mod_parallel_model, ONLY: model
USE mod_tsmp, &
#if defined CLMSA
       ONLY: tag_model_clm
#else
       ONLY: tag_model_parflow
#endif


  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step          ! Current time step
  INTEGER, INTENT(in) :: dim_obs_p     ! PE-local dimension of observation vector
  REAL, INTENT(in) :: obs_p(dim_obs_p) ! PE-local observation vector
  REAL, INTENT(out)   :: meanvar       ! Mean observation error variance

  ! local variables
  REAL :: meanvar_p                    ! PE-local Mean observation error variance 
  REAL :: sum_p                        ! PE-local sum of observation error variance
  INTEGER :: i, count
! !CALLING SEQUENCE:
! Called by: PDAF_set_forget    (as U_init_init_obs_covar)
!EOP


! *****************************
! *** Compute mean variance ***
! *****************************

  WRITE (*,*) 'TEMPLATE init_obsvar_pdaf.F90: Set mean observation variance here!'

  ! We assume that all observations have the same error.
  ! Thus, the mean variance is the error variance of each single observation.

!  meanvar = rms_obs ** 2

  ! Due to domain decomposition in our case the mean variance is computed
  ! for the full domain using the function MPI_Allreduce
#ifndef CLMSA
  if (model .eq. tag_model_parflow) then
     meanvar_p = 0
     sum_p = 0
     count = 0
     do i = 1, dim_obs_p
        if(pressure_obserr_p(i) /= 0) then
           sum_p = sum_p + pressure_obserr_p(i)
           count = count + 1 
        endif   
     enddo
     ! averaging the sum of observation errors with total no of non-zero observations
     meanvar_p = sum_p/count
     ! summing the average of observation errors and communicating it back to each rank
     call MPI_Allreduce(meanvar_p, meanvar, 1, MPI_REAL8, MPI_SUM, COMM_filter, MPIerr)
     ! to get the mean dividing the mean observation error by size of processors
     meanvar = meanvar/npes_filter 
  end if
#endif

#if defined CLMSA
  if(model .eq. tag_model_clm) then
     meanvar_p = 0
     sum_p = 0
     count = 0
     do i = 1, dim_obs_p
        if(clm_obserr_p(i) /= 0) then
           sum_p = sum_p + clm_obserr_p(i)
           count = count + 1 
        endif   
     enddo
     ! averaging the sum of observation errors with total no of non-zero observations
     meanvar_p = sum_p/count
     ! summing the average of observation errors and communicating it back to each rank
     call MPI_Allreduce(meanvar_p, meanvar, 1, MPI_REAL8, MPI_SUM, COMM_filter, MPIerr)
     ! to get the mean dividing the mean observation error by size of processors
     meanvar = meanvar/npes_filter    
  end if
#endif

!  meanvar = ?

END SUBROUTINE init_obsvar_pdaf
