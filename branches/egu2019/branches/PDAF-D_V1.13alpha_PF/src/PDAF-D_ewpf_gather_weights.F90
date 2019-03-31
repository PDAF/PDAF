! Copyright (c) 2004-2018 Lars Nerger and Paul Kirchgessner
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
! !ROUTINE: PDAF_ewpf_gather_weights --- Gather distributed weights on filter PEs
!
! !INTERFACE:
SUBROUTINE PDAF_ewpf_gather_weights(dim_ens_p, weights, screen)

! !DESCRIPTION:
! If the ensemble integration is distributed over multiple
! model tasks, this routine collects the distributed
! weight information onto the processes that perform
! the analysis step (filterpe==.true.) (should be 1 for pf).
! IMPORTANT:
! The filter revs only the -log of the weight! 
! This is important for the normalization, since before 
! the normalization one has to compute the exponent of the weights
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2014-05 - Paul Kirchgessner - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_mod_filtermpi, &
       ONLY: mype_filter, mype_couple, npes_couple, filterpe, &
       all_dim_ens_l, all_dis_ens_l, COMM_couple, &
       MPI_REALTYPE, MPIerr, MPIstatus

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_ens_p    ! Size of ensemble diffent on model_pe and filte pe
  REAL, INTENT(inout) :: weights(dim_ens_p) ! PE-local state ensemble
  INTEGER, INTENT(in) :: screen       ! Verbosity flag  
!EOP

! Local variables
  REAL :: weights_tmp(dim_ens_p)
  INTEGER :: pe_rank, col_frst, col_last  ! Counters

  
! **********************************************
! *** Gather forecast ensemble on filter PEs ***
! **********************************************

  IF (filterpe .AND. mype_filter == 0 .AND. screen > 2) &
       WRITE (*, '(a, 5x, a)') 'PDAF', '--- Gather weights on filter PEs'

  ! *** Send from model PEs that are not filter PEs ***
  subensS: IF ( (.NOT.filterpe) ) THEN
     ! Send weight to couple PEs with rank 0
     weights_tmp(1) = weights(1)
        
     CALL MPI_SEND(weights_tmp(1), 1 , MPI_REALTYPE, 0, mype_couple, &
          COMM_couple, MPIerr)

     IF ((screen>2)) WRITE (*,*) 'PDAF: gather_weights - send weight' , &
          weights(1),' from rank(couple) ',mype_couple, &
          ' in couple task ', mype_filter+1
  END IF subensS

  ! *** Receive on filter PEs ***
  subensR: IF (filterpe ) THEN
     ! Receive particle weights on filter PEs
     taskloopB: DO pe_rank = 1, npes_couple - 1
        CALL MPI_recv(weights_tmp(pe_rank+1), 1, MPI_REALTYPE, &
             pe_rank, pe_rank, COMM_couple, MPIstatus, MPIerr)
        
        IF (screen > 2) &
             WRITE (*,*) 'PDAF: gather_weights - recv weight',&
             weights(mype_couple+1),' from rank(couple): ', &
             pe_rank
     END DO taskloopB
  END IF subensR

END SUBROUTINE PDAF_ewpf_gather_weights
