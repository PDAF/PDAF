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
!enkf_clm_5.F90: Wrapper functions for CLM 5 
!-------------------------------------------------------------------------------------------
module enkf_clm_5

#include <mpif.h>

  contains

!--------------------------------------------------------------------------
! clm_init() in the TSMP-PDAF uses copied code from the start of
! cime/src/drivers/mct/main/cime_driver.F90 (tag: `cime5_6_47`) until
! before the "call cime_run()" because TSMP-PDAF proceeds with the
! initialization of other component models and of PDAF before forward
! simulation is started. Further code from `cime_driver.F90` can be
! found in the subroutines `clm_advance` and `clm_finalize`.
!
! Additionally, two major TSMP-PDAF specific changes are added
! compared to the code from CIME: (1) cime_pre_init1 is modified to
! read additional inputs: the pdaf comm and the pdaf ID and (2) at the
! end of the init we call define_clm_statevec to finish the PDAF init.
!
! All TSMP-PDAF-specific code is highlighted with comments.
! --------------------------------------------------------------------------
subroutine clm_init(finname, pdaf_id, pdaf_max, mype) bind(C,name="clm_init")
  !----------------------------------------------------------------------------
  ! share code & libs
  !----------------------------------------------------------------------------
  use shr_kind_mod,  only : r8 => SHR_KIND_R8
  use shr_kind_mod,  only : i8 => SHR_KIND_I8
  use shr_kind_mod,  only : CS => SHR_KIND_CS
  use shr_sys_mod,   only : shr_sys_irtc, shr_sys_abort
  use perf_mod,      only : t_startf, t_adj_detailf, t_stopf, t_startstop_valsf
  use ESMF,          only : ESMF_Initialize, ESMF_Finalize
  use ESMF,          only : ESMF_LogKind_Flag, ESMF_LOGKIND_NONE
  use ESMF,          only : ESMF_LOGKIND_SINGLE, ESMF_LOGKIND_MULTI
  use ESMF,          only : ESMF_LOGKIND_MULTI_ON_ERROR
  use cime_comp_mod, only : cime_pre_init1
  use cime_comp_mod, only : cime_pre_init2
  use cime_comp_mod, only : cime_init
!!>> TSMP PDAF comment out beginning
  ! use cime_comp_mod, only : cime_run
  ! use cime_comp_mod, only : cime_final
!!<< TSMP PDAF comment out end
!!>> TSMP PDAF addition beginning
  use iso_C_binding
  use enkf_clm_mod
!!<< TSMP PDAF addition end

  implicit none

!!>> TSMP PDAF addition beginning
  !--------------------------------------------------------------------------
  ! PDAF variables
  !--------------------------------------------------------------------------
  character(kind=c_char,len=1),dimension(100),intent(in) :: finname 
  integer(c_int), intent(in) :: pdaf_id
  integer(c_int), intent(in) :: pdaf_max
  integer(c_int), intent(in) :: mype
  integer(c_int) :: counter
!!<< TSMP PDAF addition end
  !--------------------------------------------------------------------------
  ! timing variables
  !--------------------------------------------------------------------------
  integer(i8) :: beg_count, end_count, irtc_rate
  real(r8)    :: cime_pre_init1_time, ESMF_Initialize_time, &
       cime_pre_init2_time, cime_init_time_adjustment

  !--------------------------------------------------------------------------
  ! For ESMF logging
  !--------------------------------------------------------------------------
  character(len=CS)       :: esmf_logfile_option
  type(ESMF_LogKind_Flag) :: esmf_logfile_kind

  !--------------------------------------------------------------------------
  ! Setup and initialize the communications and logging.
  !--------------------------------------------------------------------------
  beg_count = shr_sys_irtc(irtc_rate)

!!>> TSMP PDAF comment out beginning
  ! call cime_pre_init1(esmf_logfile_option)
!!>> TSMP PDAF addition beginning
  call cime_pre_init1(esmf_logfile_option, &
                      COMM_model_clm, &
                      pdaf_id=pdaf_id, &
                      pdaf_max=pdaf_max)
!!<< TSMP PDAF addition end

  end_count = shr_sys_irtc(irtc_rate)
  cime_pre_init1_time = real( (end_count-beg_count), r8)/real(irtc_rate, r8)

  !--------------------------------------------------------------------------
  ! Initialize ESMF.  This is done outside of the ESMF_INTERFACE ifdef
  ! because it is needed for the time manager, even if the ESMF_INTERFACE
  ! is not used.
  !--------------------------------------------------------------------------
  beg_count = shr_sys_irtc(irtc_rate)


  select case(esmf_logfile_option)
  case('ESMF_LOGKIND_SINGLE')
     esmf_logfile_kind = ESMF_LOGKIND_SINGLE
  case('ESMF_LOGKIND_MULTI')
     esmf_logfile_kind = ESMF_LOGKIND_MULTI
  case('ESMF_LOGKIND_MULTI_ON_ERROR')
     esmf_logfile_kind = ESMF_LOGKIND_MULTI_ON_ERROR
  case('ESMF_LOGKIND_NONE')
     esmf_logfile_kind = ESMF_LOGKIND_NONE
  case default
     call shr_sys_abort('CIME ERROR: invalid ESMF logfile kind '//trim(esmf_logfile_option))
  end select
!!>> TSMP PDAF addition beginning
  write(6,*) "esmf_initialize"
!!<< TSMP PDAF addition end
  call ESMF_Initialize(logkindflag=esmf_logfile_kind)

  end_count = shr_sys_irtc(irtc_rate)
  ESMF_Initialize_time = real( (end_count-beg_count), r8)/real(irtc_rate, r8)

  !--------------------------------------------------------------------------
  ! Read in the configuration information and initialize the time manager.
  !--------------------------------------------------------------------------
  ! Timer initialization has to be after determination of the maximum number
  ! of threads used across all components, so called inside of
  ! cime_pre_init2, as are t_startf and t_stopf for CPL:INIT and
  ! cime_pre_init2.
  !--------------------------------------------------------------------------
  beg_count = shr_sys_irtc(irtc_rate)

!!>> TSMP PDAF addition beginning
  write(6,*) "cime-pre-init2"
!!<< TSMP PDAF addition end
  call cime_pre_init2()

  end_count = shr_sys_irtc(irtc_rate)
  cime_pre_init2_time = real( (end_count-beg_count), r8)/real(irtc_rate, r8)

  !--------------------------------------------------------------------------
  ! Call the initialize, run and finalize routines.
  !--------------------------------------------------------------------------

  call t_startf('CPL:INIT')
  call t_adj_detailf(+1)

  call t_startstop_valsf('CPL:cime_pre_init1',  walltime=cime_pre_init1_time)
  call t_startstop_valsf('CPL:ESMF_Initialize', walltime=ESMF_Initialize_time)
  call t_startstop_valsf('CPL:cime_pre_init2',  walltime=cime_pre_init2_time)

  call cime_init()

  call t_adj_detailf(-1)
  call t_stopf('CPL:INIT')

  cime_init_time_adjustment = cime_pre_init1_time  &
       + ESMF_Initialize_time &
       + cime_pre_init2_time
  call t_startstop_valsf('CPL:INIT',  walltime=cime_init_time_adjustment, &
       callcount=0)

#if defined CLMSA
  call define_clm_statevec(mype)
#endif 


end subroutine clm_init


!--------------------------------------------------------------------------
! clm_advance is the interface for cime_run() from
! cime/src/drivers/mct/main/cime_comp_mod.F90
! However, we modify cime_run() to use the TSMP time step to make sure that
! clm_advance stops the driver loop between DA steps.
! After cime_run we call set_clm_statevec() to transfer from CLM to PDAF.
!--------------------------------------------------------------------------
subroutine clm_advance(ntstep, tstartcycle, mype) bind(C,name="clm_advance")
  use cime_comp_mod, only : cime_run
  use enkf_clm_mod, only : set_clm_statevec 
  use iso_C_binding

  implicit none
  !--------------------------------------------------------------------------
  ! PDAF variables
  !--------------------------------------------------------------------------
  integer(c_int),intent(in) :: ntstep
  integer(c_int),intent(in) :: tstartcycle
  integer(c_int),intent(in) :: mype

  ! call modified cime_run that runs for a specificied number of timesteps.
  call cime_run(ntstep)

#if defined CLMSA
  ! Calling PDAF Function to set state vector before assimiliation
  call set_clm_statevec(tstartcycle, mype)
#endif

end subroutine clm_advance
!--------------------------------------------------------------------------
! clm_finalize() calls cime_final() and ESMF_Finalize()
! Note that ESMF_Finalize() contains mpi_finalize()
! Therefor, it can cause conflicts if mpi_finalize() is called elsewhere.
!--------------------------------------------------------------------------
subroutine clm_finalize() bind(C,name="clm_finalize")
  use iso_C_binding

  ! use ESMF,          only : ESMF_Initialize, ESMF_Finalize
  use cime_comp_mod, only : cime_final

  implicit none

  call cime_final()

  !--------------------------------------------------------------------------
  ! Clean-up
  !--------------------------------------------------------------------------
  ! call ESMF_Finalize( )

end subroutine clm_finalize

end module enkf_clm_5
