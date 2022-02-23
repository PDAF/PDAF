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
!mod_parallel_model.F90: Module for control of parallel models of TerrSysMP
!-------------------------------------------------------------------------------------------

module mod_parallel_model
    use iso_c_binding
    include 'mpif.h'
save
    ! mpi related
    integer :: npes_parflow
    integer :: coupcol
#if defined PARFLOW_STAND_ALONE
! Parflow stand alone directly uses binded communicator
    integer(c_int),bind(c,name='comm_model_pdaf') :: comm_model
#else
! CLM stand alone uses comm_model directly, while TerrSysMP use this and da_comm(_clm)
    integer :: comm_model
#endif
    integer :: mype_model
    integer :: npes_model
    integer :: mype_world
    integer :: npes_world
    bind(c) :: mype_model
    bind(c) :: mype_world
    ! model input parameters
    real(c_double), bind(c) :: t_start
    !integer(c_int), bind(c) :: da_interval, model
    integer(c_int), bind(c) ::  model
    integer(c_int), bind(c, name = 'nsteps') :: total_steps
    integer :: tcycle
    interface
        subroutine read_enkfpar(parname) BIND(C, name='read_enkfpar')
            use iso_c_binding
            implicit none
            character, dimension(*), intent(in) :: parname
        end subroutine read_enkfpar
    end interface
contains
    subroutine abort_parallel
    end subroutine abort_parallel
end module mod_parallel_model
