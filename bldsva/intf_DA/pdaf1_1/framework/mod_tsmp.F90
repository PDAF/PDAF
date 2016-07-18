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
!mod_tsmp.F90: Module for binding to wrapper_tsmp.c
!-------------------------------------------------------------------------------------------

module mod_tsmp
    use iso_c_binding

    integer(c_int) , bind(c) :: enkf_subvecsize, pf_statevecsize, nprocpf, nprocclm, nproccosmo
    integer(c_int), bind(c)  :: tag_model_clm = 0
    integer(c_int), bind(c)  :: tag_model_parflow = 1
    integer(c_int), bind(c)  :: tag_model_cosmo   = 2
    type(c_ptr), bind(c)     :: pf_statevec
    real(c_double), pointer  :: pf_statevec_fortran(:)
    type(c_ptr), bind(c)     :: idx_map_subvec2state
    integer(c_int), pointer  :: idx_map_subvec2state_fortran(:)

    interface
        subroutine initialize_tsmp() bind(c)
            use iso_c_binding
            implicit none
        end subroutine initialize_tsmp
    end interface

    interface
        subroutine finalize_tsmp() bind(c)
            use iso_c_binding
            implicit none
        end subroutine finalize_tsmp
    end interface

    interface
        subroutine integrate_tsmp() bind(c)
            use iso_c_binding
            implicit none
        end subroutine integrate_tsmp
    end interface

    interface
        subroutine print_update_pfb() bind(c)
            use iso_c_binding
            implicit none
        end subroutine print_update_pfb
    end interface
    
    interface
        subroutine update_tsmp() bind(c)
            use iso_c_binding
            implicit none
        end subroutine update_tsmp
    end interface
    
end module mod_tsmp

