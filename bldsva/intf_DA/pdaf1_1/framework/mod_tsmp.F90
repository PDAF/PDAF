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
!mod_tsmp.F90: Module for binding to wrapper_tsmp.c
!-------------------------------------------------------------------------------------------

module mod_tsmp
    use iso_c_binding

    integer(c_int) , bind(c) :: enkf_subvecsize, pf_statevecsize, nprocpf, nprocclm, nproccosmo
    integer(c_int) , bind(c,name="point_obs") :: point_obs
    integer(c_int) , bind(c) :: nx_local, ny_local, nz_local, nx_glob, ny_glob, nz_glob  
    integer(c_int), bind(c)  :: tag_model_clm = 0
    integer(c_int), bind(c)  :: tag_model_parflow = 1
    integer(c_int), bind(c)  :: tag_model_cosmo   = 2
    type(c_ptr), bind(c)     :: pf_statevec
    type(c_ptr), bind(c)     :: xcoord, ycoord, zcoord
    real(c_double), pointer  :: xcoord_fortran(:), ycoord_fortran(:), zcoord_fortran(:)
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
        subroutine integrate_tsmp(tcycle) bind(c)
            use iso_c_binding
            import
            INTEGER(c_int) :: tcycle ! No. of assimilation cycle
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
   
     interface
        subroutine init_n_domains_size(n_domains_p) bind(c)
            use iso_c_binding
            import
            INTEGER(c_int) :: n_domains_p ! PE-local number of analysis domains
        end subroutine init_n_domains_size
    end interface

     interface
        subroutine init_parf_l_size(dim_l) bind(c)
            use iso_c_binding
            import
              INTEGER(c_int) :: dim_l ! Local state dimension
        end subroutine init_parf_l_size
    end interface

!!$    interface
!!$        subroutine g2l_state(domain_p, state_p, dim_l, state_l) bind(c)
!!$            use iso_c_binding
!!$            import
!!$               INTEGER(c_int) :: domain_p      ! Current local analysis domain
!!$               INTEGER(c_int) :: dim_l         ! Local state dimension
!!$               TYPE(c_ptr) :: state_p          ! PE-local full state vector 
!!$               TYPE(c_ptr) :: state_l          ! State vector on local analysis domain
!!$        end subroutine g2l_state
!!$    end interface
!!$
!!$    interface
!!$        subroutine l2g_state(domain_p, state_p, dim_l, state_l) bind(c)
!!$            use iso_c_binding
!!$            import
!!$            INTEGER(c_int) :: domain_p       ! Current local analysis domain
!!$            INTEGER(c_int) :: dim_l          ! Local state dimension
!!$            TYPE(c_ptr) :: state_l            ! State vector on local analysis domain
!!$            TYPE(c_ptr) :: state_p            ! PE-local full state vector 
!!$        end subroutine l2g_state
!!$    end interface

end module mod_tsmp

