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
!mod_clm_statistics.F90: Module for calculating CLM ensemble statistics
!-------------------------------------------------------------------------------------------

module mod_clm_statistics
  use iso_c_binding

contains
  subroutine write_clm_statistics(ts,ttot) bind(C,name="write_clm_statistics")
    integer, intent(in) :: ts,ttot

    write(*,*) "dummy sub"
  end subroutine write_clm_statistics


  character(len=256) function get_statistic_filename ()
    write(*,*) "dummy func"
    get_statistic_filename = "test"
  end function get_statistic_filename
 
end module mod_clm_statistics
