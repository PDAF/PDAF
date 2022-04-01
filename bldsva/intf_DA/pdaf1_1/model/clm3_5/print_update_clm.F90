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
!print_update_clm.F90: Module for printing updated CLM ensemble
!-------------------------------------------------------------------------------------------
! subroutine print_update_clm(ts,ttot) bind(C,name="print_update_clm")

!     use shr_kind_mod , only : r8 => shr_kind_r8
!     use clm_atmlnd   , only : clm_l2a, atm_l2a, clm_mapl2a
!     use clmtype      , only : clm3, nameg, namec
!     use subgridavemod, only : p2g, c2g
!     use domainmod    , only : latlon_type
!     use clm_varpar   , only : nlevsoi
!     use decompmod    , only : get_proc_global, get_proc_bounds, adecomp
!     use spmdgathscatmod , only : gather_data_to_master
!     use spmdmod      , only : masterproc
!     use clm_time_manager        , only : get_nstep, dtime, nelapse
!     use netcdf
!     use enkf_clm_mod, only : clmupdate_swc,clmupdate_texture,clmprint_swc

!     implicit none

!     integer, intent(in) :: ts,ttot

!     ! *** local variables ***
!     integer :: numg           ! total number of gridcells across all processors
!     integer :: numl           ! total number of landunits across all processors
!     integer :: numc           ! total number of columns across all processors
!     integer :: nump           ! total number of pfts across all processors
!     integer :: begg,endg      ! local beg/end gridcells gdc
!     integer :: begl,endl      ! local beg/end landunits
!     integer :: begc,endc      ! local beg/end columns
!     integer :: begp,endp      ! local beg/end pfts

!     integer ::   isec, info, jn, jj, ji, g1, jx    ! temporary integer
!     real(r8), pointer :: swc(:,:)
!     real(r8), pointer :: psand(:,:)
!     real(r8), pointer :: pclay(:,:)
!     real(r8), pointer :: clmstate_tmp_local(:,:)
!     real(r8), pointer :: clmstate_tmp_global(:,:)
!     real(r8), allocatable :: clmstate_out(:,:,:)
!     integer ,dimension(4) :: dimids
!     integer ,dimension(1) :: il_var_id
!     integer :: il_file_id, ncvarid(3), status
!     character(len = 300) :: update_filename
!     integer :: nerror
!     integer :: ndlon,ndlat


!     call get_proc_global(numg,numl,numc,nump)
!     call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)
!     allocate(clmstate_tmp_local(nlevsoi,begc:endc), stat=nerror)

!     ndlon = adecomp%gdc2i(numg)
!     ndlat = adecomp%gdc2j(numg)

!     if (masterproc) then
!       allocate(clmstate_tmp_global(nlevsoi,numg), stat=nerror)
!       allocate(clmstate_out(ndlon,ndlat,nlevsoi), stat=nerror)
!     end if

!     if(masterproc) then
!       call get_update_filename(update_filename)
!       if(ts.eq.1) then
!         status =  nf90_create(update_filename, NF90_CLOBBER, il_file_id)
!         status =  nf90_def_dim(il_file_id, "x", ndlon, dimids(1))
!         status =  nf90_def_dim(il_file_id, "y", ndlat, dimids(2))
!         status =  nf90_def_dim(il_file_id, "z", nlevsoi, dimids(3))
!         status =  nf90_def_dim(il_file_id, "t", ttot, dimids(4))
!         if(clmprint_swc.eq.1)     status =  nf90_def_var(il_file_id, "swc", NF90_DOUBLE, dimids, ncvarid(1))
!         if(clmupdate_texture.eq.1) status =  nf90_def_var(il_file_id, "sand", NF90_DOUBLE, dimids, ncvarid(2))
!         if(clmupdate_texture.eq.1) status =  nf90_def_var(il_file_id, "clay", NF90_DOUBLE, dimids, ncvarid(3))
!         status =  nf90_enddef(il_file_id)
!       else
!         status = nf90_open(update_filename,NF90_WRITE,il_file_id)
!       endif
!     endif


!     if(clmprint_swc.eq.1) then
!       swc  => clm3%g%l%c%cws%h2osoi_vol
!       ! swc
!       clmstate_tmp_local = transpose(swc)
!       call gather_data_to_master(clmstate_tmp_local,clmstate_tmp_global, clmlevel=nameg)

!       if(masterproc) then
!         ji = adecomp%gdc2i(numg)
!         jj = adecomp%gdc2j(numg)
!         do jn = 1, nlevsoi
!           do g1 = 1, numg
!              ji = adecomp%gdc2i(g1)
!              jj = adecomp%gdc2j(g1)
!              clmstate_out(ji,jj,jn) = clmstate_tmp_global(jn,g1)
!           end do
!         end do
!         status = nf90_inq_varid(il_file_id, "swc" , ncvarid(1))
!         status = nf90_put_var( il_file_id, ncvarid(1), clmstate_out(:,:,:), &
!                  start = (/ 1, 1, 1, ts/), count = (/ ndlon, ndlat, nlevsoi, 1 /) )
!         !status = nf90_close(il_file_id)
!       end if
!     end if

!     if(clmupdate_texture.eq.1) then
!       psand => clm3%g%l%c%cps%psand
!       pclay => clm3%g%l%c%cps%pclay
!       ! sand
!       clmstate_tmp_local = transpose(psand)
!       call gather_data_to_master(clmstate_tmp_local,clmstate_tmp_global, clmlevel=nameg)

!       if(masterproc) then
!         ji = adecomp%gdc2i(numg)
!         jj = adecomp%gdc2j(numg)
!         do jn = 1, nlevsoi
!           do g1 = 1, numg
!              ji = adecomp%gdc2i(g1)
!              jj = adecomp%gdc2j(g1)
!              clmstate_out(ji,jj,jn) = clmstate_tmp_global(jn,g1)
!           end do
!         end do
!         status = nf90_inq_varid(il_file_id, "sand" , ncvarid(2))
!         status = nf90_put_var( il_file_id, ncvarid(2), clmstate_out(:,:,:), &
!                  start = (/ 1, 1, 1, ts/), count = (/ ndlon, ndlat, nlevsoi, 1 /) )
!         !status = nf90_close(il_file_id)
!       end if

!       ! clay
!       clmstate_tmp_local = transpose(pclay)
!       call gather_data_to_master(clmstate_tmp_local,clmstate_tmp_global, clmlevel=nameg)

!       if(masterproc) then
!         ji = adecomp%gdc2i(numg)
!         jj = adecomp%gdc2j(numg)
!         do jn = 1, nlevsoi
!           do g1 = 1, numg
!              ji = adecomp%gdc2i(g1)
!              jj = adecomp%gdc2j(g1)
!              clmstate_out(ji,jj,jn) = clmstate_tmp_global(jn,g1)
!           end do
!         end do
!         status = nf90_inq_varid(il_file_id, "clay" , ncvarid(3))
!         status = nf90_put_var( il_file_id, ncvarid(3), clmstate_out(:,:,:), &
!                  start = (/ 1, 1, 1, ts/), count = (/ ndlon, ndlat, nlevsoi, 1 /) )
!         !status = nf90_close(il_file_id)
!       end if
!     end if

!     if(masterproc) then
!       status = nf90_close(il_file_id)
!       deallocate(clmstate_out)
!       deallocate(clmstate_tmp_global)
!     end if

!     deallocate(clmstate_tmp_local)

! end subroutine print_update_clm

! subroutine get_update_filename (iofile)
!     ! !USES:
!     use clm_varctl, only : caseid
!     use clm_time_manager, only : get_curr_date, get_prev_date
!     ! !ARGUMENTS:
!     implicit none
!     character(len=300),intent(inout) :: iofile
!     ! LOCAL VARIABLES:
!     character(len=256) :: cdate       !date char string
!     integer :: day                    !day (1 -> 31)
!     integer :: mon                    !month (1 -> 12)
!     integer :: yr                     !year (0 -> ...)
!     integer :: sec                    !seconds into current day
!     !-----------------------------------------------------------------------

!     call get_prev_date (yr, mon, day, sec)
!     write(cdate,'(i4.4,"-",i2.2)') yr,mon                         !other
!     call get_curr_date (yr, mon, day, sec)
!     write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
!     !iofile = trim(caseid)//".update."//trim(cdate)//".nc"
!     iofile = trim(caseid)//".update.nc"
! end subroutine get_update_filename
