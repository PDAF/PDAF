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
!print_update_clm_5.F90: Module for printing updated CLM5 ensemble
!-------------------------------------------------------------------------------------------

#if defined CLMSA
subroutine print_update_clm(ts,ttot) bind(C,name="print_update_clm")
    
    use iso_c_binding   
    use shr_kind_mod , only : r8 => shr_kind_r8
    use subgridavemod, only : p2g, c2g
    use domainMod    , only : ldomain 
    use clm_varpar   , only : nlevsoi
    use clm_varcon   , only : nameg 
    use decompmod    , only : get_proc_global, get_proc_bounds, ldecomp
    use spmdgathscatmod , only : gather_data_to_master
    use spmdmod      , only : masterproc
    use clm_time_manager        , only : get_nstep    
    use clm_instMod, only : soilstate_inst, waterstate_inst
    use netcdf
    use enkf_clm_mod, only : clmupdate_swc,clmupdate_texture,clmprint_swc

    implicit none

    integer(c_int), intent(in) :: ts,ttot

    ! *** local variables ***
    integer :: numg           ! total number of gridcells across all processors
    integer :: numl           ! total number of landunits across all processors
    integer :: numc           ! total number of columns across all processors
    integer :: nump           ! total number of pfts across all processors
    integer :: begg,endg      ! local beg/end gridcells gdc
    integer :: begl,endl      ! local beg/end landunits
    integer :: begc,endc      ! local beg/end columns
    integer :: begp,endp      ! local beg/end pfts
    
    integer ::   isec, info, jn, jj, ji, g1, jx    ! temporary integer
    real(r8), pointer :: swc(:,:)
    real(r8), pointer :: psand(:,:)
    real(r8), pointer :: pclay(:,:)
    real(r8), pointer :: porgm(:,:)
    real(r8), pointer :: clmstate_tmp_local(:)
    real(r8), pointer :: clmstate_tmp_global(:)
    real(r8), allocatable :: clmstate_out(:,:,:)
    integer ,dimension(4) :: dimids
    integer ,dimension(1) :: il_var_id
    integer :: il_file_id, ncvarid(4), status
    character(len = 300) :: update_filename
    integer :: nerror
    integer :: ndlon,ndlat


    call get_proc_global(ng=numg,nl=numl,nc=numc,np=nump)
    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)
    allocate(clmstate_tmp_local(nlevsoi*(-begc+endc)), stat=nerror)

    ndlon  = ldomain%ni
    ndlat  = ldomain%nj

    if (masterproc) then
      allocate(clmstate_tmp_global(nlevsoi*numg), stat=nerror)
      allocate(clmstate_out(ndlon,ndlat,nlevsoi), stat=nerror)
    end if

    if(masterproc) then
      call get_update_filename(update_filename)
      if(ts.eq.1) then
        status =  nf90_create(update_filename, NF90_CLOBBER, il_file_id)
        status =  nf90_def_dim(il_file_id, "x", ndlon, dimids(1))
        status =  nf90_def_dim(il_file_id, "y", ndlat, dimids(2))
        status =  nf90_def_dim(il_file_id, "z", nlevsoi, dimids(3))
        status =  nf90_def_dim(il_file_id, "t", ttot, dimids(4))

        if(clmprint_swc.eq.1) then
          status =  nf90_def_var(il_file_id, "swc", NF90_DOUBLE, dimids, ncvarid(1))
        endif
 
        if(clmupdate_texture.eq.1) then 
          status =  nf90_def_var(il_file_id, "sand", NF90_DOUBLE, dimids, ncvarid(2))
          status =  nf90_def_var(il_file_id, "clay", NF90_DOUBLE, dimids, ncvarid(3))
        endif

        ! write updates to sand, clay and organic matter
        if(clmupdate_texture.eq.2) then
          status =  nf90_def_var(il_file_id, "sand", NF90_DOUBLE, dimids, ncvarid(2))
          status =  nf90_def_var(il_file_id, "clay", NF90_DOUBLE, dimids, ncvarid(3))
          status =  nf90_def_var(il_file_id, "orgm", NF90_DOUBLE, dimids, ncvarid(4))
        endif
        status =  nf90_enddef(il_file_id)
      else
        status = nf90_open(update_filename,NF90_WRITE,il_file_id)
      endif
    endif
  
    
    if(clmprint_swc.eq.1) then
      swc  => waterstate_inst%h2osoi_vol_col
      ! swc
!      clmstate_tmp_local = pack(swc,.true.)
!      clmstate_tmp_global = clmstate_tmp_local
!      call gather_data_to_master(clmstate_tmp_local,clmstate_tmp_global,clmlevel=nameg)

      if(masterproc) then
        do jn = 1, nlevsoi
          do g1 = 1, numg
             ji = mod(ldecomp%gdc2glo(g1)-1,ldomain%ni) + 1
             jj = (ldecomp%gdc2glo(g1) - 1)/ldomain%ni + 1
             clmstate_out(ji,jj,jn) = swc(g1, jn) !clmstate_tmp_global(jn+(g1-1)*nlevsoi)
          end do
        end do
        status = nf90_inq_varid(il_file_id, "swc" , ncvarid(1))
        status = nf90_put_var( il_file_id, ncvarid(1), clmstate_out(:,:,:), &
                 start = (/ 1, 1, 1, ts/), count = (/ ndlon, ndlat, nlevsoi, 1/) )
        !status = nf90_close(il_file_id)
      end if
    end if
    
    if((clmupdate_texture.eq.1) .or. (clmupdate_texture.eq.2)) then
      psand => soilstate_inst%cellsand_col
      pclay => soilstate_inst%cellclay_col

      ! sand
      !clmstate_tmp_local = pack(psand,.true.)
      !clmstate_tmp_global = clmstate_tmp_local
!      call gather_data_to_master(clmstate_tmp_local,clmstate_tmp_global, clmlevel=nameg)
      if(masterproc) then
        do jn = 1, nlevsoi
          do g1 = 1, numg
             ji = mod(ldecomp%gdc2glo(g1)-1,ldomain%ni) + 1
             jj = (ldecomp%gdc2glo(g1) - 1)/ldomain%ni + 1
             clmstate_out(ji,jj,jn) = psand(g1,jn)!clmstate_tmp_global(jn+(g1-1)*nlevsoi)
          end do
        end do
        status = nf90_inq_varid(il_file_id, "sand" , ncvarid(2))
        status = nf90_put_var( il_file_id, ncvarid(2), clmstate_out(:,:,:), &
                 start = (/ 1, 1, 1, ts/), count = (/ ndlon, ndlat, nlevsoi, 1/) )
        !status = nf90_close(il_file_id)
      end if

      ! clay
     ! clmstate_tmp_local = pack(pclay,.true.)
     ! clmstate_tmp_global = clmstate_tmp_local
!      call gather_data_to_master(clmstate_tmp_local,clmstate_tmp_global,clmlevel=nameg)

      if(masterproc) then
        do jn = 1, nlevsoi
          do g1 = 1, numg
             ji = mod(ldecomp%gdc2glo(g1)-1,ldomain%ni) + 1
             jj = (ldecomp%gdc2glo(g1) - 1)/ldomain%ni + 1
             clmstate_out(ji,jj,jn) = pclay(g1, jn) !clmstate_tmp_global(jn+(g1-1)*nlevsoi)
          end do
        end do
        status = nf90_inq_varid(il_file_id, "clay" , ncvarid(3))
        status = nf90_put_var( il_file_id, ncvarid(3), clmstate_out(:,:,:), &
                 start = (/ 1, 1, 1, ts/), count = (/ ndlon, ndlat, nlevsoi, 1/) )
        !status = nf90_close(il_file_id)
      end if
      
      ! organic matter
      if(clmupdate_texture.eq.2) then
        porgm => soilstate_inst%cellorg_col

        !clmstate_tmp_local = pack(porgm,.true.)
        !clmstate_tmp_global = clmstate_tmp_local
!        call gather_data_to_master(clmstate_tmp_local,clmstate_tmp_global, clmlevel=nameg)

        if(masterproc) then
          do jn = 1, nlevsoi
            do g1 = 1, numg
             ji = mod(ldecomp%gdc2glo(g1)-1,ldomain%ni) + 1
             jj = (ldecomp%gdc2glo(g1) - 1)/ldomain%ni + 1
             clmstate_out(ji,jj,jn) = porgm(g1, jn) !clmstate_tmp_global(jn+(g1-1)*nlevsoi)
            end do
          end do
          status = nf90_inq_varid(il_file_id, "orgm" , ncvarid(4))
          status = nf90_put_var( il_file_id, ncvarid(4), clmstate_out(:,:,:), &
                   start = (/ 1, 1, 1, ts/), count = (/ ndlon, ndlat, nlevsoi, 1/) )
          !status = nf90_close(il_file_id)
        end if
      endif

    end if 
    
    if(masterproc) then
      status = nf90_close(il_file_id)
      deallocate(clmstate_out)
      deallocate(clmstate_tmp_global)
    end if

    deallocate(clmstate_tmp_local)

end subroutine print_update_clm
#endif

subroutine get_update_filename (iofile)
    use clm_varctl, only : caseid
    use clm_time_manager, only : get_curr_date, get_prev_date
    ! !ARGUMENTS:
    implicit none
    character(len=300),intent(inout) :: iofile
    ! LOCAL VARIABLES:
    character(len=256) :: cdate       !date char string
    integer :: day                    !day (1 -> 31)
    integer :: mon                    !month (1 -> 12)
    integer :: yr                     !year (0 -> ...)
    integer :: sec                    !seconds into current day
    !-----------------------------------------------------------------------

    call get_prev_date (yr, mon, day, sec)
    write(cdate,'(i4.4,"-",i2.2)') yr,mon 
    call get_curr_date (yr, mon, day, sec)
    !write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
    write(cdate,'(i4.4)') yr
    iofile = trim(caseid)//".update."//trim(cdate)//".nc"
    !iofile = trim(caseid)//".update.nc"
end subroutine get_update_filename
