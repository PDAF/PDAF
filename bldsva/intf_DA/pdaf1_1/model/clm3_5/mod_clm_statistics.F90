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
  !use spmdMod, only: comm_local_clm => mpicom, &
  !     npes_local_clm => npes, &
  !     mype_local_clm => iam, &
  !     is_masterproc_local_clm => masterproc
  use spmdmod, only: masterproc
  use shr_kind_mod, only: r8 => shr_kind_r8
  use mpi
  use enkf_clm_mod, only:statcomm
  use spmdgathscatmod , only : gather_data_to_master
  use subgridAveMod, only : p2g, c2g


  implicit none
  integer, parameter :: max_chars = 300

contains
  subroutine write_clm_statistics(ts,ttot) bind(C,name="write_clm_statistics")
    USE clmtype,                  ONLY : clm3,nameg,namec
    USE decompMod    , only : get_proc_global, get_proc_bounds, adecomp
    use netcdf

    integer, intent(in) :: ts,ttot
    integer :: numg           ! total number of gridcells across all processors
    integer :: numl           ! total number of landunits across all processors
    integer :: numc           ! total number of columns across all processors
    integer :: nump           ! total number of pfts across all processors
    integer :: begg,endg      ! local beg/end gridcells gdc
    integer :: begl,endl      ! local beg/end landunits
    integer :: begc,endc      ! local beg/end columns
    integer :: begp,endp      ! local beg/end pfts

    integer :: nlon,nlat

    type field_pointer
       real(r8), pointer :: p(:)
    end type field_pointer

    integer, parameter :: n_variables = 4
    character(len = max_chars) :: variable_names(n_variables)
    type(field_pointer), dimension(n_variables) :: variable_pointers
    integer :: ncvarid(n_variables)
    real(r8), pointer :: clmvar_global_g(:)
    real(r8), allocatable ::   clmvar_out(:,:)
    real(r8), pointer :: clmvar_local_g(:)

    character(len = max_chars) :: statistic_filename
    integer :: ierr,i,il_file_id,g1,nloc
    !real(r8), allocatable,target :: mm(:),sd(:)
    real(r8), pointer :: mm(:),var(:),sd(:)
    real(r8),pointer :: ptr(:)
    integer,dimension(3) :: dimids
    integer ji,jj
    integer realrank,realsize

    real(r8), pointer :: lon(:)
    real(r8), pointer :: lat(:)

    ! Get total global number of grid cells, landunits, columns and pfts
    CALL get_proc_global(numg,numl,numc,nump)
    CALL get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)


    if(masterproc) then
      lon   => clm3%g%londeg
      lat   => clm3%g%latdeg
      nlon = adecomp%gdc2i(numg)
      nlat = adecomp%gdc2j(numg)
    end if

    call mpi_comm_rank(statcomm,realrank,ierr)
    call mpi_comm_size(statcomm,realsize,ierr)

    !variable_names(1) = "FCTR"
    !variable_names(2) = "FCEV"
    !variable_names(3) = "FGEV"
    !variable_names(4) = "FSH"

    variable_names(1) = "lh_mean"
    variable_names(2) = "lh_sd"
    variable_names(3) = "sh_mean"
    variable_names(4) = "sh_sd"

    ! define netcdf output file
    if(masterproc .and. (realrank.eq.0)) then
      statistic_filename = get_statistic_filename()
      if(ts.eq.1) then
        ierr =  nf90_create(statistic_filename, NF90_CLOBBER, il_file_id)
        ierr =  nf90_def_dim(il_file_id, "lon", nlon, dimids(1))
        ierr =  nf90_def_dim(il_file_id, "lat", nlat, dimids(2))
        ierr =  nf90_def_dim(il_file_id, "t", ttot, dimids(3))
        ierr =  nf90_def_var(il_file_id, trim(variable_names(1)), NF90_DOUBLE, dimids,ncvarid(1))
        ierr =  nf90_def_var(il_file_id, trim(variable_names(2)), NF90_DOUBLE, dimids,ncvarid(2))
        ierr =  nf90_def_var(il_file_id, trim(variable_names(3)), NF90_DOUBLE, dimids,ncvarid(3))
        ierr =  nf90_def_var(il_file_id, trim(variable_names(4)), NF90_DOUBLE, dimids,ncvarid(4))
        ierr =  nf90_enddef(il_file_id)
      else
        ierr = nf90_open(statistic_filename,NF90_WRITE,il_file_id)
      endif
    end if


    ! set pointers to CLM variables
    variable_pointers(1)%p => clm3%g%l%c%p%pef%eflx_lh_vegt
    variable_pointers(2)%p => clm3%g%l%c%p%pef%eflx_lh_vege
    variable_pointers(3)%p => clm3%g%l%c%p%pef%eflx_lh_grnd
    variable_pointers(4)%p => clm3%g%l%c%p%pef%eflx_sh_tot


    ! allocate arrays
    if (masterproc) then
      allocate (clmvar_global_g(numg), stat=ierr)
      allocate (clmvar_out(nlon,nlat), stat=ierr)
    end if

    nloc = endg-begg+1
    allocate(mm(begg:endg),stat=ierr)
    allocate(var(begg:endg),stat=ierr)
    allocate(sd(begg:endg),stat=ierr)
    allocate(clmvar_local_g(begg:endg), stat=ierr)


    ! calculate statistics for sensible heat flux
    call p2g(begp, endp, begc, endc, begl, endl, begg, endg, variable_pointers(4)%p, clmvar_local_g, &
              p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    mm  = 0.0
    var = 0.0
    sd  = 0.0
    call mpi_allreduce(clmvar_local_g,mm,nloc,MPI_REAL8,MPI_SUM,statcomm,ierr)
    mm  = mm / realsize
    var = (clmvar_local_g - mm) * (clmvar_local_g - mm)
    call mpi_reduce(var,sd,nloc,MPI_REAL8,MPI_SUM,0,statcomm,ierr)
    sd = sqrt(sd/(realsize-1))


    if((realrank.eq.0)) then
      ptr => mm
      call gather_data_to_master(ptr,clmvar_global_g,clmlevel=nameg)

      if(masterproc) then
        do g1 = 1, numg
          ji = adecomp%gdc2i(g1)
          jj = adecomp%gdc2j(g1)
          clmvar_out(ji,jj) = clmvar_global_g(g1)
        end do
        ierr = nf90_inq_varid(il_file_id, trim(variable_names(3)) , ncvarid(3))
        ierr = nf90_put_var( il_file_id, ncvarid(3), clmvar_out(:,:), start = (/ 1, 1,ts /), count = (/ nlon, nlat, 1 /) )
      end if
    end if

    if((realrank.eq.0)) then
      ptr => sd
      call gather_data_to_master(ptr,clmvar_global_g,clmlevel=nameg)

      if(masterproc) then
        do g1 = 1, numg
          ji = adecomp%gdc2i(g1)
          jj = adecomp%gdc2j(g1)
          clmvar_out(ji,jj) = clmvar_global_g(g1)
        end do
        ierr = nf90_inq_varid(il_file_id, trim(variable_names(4)) , ncvarid(4))
        ierr = nf90_put_var( il_file_id, ncvarid(4), clmvar_out(:,:), start = (/ 1, 1,ts /), count = (/ nlon, nlat, 1 /) )
      end if
    end if


    ! calculate statistics for latent heat flux
    mm  = 0.0
    var = 0.0
    sd  = 0.0

    do i=1,3
      call p2g(begp, endp, begc, endc, begl, endl, begg, endg, variable_pointers(i)%p, clmvar_local_g, &
                p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
      var = var + clmvar_local_g
    end do

    call mpi_allreduce(var,mm,nloc,MPI_REAL8,MPI_SUM,statcomm,ierr)
    mm = mm / realsize
    var = (var - mm) * (var - mm)
    call mpi_reduce(var,sd,nloc,MPI_REAL8,MPI_SUM,0,statcomm,ierr)
    sd = sqrt(sd/(realsize-1))

    if((realrank.eq.0)) then
      ptr => mm
      call gather_data_to_master(ptr,clmvar_global_g,clmlevel=nameg)

      if(masterproc) then
        do g1 = 1, numg
          ji = adecomp%gdc2i(g1)
          jj = adecomp%gdc2j(g1)
          clmvar_out(ji,jj) = clmvar_global_g(g1)
        end do
        ierr = nf90_inq_varid(il_file_id, trim(variable_names(1)) , ncvarid(1))
        ierr = nf90_put_var( il_file_id, ncvarid(1), clmvar_out(:,:), start = (/ 1, 1,ts /), count = (/ nlon, nlat, 1 /) )
      end if
    end if

    if((realrank.eq.0)) then
      ptr => sd
      call gather_data_to_master(ptr,clmvar_global_g,clmlevel=nameg)

      if(masterproc) then
        do g1 = 1, numg
          ji = adecomp%gdc2i(g1)
          jj = adecomp%gdc2j(g1)
          clmvar_out(ji,jj) = clmvar_global_g(g1)
        end do
        ierr = nf90_inq_varid(il_file_id, trim(variable_names(2)) , ncvarid(2))
        ierr = nf90_put_var( il_file_id, ncvarid(2), clmvar_out(:,:), start = (/ 1, 1,ts /), count = (/ nlon, nlat, 1 /) )
      end if
    end if





    !do i=1,n_variables
    ! call p2g(begp, endp, begc, endc, begl, endl, begg, endg, variable_pointers(i)%p, clmvar_local_g, &
    !          p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    !  mm = 0.0
    !  !call mpi_reduce(clmvar_local_g,mm, nloc, mpi_float, mpi_sum, 0, statcomm, ierr)
    !  call mpi_allreduce(clmvar_local_g,mm,nloc,MPI_REAL8,MPI_SUM,statcomm,ierr)
    !  mm = mm / realsize
    !  !write(*,*) 'mm[1]=',mm(1),'tmp_local[1]=',clmvar_local_g(1)
    !  !write(*,*) 'mm[2000]=',mm(2000),'tmp_local[2000]=',clmvar_local_g(2000)

    !  !call mpi_bcast(mm, nloc, mpi_double_precision, 0, statcomm, ierr)
    !  !call mpi_bcast(mm, nloc, mpi_float, 0, statcomm, ierr)
    !  !call mpi_bcast(mm, nloc, mpi_float, 0, statcomm, ierr)
    !  !write(*,*) 'mm[1] bcast=',mm(1)

    !  ptr => mm
    !  call gather_data_to_master(ptr,clmvar_global_g,clmlevel=nameg)

    !  if(masterproc .and. (realrank.eq.0)) then
    !    do g1 = 1, numg
    !      ji = adecomp%gdc2i(g1)
    !      jj = adecomp%gdc2j(g1)
    !      clmvar_out(ji,jj) = clmvar_global_g(g1)
    !    end do
    !    ierr = nf90_inq_varid(il_file_id, trim(variable_names(i)) , ncvarid(i))
    !    ierr = nf90_put_var( il_file_id, ncvarid(i), clmvar_out(:,:), start = (/ 1, 1,ts /), count = (/ nlon, nlat, 1 /) )
    !  end if

    !end do

    ! close netcdf output file
    if(masterproc .and. (realrank.eq.0)) then
      ierr = nf90_close(il_file_id)
    end if

    ! deallocate temporary arrays
    deallocate(mm)
    deallocate(sd)
    deallocate(var)
    deallocate(clmvar_local_g)
    if (masterproc) then
      deallocate(clmvar_global_g)
      deallocate (clmvar_out)
    end if

  end subroutine write_clm_statistics



  !> @author Wolfgang Kurtz, Guowei He
  !> @date 21.03.2022
  !> @brief Return the filename for output of CLM statistics
  !> @return character(len=256) get_statistic_filename Filename
  !> @details
  !>     This functions has likely been adapted from
  !>     `set_hist_filename` from `histFileMod.F90` in CLM3.5
  character(len=256) function get_statistic_filename ()
    !
    ! !DESCRIPTION:
    !  OLD:
    ! Determine history dataset filenames.
    !  NEW:
    ! Determine statistics dataset filename.
    ! !USES:
    use iso_c_binding, only: c_f_pointer, c_loc, c_null_char
    use clm_varctl, only : caseid
    use clm_time_manager, only : get_curr_date, get_prev_date
    use enkf_clm_mod, only : outdir
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !REVISION HISTORY:
    ! Created by Mariana Vertenstein
    ! Revised by Wolfgang Kurtz, Guowei He
    !
    !EOP
    !
    ! LOCAL VARIABLES:
    ! character(len=256) :: cdate       !date char string
    ! integer :: day                    !day (1 -> 31)
    ! integer :: mon                    !month (1 -> 12)
    ! integer :: yr                     !year (0 -> ...)
    ! integer :: sec                    !seconds into current day
    character(len=100),pointer :: pchar
    integer :: s_len
    !-----------------------------------------------------------------------

    !call get_prev_date (yr, mon, day, sec)
    !write(cdate,'(i4.4,"-",i2.2)') yr,mon                         !other
    !call get_curr_date (yr, mon, day, sec)
    !write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
    !get_statistic_filename = trim(caseid)//".stat.et."//trim(cdate)//".nc"
    call c_f_pointer(c_loc(outdir),pchar)
    s_len = index(pchar,c_null_char)-1

    get_statistic_filename = trim(pchar(1:s_len))//"/"//trim(caseid)//".stat.et.nc"

  end function get_statistic_filename


end module mod_clm_statistics
