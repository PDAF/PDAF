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
  !use spmdMod, only: comm_local_clm => mpicom, &
  !     npes_local_clm => npes, &
  !     mype_local_clm => iam, &
  !     is_masterproc_local_clm => masterproc
  use spmdmod, only: masterproc
  use shr_kind_mod, only: r8 => shr_kind_r8
  use mpi
  use enkf_clm_mod, only:statcomm
  use spmdgathscatmod , only : gather_data_to_master


  implicit none
  integer, parameter :: max_chars = 300

contains
  subroutine write_clm_statistics()
    USE clmtype,                  ONLY : clm3,nameg
    USE decompMod    , only : get_proc_global, get_proc_bounds, adecomp
    use netcdf
    ! kuw: templates from send_fld_2pfl.F90
    INTEGER :: numg           ! total number of gridcells across all processors
    INTEGER :: numl           ! total number of landunits across all processors
    INTEGER :: numc           ! total number of columns across all processors
    INTEGER :: nump           ! total number of pfts across all processors
    INTEGER :: begg,endg      ! local beg/end gridcells gdc
    INTEGER :: begl,endl      ! local beg/end landunits
    INTEGER :: begc,endc      ! local beg/end columns
    INTEGER :: begp,endp      ! local beg/end pfts
    
    !INTEGER ::   isec, info, jn, jj, ji, g1, jx, i    ! temporary integer
    integer :: nlon,nlat

    type field_pointer
       real(r8), pointer :: p(:)
    end type field_pointer

    integer, parameter :: n_variables = 4
    character(len = max_chars) :: variable_names(n_variables)
    type(field_pointer), dimension(n_variables) :: variable_pointers
    integer :: ncvarid(n_variables)
    real(r8), pointer :: clmstate_tmp_global(:)
    real(r8), allocatable ::   clmstate_out(:,:)

    character(len = max_chars) :: statistic_filename
    integer :: ierr,i,il_file_id,g1,nloc
    real(r8), allocatable,target :: mm(:),sd(:)
    real(r8),pointer :: ptr(:)
    integer,dimension(2) :: dimids
    integer ji,jj
    integer realrank,realsize

    real(r8), pointer :: lon(:)
    real(r8), pointer :: lat(:)

    ! Get total global number of grid cells, landunits, columns and pfts
    write(*,*) "clmstat: get proc bounds"
    CALL get_proc_global(numg,numl,numc,nump)
    CALL get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    if(masterproc) then 
      lon   => clm3%g%londeg
      lat   => clm3%g%latdeg
      !nlon = size(lon)
      !nlat = size(lat)
      nlon = adecomp%gdc2i(numg)
      nlat = adecomp%gdc2j(numg)
      !write(*,*) "nlon=",nlon
      !write(*,*) "nlat=",nlat
    end if 

    call mpi_comm_rank(statcomm,realrank,ierr)
    call mpi_comm_size(statcomm,realsize,ierr)
    

    variable_names(1) = "FCTR"
    variable_names(2) = "FCEV"
    variable_names(3) = "FGEV"
    variable_names(4) = "FSH"

    if(masterproc .and. (realrank.eq.0)) then
      statistic_filename = get_statistic_filename()
      ierr =  nf90_create(statistic_filename, NF90_CLOBBER, il_file_id)
      ierr =  nf90_def_dim(il_file_id, "lon", nlon, dimids(1))
      ierr =  nf90_def_dim(il_file_id, "lat", nlat, dimids(2))
      ierr =  nf90_def_var(il_file_id, trim(variable_names(1)), NF90_DOUBLE, dimids,ncvarid(1))
      ierr =  nf90_def_var(il_file_id, trim(variable_names(2)), NF90_DOUBLE, dimids,ncvarid(2))
      ierr =  nf90_def_var(il_file_id, trim(variable_names(3)), NF90_DOUBLE, dimids,ncvarid(3))
      ierr =  nf90_def_var(il_file_id, trim(variable_names(4)), NF90_DOUBLE, dimids,ncvarid(4))
      ierr =  nf90_enddef(il_file_id)
      !ierr = nf90_inq_varid(il_file_id, trim(variable_names(1)) , ncvarid(1))
      !ierr = nf90_put_var( il_file_id, ncvarid(1), clmstate_out(:,:,:), &
      !         start = (/ 1, 1, 1, 1/), count = (/ ndlon, ndlat, nlevsoi, 1 /) )
      !ierr = nf90_close(il_file_id)
      !write(*,*) 'clmstat: ncdf file defined'
    end if
    

    variable_pointers(1)%p => clm3%g%l%c%p%pef%eflx_lh_vegt
    variable_pointers(2)%p => clm3%g%l%c%p%pef%eflx_lh_vege
    variable_pointers(3)%p => clm3%g%l%c%p%pef%eflx_lh_grnd
    variable_pointers(4)%p => clm3%g%l%c%p%pef%eflx_sh_tot
    
    !allocate(clmstate_tmp_local(begc:endc), stat=statur)
    if (masterproc) then
      allocate(clmstate_tmp_global(numg), stat=ierr)
      allocate (clmstate_out(nlon,nlat), stat=ierr)
    end if

    nloc = endc-begc+1
    !write(*,*)"nloc=",nloc
    !write(*,*)"size eflx_lh_vegt=",size(variable_pointers(1)%p)
    !write(*,*)"size eflx_lh_vege=",size(variable_pointers(2)%p)
    !write(*,*)"size eflx_lh_grnd=",size(variable_pointers(3)%p)
    !write(*,*)"size eflx_sh_tot=",size(variable_pointers(4)%p)
    allocate(mm(begc:endc))
    allocate(sd(begc:endc))
   
    do i=1,1!n_variables
      !write(*,*) 'mpi_reduce: variable '//trim(variable_names(i))
      !call mpi_reduce(mm, variable_pointers(i)%p, nloc, mpi_double_precision, mpi_sum, 0, statcomm, ierr)
      call mpi_reduce(variable_pointers(i)%p,mm, nloc, mpi_float, mpi_sum, 0, statcomm, ierr)
      mm = mm / realsize
      !write(*,*) 'mpi_bcast: variable '//trim(variable_names(i))
      write(*,*) 'mm[1]=',mm(1)
      !call mpi_bcast(mm, nloc, mpi_double_precision, 0, statcomm, ierr)
      call mpi_bcast(mm, nloc, mpi_float, 0, statcomm, ierr)
      write(*,*) 'mm[1] bcast=',mm(1)

      ptr => mm
      !write(*,*) 'gather_data_to_maseter: variable '//trim(variable_names(i))
      call gather_data_to_master(ptr,clmstate_tmp_global,clmlevel=nameg)

      !if(masterproc) then
      if(masterproc .and. (realrank.eq.0)) then
      write(*,*) 'clmstate_tmp_global[1]=',clmstate_tmp_global(1)
        do g1 = 1, numg
           ji = adecomp%gdc2i(g1)
           jj = adecomp%gdc2j(g1)
           clmstate_out(ji,jj) = clmstate_tmp_global(g1)
        end do
        ierr = nf90_inq_varid(il_file_id, trim(variable_names(i)) , ncvarid(i))
        ierr = nf90_put_var( il_file_id, ncvarid(i), clmstate_out(:,:), start = (/ 1, 1, 1, 1/), count = (/ nlon, nlat /) )
      end if

    end do

    !if(masterproc) then
    if(masterproc .and. (realrank.eq.0)) then
      ierr = nf90_close(il_file_id)
    end if
  
    deallocate(mm)
    deallocate(sd)
    if (masterproc) then
      deallocate(clmstate_tmp_global)
      deallocate (clmstate_out)
    end if


  end subroutine write_clm_statistics

!  subroutine compute_clm_statistics(buff, mean, variance)
!    use shr_kind_mod, only: r8 => shr_kind_r8
!    real(r8), pointer, intent(in) :: buff(:)
!    double precision, intent(out) :: mean
!    double precision, intent(out) :: variance
!
!
!    mean =  compute_mpi_mean(buff, comm_local_clm)
!    variance = compute_mpi_variance(buff, comm_local_clm, mean)
!
!  end subroutine compute_clm_statistics
!
!  subroutine write_clm_statistics_to_netcdf(filename, variable_names, means, variances)
!    use netcdf
!    implicit none
!    integer, parameter :: NDIMS = 1
!    integer, parameter :: NX = 6
!    integer, parameter :: n_statistics = 2
!    character(len = *) :: filename
!    character(len = *) :: variable_names(:)
!    double precision, intent(in) :: means(:)
!    double precision, intent(in) :: variances(:)
!    integer :: ncid, dimids
!    integer, allocatable ::  varid(:)
!    integer :: x_dimid
!
!    integer :: i
!    integer :: n_variables
!    real*8 :: value_to_netcdf

!    n_variables = size(variable_names)
!    allocate (varid(n_statistics))
!    call check( nf90_create(trim(filename), NF90_NETCDF4, ncid) )
!    
!    print *, "n_variables ", n_variables
!    call check( nf90_def_dim(ncid, "scalar", NX, x_dimid) )
!    dimids = x_dimid
!    if (.false.) then
!    do i = 1, n_variables
!!       call check( nf90_def_var(ncid, "mean"//"_"//trim(variable_names(i)),  NF90_DOUBLE, dimids, varid(1,i)) )
!       call check( nf90_def_var(ncid, "mean"//"_"//trim(variable_names(i)),  NF90_DOUBLE, varid=varid(1)) )
!!       call check( nf90_def_var(ncid, "variance"//"_"//trim(variable_names(i)),  NF90_DOUBLE, dimids, varid(2,i)) )
!       call check( nf90_def_var(ncid, "variance"//"_"//trim(variable_names(i)),  NF90_DOUBLE,varid= varid(2)) )
!    end do
!    call check( nf90_enddef(ncid) )
!    end if
!
!    call check(nf90_def_var(ncid, "mean", NF90_DOUBLE, x_dimid, varid(1)))
!    call check(nf90_def_var(ncid, "variance", NF90_DOUBLE, x_dimid, varid(2)))
!    call check( nf90_enddef(ncid) )
!    print *, "statistics field def finished"
!    if (.false.) then
!    do i = 1, n_variables
!       !            call check( nf90_put_var(ncid, varid(1,i), means(i)) )
!       !            call check( nf90_put_var(ncid, varid(2,i), variances(i)) )
!       value_to_netcdf = means(i)
!!       call check( nf90_put_var(ncid, varid(1,i), value_to_netcdf) )
!       value_to_netcdf = variances(i)
! !      call check( nf90_put_var(ncid, varid(2,i), value_to_netcdf) )
!    end do
! end if
!
!    call check(nf90_put_var(ncid, varid(1), means))
!    call check(nf90_put_var(ncid, varid(2), variances))
!    
! 
!    call check( nf90_close(ncid) )
!    print *, "statistics file closing"
!    deallocate(varid)
!  end subroutine write_clm_statistics_to_netcdf
!
!  subroutine check(ierr)
!    use netcdf
!    integer, intent ( in) :: ierr
!
!    if(ierr /= nf90_noerr) then
!       print *, trim(nf90_strerror(ierr))
!       stop "Stopped"
!    end if
!  end subroutine check
!
!
!  double precision function compute_mpi_variance(buff, comm, mean)
!    use shr_kind_mod, only: r8 => shr_kind_r8
!    real(r8), pointer, intent(in) :: buff(:)
!    integer, intent(in) :: comm
!    double precision, intent(in) :: mean
!    double precision :: sum_var_local
!    double precision :: sum_var_global
!    integer :: ierror
!    integer :: size_local
!    integer :: size_global
!
!    print *, "computing variance"
!    size_local = size(buff)
!    print *, "size of buff is ", size_local
!    call mpi_reduce(size_local, size_global, 1, mpi_integer, mpi_sum, 0, comm, ierror)
!
!    ! compute local sum of variance
!    sum_var_local = sum((buff - mean) ** 2)
!
!    ! gather and bcast because allreduce doesn't work
!    call mpi_reduce(sum_var_local, sum_var_global, 1, mpi_double_precision, mpi_sum, 0, comm, ierror)
!    if (mype_local_clm .eq. 0) compute_mpi_variance  = sum_var_global / size_global
!    call mpi_bcast(compute_mpi_variance, 1, mpi_double_precision, 0, comm, ierror)
!
!  end function compute_mpi_variance
!
!  double precision function compute_mpi_mean(buff, comm)
!    use shr_kind_mod, only: r8 => shr_kind_r8
!    real(r8), pointer, intent(in) :: buff(:)
!    integer, intent(in) :: comm
!    integer :: size_local
!    integer :: size_global
!    integer :: ierror
!    double precision :: sum_local
!    integer :: npes_local
!
!    call mpi_comm_size(comm, npes_local, ierror)
!    print *, npes_local, " in this CLM component"
!    print *, "computing mean"
!    size_local = size(buff)
!    print *, "size of buff is ", size_local
!    call mpi_reduce(size_local, size_global, 1, mpi_integer, mpi_sum, 0, comm, ierror)
!    call mpi_bcast(size_global, 1, mpi_integer, 0, comm, ierror)
!    call mpi_allreduce(size_local, size_global, 1, mpi_integer, mpi_sum, 0, comm, ierror)
!    print *, "global size of buff is ", size_global
!
!    sum_local = sum(buff)
!    print *, "local sum is ",sum_local
!    call mpi_reduce(sum_local, compute_mpi_mean, 1, mpi_double_precision, mpi_sum, 0, comm, ierror)
!    compute_mpi_mean = compute_mpi_mean / size_global
!    call mpi_bcast(compute_mpi_mean, 1, mpi_double_precision, 0, comm, ierror)
!
!    !call mpi_allreduce(sum_local, compute_mpi_mean, 1, mpi_double_precision, mpi_sum, 0, comm, ierror)
!  end function compute_mpi_mean
!
  character(len=256) function get_statistic_filename ()
    !
    ! !DESCRIPTION:
    ! Determine history dataset filenames.
    !
    ! !USES:
    use clm_varctl, only : caseid
    use clm_time_manager, only : get_curr_date, get_prev_date
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !REVISION HISTORY:
    ! Created by Mariana Vertenstein
    !
    !EOP
    !
    ! LOCAL VARIABLES:
    character(len=256) :: cdate       !date char string
    integer :: day                    !day (1 -> 31)
    integer :: mon                    !month (1 -> 12)
    integer :: yr                     !year (0 -> ...)
    integer :: sec                    !seconds into current day
    !-----------------------------------------------------------------------

    call get_prev_date (yr, mon, day, sec)
    write(cdate,'(i4.4,"-",i2.2)') yr,mon                         !other
    call get_curr_date (yr, mon, day, sec)
    write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
    get_statistic_filename = trim(caseid)//".stat.et."//trim(cdate)//".nc"

  end function get_statistic_filename
 
 
end module mod_clm_statistics
