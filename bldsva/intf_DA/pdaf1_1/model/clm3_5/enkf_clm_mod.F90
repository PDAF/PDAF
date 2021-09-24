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
!enkf_clm_mod.F90: Module for CLM
!-------------------------------------------------------------------------------------------

module enkf_clm_mod
 
  use iso_c_binding

! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8, SHR_KIND_CL
  use shr_orb_mod          
  use clm_varorb      , only : eccen, mvelpp, lambm0, obliqr, obliq, &
                               iyear_AD, nmvelp
  use clm_comp        , only : clm_init0, clm_init1, clm_init2, clm_run1, clm_run2
  use clm_time_manager, only : is_last_step, advance_timestep, get_nstep
  use atmdrvMod       , only : atmdrv, atmdrv_init
  use abortutils      , only : endrun
  use controlMod      , only : control_setNL
  use clm_mct_mod
  use spmdMod  
  use ESMF_Mod
  use perf_mod

#if ((defined COUP_OAS_COS || defined COUP_OAS_PFL) && (!defined CLMSA))
  use oas_clm_vardef , only : kl_comm
#endif
!
! !ARGUMENTS:
    implicit none

#if (defined CLMSA)
  integer :: da_comm
  integer :: clm_statevecsize
  integer :: clm_varsize
  integer :: clm_begg,clm_endg
  real(r8),allocatable :: clm_statevec(:)
  integer(c_int),bind(C,name="clmupdate_swc")     :: clmupdate_swc
  integer(c_int),bind(C,name="clmupdate_texture") :: clmupdate_texture
  integer(c_int),bind(C,name="clmprint_swc")      :: clmprint_swc
#endif
  integer(c_int),bind(C,name="clmprint_et")       :: clmprint_et

  integer  :: nstep     ! time step index
  real(r8) :: dtime     ! time step increment (sec)
  integer  :: ier       ! error code

  character(c_char),bind(C,name="outdir"),target :: outdir

  logical  :: log_print    ! true=> print diagnostics
  real(r8) :: eccf         ! earth orbit eccentricity factor
  logical  :: mpi_running  ! true => MPI is initialized 
  integer  :: mpicom_glob  ! MPI communicator

  character(len=SHR_KIND_CL) :: nlfilename = " "
  integer :: ierror, lengths_of_types, i
  logical :: flag
  integer(c_int),bind(C,name="clmprefixlen") :: clmprefixlen
  integer :: statcomm

  contains

#if defined CLMSA
  subroutine define_clm_statevec()
    use shr_kind_mod, only: r8 => shr_kind_r8
    use decompMod , only : get_proc_bounds
    use clm_varpar   , only : nlevsoi

    implicit none

    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    clm_begg     = begg
    clm_endg     = endg  

    if(clmupdate_swc.eq.1) then
      clm_varsize      =  (endg-begg+1) * nlevsoi
      clm_statevecsize =  (endg-begg+1) * nlevsoi
    endif

    if(clmupdate_swc.eq.2) then
      clm_varsize      =  (endg-begg+1) * nlevsoi
      clm_statevecsize =  (endg-begg+1) * (nlevsoi+1)
    endif

    if(clmupdate_texture.eq.1) then
        clm_statevecsize = clm_statevecsize + 2*((endg-begg+1)*nlevsoi)
    endif

    !write(*,*) 'clm_statevecsize is ',clm_statevecsize
    IF (allocated(clm_statevec)) deallocate(clm_statevec)
    allocate(clm_statevec(clm_statevecsize))
  end subroutine

  subroutine set_clm_statevec()
    USE clmtype      , only : clm3
    USE clm_varpar   , only : nlevsoi
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    real(r8), pointer :: swc(:,:)
    real(r8), pointer :: psand(:,:)
    real(r8), pointer :: pclay(:,:)
    integer :: i,j,cc=1,offset=0
    real(r8) :: swcave
    real(r8) :: swctmp(10)

    swc   => clm3%g%l%c%cws%h2osoi_vol
    psand => clm3%g%l%c%cps%psand
    pclay => clm3%g%l%c%cps%pclay

    ! calculate shift when CRP data are assimilated
    if(clmupdate_swc.eq.2) then
      offset = clm_endg-clm_begg+1
    endif

    ! write swc values to state vector
    cc = 1
    do i=1,nlevsoi
      do j=clm_begg,clm_endg
        clm_statevec(cc+offset) = swc(j,i)
        cc = cc + 1
      end do
    end do

    ! write average swc to state vector (CRP assimilation)
    if(clmupdate_swc.eq.2) then
      cc = 1
      do j=clm_begg,clm_endg
        do i=1,nlevsoi
          swctmp(i) = swc(j,i)
        end do
        call average_swc_crp(swctmp,swcave)
        clm_statevec(cc) = swcave
        cc = cc + 1
      end do

    endif

    ! write texture values to state vector (if desired)
    if(clmupdate_texture.eq.1) then
      cc = 1
      do i=1,nlevsoi
        do j=clm_begg,clm_endg
          clm_statevec(cc+1*clm_varsize+offset) = psand(j,i)
          clm_statevec(cc+2*clm_varsize+offset) = pclay(j,i)
          cc = cc + 1
        end do
      end do
    endif
  end subroutine 

  subroutine update_clm()
    USE clmtype      , only : clm3
    USE clm_varpar   , only : nlevsoi
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clm_varcon   , only : denh2o,denice

    implicit none
    real(r8), pointer :: swc(:,:)
    real(r8), pointer :: watsat(:,:)
    real(r8), pointer :: psand(:,:)
    real(r8), pointer :: pclay(:,:)

    real(r8), pointer :: dz(:,:)          ! layer thickness depth (m)
    real(r8), pointer :: h2osoi_liq(:,:)  ! liquid water (kg/m2)
    real(r8), pointer :: h2osoi_ice(:,:)  
    real(r8)  :: rliq,rice

    integer :: i,j,cc=1,offset=0

    swc   => clm3%g%l%c%cws%h2osoi_vol
    watsat => clm3%g%l%c%cps%watsat
    psand => clm3%g%l%c%cps%psand
    pclay => clm3%g%l%c%cps%pclay

    dz            => clm3%g%l%c%cps%dz
    h2osoi_liq    => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_ice    => clm3%g%l%c%cws%h2osoi_ice

    ! calculate shift when CRP data are assimilated
    if(clmupdate_swc.eq.2) then
      offset = clm_endg-clm_begg+1
    endif

    ! write updated swc back to CLM
    cc = 1
    do i=1,nlevsoi
      do j=clm_begg,clm_endg
        rliq = h2osoi_liq(j,i)/(dz(j,i)*denh2o*swc(j,i))
        rice = h2osoi_ice(j,i)/(dz(j,i)*denice*swc(j,i))
        !h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)
        if(clm_statevec(cc+offset).le.0.0) then
          swc(j,i)   = 0.05
        else if(clm_statevec(cc+offset).ge.watsat(j,i)) then
          swc(j,i) = watsat(j,i)
        else
          swc(j,i)   = clm_statevec(cc+offset)
        endif
        ! update liquid water content
        h2osoi_liq(j,i) = swc(j,i) * dz(j,i)*denh2o*rliq
        ! update ice content
        h2osoi_ice(j,i) = swc(j,i) * dz(j,i)*denice*rice
        cc = cc + 1
      end do
    end do

    !! update liquid water content
    !do j=clm_begg,clm_endg
    !  do i=1,nlevsoi
    !    h2osoi_liq(j,i) = swc(j,i) * dz(j,i)*denh2o
    !  end do
    !end do

    ! write updated texture back to CLM
    if(clmupdate_texture.eq.1) then
      cc = 1
      do i=1,nlevsoi
        do j=clm_begg,clm_endg
          psand(j,i) = clm_statevec(cc+1*clm_varsize+offset)
          pclay(j,i) = clm_statevec(cc+2*clm_varsize+offset)
          cc = cc + 1
        end do
      end do
      call clm_correct_texture
      call clm_texture_to_parameters
    endif
  end subroutine 

  subroutine clm_correct_texture
    use clmtype
    use clm_varpar  , only : nlevsoi
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    integer :: c,lev
    real(r8) :: clay,sand,ttot
    real(r8), pointer :: psand(:,:)
    real(r8), pointer :: pclay(:,:)

    psand => clm3%g%l%c%cps%psand
    pclay => clm3%g%l%c%cps%pclay

    do c = clm_begg,clm_endg
      do lev = 1,nlevsoi
         clay = pclay(c,lev)
         sand = psand(c,lev)

         if(sand.le.0.0) sand = 1.0
         if(clay.le.0.0) clay = 1.0
         ttot = sand + clay
         if(ttot.gt.100) then
             sand = sand/ttot * 100
             clay = clay/ttot * 100
         end if
          
         pclay(c,lev) = clay
         psand(c,lev) = sand
      end do
    end do
  end subroutine

  subroutine clm_texture_to_parameters()
    use clmtype
    use clm_varpar  , only : nlevsoi
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    integer :: c,lev
    
    real(r8) :: clay,sand        ! temporaries
    real(r8) :: bd               ! bulk density of dry soil material [kg/m^3]
    real(r8) :: xksat            ! maximum hydraulic conductivity of soil [mm/s]
    real(r8) :: tkm              ! mineral conductivity
    real(r8), pointer :: watsat(:,:)        ! volumetric soil water at saturation (porosity) (nlevsoi) 
    real(r8), pointer :: bsw(:,:)           ! Clapp and Hornberger "b" (nlevsoi)  
    real(r8), pointer :: bsw2(:,:)          ! Clapp and Hornberger "b" for CN code
    real(r8), pointer :: psisat(:,:)        ! soil water potential at saturation for CN code (MPa)
    real(r8), pointer :: vwcsat(:,:)        ! volumetric water content at saturation for CN code (m3/m3)
    real(r8), pointer :: hksat(:,:)         ! hydraulic conductivity at saturation (mm H2O /s) (nlevsoi) 
    real(r8), pointer :: sucsat(:,:)        ! minimum soil suction (mm) (nlevsoi) 
    real(r8), pointer :: tkmg(:,:)          ! thermal conductivity, soil minerals  [W/m-K] (new) (nlevsoi) 
    real(r8), pointer :: tksatu(:,:)        ! thermal conductivity, saturated soil [W/m-K] (new) (nlevsoi) 
    real(r8), pointer :: tkdry(:,:)         ! thermal conductivity, dry soil (W/m/Kelvin) (nlevsoi) 
    real(r8), pointer :: csol(:,:)          ! heat capacity, soil solids (J/m**3/Kelvin) (nlevsoi) 
    real(r8), pointer :: watdry(:,:)        ! btran parameter for btran=0
    real(r8), pointer :: watopt(:,:)        ! btran parameter for btran = 1
    real(r8), pointer :: psand(:,:)
    real(r8), pointer :: pclay(:,:)
    watsat          => clm3%g%l%c%cps%watsat
    tkmg            => clm3%g%l%c%cps%tkmg
    bsw             => clm3%g%l%c%cps%bsw
    bsw2            => clm3%g%l%c%cps%bsw2
    psisat          => clm3%g%l%c%cps%psisat
    vwcsat          => clm3%g%l%c%cps%vwcsat
    hksat           => clm3%g%l%c%cps%hksat
    sucsat          => clm3%g%l%c%cps%sucsat
    tkmg            => clm3%g%l%c%cps%tkmg
    tksatu          => clm3%g%l%c%cps%tksatu
    tkdry           => clm3%g%l%c%cps%tkdry
    csol            => clm3%g%l%c%cps%csol
    watdry          => clm3%g%l%c%cps%watdry  
    watopt          => clm3%g%l%c%cps%watopt  
    psand            => clm3%g%l%c%cps%psand
    pclay            => clm3%g%l%c%cps%pclay

    do c = clm_begg,clm_endg
      do lev = 1,nlevsoi
         clay = pclay(c,lev)
         sand = psand(c,lev)
         watsat(c,lev) = 0.489_r8 - 0.00126_r8*sand
         bd = (1._r8-watsat(c,lev))*2.7e3_r8
         xksat = 0.0070556_r8 *( 10._r8**(-0.884_r8+0.0153_r8*sand) ) ! mm/s
         tkm = (8.80_r8*sand+2.92_r8*clay)/(sand+clay)          ! W/(m K)

         bsw(c,lev) = 2.91_r8 + 0.159_r8*clay
         bsw2(c,lev) = -(3.10_r8 + 0.157_r8*clay - 0.003_r8*sand)
         psisat(c,lev) = -(exp((1.54_r8 - 0.0095_r8*sand + 0.0063_r8*(100.0_r8-sand-clay))*log(10.0_r8))*9.8e-5_r8)
         vwcsat(c,lev) = (50.5_r8 - 0.142_r8*sand - 0.037_r8*clay)/100.0_r8
         hksat(c,lev) = xksat
         sucsat(c,lev) = 10._r8 * ( 10._r8**(1.88_r8-0.0131_r8*sand) )
         tkmg(c,lev) = tkm ** (1._r8- watsat(c,lev))
         tksatu(c,lev) = tkmg(c,lev)*0.57_r8**watsat(c,lev)
         tkdry(c,lev) = (0.135_r8*bd + 64.7_r8) / (2.7e3_r8 - 0.947_r8*bd)
         csol(c,lev) = (2.128_r8*sand+2.385_r8*clay) / (sand+clay)*1.e6_r8  ! J/(m3 K)
         watdry(c,lev) = watsat(c,lev) * (316230._r8/sucsat(c,lev)) ** (-1._r8/bsw(c,lev)) 
         watopt(c,lev) = watsat(c,lev) * (158490._r8/sucsat(c,lev)) ** (-1._r8/bsw(c,lev)) 
      end do
    end do
  end subroutine
 
  subroutine  average_swc_crp(profdat,profave)
    use clm_varcon  , only : zsoi

    implicit none

    real(r8),intent(in)  :: profdat(10)
    real(r8),intent(out) :: profave

    real(r8) :: w(10)
    real(r8) :: mold=0,mnew=0,delta=1,g,wtmp=0
    integer  :: iter=0,i,j

    w = 1/10
    do while(delta>0.0001)

      ! calculate initial weighted mean
      mold = 0
      do i=1,10
        mold = mold + w(i)*profdat(i)
      end do

      ! calculate g factor
      g = -5.8/(log(0.14)*(mold+0.0829))

      ! calculate new weights
      w(1) =  1-exp(-zsoi(1)*100/g)
      do i=2,10
        ! calculate sum of previous weights
        wtmp = 0
        do j=1,(i-1)
          wtmp = wtmp + w(j)
        end do
        w(i) = 1-exp(-zsoi(i)*100/g)-wtmp
      end do

      ! calculate new weighted mean
      mnew = 0
      do i=1,10
        mnew = mnew + w(i)*profdat(i)
      end do
    
      ! compare old and new weighted mean
      iter=iter+1
      delta = abs(mnew-mold)
    end do
    
    profave = mnew 

  end subroutine average_swc_crp
#endif

  subroutine domain_def_clm(lon_clmobs, lat_clmobs, dim_obs, longxy, latixy, longxy_obs, latixy_obs)
!#if defined CLMSA
    USE domainMod, ONLY: alatlon
    USE decompMod, ONLY: get_proc_total, get_proc_bounds_atm, adecomp
    USE spmdMod,   ONLY: npes, iam ! number of processors for clm and processor number
   ! USE mod_read_obs, ONLY: lon_clmobs, lat_clmobs
    USE clmtype,    ONLY : clm3
   ! USE mod_assimilation, ONLY: dim_obs
!#endif

    implicit none
    real, intent(in) :: lon_clmobs(:)
    real, intent(in) :: lat_clmobs(:)
    integer, intent(in) :: dim_obs
    integer, allocatable, intent(inout) :: longxy(:)
    integer, allocatable, intent(inout) :: latixy(:)
    integer, allocatable, intent(inout) :: longxy_obs(:)
    integer, allocatable, intent(inout) :: latixy_obs(:)
    integer :: ni, nj, ii, jj, kk, cid, ier, ncells, nlunits, &
               ncols, npfts, counter
    real :: minlon, minlat, maxlon, maxlat
    real(r8), pointer :: lon(:)
    real(r8), pointer :: lat(:)
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices

    lon   => alatlon%lonc
    lat   => alatlon%latc
    ni  = alatlon%ni
    nj  = alatlon%nj

    ! get total number of gridcells, landunits,
    ! columns and pfts for any processor
    call get_proc_total(iam, ncells, nlunits, ncols, npfts)

    ! beg and end gridcell for atm
    call get_proc_bounds_atm(begg, endg)

   !print *,'ni, nj ', ni, nj
   !print *,'cells per processor ', ncells
   !print *,'begg, endg ', begg, endg
 
    ! allocate vector with size of elements in x directions * size of elements in y directions 
    if(allocated(longxy)) deallocate(longxy)
    allocate(longxy(ncells), stat=ier)
    if(allocated(latixy)) deallocate(latixy)
    allocate(latixy(ncells), stat=ier)
    ! initialize vector with zero values
    longxy(:) = 0
    latixy(:) = 0
    counter = 1
    do ii = 1, nj
      do jj = 1, ni
         cid = (ii-1)*ni + jj
         do kk = begg, endg  
            if(cid == adecomp%gdc2glo(kk)) then  
               latixy(counter) = ii
               longxy(counter) = jj
               counter = counter + 1
            endif   
         enddo  
      end do
    end do

    ! set intial values for max/min of lon/lat
    minlon = 999
    minlat = 999
    maxlon = -999
    maxlat = -999

    ! looping over all cell centers to get min/max longitude and latitude
    minlon = MINVAL(lon(:) + 180)
    maxlon = MAXVAL(lon(:) + 180)
    minlat = MINVAL(lat(:) + 90) 
    maxlat = MAXVAL(lat(:) + 90)

    if(allocated(longxy_obs)) deallocate(longxy_obs)
    allocate(longxy_obs(dim_obs), stat=ier) 
    if(allocated(latixy_obs)) deallocate(latixy_obs)
    allocate(latixy_obs(dim_obs), stat=ier)
    do i = 1, dim_obs
       if(((lon_clmobs(i) + 180) - minlon) /= 0 .and. ((lat_clmobs(i) + 90) - minlat) /= 0) then
          longxy_obs(i) = ceiling(((lon_clmobs(i) + 180) - minlon) * ni / (maxlon - minlon)) !+ 1
          latixy_obs(i) = ceiling(((lat_clmobs(i) + 90) - minlat) * nj / (maxlat - minlat)) !+ 1
          !print *,'longxy_obs(i) , latixy_obs(i) ', longxy_obs(i) , latixy_obs(i)
       else if(((lon_clmobs(i) + 180) - minlon) == 0 .and. ((lat_clmobs(i) + 90) - minlat) == 0) then
          longxy_obs(i) = 1
          latixy_obs(i) = 1
       else if(((lon_clmobs(i) + 180) - minlon) == 0) then
          longxy_obs(i) = 1
          latixy_obs(i) = ceiling(((lat_clmobs(i) + 90) - minlat) * nj / (maxlat - minlat))
       else if(((lat_clmobs(i) + 90) - minlat) == 0) then
          longxy_obs(i) = ceiling(((lon_clmobs(i) + 180) - minlon) * ni / (maxlon - minlon))
          latixy_obs(i) = 1
       endif   
    end do   
    ! deallocate temporary arrays
    !deallocate(longxy)
    !deallocate(latixy)
    !deallocate(longxy_obs)
    !deallocate(latixy_obs)    

  end subroutine  

#if defined CLMSA
  subroutine init_clm_l_size(dim_l)
    use clm_varpar   , only : nlevsoi

    implicit none
    integer, intent(out) :: dim_l   
    integer              :: nshift
      
    if(clmupdate_swc.eq.1) then
      dim_l = nlevsoi
      nshift = nlevsoi
    endif

    if(clmupdate_swc.eq.2) then
      dim_l = nlevsoi + 1
      nshift = nlevsoi + 1
    endif

    if(clmupdate_texture.eq.1) then
      dim_l = 2*nlevsoi + nshift
    endif    

  end subroutine    
#endif

end module enkf_clm_mod

