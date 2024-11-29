!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz, Guowei He and Mukund Pondkule (Forschungszentrum Juelich GmbH)
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
!enkf_clm_mod.F90: Module for CLM
!-------------------------------------------------------------------------------------------

module enkf_clm_mod

  use iso_c_binding

! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8, SHR_KIND_CL

! !ARGUMENTS:
    implicit none

#if (defined CLMSA)
  integer :: COMM_model_clm
  integer :: clm_statevecsize
  integer :: clm_paramsize !hcp: Size of CLM parameter vector (f.e. LAI)
  integer :: clm_varsize
  integer :: clm_begg,clm_endg
  integer :: clm_begc,clm_endc
  integer :: clm_begp,clm_endp
  real(r8),allocatable :: clm_statevec(:)
  integer,allocatable :: state_pdaf2clm_c_p(:)
  integer,allocatable :: state_pdaf2clm_j_p(:)
  ! clm_paramarr: Contains LAI used in obs_op_pdaf for computing model
  ! LST in LST assimilation (clmupdate_T)
  real(r8),allocatable :: clm_paramarr(:)  !hcp CLM parameter vector (f.e. LAI)
  integer, allocatable :: state_clm2pdaf_p(:,:) !Index of column in hydraulic active state vector (nlevsoi,endc-begc+1)
  integer(c_int),bind(C,name="clmupdate_swc")     :: clmupdate_swc
  integer(c_int),bind(C,name="clmupdate_T")     :: clmupdate_T  ! by hcp
  integer(c_int),bind(C,name="clmupdate_texture") :: clmupdate_texture
  integer(c_int),bind(C,name="clmprint_swc")      :: clmprint_swc
#endif
  integer(c_int),bind(C,name="clmprint_et")       :: clmprint_et
  integer(c_int),bind(C,name="clmstatevec_allcol")       :: clmstatevec_allcol
  integer(c_int),bind(C,name="clmstatevec_only_active")  :: clmstatevec_only_active
  integer(c_int),bind(C,name="clmstatevec_max_layer")  :: clmstatevec_max_layer
  integer(c_int),bind(C,name="clmt_printensemble")       :: clmt_printensemble
  integer(c_int),bind(C,name="clmwatmin_switch")         :: clmwatmin_switch

  integer  :: nstep     ! time step index
  real(r8) :: dtime     ! time step increment (sec)
  integer  :: ier       ! error code

  character(kind=c_char,len=100),bind(C,name="outdir"),target :: outdir

  logical  :: log_print    ! true=> print diagnostics
  real(r8) :: eccf         ! earth orbit eccentricity factor
  logical  :: mpi_running  ! true => MPI is initialized
  integer  :: mpicom_glob  ! MPI communicator

  character(len=SHR_KIND_CL) :: nlfilename = " "
  integer :: ierror, lengths_of_types, i
  logical :: flag
  integer(c_int),bind(C,name="clmprefixlen") :: clmprefixlen
  integer :: COMM_couple_clm    ! CLM-version of COMM_couple
                                ! (currently not used for clm5_0)
  logical :: newgridcell        !only clm5_0

  contains

#if defined CLMSA
  subroutine define_clm_statevec(mype)
    use shr_kind_mod, only: r8 => shr_kind_r8
    use decompMod , only : get_proc_bounds
    use clm_varpar   , only : nlevsoi
    use clm_varcon , only : ispval
    use ColumnType , only : col

    implicit none

    integer,intent(in) :: mype

    integer :: i
    integer :: j
    integer :: jj
    integer :: c
    integer :: g
    integer :: cc
    integer :: cccheck

    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices


    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

#ifdef PDAF_DEBUG
    WRITE(*,"(a,i5,a,i10,a,i10,a,i10,a,i10,a,i10,a,i10,a,i10,a,i10,a)") &
      "TSMP-PDAF mype(w)=", mype, " define_clm_statevec, CLM5-bounds (g,l,c,p)----",&
      begg,",",endg,",",begl,",",endl,",",begc,",",endc,",",begp,",",endp," -------"
#endif

    clm_begg     = begg
    clm_endg     = endg
    clm_begc     = begc
    clm_endc     = endc
    clm_begp     = begp
    clm_endp     = endp

    if(clmupdate_swc.eq.1) then
      if(clmstatevec_allcol.eq.1) then

        if(clmstatevec_only_active .eq. 1) then

          IF (allocated(state_clm2pdaf_p)) deallocate(state_clm2pdaf_p)
          allocate(state_clm2pdaf_p(begc:endc,nlevsoi))

          cc = 0

          do i=1,nlevsoi
            do c=clm_begc,clm_endc
              ! Only take into account layers above input maximum layer
              if(i<=clmstatevec_max_layer) then
                ! Only take into account hydrologically active columns
                ! and layers above bedrock
                if(col%hydrologically_active(c) .and. i<=col%nbedrock(c)) then
                  cc = cc + 1
                  state_clm2pdaf_p(c,i) = cc
                else
                  state_clm2pdaf_p(c,i) = ispval
                end if
              else
                state_clm2pdaf_p(c,i) = ispval
              end if
            end do
          end do

          ! Set `clm_varsize`, even though it is currently not used
          ! for `clmupdate_swc.eq.1`
          clm_varsize = clm_statevecsize
          clm_statevecsize = cc

          IF (allocated(state_pdaf2clm_c_p)) deallocate(state_pdaf2clm_c_p)
          allocate(state_pdaf2clm_c_p(clm_statevecsize))
          IF (allocated(state_pdaf2clm_j_p)) deallocate(state_pdaf2clm_j_p)
          allocate(state_pdaf2clm_j_p(clm_statevecsize))

          cc = 0

          do i=1,nlevsoi
            do c=clm_begc,clm_endc
              ! Only take into account layers above input maximum layer
              if(i<=clmstatevec_max_layer) then
                ! Only take into account hydrologically active columns
                ! and layers above bedrock
                if(col%hydrologically_active(c) .and. i<=col%nbedrock(c)) then
                  cc = cc + 1
                  state_pdaf2clm_c_p(cc) = c
                  state_pdaf2clm_j_p(cc) = i
                end if
              end if
            end do
          end do

        else

          IF (allocated(state_clm2pdaf_p)) deallocate(state_clm2pdaf_p)
          allocate(state_clm2pdaf_p(begc:endc,nlevsoi))

          do i=1,nlevsoi
            do c=clm_begc,clm_endc
              state_clm2pdaf_p(c,i) = (c - clm_begc + 1) + (i - 1)*(clm_endc - clm_begc + 1)
            end do
          end do

          ! #cols values per grid-cell
          clm_varsize      =  (endc-begc+1) * nlevsoi
          clm_statevecsize =  (endc-begc+1) * nlevsoi

          IF (allocated(state_pdaf2clm_c_p)) deallocate(state_pdaf2clm_c_p)
          allocate(state_pdaf2clm_c_p(clm_statevecsize))
          IF (allocated(state_pdaf2clm_j_p)) deallocate(state_pdaf2clm_j_p)
          allocate(state_pdaf2clm_j_p(clm_statevecsize))

          cc = 0

          do i=1,nlevsoi
            do c=clm_begc,clm_endc
              cc = cc + 1
              state_pdaf2clm_c_p(cc) = c
              state_pdaf2clm_j_p(cc) = i
            end do
          end do

        end if

      else

        IF (allocated(state_clm2pdaf_p)) deallocate(state_clm2pdaf_p)
        allocate(state_clm2pdaf_p(begc:endc,nlevsoi))

        do i=1,nlevsoi
          do c=clm_begc,clm_endc
            ! All columns in a gridcell are assigned the updated
            ! gridcell-SWC
            state_clm2pdaf_p(c,i) = (col%gridcell(c) - clm_begg + 1) + (i - 1)*(clm_endg - clm_begg + 1)
          end do
        end do

        ! One value per grid-cell
        clm_varsize      =  (endg-begg+1) * nlevsoi
        clm_statevecsize =  (endg-begg+1) * nlevsoi

        IF (allocated(state_pdaf2clm_c_p)) deallocate(state_pdaf2clm_c_p)
        allocate(state_pdaf2clm_c_p(clm_statevecsize))
        IF (allocated(state_pdaf2clm_j_p)) deallocate(state_pdaf2clm_j_p)
        allocate(state_pdaf2clm_j_p(clm_statevecsize))

        cc = 0

        do i=1,nlevsoi
          do j=clm_begg,clm_endg

            ! SWC from the first column of each gridcell
            newgridcell = .true.
            do jj=clm_begc,clm_endc
              g = col%gridcell(jj)
              if (g .eq. j) then
                if (newgridcell) then
                  newgridcell = .false.
                  ! Possibliy: Add state_pdaf2clm_g_p
                  state_pdaf2clm_c_p(cc) = jj
                  state_pdaf2clm_j_p(cc) = i
                end if
              end if
            end do

            cc = cc + 1
          end do
        end do



      end if
    endif

    if(clmupdate_swc.eq.2) then
      error stop "Not implemented: clmupdate_swc.eq.2"
    endif

    if(clmupdate_texture.eq.1) then
        clm_statevecsize = clm_statevecsize + 2*((endg-begg+1)*nlevsoi)
    endif

    if(clmupdate_texture.eq.2) then
        clm_statevecsize = clm_statevecsize + 3*((endg-begg+1)*nlevsoi)
    endif

    !hcp LST DA
    if(clmupdate_T.eq.1) then
      error stop "Not implemented: clmupdate_T.eq.1"
    endif
    !end hcp

#ifdef PDAF_DEBUG
    ! Debug output of clm_statevecsize
    WRITE(*, '(a,x,a,i5,x,a,i10)') "TSMP-PDAF-debug", "mype(w)=", mype, "define_clm_statevec: clm_statevecsize=", clm_statevecsize
#endif

    !write(*,*) 'clm_statevecsize is ',clm_statevecsize
    IF (allocated(clm_statevec)) deallocate(clm_statevec)
    if ((clmupdate_swc.ne.0) .or. (clmupdate_T.ne.0) .or. (clmupdate_texture.ne.0)) then
      !hcp added condition
      allocate(clm_statevec(clm_statevecsize))
    end if

    !write(*,*) 'clm_paramsize is ',clm_paramsize
    if (allocated(clm_paramarr)) deallocate(clm_paramarr)         !hcp
    if ((clmupdate_T.ne.0)) then  !hcp
      error stop "Not implemented clmupdate_T.NE.0"
    end if

  end subroutine define_clm_statevec

  subroutine cleanup_clm_statevec()

    implicit none

    ! Deallocate arrays from `define_clm_statevec`
    IF (allocated(clm_statevec)) deallocate(clm_statevec)
    IF (allocated(state_pdaf2clm_c_p)) deallocate(state_pdaf2clm_c_p)
    IF (allocated(state_pdaf2clm_j_p)) deallocate(state_pdaf2clm_j_p)
    IF (allocated(state_clm2pdaf_p)) deallocate(state_clm2pdaf_p)

  end subroutine cleanup_clm_statevec

  subroutine set_clm_statevec(tstartcycle, mype)
    use clm_instMod, only : soilstate_inst, waterstate_inst
    use clm_varpar   , only : nlevsoi
    ! use clm_varcon, only: nameg, namec
    ! use GetGlobalValuesMod, only: GetGlobalWrite
    use ColumnType , only : col
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    integer,intent(in) :: tstartcycle
    integer,intent(in) :: mype
    real(r8), pointer :: swc(:,:)
    real(r8), pointer :: psand(:,:)
    real(r8), pointer :: pclay(:,:)
    real(r8), pointer :: porgm(:,:)
    integer :: i,j,jj,g,cc=0,offset=0
    character (len = 34) :: fn    !TSMP-PDAF: function name for state vector output
    character (len = 34) :: fn2    !TSMP-PDAF: function name for swc output

    swc   => waterstate_inst%h2osoi_vol_col
    psand => soilstate_inst%cellsand_col
    pclay => soilstate_inst%cellclay_col
    porgm => soilstate_inst%cellorg_col

#ifdef PDAF_DEBUG
    IF(clmt_printensemble == tstartcycle + 1 .OR. clmt_printensemble < 0) THEN
      ! TSMP-PDAF: Debug output of CLM swc
      WRITE(fn2, "(a,i5.5,a,i5.5,a)") "swcstate_", mype, ".integrate.", tstartcycle + 1, ".txt"
      OPEN(unit=71, file=fn2, action="write")
      WRITE (71,"(es22.15)") swc(:,:)
      CLOSE(71)
    END IF
#endif

    ! calculate shift when CRP data are assimilated
    if(clmupdate_swc.eq.2) then
      error stop "Not implemented clmupdate_swc.eq.2"
    endif

    if(clmupdate_swc.ne.0) then
      ! write swc values to state vector
      do cc = 1, clm_statevecsize
        clm_statevec(cc) = swc(state_pdaf2clm_c_p(cc), state_pdaf2clm_j_p(cc))
      end do
    endif

    !hcp  LAI
    if(clmupdate_T.eq.1) then
      error stop "Not implemented: clmupdate_T.eq.1"
    endif
    !end hcp  LAI

    ! write average swc to state vector (CRP assimilation)
    if(clmupdate_swc.eq.2) then
      error stop "Not implemented: clmupdate_swc.eq.2"
    endif

    ! write texture values to state vector (if desired)
    if(clmupdate_texture.ne.0) then
      cc = 1
      do i=1,nlevsoi
        do j=clm_begg,clm_endg
          clm_statevec(cc+1*clm_varsize+offset) = psand(j,i)
          clm_statevec(cc+2*clm_varsize+offset) = pclay(j,i)
          if(clmupdate_texture.eq.2) then
            !incl. organic matter values
            clm_statevec(cc+3*clm_varsize+offset) = porgm(j,i)
          end if
          cc = cc + 1
        end do
      end do
    endif

#ifdef PDAF_DEBUG
    IF(clmt_printensemble == tstartcycle + 1 .OR. clmt_printensemble < 0) THEN
      ! TSMP-PDAF: For debug runs, output the state vector in files
      WRITE(fn, "(a,i5.5,a,i5.5,a)") "clmstate_", mype, ".integrate.", tstartcycle + 1, ".txt"
      OPEN(unit=71, file=fn, action="write")
      DO i = 1, clm_statevecsize
        WRITE (71,"(es22.15)") clm_statevec(i)
      END DO
      CLOSE(71)
    END IF
#endif

  end subroutine set_clm_statevec

  subroutine update_clm(tstartcycle, mype) bind(C,name="update_clm")
    use clm_varpar   , only : nlevsoi
    use clm_time_manager  , only : update_DA_nstep
    use shr_kind_mod , only : r8 => shr_kind_r8
    use ColumnType , only : col
    use clm_instMod, only : soilstate_inst, waterstate_inst
    use clm_varcon      , only : denh2o, denice, watmin
    use clm_varcon      , only : ispval
    use clm_varcon      , only : spval

    implicit none

    integer,intent(in) :: tstartcycle
    integer,intent(in) :: mype

    real(r8), pointer :: swc(:,:)
    real(r8), pointer :: watsat(:,:)
    real(r8), pointer :: psand(:,:)
    real(r8), pointer :: pclay(:,:)
    real(r8), pointer :: porgm(:,:)

    real(r8), pointer :: dz(:,:)          ! layer thickness depth (m)
    real(r8), pointer :: h2osoi_liq(:,:)  ! liquid water (kg/m2)
    real(r8), pointer :: h2osoi_ice(:,:)
    real(r8)  :: rliq,rice
    real(r8)  :: watmin_check      ! minimum soil moisture for checking clm_statevec (mm)
    real(r8)  :: watmin_set        ! minimum soil moisture for setting swc (mm)
    real(r8)  :: swc_update        ! updated SWC in loop

    integer :: i,j,jj,g,cc=0,offset=0
    character (len = 31) :: fn    !TSMP-PDAF: function name for state vector outpu
    character (len = 31) :: fn2    !TSMP-PDAF: function name for state vector outpu
    character (len = 32) :: fn3    !TSMP-PDAF: function name for state vector outpu
    character (len = 32) :: fn4    !TSMP-PDAF: function name for state vector outpu
    character (len = 32) :: fn5    !TSMP-PDAF: function name for state vector outpu
    character (len = 32) :: fn6    !TSMP-PDAF: function name for state vector outpu

    logical :: swc_zero_before_update = .false.

#ifdef PDAF_DEBUG
    IF(clmt_printensemble == tstartcycle .OR. clmt_printensemble < 0) THEN
      ! TSMP-PDAF: For debug runs, output the state vector in files
      WRITE(fn, "(a,i5.5,a,i5.5,a)") "clmstate_", mype, ".update.", tstartcycle, ".txt"
      OPEN(unit=71, file=fn, action="write")
      DO i = 1, clm_statevecsize
        WRITE (71,"(es22.15)") clm_statevec(i)
      END DO
      CLOSE(71)
    END IF
#endif

    swc   => waterstate_inst%h2osoi_vol_col
    watsat => soilstate_inst%watsat_col
    psand => soilstate_inst%cellsand_col
    pclay => soilstate_inst%cellclay_col
    porgm => soilstate_inst%cellorg_col

    dz            => col%dz
    h2osoi_liq    => waterstate_inst%h2osoi_liq_col
    h2osoi_ice    => waterstate_inst%h2osoi_ice_col

#ifdef PDAF_DEBUG
    IF(clmt_printensemble == tstartcycle .OR. clmt_printensemble < 0) THEN
      ! TSMP-PDAF: For debug runs, output the state vector in files
      WRITE(fn5, "(a,i5.5,a,i5.5,a)") "h2osoi_liq", mype, ".bef_up.", tstartcycle, ".txt"
      OPEN(unit=71, file=fn5, action="write")
      WRITE (71,"(es22.15)") h2osoi_liq(:,:)
      CLOSE(71)
    END IF
#endif
#ifdef PDAF_DEBUG
    IF(clmt_printensemble == tstartcycle .OR. clmt_printensemble < 0) THEN
      ! TSMP-PDAF: For debug runs, output the state vector in files
      WRITE(fn6, "(a,i5.5,a,i5.5,a)") "h2osoi_ice", mype, ".bef_up.", tstartcycle, ".txt"
      OPEN(unit=71, file=fn6, action="write")
      WRITE (71,"(es22.15)") h2osoi_ice(:,:)
      CLOSE(71)
    END IF
#endif

    ! calculate shift when CRP data are assimilated
    if(clmupdate_swc.eq.2) then
      error stop "Not implemented: clmupdate_swc.eq.2"
    endif

    ! CLM5: Update the Data Assimulation time-step to the current time
    ! step, since DA has been done. Used by CLM5 to skip BalanceChecks
    ! directly after the DA step.
    call update_DA_nstep()

    ! write updated swc back to CLM
    if(clmupdate_swc.ne.0) then

        ! Set minimum soil moisture for checking the state vector and
        ! for setting minimum swc for CLM
        if(clmwatmin_switch.eq.3) then
          ! CLM3.5 type watmin
          watmin_check = 0.00
          watmin_set = 0.05
        else if(clmwatmin_switch.eq.5) then
          ! CLM5.0 type watmin
          watmin_check = watmin
          watmin_set = watmin
        else
          ! Default
          watmin_check = 0.0
          watmin_set = 0.0
        end if

        ! cc = 0
        do i=1,nlevsoi
          ! CLM3.5: iterate over grid cells
          ! CLM5.0: iterate over columns
          ! do j=clm_begg,clm_endg
            do j=clm_begc,clm_endc

              ! Update only those SWCs that are not excluded by ispval
              if(state_clm2pdaf_p(j,i) .ne. ispval) then

                if(swc(j,i).eq.0.0) then
                  swc_zero_before_update = .true.

                  ! Zero-SWC leads to zero denominator in computation of
                  ! rliq/rice, therefore setting rliq/rice to special
                  ! value
                  rliq = spval
                  rice = spval
                else
                  swc_zero_before_update = .false.

                  rliq = h2osoi_liq(j,i)/(dz(j,i)*denh2o*swc(j,i))
                  rice = h2osoi_ice(j,i)/(dz(j,i)*denice*swc(j,i))
                  !h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)
                end if

                swc_update = clm_statevec(state_clm2pdaf_p(j,i))

                if(swc_update.le.watmin_check) then
                  swc(j,i) = watmin_set
                else if(swc_update.ge.watsat(j,i)) then
                  swc(j,i) = watsat(j,i)
                else
                  swc(j,i)   = swc_update
                endif

                if (isnan(swc(j,i))) then
                  swc(j,i) = watmin_set
                  print *, "WARNING: swc at j,i is nan: ", j, i
                endif

                if(swc_zero_before_update) then
                  ! This case should not appear for hydrologically
                  ! active columns/layers, where always: swc > watmin
                  !
                  ! If you want to make sure that no zero SWCs appear in
                  ! the code, comment out the error stop

#ifdef PDAF_DEBUG
                  ! error stop "ERROR: Update of zero-swc"
                  print *, "WARNING: Update of zero-swc"
                  print *, "WARNING: Any new H2O added to h2osoi_liq(j,i) with j,i = ", j, i
#endif
                  h2osoi_liq(j,i) = swc(j,i) * dz(j,i)*denh2o
                  h2osoi_ice(j,i) = 0.0
                else
                  ! update liquid water content
                  h2osoi_liq(j,i) = swc(j,i) * dz(j,i)*denh2o*rliq
                  ! update ice content
                  h2osoi_ice(j,i) = swc(j,i) * dz(j,i)*denice*rice
                end if

              end if
              ! cc = cc + 1
            end do
        end do

#ifdef PDAF_DEBUG
        IF(clmt_printensemble == tstartcycle .OR. clmt_printensemble < 0) THEN
          ! TSMP-PDAF: For debug runs, output the state vector in files
          WRITE(fn3, "(a,i5.5,a,i5.5,a)") "h2osoi_liq", mype, ".update.", tstartcycle, ".txt"
          OPEN(unit=71, file=fn3, action="write")
          WRITE (71,"(es22.15)") h2osoi_liq(:,:)
          CLOSE(71)
        END IF
#endif
#ifdef PDAF_DEBUG
        IF(clmt_printensemble == tstartcycle .OR. clmt_printensemble < 0) THEN
          ! TSMP-PDAF: For debug runs, output the state vector in files
          WRITE(fn4, "(a,i5.5,a,i5.5,a)") "h2osoi_ice", mype, ".update.", tstartcycle, ".txt"
          OPEN(unit=71, file=fn4, action="write")
          WRITE (71,"(es22.15)") h2osoi_ice(:,:)
          CLOSE(71)
        END IF
#endif
#ifdef PDAF_DEBUG
        IF(clmt_printensemble == tstartcycle .OR. clmt_printensemble < 0) THEN
          ! TSMP-PDAF: For debug runs, output the state vector in files
          WRITE(fn2, "(a,i5.5,a,i5.5,a)") "swcstate_", mype, ".update.", tstartcycle, ".txt"
          OPEN(unit=71, file=fn2, action="write")
          WRITE (71,"(es22.15)") swc(:,:)
          CLOSE(71)
        END IF
#endif

    endif

    !hcp: TG, TV
    if(clmupdate_T.EQ.1) then
      error stop "Not implemented: clmupdate_T.eq.1"
    endif
    ! end hcp TG, TV

    !! update liquid water content
    !do j=clm_begg,clm_endg
    !  do i=1,nlevsoi
    !    h2osoi_liq(j,i) = swc(j,i) * dz(j,i)*denh2o
    !  end do
    !end do

    ! write updated texture back to CLM
    if(clmupdate_texture.ne.0) then
      cc = 1
      do i=1,nlevsoi
        do j=clm_begg,clm_endg
          psand(j,i) = clm_statevec(cc+1*clm_varsize+offset)
          pclay(j,i) = clm_statevec(cc+2*clm_varsize+offset)
          if(clmupdate_texture.eq.2) then
            ! incl. organic matter
            porgm(j,i) = clm_statevec(cc+3*clm_varsize+offset)
          end if
          cc = cc + 1
        end do
      end do
      call clm_correct_texture
      call clm_texture_to_parameters
    endif

  end subroutine update_clm

  subroutine clm_correct_texture()

    use clm_varpar   , only : nlevsoi
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clm_instMod, only : soilstate_inst
    use CNSharedParamsMod, only : CNParamsShareInst

    implicit none

    integer :: c,lev
    real(r8) :: clay,sand,ttot
    real(r8) :: orgm
    real(r8), pointer :: psand(:,:)
    real(r8), pointer :: pclay(:,:)
    real(r8), pointer :: porgm(:,:)

    psand => soilstate_inst%cellsand_col
    pclay => soilstate_inst%cellclay_col
    porgm => soilstate_inst%cellorg_col

    do c = clm_begg,clm_endg
      do lev = 1,nlevsoi
         clay = pclay(c,lev)
         sand = psand(c,lev)

         if(sand.le.0.0) sand = 1.0
         if(clay.le.0.0) clay = 1.0

         ttot = sand + clay
         if(ttot.gt.100) then
             sand = sand/ttot * 100.0
             clay = clay/ttot * 100.0
         end if

         pclay(c,lev) = clay
         psand(c,lev) = sand

         if(clmupdate_texture.eq.2) then
           orgm = (porgm(c,lev) / CNParamsShareInst%organic_max) * 100.0
           if(orgm.le.0.0) orgm = 0.0
           if(orgm.ge.100.0) orgm = 100.0
           porgm(c,lev) = (orgm / 100.0)* CNParamsShareInst%organic_max
         end if

       end do
    end do
  end subroutine clm_correct_texture

  subroutine clm_texture_to_parameters()
    use clm_varpar   , only : nlevsoi
    use clm_varcon   , only : zsoi, secspday
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clm_instMod, only : soilstate_inst
    use FuncPedotransferMod , only : pedotransf, get_ipedof
    use CNSharedParamsMod, only : CNParamsShareInst
    implicit none

    real(r8)           :: om_tkm         = 0.25_r8      
    ! thermal conductivity of organic soil (Farouki, 1986) [W/m/K]
    real(r8)           :: om_watsat_lake = 0.9_r8       
    ! porosity of organic soil
    real(r8)           :: om_hksat_lake  = 0.1_r8       
    ! saturated hydraulic conductivity of organic soil [mm/s]
    real(r8)           :: om_sucsat_lake = 10.3_r8      
    ! saturated suction for organic matter (Letts, 2000)
    real(r8)           :: om_b_lake      = 2.7_r8       
    ! Clapp Hornberger paramater for oragnic soil (Letts, 2000) (lake)
    real(r8)           :: om_watsat                     
    ! porosity of organic soil
    real(r8)           :: om_hksat                      
    ! saturated hydraulic conductivity of organic soil [mm/s]
    real(r8)           :: om_sucsat                     
    ! saturated suction for organic matter (mm)(Letts, 2000)
    real(r8)           :: om_csol        = 2.5_r8       
    ! heat capacity of peat soil *10^6 (J/K m3) (Farouki, 1986)
    real(r8)           :: om_tkd         = 0.05_r8      
    ! thermal conductivity of dry organic soil (Farouki, 1981)
    real(r8)           :: om_b                          
    ! Clapp Hornberger paramater for oragnic soil (Letts, 2000)
    real(r8)           :: zsapric        = 0.5_r8
    ! depth (m) that organic matter takes on characteristics of sapric peat
    real(r8)           :: pcalpha        = 0.5_r8       ! percolation threshold
    real(r8)           :: pcbeta         = 0.139_r8     ! percolation exponent
    real(r8)           :: perc_frac    ! "percolating" fraction of organic soil
    real(r8)           :: perc_norm    ! normalize to 1 when 100% organic soil
    real(r8)           :: uncon_hksat  ! series conductivity of mineral/organic soil
    real(r8)           :: uncon_frac   ! fraction of "unconnected" soil
    real(r8)           :: bd           ! bulk density of dry soil material [kg/m^3]
    real(r8)           :: tkm          ! mineral conductivity 
    real(r8)           :: xksat        ! maximum hydraulic conductivity of soil [mm/s]
    
    integer  :: ipedof,c,lev
    real(r8) :: clay,sand,om_frac
    real(r8), pointer :: psand(:,:)
    real(r8), pointer :: pclay(:,:)
    real(r8), pointer :: porgm(:,:)

    psand => soilstate_inst%cellsand_col
    pclay => soilstate_inst%cellclay_col
    porgm   => soilstate_inst%cellorg_col

    do c = clm_begg,clm_endg
      do lev = 1,nlevsoi
         clay = pclay(c,lev)
         sand = psand(c,lev)
         om_frac = porgm(c,lev) / CNParamsShareInst%organic_max

         ipedof=get_ipedof(0)
         call pedotransf(ipedof, sand, clay, &
                         soilstate_inst%watsat_col(c,lev), &
                         soilstate_inst%bsw_col(c,lev), &
                         soilstate_inst%sucsat_col(c,lev), &
                         xksat)

         om_watsat         = max(0.93_r8 - 0.1_r8   *(zsoi(lev)/zsapric), 0.83_r8)
         om_b              = min(2.7_r8  + 9.3_r8   *(zsoi(lev)/zsapric), 12.0_r8)
         om_sucsat         = min(10.3_r8 - 0.2_r8   *(zsoi(lev)/zsapric), 10.1_r8)
         om_hksat          = max(0.28_r8 - 0.2799_r8*(zsoi(lev)/zsapric), xksat)

         soilstate_inst%bd_col(c,lev) = &
              (1._r8 - soilstate_inst%watsat_col(c,lev))*2.7e3_r8
         
         soilstate_inst%watsat_col(c,lev) = &
              (1._r8 - om_frac) * soilstate_inst%watsat_col(c,lev) + om_watsat*om_frac
         
         tkm = (1._r8 - om_frac) * (8.80_r8*sand+2.92_r8*clay)/(sand+clay)+om_tkm*om_frac 
         ! W/(m K)
         soilstate_inst%bsw_col(c,lev) = &
              (1._r8-om_frac) * (2.91_r8 + 0.159_r8*clay) + om_frac*om_b   
         
         soilstate_inst%sucsat_col(c,lev) = &
              (1._r8-om_frac) * soilstate_inst%sucsat_col(c,lev) + om_sucsat*om_frac  
         
         soilstate_inst%hksat_min_col(c,lev) = xksat

         ! perc_frac is zero unless perf_frac greater than percolation threshold
         if (om_frac > pcalpha) then
            perc_norm = (1._r8 - pcalpha)**(-pcbeta)
            perc_frac = perc_norm*(om_frac - pcalpha)**pcbeta
         else
            perc_frac = 0._r8
         endif

         ! uncon_frac is fraction of mineral soil plus fraction of
         ! "nonpercolating" organic soil
         uncon_frac = (1._r8-om_frac)+(1._r8-perc_frac)*om_frac

         ! uncon_hksat is series addition of mineral/organic
         ! conductivites
         if (om_frac < 1._r8) then
            uncon_hksat = uncon_frac/((1._r8-om_frac)/xksat &
                 +((1._r8-perc_frac)*om_frac)/om_hksat)
         else
            uncon_hksat = 0._r8
         end if

         soilstate_inst%hksat_col(c,lev)  = uncon_frac*uncon_hksat + &
             (perc_frac*om_frac)*om_hksat

         soilstate_inst%tkmg_col(c,lev)   = tkm ** (1._r8 - &
             soilstate_inst%watsat_col(c,lev))           

         soilstate_inst%tksatu_col(c,lev) = &
             soilstate_inst%tkmg_col(c,lev)*0.57_r8**soilstate_inst%watsat_col(c,lev)

         soilstate_inst%tkdry_col(c,lev)  = &
             ((0.135_r8*soilstate_inst%bd_col(c,lev) + 64.7_r8) / &
             (2.7e3_r8 - 0.947_r8*soilstate_inst%bd_col(c,lev)) &
             )*(1._r8-om_frac) + om_tkd*om_frac  

         soilstate_inst%csol_col(c,lev)   = &
             ((1._r8-om_frac)*(2.128_r8*sand+2.385_r8*clay) / (sand+clay) + &
             om_csol*om_frac)*1.e6_r8

         soilstate_inst%watdry_col(c,lev) = &
             soilstate_inst%watsat_col(c,lev) * &
             (316230._r8/soilstate_inst%sucsat_col(c,lev)  &
             ) ** (-1._r8/soilstate_inst%bsw_col(c,lev)) 

         soilstate_inst%watopt_col(c,lev) = & 
             soilstate_inst%watsat_col(c,lev) * &
             (158490._r8/soilstate_inst%sucsat_col(c,lev)  &
             ) ** (-1._r8/soilstate_inst%bsw_col(c,lev)) 

         ! secspday (day/sec)
         soilstate_inst%watfc_col(c,lev) = &
             soilstate_inst%watsat_col(c,lev) * &
             (0.1_r8 / (soilstate_inst%hksat_col(c,lev)*secspday) &
             )**(1._r8/(2._r8*soilstate_inst%bsw_col(c,lev)+3._r8))

      end do
    end do
  end subroutine clm_texture_to_parameters

  ! subroutine  average_swc_crp(profdat,profave)
  !   use clm_varcon  , only : zsoi

  !   implicit none

  !   real(r8),intent(in)  :: profdat(10)
  !   real(r8),intent(out) :: profave

  !   error stop "Not implemented average_swc_crp"
  ! end subroutine average_swc_crp
#endif

  subroutine domain_def_clm(lon_clmobs, lat_clmobs, dim_obs, &
                            longxy, latixy, longxy_obs, latixy_obs)
    use spmdMod,   only : npes, iam
    use domainMod, only : ldomain
    use decompMod, only : get_proc_total, get_proc_bounds, ldecomp

    implicit none
    real, intent(in) :: lon_clmobs(:)
    real, intent(in) :: lat_clmobs(:)
    integer, intent(in) :: dim_obs
    integer, allocatable, intent(inout) :: longxy(:)
    integer, allocatable, intent(inout) :: latixy(:)
    integer, allocatable, intent(inout) :: longxy_obs(:)
    integer, allocatable, intent(inout) :: latixy_obs(:)
    integer :: ni, nj, ii, jj, kk, cid, ier, ncells, nlunits
    integer :: ncols, counter
    integer :: npatches, ncohorts
    real :: minlon, minlat, maxlon, maxlat
    real(r8), pointer :: lon(:)
    real(r8), pointer :: lat(:)
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices


    lon => ldomain%lonc
    lat => ldomain%latc
    ni = ldomain%ni
    nj = ldomain%nj

    ! get total number of gridcells, landunits,
    ! columns, patches and cohorts on processor

    call get_proc_total(iam, ncells, nlunits, ncols, npatches, ncohorts)

    ! beg and end gridcell
    call get_proc_bounds(begg=begg, endg=endg)
    
   !print *,'ni, nj ', ni, nj
   !print *,'cells per processor ', ncells
   !print *,'begg, endg ', begg, endg

    ! allocate vector with size of elements in x directions * size of elements in y directions
    allocate(longxy(ncells), stat=ier)
    allocate(latixy(ncells), stat=ier)

    ! initialize vector with zero values
    longxy(:) = 0
    latixy(:) = 0
  
    ! fill vector with index values
    counter = 1
    do ii = 1, nj
      do jj = 1, ni
        cid = (ii-1)*ni + jj
        do kk = begg, endg
          if(cid == ldecomp%gdc2glo(kk)) then
            latixy(counter) = ii
            longxy(counter) = jj
            counter = counter + 1
          end if
        end do
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

    allocate(longxy_obs(dim_obs), stat=ier)
    allocate(latixy_obs(dim_obs), stat=ier)

    do i = 1, dim_obs
       if(((lon_clmobs(i) + 180) - minlon) /= 0 .and. &
         ((lat_clmobs(i) + 90) - minlat) /= 0) then
          longxy_obs(i) = ceiling(((lon_clmobs(i) + 180) - minlon) * ni / (maxlon - minlon)) !+ 1
          latixy_obs(i) = ceiling(((lat_clmobs(i) + 90) - minlat) * nj / (maxlat - minlat)) !+ 1
          !print *,'longxy_obs(i) , latixy_obs(i) ', longxy_obs(i) , latixy_obs(i)
        else if(((lon_clmobs(i) + 180) - minlon) == 0 .and. &
                ((lat_clmobs(i) + 90) - minlat) == 0) then
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

  end subroutine domain_def_clm

  !> @author  Mukund Pondkule, Johannes Keller
  !> @date    27.03.2023
  !> @brief   Set indices of grid cells with lon/lat smaller than observation locations
  !> @details
  !>    This routine sets the indices of grid cells with lon/lat
  !>    smaller than observation locations.
  subroutine get_interp_idx(lon_clmobs, lat_clmobs, dim_obs, longxy_obs_floor, latixy_obs_floor)

    USE domainMod, ONLY: ldomain
    ! USE decompMod, ONLY: get_proc_total, get_proc_bounds_atm, adecomp
    ! USE spmdMod,   ONLY: npes, iam ! number of processors for clm and processor number
   ! USE mod_read_obs, ONLY: lon_clmobs, lat_clmobs
    ! USE clmtype,    ONLY : clm3
   ! USE mod_assimilation, ONLY: dim_obs
!#endif

    implicit none
    real, intent(in) :: lon_clmobs(:)
    real, intent(in) :: lat_clmobs(:)
    integer, intent(in) :: dim_obs
    integer, allocatable, intent(inout) :: longxy_obs_floor(:)
    integer, allocatable, intent(inout) :: latixy_obs_floor(:)
    integer :: i

    integer :: ni, nj
    ! integer :: ii, jj, kk, cid
    integer :: ier
    ! integer :: ncells, nlunits,ncols, npfts
    ! integer :: counter

    real :: minlon, minlat, maxlon, maxlat
    real(r8), pointer :: lon(:)
    real(r8), pointer :: lat(:)
    ! integer :: begg, endg   ! per-proc gridcell ending gridcell indices

    lon => ldomain%lonc
    lat => ldomain%latc
    ni = ldomain%ni
    nj = ldomain%nj

    ! ! get total number of gridcells, landunits,
    ! ! columns and pfts for any processor
    ! call get_proc_total(iam, ncells, nlunits, ncols, npfts)

    ! ! beg and end gridcell for atm
    ! call get_proc_bounds_atm(begg, endg)

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

    if(allocated(longxy_obs_floor)) deallocate(longxy_obs_floor)
    allocate(longxy_obs_floor(dim_obs), stat=ier)
    if(allocated(latixy_obs_floor)) deallocate(latixy_obs_floor)
    allocate(latixy_obs_floor(dim_obs), stat=ier)
    do i = 1, dim_obs
       if(((lon_clmobs(i) + 180) - minlon) /= 0 .and. ((lat_clmobs(i) + 90) - minlat) /= 0) then
          longxy_obs_floor(i) = floor(((lon_clmobs(i) + 180) - minlon) * ni / (maxlon - minlon)) !+ 1
          latixy_obs_floor(i) = floor(((lat_clmobs(i) + 90) - minlat) * nj / (maxlat - minlat)) !+ 1
          !print *,'longxy_obs(i) , latixy_obs(i) ', longxy_obs(i) , latixy_obs(i)
       else if(((lon_clmobs(i) + 180) - minlon) == 0 .and. ((lat_clmobs(i) + 90) - minlat) == 0) then
          longxy_obs_floor(i) = 1
          latixy_obs_floor(i) = 1
       else if(((lon_clmobs(i) + 180) - minlon) == 0) then
          longxy_obs_floor(i) = 1
          latixy_obs_floor(i) = floor(((lat_clmobs(i) + 90) - minlat) * nj / (maxlat - minlat))
       else if(((lat_clmobs(i) + 90) - minlat) == 0) then
          longxy_obs_floor(i) = floor(((lon_clmobs(i) + 180) - minlon) * ni / (maxlon - minlon))
          latixy_obs_floor(i) = 1
       endif
    end do

  end subroutine get_interp_idx

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
      error stop "Not implemented: clmupdate_swc.eq.2"
      ! dim_l = nlevsoi + 1
      ! nshift = nlevsoi + 1
    endif

    if(clmupdate_texture.eq.1) then
      dim_l = 2*nlevsoi + nshift
    endif

    if(clmupdate_texture.eq.2) then
      dim_l = 3*nlevsoi + nshift
    endif

  end subroutine init_clm_l_size
#endif

end module enkf_clm_mod

