! !Program: rewrite_aux3d --- Rewrite mesh file aux3d.out into netcdf
!
! !INTERFACE:

program rewrite_aux3d

! This program reads the text file aux3d and rewrites it in a NetCDF format.

  implicit none
#include "netcdf.inc"

  character(len=160) :: meshpath
  character(len=120) :: outpath, outfile, ncfile_out
  character(len=150) :: attstr
  integer :: i, node, k, s
  integer :: nod2d, nod3d, nlayers
  integer, allocatable :: aux3d(:,:)
  integer, allocatable :: nlayer(:)
  integer, allocatable :: nod2d_above_nod3d(:)
  integer :: stat(100)
  integer :: dimids(2)
  integer :: ncid_out, dimid_n2d, dimid_n3d, dimid_nlayers
  integer :: id_aux3d, id_nlayer, id_n3d_to_n2d



  ! Path to the fesom mesh
  meshpath = '/scratch/usr/hbkqtang/old_work2/input/CORE2_final/'

  outpath = '/scratch/usr/hzfblner/tmpSST/'
  outpath = './'
  outfile = 'aux3d.nc'


  ! Read mesh information

  open (20,file=trim(MeshPath)//'nod2d.out',  status='old')
  read(20,*) nod2D 
  close(20)

  open (20,file=trim(MeshPath)//'nod3d.out',  status='old')
  read(20,*) nod3D 
  close(20)

  open (21,file=trim(MeshPath)//'aux3d.out',  status='old')
  read(21,*) nlayers

  write(*,*) 'Mesh dimensions: nod2d, nlayers', nod2d, nlayers

  ! *** Read aux3d into array
  allocate(aux3d(nlayers,nod2d))
  do node=1, nod2d
     do k=1, nlayers
        read(21,*) aux3d(k, node)
     end do
  end do
  close(21)
     
  ! *** Initilize array holding number of layers
  allocate(nlayer(nod2d))
  nlayer=0

  do node=1, nod2d
     layerloop: do k=1, nlayers
        if(aux3d(k, node)>0) then
           nlayer(node) = k
        else
           exit layerloop
        end if
     end do layerloop
  end do


  ! *** Initialize array pointing from node3d to node2d
  allocate(nod2d_above_nod3d(nod3d))
  do node=1, nod2d
     do k=1, nlayer(node)
        nod2d_above_nod3d(aux3d(k,node)) = aux3d(1,node)
     end do
  end do


  ! *** Write aux3d information into netcdf

  ncfile_out = trim(outpath)//trim(outfile)
  write (*,*) 'Write aux3d information to file: ',trim(ncfile_out)

  s = 1
  stat(s) = NF_CREATE(ncfile_out, 0, ncid_out)

  attstr  = 'aux3d mesh information in NetCDF format'
  s = s + 1
  stat(s) = NF_PUT_ATT_TEXT(ncid_out, NF_GLOBAL, 'title', len_trim(attstr), &
       trim(attstr))

  ! Define dimensions
  s = s + 1
  stat(s) = NF_DEF_DIM(ncid_out, 'nodes_2D', nod2D, dimid_n2d)
  s = s + 1
  stat(s) = NF_DEF_DIM(ncid_out, 'nodes_3D', nod3D, dimid_n3d)
  s = s + 1
  stat(s) = NF_DEF_DIM(ncid_out, 'nlayers', nlayers, dimid_nlayers)

  dimids(1) = dimid_nlayers
  dimids(2) = dimid_n2D
  s = s + 1
  stat(s) = NF_DEF_VAR(ncid_out, 'aux3d', NF_INT, 2, dimids, id_aux3d)

  dimids(1) = dimid_n2D
  s = s + 1
  stat(s) = NF_DEF_VAR(ncid_out, 'nlayer', NF_INT, 1, dimids(1), id_nlayer)

  dimids(1) = dimid_n3D
  s = s + 1
  stat(s) = NF_DEF_VAR(ncid_out, 'nod2d_above_nod3d', NF_INT, 1, dimids(1), id_n3d_to_n2d)

  s = s + 1
  stat(s) = NF_ENDDEF(ncid_out)

  ! Write data

  s = s + 1
  stat(s) = NF_PUT_VAR_INT(ncid_out, id_aux3d, aux3d)
  s = s+1
  stat(s) = NF_PUT_VAR_INT(ncid_out, id_nlayer, nlayer)
  s = s+1
  stat(s) = NF_PUT_VAR_INT(ncid_out, id_n3d_to_n2d, nod2d_above_nod3d)


  ! Close file
  s = s+1
  stat(s) = NF_CLOSE(ncid_out)

  do i = 1,  s
     if (stat(i) /= NF_NOERR) &
          write(*, *) 'NetCDF error in writing of output file, no.', i
  end do

  open (20,file='nod2d_above_nod3d.out',  status='replace')
  do node = 1, nod3d
     write (20, *) nod2d_above_nod3d(node)
  end do
  close(20)

  deallocate(aux3d, nlayer, nod2d_above_nod3d)

end program rewrite_aux3d
