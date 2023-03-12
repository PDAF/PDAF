#!/bin/tcsh

# ARCH and ARCH_MPI specify PDAF_ARCH without and with PDAF
setenv ARCH osx_gfortran
setenv ARCH_MPI osx_gfortran_openmpi
setenv DA_SPECS "-filtertype 7"

echo "------------------ COMPILING ----------------"

setenv PDAF_ARCH $ARCH
echo -------------- PDAF_ARCH: $PDAF_ARCH

echo "------------ offline_2D_serial --------------"
cd offline_2D_serial
make cleanall
make
cd ..

echo "------------ offline_2D_openmp ---------------"
cd offline_2D_openmp
make clean
make cleandata
make
cd ..

echo "------------ offline_2D_parallel ---------------"
setenv PDAF_ARCH $ARCH_MPI
echo PDAF_ARCH: $PDAF_ARCH
cd offline_2D_parallel
make cleanall
make
cd ..

echo "------------ online_2D_serialmodel ---------------"
setenv PDAF_ARCH $ARCH_MPI
echo PDAF_ARCH: $PDAF_ARCH
cd online_2D_serialmodel
make clean
make cleandata
make model
make model_pdaf
cd ..

echo "------------ online_2D_serialmodel_openmp ---------------"
setenv PDAF_ARCH $ARCH_MPI
echo PDAF_ARCH: $PDAF_ARCH
cd online_2D_serialmodel_openmp
make clean
make cleandata
make model
make model_pdaf
cd ..

echo "------------ online_2D_parallelmodel ---------------"
setenv PDAF_ARCH $ARCH_MPI
echo PDAF_ARCH: $PDAF_ARCH
cd online_2D_parallelmodel
make clean
make cleandata
make model
make model_pdaf
cd ..

echo "------------ online_2D_parallelmodel_fullpar ---------------"
setenv PDAF_ARCH $ARCH_MPI
echo PDAF_ARCH: $PDAF_ARCH
cd online_2D_parallelmodel_fullpar
make clean
make cleandata
make model
make model_pdaf
cd ..

echo "------------ online_2D_parallelmodel_fullpar_1fpe ---------------"
setenv PDAF_ARCH $ARCH_MPI
echo PDAF_ARCH: $PDAF_ARCH
cd online_2D_parallelmodel_fullpar_1fpe
make clean
make cleandata
make model
make model_pdaf
cd ..


echo  " "
echo "-------------------- RUNNING ----------------"

echo "------------ offline_2D_serial --------------"
setenv OMP_NUM_THREADS 1
cd offline_2D_serial
./PDAF_offline $DA_SPECS > ../out.offline_2D_serial
cd ..

echo "------------ offline_2D_openmp ---------------"
setenv OMP_NUM_THREADS 4
cd offline_2D_openmp
./PDAF_offline $DA_SPECS > ../out.offline_2D_openmp
cd ..

echo "------------ offline_2D_parallel ---------------"
setenv OMP_NUM_THREADS 1
cd offline_2D_parallel
mpirun -np 4 ./PDAF_offline $DA_SPECS > ../out.offline_2D_parallel
cd ..


echo "------------ online_2D_serialmodel ---------------"
setenv OMP_NUM_THREADS 1
cd online_2D_serialmodel
mpirun -np 9 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_serialmodel
cd ..

echo "------------ online_2D_serialmodel_openmp ---------------"
setenv OMP_NUM_THREADS 2
cd online_2D_serialmodel_openmp
mpirun -np 9 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_serialmodel_openmp
cd ..

echo "------------ online_2D_parallelmodel ---------------"
setenv OMP_NUM_THREADS 1
cd online_2D_parallelmodel
mpirun -np 18 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_parallelmodel
cd ..

echo "------------ online_2D_parallelmodel_fullpar ---------------"
setenv OMP_NUM_THREADS 1
cd online_2D_parallelmodel_fullpar
mpirun -np 20 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_parallelmodel_fullpar
cd ..

echo "------------ online_2D_parallelmodel ---------------"
setenv OMP_NUM_THREADS 1
cd online_2D_parallelmodel_fullpar_1fpe
mpirun -np 19 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_parallelmodel_fullpar_1fpe
cd ..

