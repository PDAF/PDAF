#!/bin/tcsh

# ARCH and ARCH_MPI specify PDAF_ARCH without and with PDAF
setenv ARCH osx_gfortran_openmpi #linux_gfortran
setenv ARCH_MPI osx_gfortran_openmpi
setenv DA_SPECS "-filtertype 7"
setenv DA_SPECS2 "-filtertype 6"

echo "------------------ COMPILING ----------------"

setenv PDAF_ARCH $ARCH
echo -------------- PDAF_ARCH: $PDAF_ARCH

echo "------------ offline_2D_serial --------------"
cd offline_2D_serial
make cleanall
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

echo "------------ online_2D_parallelmodel ---------------"
setenv PDAF_ARCH $ARCH_MPI
echo PDAF_ARCH: $PDAF_ARCH
cd online_2D_parallelmodel
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
python verification/check_offline.py offline_2D_serial

echo "------------ offline_2D_openmp ---------------"
setenv OMP_NUM_THREADS 4
cd offline_2D_serial
./PDAF_offline $DA_SPECS > ../out.offline_2D_openmp
cd ..
python verification/check_offline.py offline_2D_serial

echo "------------ offline_2D_parallel ---------------"
setenv OMP_NUM_THREADS 1
cd offline_2D_parallel
mpirun -np 4 ./PDAF_offline $DA_SPECS > ../out.offline_2D_parallel
cd ..
python verification/check_offline.py offline_2D_parallel

echo "------------ online_2D_serialmodel LESTKF ---------------"
setenv OMP_NUM_THREADS 1
cd online_2D_serialmodel
make cleandata
mpirun --oversubscribe -np 9 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_serialmodel
cd ..
python verification/check_online2.py online_2D_serialmodel online_2D_serialmodel

echo "------------ online_2D_parallelmodel LESTKF ---------------"
setenv OMP_NUM_THREADS 1
cd online_2D_parallelmodel
make cleandata
mpirun --oversubscribe -np 18 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_parallelmodel
cd ..
python verification/check_online2.py online_2D_parallelmodel online_2D_parallelmodel


echo "------------ online_2D_serialmodel ESTKF ---------------"
setenv OMP_NUM_THREADS 1
cd online_2D_serialmodel
make cleandata
mpirun --oversubscribe -np 9 ./model_pdaf -dim_ens 9 $DA_SPECS2 > ../out.online_2D_serialmodel_ESTKF
cd ..
python verification/check_online2.py online_2D_serialmodel online_2D_serialmodel_ESTKF

echo "------------ online_2D_parallelmodel ESTKF ---------------"
setenv OMP_NUM_THREADS 1
cd online_2D_parallelmodel
make cleandata
mpirun --oversubscribe -np 18 ./model_pdaf -dim_ens 9 $DA_SPECS2 > ../out.online_2D_parallelmodel_ESTKF
cd ..
python verification/check_online2.py online_2D_parallelmodel online_2D_parallelmodel_ESTKF
