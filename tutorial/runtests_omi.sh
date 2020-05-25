#!/bin/tcsh

# ARCH and ARCH_MPI specify PDAF_ARCH without and with PDAF
setenv ARCH osx_gfortran_openmpi #linux_gfortran
setenv ARCH_MPI osx_gfortran_openmpi
setenv DA_SPECS "-filtertype 7"
setenv DA_SPECS2 "-filtertype 6"

echo "------------------ COMPILING ----------------"

setenv PDAF_ARCH $ARCH
echo -------------- PDAF_ARCH: $PDAF_ARCH

echo "------------ online_2D_serialmodel_omi ---------------"
setenv PDAF_ARCH $ARCH_MPI
echo PDAF_ARCH: $PDAF_ARCH
cd online_2D_serialmodel_omi
make clean
make cleandata
make model
make model_pdaf
cd ..

echo "------------ online_2D_serialmodel_openmp_omi ---------------"
setenv PDAF_ARCH $ARCH_MPI
echo PDAF_ARCH: $PDAF_ARCH
cd online_2D_serialmodel_openmp_omi
make clean
make cleandata
make model
make model_pdaf
cd ..

echo "------------ online_2D_parallelmodel_omi ---------------"
setenv PDAF_ARCH $ARCH_MPI
echo PDAF_ARCH: $PDAF_ARCH
cd online_2D_parallelmodel_omi
make clean
make cleandata
make model
make model_pdaf
cd ..


echo  " "
echo "-------------------- RUNNING ----------------"


echo "------------ online_2D_serialmodel_omi LESTKF ---------------"
setenv OMP_NUM_THREADS 1
cd online_2D_serialmodel_omi
make cleandata
mpirun --oversubscribe -np 9 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_serialmodel_omi
cd ..
python verification/check_online2.py online_2D_serialmodel_omi online_2D_serialmodel

echo "------------ online_2D_serialmodel_openmp_omi LESTKF ---------------"
setenv OMP_NUM_THREADS 2
cd online_2D_serialmodel_openmp_omi
make cleandata
mpirun --oversubscribe -np 9 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_serialmodel_openmp_omi
cd ..
python verification/check_online2.py online_2D_serialmodel_openmp_omi online_2D_serialmodel_openmp

echo "------------ online_2D_parallelmodel_omi LESTKF ---------------"
setenv OMP_NUM_THREADS 1
cd online_2D_parallelmodel_omi
make cleandata
mpirun --oversubscribe -np 18 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_parallelmodel_omi
cd ..
python verification/check_online2.py online_2D_parallelmodel_omi online_2D_parallelmodel


echo "------------ online_2D_serialmodel_omi ESTKF ---------------"
setenv OMP_NUM_THREADS 1
cd online_2D_serialmodel_omi
make cleandata
mpirun --oversubscribe -np 9 ./model_pdaf -dim_ens 9 $DA_SPECS2 > ../out.online_2D_serialmodel_omi_ESTKF
cd ..
python verification/check_online2.py online_2D_serialmodel_omi online_2D_serialmodel_ESTKF

echo "------------ online_2D_parallelmodel_omi ESTKF ---------------"
setenv OMP_NUM_THREADS 1
cd online_2D_parallelmodel_omi
make cleandata
mpirun --oversubscribe -np 18 ./model_pdaf -dim_ens 9 $DA_SPECS2 > ../out.online_2D_parallelmodel_omi_ESTKF
cd ..
python verification/check_online2.py online_2D_parallelmodel_omi online_2D_parallelmodel_ESTKF
