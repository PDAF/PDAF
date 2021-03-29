#!/bin/tcsh

# ARCH and ARCH_MPI specify PDAF_ARCH without and with PDAF
setenv ARCH linux_gfortran
setenv ARCH_MPI linux_gfortran_openmpi
setenv DA_SPECS "-filtertype 7"
setenv DA_SPECS2 "-filtertype 6"
setenv DA_SPECS3 "-filtertype 7 -assim_A .false. -assim_B .true"

echo "------------------ COMPILING ----------------"

setenv PDAF_ARCH $ARCH
echo -------------- PDAF_ARCH: $PDAF_ARCH

echo "------------ online_2D_serialmodel ---------------"
setenv PDAF_ARCH $ARCH_MPI
echo PDAF_ARCH: $PDAF_ARCH
cd online_2D_serialmodel
make clean
make cleandata
make model
make model_pdaf
cd ..

echo  " "
echo "-------------------- RUNNING ----------------"


echo "------------ online_2D_serialmodel LESTKF ---------------"
setenv OMP_NUM_THREADS 1
cd online_2D_serialmodel
make cleandata
mpirun --oversubscribe -np 9 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_serialmodel
cd ..
python verification/check_online2.py online_2D_serialmodel online_2D_serialmodel

echo "------------ online_2D_serialmodel_openmp LESTKF ---------------"
setenv OMP_NUM_THREADS 2
cd online_2D_serialmodel
make cleandata
mpirun --oversubscribe -np 9 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_serialmodel_openmp
cd ..
python verification/check_online.py online_2D_serialmodel

echo "------------ online_2D_serialmodel ESTKF ---------------"
setenv OMP_NUM_THREADS 1
cd online_2D_serialmodel
make cleandata
mpirun --oversubscribe -np 9 ./model_pdaf -dim_ens 9 $DA_SPECS2 > ../out.online_2D_serialmodel_ESTKF
cd ..
python verification/check_online2.py online_2D_serialmodel online_2D_serialmodel_ESTKF

