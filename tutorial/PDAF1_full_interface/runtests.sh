#!/bin/bash

# ARCH and ARCH_MPI specify PDAF_ARCH without and with PDAF
export ARCH=linux_gfortran_openmpi
export DA_SPECS="-filtertype 7 -screen 1"
export DA_SPECS2="-filtertype 6 -screen 1"

echo "------------------ COMPILING ----------------"

COMPILEPDAF=0
COMPILE=1

if [ $COMPILE -eq 1 ]
then

    export PDAF_ARCH=$ARCH
    echo -------------- PDAF_ARCH: $PDAF_ARCH

    echo "------------ offline_2D_serial --------------"
    cd offline_2D_serial
    if [ $COMPILEPDAF -eq 1]
    then
        make cleanall
    else
        make clean
        make cleandataq
    fi
    make
    cd ..

    echo "------------ offline_2D_parallel ---------------"
    cd offline_2D_parallel
    make clean
    make cleandata
    make
    cd ..

    echo "------------ online_2D_serialmodel ---------------"
    export PDAF_ARCH=$ARCH
    echo PDAF_ARCH: $PDAF_ARCH
    cd online_2D_serialmodel
    make clean
    make cleandata
    make model
    make model_pdaf
    cd ..

    echo "------------ online_2D_parallelmodel ---------------"
    export PDAF_ARCH=$ARCH
    echo PDAF_ARCH: $PDAF_ARCH
    cd online_2D_parallelmodel
    make clean
    make cleandata
    make model
    make model_pdaf
    cd ..

else
    echo "Compilation is deactivated!"
fi


echo  " "
echo "-------------------- RUNNING ----------------"


echo "------------ offline_2D_serial ----------------------------------"
export OMP_NUM_THREADS=1
cd offline_2D_serial
make cleandataq
./PDAF_offline $DA_SPECS > ../out.offline_2D_serial
cd ..
python ../verification/check_offline_cl.py offline_2D_serial

echo "------------ offline_2D_openmp ----------------------------------"
export OMP_NUM_THREADS=4
cd offline_2D_serial
make cleandataq
./PDAF_offline $DA_SPECS > ../out.offline_2D_openmp
cd ..
python ../verification/check_offline_cl.py offline_2D_serial

echo "------------ offline_2D_parallel --------------------------------"
export OMP_NUM_THREADS=1
cd offline_2D_parallel
make cleandataq
mpirun -np 4 ./PDAF_offline $DA_SPECS > ../out.offline_2D_parallel
cd ..
python ../verification/check_offline_cl.py offline_2D_parallel

echo "------------ online_2D_serialmodel LESTKF -----------------------"
export OMP_NUM_THREADS=1
cd online_2D_serialmodel
make cleandataq
mpirun --oversubscribe -np 9 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_serialmodel
cd ..
python ../verification/check_online_cl.py online_2D_serialmodel

echo "------------ online_2D_serialmodel_openmp LESTKF ----------------"
export OMP_NUM_THREADS=2
cd online_2D_serialmodel
make cleandataq
mpirun --oversubscribe -np 9 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_serialmodel_openmp
cd ..
python ../verification/check_online_cl.py online_2D_serialmodel

echo "------------ online_2D_parallelmodel LESTKF ---------------------"
export OMP_NUM_THREADS=1
cd online_2D_parallelmodel
make cleandataq
mpirun --oversubscribe -np 18 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_parallelmodel
cd ..
python ../verification/check_online_cl.py online_2D_parallelmodel

echo "------------ online_2D_serialmodel ESTKF ------------------------"
export OMP_NUM_THREADS=1
cd online_2D_serialmodel
make cleandataq
mpirun --oversubscribe -np 9 ./model_pdaf -dim_ens 9 $DA_SPECS2 > ../out.online_2D_serialmodel_ESTKF
cd ..
python ../verification/check_online2_cl.py online_2D_serialmodel online_2D_serialmodel_ESTKF

echo "------------ online_2D_parallelmodel ESTKF ----------------------"
export OMP_NUM_THREADS=1
cd online_2D_parallelmodel
make cleandataq
mpirun --oversubscribe -np 18 ./model_pdaf -dim_ens 9 $DA_SPECS2 > ../out.online_2D_parallelmodel_ESTKF
cd ..
python ../verification/check_online2_cl.py online_2D_parallelmodel online_2D_parallelmodel_ESTKF
