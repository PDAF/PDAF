#!/bin/bash

# ARCH specifies PDAF_ARCH without and with PDAF
export ARCH=linux_gfortran_openmpi
export DA_SPECS="-filtertype 7 -screen 1"
export DA_SPECS2="-filtertype 6 -screen 1"
export DA_SPECS3="-filtertype 7 screen 1 -assim_A .false. -assim_B .true"

COMPILE=1
RUN_OFFLINE=1
RUN_ONLINE_SERIAL=1
RUN_ONLINE_PARALLEL=1

echo "------------------ COMPILING ----------------"

if [ $COMPILE -eq 1 ]
then

    export PDAF_ARCH=$ARCH
    echo -------------- PDAF_ARCH: $PDAF_ARCH

    echo "------------ offline_2D_serial --------------"
    cd offline_2D_serial
    make cleanall
    make
    cd ..

    echo "------------ offline_2D_parallel ---------------"
    export PDAF_ARCH=$ARCH
    echo PDAF_ARCH: $PDAF_ARCH
    cd offline_2D_parallel
    make cleanall
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

    echo "------------ online_2D_serialmodel_2fields ---------------"
    export PDAF_ARCH=$ARCH
    echo PDAF_ARCH: $PDAF_ARCH
    cd online_2D_serialmodel_2fields
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

    echo "------------ online_2D_parallelmodel_fullpar ---------------"
    export PDAF_ARCH=$ARCH
    echo PDAF_ARCH: $PDAF_ARCH
    cd online_2D_parallelmodel_fullpar
    make clean
    make cleandata
    make model_pdaf
    cd ..

    echo "------------ online_2D_parallelmodel_fullpar_1fpe ---------------"
    export PDAF_ARCH=$ARCH
    echo PDAF_ARCH: $PDAF_ARCH
    cd online_2D_parallelmodel_fullpar_1fpe
    make clean
    make cleandata
    make model_pdaf
    cd ..

else
    echo "Compilation is deactivated!"
fi

echo  " "
echo "-------------------- RUNNING ----------------"

#--------- OFFLINE -------------

if [ $RUN_OFFLINE -eq 1 ]
then

    echo "------------ offline_2D_serial ---------------------------------------------"
    export OMP_NUM_THREADS=1
    cd offline_2D_serial
    make cleandataq
    ./PDAF_offline $DA_SPECS > ../out.offline_2D_serial
    cd ..
    python verification/check_offline.py offline_2D_serial

    echo "------------ offline_2D_serial ---------------------------------------------"
    export OMP_NUM_THREADS=1
    cd offline_2D_serial
    make cleandataq
    ./PDAF_offline $DA_SPECS2 > ../out.offline_2D_serial_ESTKF
    cd ..
    python verification/check_offline2.py offline_2D_serial offline_2D_serial_ESTKF

    echo "------------ offline_2D_openmp ---------------------------------------------"
    export OMP_NUM_THREADS=4
    cd offline_2D_serial
    make cleandataq
    ./PDAF_offline $DA_SPECS > ../out.offline_2D_openmp
    cd ..
    python verification/check_offline.py offline_2D_serial

    echo "------------ offline_2D_parallel -------------------------------------------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    mpirun -np 4 ./PDAF_offline $DA_SPECS > ../out.offline_2D_parallel
    cd ..
    python verification/check_offline.py offline_2D_parallel
fi

#--------- ONLINE SERIAL -------------

if [ $RUN_ONLINE_SERIAL -eq 1 ]
then

    echo "------------ online_2D_serialmodel LESTKF ----------------------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_serialmodel
    make cleandataq
    mpirun --oversubscribe -np 9 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_serialmodel
    cd ..
    python verification/check_online2.py online_2D_serialmodel online_2D_serialmodel

    echo "------------ online_2D_serialmodel_openmp LESTKF ---------------------------"
    export OMP_NUM_THREADS=2
    cd online_2D_serialmodel
    make cleandataq
    mpirun --oversubscribe -np 9 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_serialmodel_openmp
    cd ..
    python verification/check_online.py online_2D_serialmodel

    echo "------------ online_2D_serialmodel_2fields LESTKF --------------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_serialmodel_2fields
    make cleandataq
    mpirun --oversubscribe -np 9 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_serialmodel_2fields
    cd ..
    python verification/check_online2.py online_2D_serialmodel_2fields online_2D_serialmodel

    echo "------------ online_2D_serialmodel_2fields LESTKF obs-type B ---------------"
    export OMP_NUM_THREADS=1
    cd online_2D_serialmodel_2fields
    make cleandataq
    mpirun --oversubscribe -np 9 ./model_pdaf -dim_ens 9 $DA_SPECS3 > ../out.online_2D_serialmodel_2fields_obsB
    cd ..
    python verification/check_online3.py online_2D_serialmodel_2fields online_2D_serialmodel_2fields_obsB

    echo "------------ online_2D_serialmodel ESTKF -----------------------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_serialmodel
    make cleandataq
    mpirun --oversubscribe -np 9 ./model_pdaf -dim_ens 9 $DA_SPECS2 > ../out.online_2D_serialmodel_ESTKF
    cd ..
    python verification/check_online2.py online_2D_serialmodel online_2D_serialmodel_ESTKF

    echo "------------ online_2D_serialmodel_2fields ESTKF ---------------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_serialmodel_2fields
    make cleandataq
    mpirun --oversubscribe -np 9 ./model_pdaf -dim_ens 9 $DA_SPECS2 > ../out.online_2D_serialmodel_2fields_ESTKF
    cd ..
    python verification/check_online2.py online_2D_serialmodel_2fields online_2D_serialmodel_ESTKF
fi

#--------- ONLINE PARALLEL -------------

if [ $RUN_ONLINE_PARALLEL -eq 1 ]
then

    echo "------------ online_2D_parallelmodel LESTKF --------------------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel
    make cleandataq
    mpirun --oversubscribe -np 18 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_parallelmodel
    cd ..
    python verification/check_online.py online_2D_parallelmodel

    echo "------------ online_2D_parallelmodel ESTKF ---------------------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel
    make cleandataq
    mpirun --oversubscribe -np 18 ./model_pdaf -dim_ens 9 $DA_SPECS2 > ../out.online_2D_parallelmodel_ESTKF
    cd ..
    python verification/check_online2.py online_2D_parallelmodel online_2D_parallelmodel_ESTKF

    echo "------------ online_2D_parallelmodel_fullpar LESTKF ------------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel_fullpar
    make cleandataq
    mpirun --oversubscribe -np 20 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_parallelmodel_fullpar
    cd ..
    python verification/check_online.py online_2D_parallelmodel_fullpar

    echo "------------ online_2D_parallelmodel_fullpar_1fpe LESTKF -------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel_fullpar_1fpe
    make cleandataq
    mpirun --oversubscribe -np 19 ./model_pdaf -dim_ens 9 $DA_SPECS > ../out.online_2D_parallelmodel_fullpar_1fpe
    cd ..
    python verification/check_online.py online_2D_parallelmodel_fullpar_1fpe
fi
