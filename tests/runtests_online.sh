#!/bin/bash

# ARCH specifies PDAF_ARCH without and with PDAF
export ARCH=linux_gfortran_openmpi
export DA_SPECS="-filtertype 7 -screen 1"
export DA_SPECS2="-filtertype 6 -screen 1"
export DA_SPECS3="-filtertype 7 screen 1 -assim_A .false. -assim_B .true"

COMPILEPDAF=0
COMPILE=0
RUN_ONLINE_SERIAL=1
RUN_ONLINE_PARALLEL=1

echo "------------------ COMPILING ----------------"

if [ $COMPILE -eq 1 ]
then

    export PDAF_ARCH=$ARCH
    echo -------------- PDAF_ARCH: $PDAF_ARCH

    echo "------------ online_2D_parallelmodel ---------------"
    export PDAF_ARCH=$ARCH
    echo PDAF_ARCH: $PDAF_ARCH
    cd online_2D_parallelmodel
    if [ $COMPILEPDAF -eq 1 ]
    then
	make cleanall
    else
	make clean
	make cleandata
    fi
    make
    cd ..

else
    echo "Compilation is deactivated!"
fi

echo  " "
echo "-------------------- RUNNING ----------------"


#--------- ONLINE SERIAL -------------

if [ $RUN_ONLINE_SERIAL -eq 1 ]
then

    echo "------------ online_2D_serialmodel LESTKF ----------------------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel
    make cleandataq
    mpirun --oversubscribe -np 9 ./PDAF_online -dim_ens 9 $DA_SPECS > ../out.online_2D_serialmodel
    cd ..
    python verification/check_online2.py online_2D_parallelmodel online_2D_serialmodel

    echo "------------ online_2D_serialmodel_openmp LESTKF ---------------------------"
    export OMP_NUM_THREADS=6
    cd online_2D_parallelmodel
    make cleandataq
    mpirun --oversubscribe -np 9 ./PDAF_online -dim_ens 9 $DA_SPECS > ../out.online_2D_serialmodel_openmp
    cd ..
    python verification/check_online2.py online_2D_parallelmodel online_2D_serialmodel

    echo "------------ online_2D_serialmodel ESTKF -----------------------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel
    make cleandataq
    mpirun --oversubscribe -np 9 ./PDAF_online -dim_ens 9 $DA_SPECS2 > ../out.online_2D_serialmodel_ESTKF
    cd ..
    python verification/check_online2.py online_2D_parallelmodel online_2D_serialmodel_ESTKF
fi

#--------- ONLINE PARALLEL -------------

if [ $RUN_ONLINE_PARALLEL -eq 1 ]
then

    echo "------------ online_2D_parallelmodel LESTKF --------------------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel
    make cleandataq
    mpirun --oversubscribe -np 18 ./PDAF_online -dim_ens 9 $DA_SPECS > ../out.online_2D_parallelmodel
    cd ..
    python verification/check_online2.py online_2D_parallelmodel online_2D_serialmodel

    echo "------------ online_2D_parallelmodel ESTKF ---------------------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel
    make cleandataq
    mpirun --oversubscribe -np 18 ./PDAF_online -dim_ens 9 $DA_SPECS2 > ../out.online_2D_parallelmodel_ESTKF
    cd ..
    python verification/check_online2.py online_2D_parallelmodel online_2D_serialmodel_ESTKF
fi
