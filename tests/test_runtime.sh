#!/bin/bash

# ARCH specifies PDAF_ARCH without and with PDAF
export ARCH=linux_gfortran_openmpi
DA_SPECS=" -dim_ens 8 -forget 0.8 -screen 1 -cradius 500.0 -gridsize 3"
RUNSTR="mpirun -np 1 ./PDAF_offline"
RUNPAR="mpirun -np 4 ./PDAF_offline"

COMPILEPDAF=1
COMPILE=1
GENERATE_INPUTS=1
TEST_PAR=1
TEST_ENSRF=1
DO_CHECK=0

echo "------------------ GENERATE INPUT FILES ----------------"

if [ $GENERATE_INPUTS -eq 1 ]
then
     cd inputs_offline.512x2048
     python generate_test.py
     cd ..
else
    echo "Generation of input files is deactivated! Please ensure that the files exist."
fi

echo "------------------ COMPILING ----------------"

if [ $COMPILE -eq 1 ]
then

    export PDAF_ARCH=$ARCH
    echo -------------- PDAF_ARCH: $PDAF_ARCH

    echo "------------ offline_2D_parallel ---------------"
    export PDAF_ARCH=$ARCH
    echo PDAF_ARCH: $PDAF_ARCH
    cd offline_2D_parallel
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

#--------- OFFLINE -------------

if [ $TEST_PAR -eq 1 ]
then

    # ESTKF ##############

    echo "     +++++++++++++ ESTKF offline large grid +++++++++++++"

    FTYPE=6
    STYPE=0
    echo "-------offline_2D, serial, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_grid2_filter${FTYPE}s${STYPE}
    cd ..
    if [ $DO_CHECK -eq 1 ]
    then
	python verification/check_offline2.py offline_2D_parallel offline_2D_grid2_ftype${FTYPE}s${STYPE}
    fi

    FTYPE=6
    STYPE=0
    echo "-------offline_2D, MPI 4 tasks, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNPAR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNPAR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_grid2_MPI4_filter${FTYPE}s${STYPE}
    cd ..
    if [ $DO_CHECK -eq 1 ]
    then
	python verification/check_offline2.py offline_2D_parallel offline_2D_grid2_ftype${FTYPE}s${STYPE}
    fi

    # LESTKF ##############

    echo "     +++++++++++++ LESTKF offline  large grid +++++++++++++"

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}
    cd ..
    if [ $DO_CHECK -eq 1 ]
    then
	python verification/check_offline2.py offline_2D_parallel offline_2D_grid2_ftype${FTYPE}s${STYPE}
    fi
    
    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=10, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=10
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_grid2_OMP10_filter${FTYPE}s${STYPE}
    cd ..
    if [ $DO_CHECK -eq 1 ]
    then
	python verification/check_offline2.py offline_2D_parallel offline_2D_grid2_ftype${FTYPE}s${STYPE}
    fi

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, MPInp=4, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNPAR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNPAR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_grid2_MPI4_filter${FTYPE}s${STYPE}
    cd ..
    if [ $DO_CHECK -eq 1 ]
    then
	python verification/check_offline2.py offline_2D_parallel offline_2D_grid2_ftype${FTYPE}s${STYPE}
    fi
fi

if [ $TEST_ENSRF -eq 1 ]
then
    
    # ENSRF ##############

    echo "     +++++++++++++ ENSRF offline large grid +++++++++++++"

#     FTYPE=13
#     STYPE=1
#     echo "-------offline_2D, serial, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
#     export OMP_NUM_THREADS=1
#     cd offline_2D_parallel
#     make cleandataq
#     echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
#     $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE #> ../out.offline_2D_grid2_filter${FTYPE}s${STYPE}
#     cd ..
#     if [ $DO_CHECK -eq 1 ]
#     then
# 	python verification/check_offline2.py offline_2D_parallel offline_2D_grid2_ftype${FTYPE}s${STYPE}
#     fi

    FTYPE=13
    STYPE=1
    echo "-------offline_2D, MPI 4 tasks, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNPAR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNPAR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_grid2_MPI4_filter${FTYPE}s${STYPE}
    cd ..
    if [ $DO_CHECK -eq 1 ]
    then
	python verification/check_offline2.py offline_2D_parallel offline_2D_grid2_ftype${FTYPE}s${STYPE}
    fi

fi
