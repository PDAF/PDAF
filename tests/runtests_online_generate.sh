#!/bin/bash

# ARCH specifies PDAF_ARCH without and with PDAF
export ARCH=linux_gfortran_openmpi
DA_SPECS=" -dim_ens 4 -forget 0.8 -screen 1 -cradius 5.0 -locweight 2 -delt_obs 9"
DA_SPECS_IAU=" -dim_ens 4 -forget 0.8 -screen 1 -cradius 5.0 -locweight 2 -delt_obs 9 -type_iau 1 -steps_iau 5"
RUNSTR="mpirun -np 4 ./PDAF_online"
RUNPAR="mpirun -np 8 ./PDAF_online"


COMPILEPDAF=0
COMPILE=0
TEST_SUBTYPES=1
TEST_IAU=1
TEST_PARALLEL=1

echo "------------------ COMPILING ----------------"

if [ $COMPILE -eq 1 ]
then

    export PDAF_ARCH=$ARCH
    echo -------------- PDAF_ARCH: $PDAF_ARCH

    echo "------------ online_2d_parallelmodel ---------------"
    export PDAF_ARCH=$ARCH
    echo PDAF_ARCH: $PDAF_ARCH
    cd online_2d_parallelmodel
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

#--------- ONLINE SUBTYPE -------------

if [ $TEST_SUBTYPES -eq 1 ]
then

    # SEIK ##############

    echo "     +++++++++++++ SEIK online 1 task per model +++++++++++++"

    FTYPE=1
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}

    FTYPE=1
    STYPE=1
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}

    FTYPE=1
    STYPE=4
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}

    # LSEIK ##############

    echo "     +++++++++++++ LSEIK online 1 task per model +++++++++++++"

    FTYPE=3
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}

    # EnKF ##############

    echo "     +++++++++++++ EnKF online 1 task per model +++++++++++++"

    FTYPE=2
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}

    FTYPE=2
    STYPE=1
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}


    # LEnKF ##############

    echo "     +++++++++++++ LEnKF online 1 task per model +++++++++++++"

    FTYPE=8
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}


    # ETKF ##############

    echo "     +++++++++++++ ETKF online 1 task per model +++++++++++++"

    FTYPE=4
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}

    FTYPE=4
    STYPE=1
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}


    # LETKF ##############

    echo "     +++++++++++++ LETKF online 1 task per model +++++++++++++"

    FTYPE=5
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}


    # ESTKF ##############

    echo "     +++++++++++++ ESTKF online 1 task per model +++++++++++++"

    FTYPE=6
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}


    # LESTKF ##############

    echo "     +++++++++++++ LESTKF online 1 task per model +++++++++++++"

    FTYPE=7
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}


    # NETF ##############

    echo "     +++++++++++++ NETF online 1 task per model +++++++++++++"

    FTYPE=9
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}


    # LNETF ##############

    echo "     +++++++++++++ LNETF online 1 task per model +++++++++++++"

    FTYPE=10
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}


    # LKNETF ##############

    echo "     +++++++++++++ LKNETF online 1 task per model +++++++++++++"

    FTYPE=11
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}

    FTYPE=11
    STYPE=1
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}

    FTYPE=11
    STYPE=4
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}


    # PF ##############

    echo "     +++++++++++++ PF online 1 task per model +++++++++++++"

    FTYPE=12
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}


    # ENSRF ##############

    echo "     +++++++++++++ ENSRF online 1 task per model +++++++++++++"

    FTYPE=13
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype13s${STYPE}

    FTYPE=13
    STYPE=1
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype13s${STYPE}

fi

#--------- ONLINE SUBTYPE IAU -------------

if [ $TEST_IAU -eq 1 ]
then

    # SEIK ##############

    echo "     +++++++++++++ SEIK online IAU 1 task per model +++++++++++++"

    FTYPE=1
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS_IAU -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS_IAU  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_IAU_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_IAU_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_IAU_ftype${FTYPE}s${STYPE}

    # LSEIK ##############

    echo "     +++++++++++++ LSEIK online IAU 1 task per model +++++++++++++"

    FTYPE=3
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS_IAU -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS_IAU  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_IAU_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_IAU_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_IAU_ftype${FTYPE}s${STYPE}

    # EnKF ##############

    echo "     +++++++++++++ EnKF online IAU 1 task per model +++++++++++++"

    FTYPE=2
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS_IAU -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS_IAU  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_IAU_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_IAU_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_IAU_ftype${FTYPE}s${STYPE}

    # LEnKF ##############

    echo "     +++++++++++++ LEnKF online IAU 1 task per model +++++++++++++"

    FTYPE=8
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS_IAU -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS_IAU  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_IAU_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_IAU_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_IAU_ftype${FTYPE}s${STYPE}


    # ETKF ##############

    echo "     +++++++++++++ ETKF online IAU 1 task per model +++++++++++++"

    FTYPE=4
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS_IAU -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS_IAU  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_IAU_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_IAU_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_IAU_ftype${FTYPE}s${STYPE}


    # LETKF ##############

    echo "     +++++++++++++ LETKF online IAU 1 task per model +++++++++++++"

    FTYPE=5
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS_IAU -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS_IAU  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_IAU_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_IAU_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_IAU_ftype${FTYPE}s${STYPE}


    # ESTKF ##############

    echo "     +++++++++++++ ESTKF online IAU 1 task per model +++++++++++++"

    FTYPE=6
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS_IAU -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS_IAU  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_IAU_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_IAU_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_IAU_ftype${FTYPE}s${STYPE}


    # LESTKF ##############

    echo "     +++++++++++++ LESTKF online IAU 1 task per model +++++++++++++"

    FTYPE=7
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS_IAU -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS_IAU  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_IAU_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_IAU_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_IAU_ftype${FTYPE}s${STYPE}


    # NETF ##############

    echo "     +++++++++++++ NETF online IAU 1 task per model +++++++++++++"

    FTYPE=9
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS_IAU -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS_IAU  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_IAU_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_IAU_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_IAU_ftype${FTYPE}s${STYPE}


    # LNETF ##############

    echo "     +++++++++++++ LNETF online IAU 1 task per model +++++++++++++"

    FTYPE=10
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS_IAU -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS_IAU  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_IAU_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_IAU_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_IAU_ftype${FTYPE}s${STYPE}


    # LKNETF ##############

    echo "     +++++++++++++ LKNETF online IAU 1 task per model +++++++++++++"

    FTYPE=11
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS_IAU -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS_IAU  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_IAU_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_IAU_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_IAU_ftype${FTYPE}s${STYPE}


    # PF ##############

    echo "     +++++++++++++ PF online IAU 1 task per model +++++++++++++"

    FTYPE=12
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS_IAU -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS_IAU  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_IAU_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_IAU_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_IAU_ftype${FTYPE}s${STYPE}


    # ENSRF ##############

    echo "     +++++++++++++ ENSRF online IAU 1 task per model +++++++++++++"

    FTYPE=13
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNSTR $DA_SPECS_IAU -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS_IAU  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_IAU_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_IAU_ftype${FTYPE}s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_IAU_ftype13s${STYPE}

fi


if [ $TEST_PARALLEL -eq 1 ]
then

    echo "     +++++++++++++ PARALLEL TESTS +++++++++++++"

    # We run with 3 MPI tasks since then task=1 has dim_obs_p=0, i.e. there are no observation
    # on the middle process subdomain. This is a special case that has to be correctly handled by OMI.
    # The results have to be identical to those of the serial runs.


    # EnKF ##############

    echo "     +++++++++++++ EnKF online parallel +++++++++++++"

    FTYPE=2
    STYPE=0
    echo "-------online_2D, parallel, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNPAR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNPAR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_par_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}

    FTYPE=2
    STYPE=1
    echo "-------online_2D, parallel, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNPAR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNPAR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_par_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}


    # LEnKF ##############

    echo "     +++++++++++++ LEnKF online parallel +++++++++++++"

    FTYPE=8
    STYPE=0
    echo "-------online_2D, parallel, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNPAR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNPAR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_par_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}


    # ESTKF ##############

    echo "     +++++++++++++ ESTKF online parallel +++++++++++++"

    FTYPE=6
    STYPE=0
    echo "-------online_2D, parallel, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNPAR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNPAR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_par_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}


    # LESTKF ##############

    echo "     +++++++++++++ LESTKF online parallel +++++++++++++"

    FTYPE=7
    STYPE=0
    echo "-------online_2D, parallel, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNPAR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNPAR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_par_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}


    # NETF ##############

    echo "     +++++++++++++ NETF online parallel +++++++++++++"

    FTYPE=9
    STYPE=0
    echo "-------online_2D, parallel, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNPAR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNPAR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_par_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}


    # LNETF ##############

    echo "     +++++++++++++ LNETF online parallel +++++++++++++"

    FTYPE=10
    STYPE=0
    echo "-------online_2D, parallel, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNPAR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNPAR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_par_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype${FTYPE}s${STYPE}


    # ENSRF ##############

    echo "     +++++++++++++ ENSRF online 1 task per model +++++++++++++"

    FTYPE=13
    STYPE=0
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNPAR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype13s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype13s${STYPE}

    FTYPE=13
    STYPE=1
    echo "-------online_2D, 1 task per model, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd online_2d_parallelmodel
    make cleandataq
    echo $RUNPAR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.online_2D_filter${FTYPE}s${STYPE}
    mv state*txt ens*txt ../verification/online_2D_ftype13s${STYPE}/
    cd ..
#    python verification/check_online2.py online_2d_parallelmodel online_2D_ftype13s${STYPE}

fi
