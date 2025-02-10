#!/bin/bash

# ARCH specifies PDAF_ARCH without and with PDAF
export ARCH=linux_gfortran_openmpi
DA_SPECS=" -dim_ens 8 -forget 0.8 -screen 1 -cradius 5.0"
RUNSTR="mpirun -np 1 ./PDAF_offline"
RUNPAR="mpirun -np 3 ./PDAF_offline"


COMPILEPDAF=0
COMPILE=1
TEST_SUBTYPES=1
TEST_OPTIONS=1
TEST_PARALLEL=1

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

if [ $TEST_SUBTYPES -eq 1 ]
then

    # SEIK ##############

    echo "     +++++++++++++ SEIK offline serial +++++++++++++"

    FTYPE=1
    STYPE=0
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}

    FTYPE=1
    STYPE=1
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}

    FTYPE=1
    STYPE=4
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}

    # LSEIK ##############

    echo "     +++++++++++++ LSEIK offline serial +++++++++++++"

    FTYPE=3
    STYPE=0
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}

    # EnKF ##############

    echo "     +++++++++++++ EnKF offline serial +++++++++++++"

    FTYPE=2
    STYPE=0
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}

    FTYPE=2
    STYPE=1
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}


    # LEnKF ##############

    echo "     +++++++++++++ LEnKF offline serial +++++++++++++"

    FTYPE=8
    STYPE=0
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}


    # ETKF ##############

    echo "     +++++++++++++ ETKF offline serial +++++++++++++"

    FTYPE=4
    STYPE=0
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}

    FTYPE=4
    STYPE=1
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}


    # LETKF ##############

    echo "     +++++++++++++ LETKF offline serial +++++++++++++"

    FTYPE=5
    STYPE=0
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}


    # ESTKF ##############

    echo "     +++++++++++++ ESTKF offline serial +++++++++++++"

    FTYPE=6
    STYPE=0
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}


    # LESTKF ##############

    echo "     +++++++++++++ LESTKF offline serial +++++++++++++"

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}


    # NETF ##############

    echo "     +++++++++++++ NETF offline serial +++++++++++++"

    FTYPE=9
    STYPE=0
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}


    # LNETF ##############

    echo "     +++++++++++++ LNETF offline serial +++++++++++++"

    FTYPE=10
    STYPE=0
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}


    # LKNETF ##############

    echo "     +++++++++++++ LKNETF offline serial +++++++++++++"

    FTYPE=11
    STYPE=0
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}

    FTYPE=11
    STYPE=1
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}

    FTYPE=11
    STYPE=4
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}


    # PF ##############

    echo "     +++++++++++++ PF offline serial +++++++++++++"

    FTYPE=12
    STYPE=0
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}


fi

if [ $TEST_OPTIONS -eq 1 ]
then
    # Test Special settings

    echo " "
    echo "     +++++++++++++ Test type_forget in LESTKF +++++++++++++"

    # TYPE_FORGET ##############

    FTYPE=7
    STYPE=0
    TFORGET=1
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, type_forget="$TFORGET " -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -type_forget $TFORGET
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -type_forget $TFORGET > ../out.offline_2D_filter${FTYPE}s${STYPE}tf${TFORGET}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}tf${TFORGET}

    FTYPE=7
    STYPE=0
    TFORGET=2
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, type_forget="$TFORGET " -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -type_forget $TFORGET
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -type_forget $TFORGET > ../out.offline_2D_filter${FTYPE}s${STYPE}tf${TFORGET}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}tf${TFORGET}


    # TYPE_TRANS ##############

    echo " "
    echo "     +++++++++++++ Test type_trans in LESTKF +++++++++++++"

    FTYPE=7
    STYPE=0
    TTRANS=1
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, type_trans="$TTRANS " -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -type_trans $TTRANS
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -type_trans $TTRANS > ../out.offline_2D_filter${FTYPE}s${STYPE}ttrans${TTRANS}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}ttrans${TTRANS}

    FTYPE=7
    STYPE=0
    TTRANS=2
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, type_trans="$TTRANS " -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -type_trans $TTRANS
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -type_trans $TTRANS > ../out.offline_2D_filter${FTYPE}s${STYPE}ttrans${TTRANS}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}ttrans${TTRANS}


    # TYPE_SQRT ##############

    echo " "
    echo "     +++++++++++++ Test type_sqrt in LSEIK, LESTKF and SEIK, ESTKF +++++++++++++"

    FTYPE=3
    STYPE=0
    TSQRT=1
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, type_sqrt="$TSQRT " -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -type_sqrt $TSQRT
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -type_sqrt $TSQRT > ../out.offline_2D_filter${FTYPE}s${STYPE}tsqrt${TSQRT}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}tsqrt${TSQRT}

    FTYPE=7
    STYPE=0
    TSQRT=1
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, type_sqrt="$TSQRT " -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -type_sqrt $TSQRT
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -type_sqrt $TSQRT > ../out.offline_2D_filter${FTYPE}s${STYPE}tsqrt${TSQRT}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}tsqrt${TSQRT}


    FTYPE=1
    STYPE=0
    TSQRT=1
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, type_sqrt="$TSQRT " -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -type_sqrt $TSQRT
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -type_sqrt $TSQRT > ../out.offline_2D_filter${FTYPE}s${STYPE}tsqrt${TSQRT}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}tsqrt${TSQRT}

    FTYPE=6
    STYPE=0
    TSQRT=1
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, type_sqrt="$TSQRT " -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -type_sqrt $TSQRT
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -type_sqrt $TSQRT > ../out.offline_2D_filter${FTYPE}s${STYPE}tsqrt${TSQRT}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}tsqrt${TSQRT}


    # ASSIM_A, ASSIM_B ##############

    echo " "
    echo "     +++++++++++++ Test selecting observation types in ESTKF and LESTKF +++++++++++++"

    FTYPE=6
    STYPE=0
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, obs=B -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -assim_A .false. -assim_B .true. 
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -assim_A .false. -assim_B .true. > ../out.offline_2D_filter${FTYPE}s${STYPE}obsB
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}obsB

    FTYPE=6
    STYPE=0
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, obs=A+B -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -assim_A .true. -assim_B .true. 
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -assim_A .true. -assim_B .true. > ../out.offline_2D_filter${FTYPE}s${STYPE}obsAB
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}obsAB

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, obs=B -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -assim_A .false. -assim_B .true. 
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -assim_A .false. -assim_B .true. > ../out.offline_2D_filter${FTYPE}s${STYPE}obsB
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}obsB

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, serial, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, obs=A+B -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -assim_A .true. -assim_B .true. 
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -assim_A .true. -assim_B .true. > ../out.offline_2D_filter${FTYPE}s${STYPE}obsAB
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}obsAB

fi

if [ $TEST_PARALLEL -eq 1 ]
then

    echo "     +++++++++++++ PARALLEL TESTS +++++++++++++"

    # We run with 3 MPI tasks since then task=1 has dim_obs_p=0, i.e. there are no observation
    # on the middle process subdomain. This is a special case that has to be correctly handled by OMI.
    # The results have to be identical to those of the serial runs.


    # EnKF ##############

    echo "     +++++++++++++ EnKF offline parallel +++++++++++++"

    FTYPE=2
    STYPE=0
    echo "-------offline_2D, parallel, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNPAR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNPAR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_par_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}

    FTYPE=2
    STYPE=1
    echo "-------offline_2D, parallel, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNPAR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNPAR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_par_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}


    # LEnKF ##############

    echo "     +++++++++++++ LEnKF offline parallel +++++++++++++"

    FTYPE=8
    STYPE=0
    echo "-------offline_2D, parallel, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNPAR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNPAR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_par_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}


    # ESTKF ##############

    echo "     +++++++++++++ ESTKF offline parallel +++++++++++++"

    FTYPE=6
    STYPE=0
    echo "-------offline_2D, parallel, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNPAR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNPAR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_par_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}


    # LESTKF ##############

    echo "     +++++++++++++ LESTKF offline parallel +++++++++++++"

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, parallel, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNPAR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNPAR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_par_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}


    # NETF ##############

    echo "     +++++++++++++ NETF offline parallel +++++++++++++"

    FTYPE=9
    STYPE=0
    echo "-------offline_2D, parallel, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNPAR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNPAR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_par_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}


    # LNETF ##############

    echo "     +++++++++++++ LNETF offline parallel +++++++++++++"

    FTYPE=10
    STYPE=0
    echo "-------offline_2D, parallel, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=4
    cd offline_2D_parallel
    make cleandataq
    echo $RUNPAR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNPAR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_par_filter${FTYPE}s${STYPE}
    cd ..
    python verification/check_offline2.py offline_2D_parallel offline_2D_ftype${FTYPE}s${STYPE}

fi
