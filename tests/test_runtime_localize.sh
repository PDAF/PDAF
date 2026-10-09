#!/bin/bash

# ARCH specifies PDAF_ARCH without and with PDAF
export ARCH=linux_gfortran_openmpi
DA_SPECS=" -dim_ens 8 -forget 0.8 -screen 1 -cradius 500.0 -gridsize 3"
DA_SPECS_GEO=" -dim_ens 8 -forget 0.8 -screen 1 -gridsize 3 -cradius 1.0e4"
RUNSTR="mpirun -np 1 ./PDAF_offline"
RUNPAR="mpirun -np 4 ./PDAF_offline"

COMPILEPDAF=0
COMPILE=1
GENERATE_INPUTS=1
TEST_PAR=1         # Cartesian coordinates (type_coords=0)
TEST_PAR_GEO=1     # Geographic coordinates, simplified distance (type_coords=2)
TEST_PAR_GC=1      # Geographic coordinates, great circle distance (type_coords=3)
TEST_PAR_GEO_NONISO=1
TEST_PAR_GC_NONISO=1
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


    # LESTKF ##############

    echo "     +++++++++++++ LESTKF offline  large grid +++++++++++++"

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 0 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -omi_search_type 0
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -omi_search_type 0 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}search0
    cd ..
    if [ $DO_CHECK -eq 1 ]
    then
	python verification/check_offline2.py offline_2D_parallel offline_2D_grid2_ftype${FTYPE}s${STYPE}
    fi

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -omi_search_type 1
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -omi_search_type 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}search1
    cd ..
    if [ $DO_CHECK -eq 1 ]
    then
	python verification/check_offline2.py offline_2D_parallel offline_2D_grid2_ftype${FTYPE}s${STYPE}
    fi

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}search2
    cd ..
    if [ $DO_CHECK -eq 1 ]
    then
	python verification/check_offline2.py offline_2D_parallel offline_2D_grid2_ftype${FTYPE}s${STYPE}
    fi

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 11, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -omi_search_type 11 -omi_sort_dir 1
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -omi_search_type 11 -omi_sort_dir 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}search11dir1
    cd ..
    if [ $DO_CHECK -eq 1 ]
    then
	python verification/check_offline2.py offline_2D_parallel offline_2D_grid2_ftype${FTYPE}s${STYPE}
    fi

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 11, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -omi_search_type 11 -omi_sort_dir 2
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -omi_search_type 11 -omi_sort_dir 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}search11dir2
    cd ..
    if [ $DO_CHECK -eq 1 ]
    then
	python verification/check_offline2.py offline_2D_parallel offline_2D_grid2_ftype${FTYPE}s${STYPE}
    fi

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 12, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -omi_search_type 12 -omi_sort_dir 1
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -omi_search_type 12 -omi_sort_dir 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}search12dir1
    cd ..
    if [ $DO_CHECK -eq 1 ]
    then
	python verification/check_offline2.py offline_2D_parallel offline_2D_grid2_ftype${FTYPE}s${STYPE}
    fi

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 12, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -omi_search_type 12 -omi_sort_dir 2
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -omi_search_type 12 -omi_sort_dir 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}search12dir2
    cd ..
    if [ $DO_CHECK -eq 1 ]
    then
	python verification/check_offline2.py offline_2D_parallel offline_2D_grid2_ftype${FTYPE}s${STYPE}
    fi

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 21, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -omi_search_type 21 -omi_sort_dir 1
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -omi_search_type 21 -omi_sort_dir 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}search21dir1
    cd ..
    if [ $DO_CHECK -eq 1 ]
    then
	python verification/check_offline2.py offline_2D_parallel offline_2D_grid2_ftype${FTYPE}s${STYPE}
    fi

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 21, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -omi_search_type 21 -omi_sort_dir 2
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -omi_search_type 21 -omi_sort_dir 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}search21dir2
    cd ..
    if [ $DO_CHECK -eq 1 ]
    then
	python verification/check_offline2.py offline_2D_parallel offline_2D_grid2_ftype${FTYPE}s${STYPE}
    fi

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 22, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -omi_search_type 22 -omi_sort_dir 1
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -omi_search_type 22 -omi_sort_dir 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}search22dir1
    cd ..
    if [ $DO_CHECK -eq 1 ]
    then
	python verification/check_offline2.py offline_2D_parallel offline_2D_grid2_ftype${FTYPE}s${STYPE}
    fi

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 22, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS -filtertype $FTYPE -subtype $STYPE -omi_search_type 22 -omi_sort_dir 2
    $RUNSTR $DA_SPECS  -filtertype $FTYPE -subtype $STYPE -omi_search_type 22 -omi_sort_dir 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}search22dir2
    cd ..
    if [ $DO_CHECK -eq 1 ]
    then
	python verification/check_offline2.py offline_2D_parallel offline_2D_grid2_ftype${FTYPE}s${STYPE}
    fi
fi

if [ $TEST_PAR_GEO -eq 1 ]
then

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 0 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 0
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 0 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords2search0
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords2search1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 0 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords2search2
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 11, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 11 -omi_sort_dir 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 11 -omi_sort_dir 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords2search11dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 11, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 11 -omi_sort_dir 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 11 -omi_sort_dir 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords2search11dir2
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 12, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 12 -omi_sort_dir 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 12 -omi_sort_dir 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords2search12dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 12, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 12 -omi_sort_dir 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 12 -omi_sort_dir 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords2search12dir2
    cd ..


    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 21, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 21 -omi_sort_dir 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 21 -omi_sort_dir 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords2search21dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 21, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 21 -omi_sort_dir 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 21 -omi_sort_dir 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords2search21dir2
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 22, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 22 -omi_sort_dir 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 22 -omi_sort_dir 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords2search22dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 22, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 22 -omi_sort_dir 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 22 -omi_sort_dir 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords2search22dir2
    cd ..
fi



if [ $TEST_PAR_GEO_NONISO -eq 1 ]
then

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 0 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 0 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 0 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords2search0
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 1 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 1 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords2search1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 0 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 2 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 2 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords2search2
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 11, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 11 -omi_sort_dir 1 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 11 -loc_noniso 1 -omi_sort_dir 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords2search11dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 11, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 11 -omi_sort_dir 2 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 11 -omi_sort_dir 2 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords2search11dir2
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 12, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 12 -omi_sort_dir 1 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 12 -omi_sort_dir 1 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords2search12dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 12, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 12 -omi_sort_dir 2 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 12 -omi_sort_dir 2 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords2search12dir2
    cd ..


    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 21, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 21 -omi_sort_dir 1 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 21 -omi_sort_dir 1 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords2search21dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 21, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 21 -omi_sort_dir 2 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 21 -omi_sort_dir 2 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords2search21dir2
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 22, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 22 -omi_sort_dir 1 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 22 -omi_sort_dir 1 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords2search22dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 22, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 22 -omi_sort_dir 2 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 22 -omi_sort_dir 2 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords2search22dir2
    cd ..

    echo " +++++ Test noniso 2 +++++"

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 0 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 0 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 0 -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords2search0
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 1 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 1 -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords2search1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 0 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 2 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 2 -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords2search2
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 11, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 11 -omi_sort_dir 1 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 11 -loc_noniso 2 -omi_sort_dir 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords2search11dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 11, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 11 -omi_sort_dir 2 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 11 -omi_sort_dir 2 -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords2search11dir2
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 12, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 12 -omi_sort_dir 1 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 12 -omi_sort_dir 1 -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords2search12dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 12, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 12 -omi_sort_dir 2 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 12 -omi_sort_dir 2 -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords2search12dir2
    cd ..


    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 21, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 21 -omi_sort_dir 1 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 21 -omi_sort_dir 1 -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords2search21dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 21, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 21 -omi_sort_dir 2 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 21 -omi_sort_dir 2 -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords2search21dir2
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 22, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 22 -omi_sort_dir 1 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 22 -omi_sort_dir -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords2search22dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 22, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 22 -omi_sort_dir 2 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 2 -omi_search_type 22 -omi_sort_dir 2 -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords2search22dir2
    cd ..

fi

if [ $TEST_PAR_GC -eq 1 ]
then

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 0 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 0
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 0 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords3search0
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords3search1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 0 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords3search2
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 11, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 11 -omi_sort_dir 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 11 -omi_sort_dir 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords3search11dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 11, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 11 -omi_sort_dir 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 11 -omi_sort_dir 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords3search11dir2
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 12, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 12 -omi_sort_dir 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 12 -omi_sort_dir 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords3search12dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 12, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 12 -omi_sort_dir 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 12 -omi_sort_dir 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords3search12dir2
    cd ..


    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 21, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 21 -omi_sort_dir 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 21 -omi_sort_dir 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords3search21dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 21, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 21 -omi_sort_dir 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 21 -omi_sort_dir 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords3search21dir2
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 22, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 22 -omi_sort_dir 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 22 -omi_sort_dir 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords3search22dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 22, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 22 -omi_sort_dir 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 22 -omi_sort_dir 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}coords3search22dir2
    cd ..

fi



if [ $TEST_PAR_GC_NONISO -eq 1 ]
then

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 0 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 0 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 0 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords3search0
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 1 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 1 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords3search1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 0 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 2 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 2 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords3search2
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 11, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 11 -omi_sort_dir 1 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 11 -omi_sort_dir 1 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords3search11dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 11, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 11 -omi_sort_dir 2 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 11 -omi_sort_dir 2 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords3search11dir2
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 12, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 12 -omi_sort_dir 1 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 12 -omi_sort_dir 1 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords3search12dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 12, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 12 -omi_sort_dir 2 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 12 -omi_sort_dir 2 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords3search12dir2
    cd ..


    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 21, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 21 -omi_sort_dir 1 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 21 -omi_sort_dir 1 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords3search21dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 21, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 21 -omi_sort_dir 2 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 21 -omi_sort_dir 2 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords3search21dir2
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 22, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 22 -omi_sort_dir 1 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 22 -omi_sort_dir 1 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords3search22dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 22, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 22 -omi_sort_dir 2 -loc_noniso 1
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 22 -omi_sort_dir 2 -loc_noniso 1 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso1coords3search22dir2
    cd ..

    echo " +++++ Test noniso 2 +++++"



    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 0 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 0 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 0 -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords3search0
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 1 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 1 -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords3search1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 0 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 2 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 2 -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords3search2
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 11, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 11 -omi_sort_dir 1 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 11 -omi_sort_dir 1 -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords3search11dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 11, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 11 -omi_sort_dir 2 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 11 -omi_sort_dir 2 -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords3search11dir2
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 12, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 12 -omi_sort_dir 1 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 12 -omi_sort_dir 1 -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords3search12dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 12, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 12 -omi_sort_dir 2 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 12 -omi_sort_dir 2 -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords3search12dir2
    cd ..


    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 21, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 21 -omi_sort_dir 1 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 21 -omi_sort_dir 1 -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords3search21dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 21, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 21 -omi_sort_dir 2 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 21 -omi_sort_dir 2 -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords3search21dir2
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 22, dir 1 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 22 -omi_sort_dir 1 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 22 -omi_sort_dir 1 -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords3search22dir1
    cd ..

    FTYPE=7
    STYPE=0
    echo "-------offline_2D, threads=1, filtertype="$FTYPE ", subtype="$STYPE ", forget 0.8, search type 22, dir 2 -----------"
    export OMP_NUM_THREADS=1
    cd offline_2D_parallel
    make cleandataq
    echo $RUNSTR $DA_SPECS_GEO -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 22 -omi_sort_dir 2 -loc_noniso 2
    $RUNSTR $DA_SPECS_GEO  -filtertype $FTYPE -subtype $STYPE -type_coords 3 -omi_search_type 22 -omi_sort_dir 2 -loc_noniso 2 > ../out.offline_2D_grid2_OMP1_filter${FTYPE}s${STYPE}noniso2coords3search22dir2
    cd ..

fi