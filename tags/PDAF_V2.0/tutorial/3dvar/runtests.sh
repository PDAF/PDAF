#!/bin/bash

# ARCH and ARCH_MPI specify PDAF_ARCH without and with PDAF
export ARCH=linux_gfortran_openmpi
export DA_SPECS_3DVar="-dim_ens 1 -dim_cvec 9 -filtertype 200 -type_3dvar 0"
export DA_SPECS_3DEnVar="-dim_ens 9 -filtertype 200"
export DA_SPECS_hyb3DVar="-dim_ens 9 -filtertype 200"

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

#--------- OFFLINE -------------

if [ $RUN_OFFLINE -eq 1 ]
then

    echo " "
    echo "_______ offline_2D_serial ________"
    echo " "

    echo "------------ offline_2D_serial 3D-Var LFBGS -----------------------------"
    export OMP_NUM_THREADS=1
    cd offline_2D_serial
    make cleandataq
    ./PDAF_offline $DA_SPECS_3DVar  -type_opt 1 > ../out.offline_2D_serial_3dv_opt1
    cd ..
    python ../verification/check_offline_var.py offline_2D_serial offline_2D_serial_var_bfgs

    echo "------------ offline_2D_serial 3D-Var CG+ -------------------------------"
    export OMP_NUM_THREADS=1
    cd offline_2D_serial
    make cleandataq
    ./PDAF_offline $DA_SPECS_3DVar  -type_opt 2 > ../out.offline_2D_serial_3dv_opt2
    cd ..
    python ../verification/check_offline_var.py offline_2D_serial offline_2D_serial_var_cg

    echo "------------ offline_2D_serial 3D-Var plain CG --------------------------"
    export OMP_NUM_THREADS=1
    cd offline_2D_serial
    make cleandataq
    ./PDAF_offline $DA_SPECS_3DVar  -type_opt 3 > ../out.offline_2D_serial_3dv_opt3
    cd ..
    python ../verification/check_offline_var.py offline_2D_serial offline_2D_serial_var_cg

    echo "------------ offline_2D_serial 3D-EnVar/LESTKF CG+ ----------------------"
    export OMP_NUM_THREADS=1
    cd offline_2D_serial
    make cleandataq
    ./PDAF_offline $DA_SPECS_3DEnVar  -type_3dvar 1 -type_opt 2 > ../out.offline_2D_serial_3dlenvar_opt2
    cd ..
    python ../verification/check_offline_envar.py offline_2D_serial offline_2D_serial_lenvar_cg

    echo "------------ offline_2D_serial 3D-EnVar/ESTKF CG+ -----------------------"
    echo "   -> compare to ESTKF case offline_2D_serial_ESTKF, expect <1.e-9"
    export OMP_NUM_THREADS=1
    cd offline_2D_serial
    make cleandataq
    ./PDAF_offline $DA_SPECS_3DEnVar  -type_3dvar 4 -type_opt 2 > ../out.offline_2D_serial_3denvar_opt2
    cd ..
    python ../verification/check_offline_envar.py offline_2D_serial offline_2D_serial_ESTKF

    echo "------------ offline_2D_serial hybrid 3D-Var/LESTKF CG+ -----------------"
    export OMP_NUM_THREADS=1
    cd offline_2D_serial
    make cleandataq
    ./PDAF_offline $DA_SPECS_hyb3DVar -type_3dvar 6 -type_opt 2 -dim_cvec 9 > \
	../out.offline_2D_serial_3dlhybvar_opt2
    cd ..
    python ../verification/check_offline_envar.py offline_2D_serial offline_2D_serial_lenvar_cg

    echo "------------ offline_2D_serial hybrid 3D-Var/ESTKF CG+ ------------------"
    echo "   -> compare to ESTKF case offline_2D_serial_ESTKF, expect <1.e-9"
    export OMP_NUM_THREADS=1
    cd offline_2D_serial
    make cleandataq
    ./PDAF_offline $DA_SPECS_hyb3DVar -type_3dvar 7 -type_opt 2 -dim_cvec 9 > \
	../out.offline_2D_serial_3dhybvar_opt2
    cd ..
    python ../verification/check_offline_envar.py offline_2D_serial offline_2D_serial_ESTKF

fi

#--------- ONLINE SERIAL -------------

if [ $RUN_ONLINE_SERIAL -eq 1 ]
then

    echo " "
    echo "_______ online_2D_serialmodel ________"
    echo " "
    echo "------------ online_2D_serialmodel 3D-Var BFGS --------------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_serialmodel
    make cleandataq
    ./model_pdaf $DA_SPECS_3DVar  -type_opt 1 > ../out.online_2D_serial_3dv_opt1
    cd ..
    python ../verification/check_online_var2.py online_2D_serialmodel online_2D_serialmodel_var_bfgs

    echo "------------ online_2D_serialmodel 3D-Var CG+ ----------------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_serialmodel
    make cleandataq
    ./model_pdaf $DA_SPECS_3DVar  -type_opt 2 > ../out.online_2D_serial_3dv_opt2
    cd ..
    python ../verification/check_online_var2.py online_2D_serialmodel online_2D_serialmodel_var_cg

    echo "------------ online_2D_serialmodel 3D-Var plain CG -----------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_serialmodel
    make cleandataq
    ./model_pdaf $DA_SPECS_3DVar  -type_opt 3 > ../out.online_2D_serial_3dv_opt3
    cd ..
    python ../verification/check_online_var2.py online_2D_serialmodel online_2D_serialmodel_var_cg

    echo "------------ online_2D_serialmodel 3D-EnVar/LESTKF CG+ -------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_serialmodel
    make cleandataq
    mpirun --oversubscribe -np 9 ./model_pdaf $DA_SPECS_3DEnVar -subtype 1 -type_opt 2 > \
	../out.online_2D_serialmodel_3dlenvar_opt2
    cd ..
    python ../verification/check_online_envar2.py 1.e-9 online_2D_serialmodel online_2D_serialmodel_lenvar_cg

    echo "------------ online_2D_serialmodel 3D-EnVar/ESTKF CG+ --------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_serialmodel
    make cleandataq
    mpirun --oversubscribe -np 9 ./model_pdaf $DA_SPECS_3DEnVar -subtype 4 -type_opt 2 > \
	../out.online_2D_serialmodel_3denvar_opt2
    cd ..
#echo "   -> compare to ESTKF case online_2D_serialmodel_ESTKF - expect accuracy < 1.e-5"
#python ../verification/check_online_envar2.py 1.e-5 online_2D_serialmodel online_2D_serialmodel_ESTKF
#echo "   -> compare to ESTKF case online_2D_serialmodel_envar_cg_ESTKF - expect accuracy < 1.e-9"
    python ../verification/check_online_envar2.py 1.e-9 online_2D_serialmodel online_2D_serialmodel_envar_cg_ESTKF

    echo "------------ online_2D_serialmodel 3D-EnVar/LESTKF plain CG --------------"
    export OMP_NUM_THREADS=1
    cd online_2D_serialmodel
    make cleandataq
    mpirun --oversubscribe -np 9 ./model_pdaf $DA_SPECS_3DEnVar -subtype 1 -type_opt 3 > \
	../out.online_2D_serialmodel_3dlenvar_opt3
    cd ..
    python ../verification/check_online_envar2.py 1.e-9 online_2D_serialmodel online_2D_serialmodel_lenvar_cg

    echo "------------ online_2D_serialmodel hybrid 3D-Var/LESTKF CG+ --------------"
    export OMP_NUM_THREADS=1
    cd online_2D_serialmodel
    make cleandataq
    mpirun --oversubscribe -np 9 ./model_pdaf $DA_SPECS_hyb3DVar -subtype 6 -type_opt 2 -dim_cvec 9 > ../out.online_2D_serialmodel_3dlhybvar_opt2
    cd ..
    python ../verification/check_online_envar2.py 1.e-9 online_2D_serialmodel online_2D_serialmodel_lhybvar_cg

    echo "------------ online_2D_serialmodel hybrid 3D-Var/ESTKF CG+ ---------------"
    export OMP_NUM_THREADS=1
    cd online_2D_serialmodel
    make cleandataq
    mpirun --oversubscribe -np 9 ./model_pdaf $DA_SPECS_hyb3DVar -subtype 7 -type_opt 2 -dim_cvec 9 > ../out.online_2D_serialmodel_3dhybvar_opt2
    cd ..
    python ../verification/check_online_envar2.py 1.e-9 online_2D_serialmodel online_2D_serialmodel_hybvar_cg_ESTKF

    echo "------------ online_2D_serialmodel 3D-EnVar/ESTKF plain CG ---------------"
    export OMP_NUM_THREADS=1
    cd online_2D_serialmodel
    make cleandataq
    mpirun --oversubscribe -np 9 ./model_pdaf $DA_SPECS_3DEnVar -subtype 4 -type_opt 3 > ../out.online_2D_serialmodel_3denvar_opt3
    cd ..
    python ../verification/check_online_envar2.py 1.e-5 online_2D_serialmodel online_2D_serialmodel_envar_plaincg_ESTKF
#echo "   -> compare to ESTKF case online_2D_serialmodel_ESTKF - expect accuracy < 1.e-9"
#python ../verification/check_online_envar2.py 1.e-9 online_2D_serialmodel online_2D_serialmodel_ESTKF
fi


#--------- ONLINE PARALLEL -------------

if [ $RUN_ONLINE_PARALLEL -eq 1 ]
then

    echo " "
    echo "_______ online_2D_parallelmodel ________"
    echo " "
    echo "------------ online_2D_parallelmodel 3D-Var BFGS -----------------------------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel
    make cleandataq
    ./model_pdaf $DA_SPECS_3DVar  -type_opt 1 > ../out.online_2D_parallelmodel_3dv_opt1
    cd ..
    python ../verification/check_online_var2.py online_2D_parallelmodel online_2D_serialmodel_var_bfgs

    echo "------------ online_2D_parallelmodel 3D-Var CG+ parallel ---------------------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel
    make cleandataq
    mpirun -np 2 ./model_pdaf $DA_SPECS_3DVar  -type_opt 12 > ../out.online_2D_parallelmodel_3dv_opt12
    cd ..
    python ../verification/check_online_var2.py online_2D_parallelmodel online_2D_serialmodel_var_cg

    echo "------------ online_2D_parallelmodel 3D-Var plain CG parallel ----------------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel
    make cleandataq
    mpirun -np 2 ./model_pdaf $DA_SPECS_3DVar  -type_opt 13 > ../out.online_2D_parallelmodel_3dv_opt13
    cd ..
    python ../verification/check_online_var2.py online_2D_parallelmodel online_2D_serialmodel_var_cg

    echo "------------ online_2D_parallelmodel 3D-EnVar/LESTKF CG+ parallel ------------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel
    make cleandataq
    mpirun --oversubscribe -np 18 ./model_pdaf $DA_SPECS_3DEnVar -subtype 1 -type_opt 12 > ../out.online_2D_parallelmodel_3dlenvar_opt12
    cd ..
    python ../verification/check_online_envar2.py 1.e-9 online_2D_parallelmodel online_2D_serialmodel_lenvar_cg

    echo "------------ online_2D_parallelmodel 3D-EnVar/ESTKF CG+ parallel -------------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel
    make cleandataq
    mpirun --oversubscribe -np 18 ./model_pdaf $DA_SPECS_3DEnVar -subtype 4 -type_opt 12 > ../out.online_2D_parallelmodel_3denvar_opt12
    cd ..
#echo "   -> compare to ESTKF case online_2D_serialmodel_ESTKF - expect accuracy < 1.e-5"
#python ../verification/check_online_envar2.py 1.e-5 online_2D_parallelmodel online_2D_serialmodel_ESTKF
#echo "   -> compare to ESTKF case offline_2D_parallel_ESTKF - expect accuracy < 1.e-9"
    python ../verification/check_online_envar2.py 1.e-9 online_2D_parallelmodel online_2D_serialmodel_envar_cg_ESTKF

    echo "------------ online_2D_parallelmodel 3D-EnVar/LESTKF CG+ parallel ------------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel
    mpirun --oversubscribe -np 18 ./model_pdaf $DA_SPECS_3DEnVar -subtype 1 -type_opt 12 > ../out.online_2D_parallelmodel_3dlenvar_opt12
    cd ..
    python ../verification/check_online_envar2.py 1.e-9 online_2D_parallelmodel online_2D_serialmodel_lenvar_cg

    echo "------------ online_2D_parallelmodel hybrid 3D-Var/LESTKF CG+ parallel -------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel
    make cleandataq
    mpirun --oversubscribe -np 18 ./model_pdaf $DA_SPECS_hyb3DVar -subtype 6 -type_opt 12 -dim_cvec 9 > ../out.online_2D_parallelmodel_3dlhybvar_opt12
    cd ..
    python ../verification/check_online_envar2.py 1.e-9 online_2D_parallelmodel online_2D_serialmodel_lhybvar_cg

    echo "------------ online_2D_parallelmodel hybrid 3D-Var/ESTKF CG+ parallel --------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel
    make cleandataq
    mpirun --oversubscribe -np 18 ./model_pdaf $DA_SPECS_hyb3DVar -subtype 7 -type_opt 12 -dim_cvec 9 > ../out.online_2D_parallelmodel_3dhybvar_opt12
    cd ..
    python ../verification/check_online_envar2.py 1.e-9 online_2D_parallelmodel online_2D_serialmodel_hybvar_cg_ESTKF
    

    echo "------------ online_2D_parallelmodel 3D-EnVar/LESTKF plain CG parallel -------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel
    make cleandataq
    mpirun --oversubscribe -np 18 ./model_pdaf $DA_SPECS_3DEnVar -subtype 1 -type_opt 13 > ../out.online_2D_parallelmodel_3dlenvar_opt13
    cd ..
    python ../verification/check_online_envar2.py 1.e-9 online_2D_parallelmodel online_2D_serialmodel_lenvar_cg

    echo "------------ online_2D_parallelmodel 3D-EnVar/ESTKF plain CG parallel --------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel
    make cleandataq
    mpirun --oversubscribe -np 18 ./model_pdaf $DA_SPECS_3DEnVar -subtype 4 -type_opt 13 > ../out.online_2D_parallelmodel_3denvar_opt13
    cd ..
    python ../verification/check_online_envar2.py 1.e-5 online_2D_parallelmodel  online_2D_serialmodel_envar_plaincg_ESTKF
#echo "   -> compare to ESTKF case online_2D_serialmodel_ESTKF - expect accuracy < 1.e-9"
#python ../verification/check_online_envar2.py 1.e-9 online_2D_parallelmodel online_2D_serialmodel_ESTKF

    echo "------------ online_2D_parallelmodel 3D-EnVar/LESTKF plain CG parallel -------------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel
    make cleandataq
    mpirun --oversubscribe -np 18 ./model_pdaf $DA_SPECS_3DEnVar -subtype 1 -type_opt 13 > ../out.online_2D_parallelmodel_3dlenvar_opt13
    cd ..
    python ../verification/check_online_envar2.py 1.e-9 online_2D_parallelmodel online_2D_serialmodel_lenvar_cg

    echo "------------ online_2D_parallelmodel hybrid 3D-Var/LESTKF plain CG parallel --------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel
    make cleandataq
    mpirun --oversubscribe -np 18 ./model_pdaf $DA_SPECS_hyb3DVar -subtype 6 -type_opt 13 -dim_cvec 9 > ../out.online_2D_parallelmodel_3dlhybvar_opt13
    cd ..
    python ../verification/check_online_envar2.py 1.e-9 online_2D_parallelmodel online_2D_serialmodel_lhybvar_plaincg

    echo "------------ online_2D_parallelmodel hybrid 3D-Var/ESTKF plain CG parallel ---------------"
    export OMP_NUM_THREADS=1
    cd online_2D_parallelmodel
    make cleandataq
    mpirun --oversubscribe -np 18 ./model_pdaf $DA_SPECS_hyb3DVar -subtype 7 -type_opt 13 -dim_cvec 9 > ../out.online_2D_parallelmodel_3dhybvar_opt13
    cd ..
    python ../verification/check_online_envar2.py 1.e-9 online_2D_parallelmodel online_2D_serialmodel_hybvar_cg_ESTKF

fi