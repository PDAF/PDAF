#!/bin/bash  -x
#SBATCH --nodes=6
#SBATCH --ntasks=576
#SBATCH --ntasks-per-node=96
#SBATCH --time=00:30:00
#SBATCH --partition=standard96

# *** start of job script ***
ulimit -s unlimited

export NOPP=1

echo $NCPUS
source ~/fesom_echam6_oasis3-mct/env.sh

module load intel/18.0.6
module load impi/2018.5

#export I_MPI_STATS=1-20

export OMP_NUM_THREADS=1

# for mpirun 
export SLURM_CPU_BIND=none
export I_MPI_HYDRA_BRANCH_COUNT=-1
export I_MPI_HYDRA_TOPOLIB=ipl
export I_MPI_FABRICS=shm:ofi
export I_MPI_OFI_PROVIDER=psm2

# Ensmeble size
export NENS=2

# Whether pairs are used
export pair=T

###
read -r -d '' AWI_FESOM_YAML << 'EOF_AWI_FESOM_YAML'
---
output_schedules:
  - vars: [evs,fsitherm,hfds,opottemptend,pbo,prlq,prsn,rsdo,siarean,siareas,sidmassevapsubl,sidmasssi,sidmassth,sidmasstranx,sidmasstrany,siextentn,siextents,sifllatstop,sisnconc,sisnmass,sisnthick,sispeed,sivol,sivoln,sivols,soga,thetaoga,u2o,uo,uso,uto,v2o,vo,volo,vso,vto,w2o,wfo,wo,wso,wto,zossq,rho,uhice,uhsnow,urho,uv,vhice,vhsnow,virtual_salt,vrho,mlotst,omldamax,sic,sistrxdtop,sistrxubot,sistrydtop,sistryubot,sithick,sitimefrac,siu,siv,so,sos,tauuo,tauvo,thetao,tos,zos,flice,wnet,evap,runoff,thdgr,thdgrsn,tso]
    unit: d
    rate: 1
  - vars: [restart]
    unit: m
    first: 12
    rate: 1
  - vars: [lwrd,olat,olwout,osen,relax_salt,shum,tair,uwind,vwind]
    unit: y
    first: 1000
    rate: 1
EOF_AWI_FESOM_YAML
export AWI_FESOM_YAML
###


# ======================================= USER SECTION ========================================
# ***** time extents *****
export BEGIN_YEAR=2016
export BEGIN_MONTH=01
export BEGIN_DAY=01

export END_YEAR=2700
export END_MONTH=01
export END_DAY=01
# only one of the following three should not be zero
export RESTART_YEAR=0
export RESTART_MONTH=0
export RESTART_DAY=10

# ***** Platform juwels  *****
export EXPID=SCEN
# *****  Base Directories **********
export BASE_DIR=/home/${USER}/fesom_echam6_oasis3-mct.GMD   # root path to the model directories
export WORK_DIR=/scratch/usr/${USER}/AWI-CM-PDAF        # root path to the output directories
export ENS_DIR=${WORK_DIR}/gmd_scaling2                      # output for all erstarts and model logs

mkdir $ENS_DIR
for((ENS=1;ENS<=$NENS;ENS++))
do
  ENSstr=`printf %02d $ENS`
  echo $ENSstr
  export RUN_DIR=${ENS_DIR}/${ENSstr}'/'                # output for all erstarts and model logs
  mkdir $RUN_DIR
  #ECHAM resolution
  export RES=T63
  export LEVELS=L47		# number of levels
  export OCERES=CORE2		# GR15 or TR01
                                # used for lsm
  #ECHAM relevant directories
  export CEXPER=TST
  #FESOM relevant directories
  export FESOM_meshpath=/scratch/usr/hbkqtang/old_work2/input/CORE2_final/
  export FESOM_opbndpath=
  export FESOM_climatedatapath=/scratch/usr/hbkqtang/old_work2/input/hydrography/
  export FESOM_forcingdatapath=
  export FESOM_restartflag=last
  export FESOM_tideforcingpath=
  export FESOM_resultpath=${RUN_DIR}
  # use FESOM spinup to initialize the start
  export FESOM_spinup=T
  export FESOM_spinup_path=/scratch/usr/hbkqtang/old_work2/fesom_spinup_2016/
  export FESOM_spinup_year=2016

  # use ECHAM spinup to initialize the start
  export ECHAM_spinup=F
  export ECHAM_spinup_path=/scratch/usr/hbkqtang/old_work2/echam_spinup_2016/
  export fdate=2000

  # use OASIS spinup to initialize the start
  export OASIS_restart=T
  export OASIS_restart_path=/scratch/usr/hbkqtang/old_work2/input/

  # POSTPROCESSING
  # if the postprocessing to be applied to the output 
  # 1. set POSTPROCESS to T
  # 2. modify the part at the end of the script if needed
  export POSTPROCESS=F
  export POSTPROCESS_PATH=${WORK_DIR}/test2_output/

  # ***** ECHAM process topology (NPROCA > NPROCB) *****
  export NPROCA=12 # num. cpus for lat-dimension (must be <= nlat/2)
  #export NPROCA=12 # num. cpus for lat-dimension (must be <= nlat/2)
  export NPROCB=8  # num. cpus for lon-dimension

  # ***** FESOM process number *****
  export NCPUOC=192

  # ***** Time step FESOM *****
  export OCESTP=900

  # ***** Time step ECHAM *****
  export ATMSTP=450

  # ***** Coupling Time step *****
  export CPL_FREQ=3600 #1800 #3600

  # ***** Define grid dimensions for ATM and OCE
  if test "${RES}" = T63
  then
      nxa=192
      nya=96
  fi
  if test "${RES}" = T127
  then
      nxa=384
      nya=192
  fi
  if test "${OCERES}" = CORE2
  then
      nxo=126859
      nyo=1
  fi
  if test "${OCERES}" = AGUV
  then
      nxo=810471
      nyo=1
  fi
  if test "${OCERES}" = GLOB
  then
      nxo=830305
      nyo=1
  fi
  if test "${OCERES}" = REF87K
  then
      nxo=86803
      nyo=1
  fi
  if test "${OCERES}" = BOLD
  then
      nxo=1306775
      nyo=1
  fi

  # RESTART (comment all if initial run)
  #if test "x${FIRSTRUN}" != "xFALSE"; then
  #export FIRSTRUN="FALSE"
  #echo $FIRSTRUN 
  ## pick the following data from any clock up (fesom.clock or echam.namelist)
  #export SCRIPT_END_YEAR=2016
  #export SCRIPT_END_MONTH=01
  #export SCRIPT_END_DAY=01
  #fi

  export INI_DATA=/scratch/usr/hbkqtang/old_work2/input/pool-data/
  export INIECH=${INI_DATA}/ECHAM6/input/r0004
  export INIJSB=${INI_DATA}/JSBACH/input/r0007
  export INIECHFES=$FESOM_meshpath/tarfiles${RES}/input/echam6
  export INIJSBFES=$FESOM_meshpath/tarfiles${RES}/input/jsbach

  export CDO_PATH=
  export NCO_PATH=
  #####################################################################
  #
  # ==================================== END OF USER SECTION ====================================
  
  source ${BASE_DIR}/run_scripts/include_oasis3mct/cpl_settings.inc
  if test "${ENS}" = 1; then
        source ${BASE_DIR}/run_scripts/include_oasis3mct/cpl_timeop.inc
  fi

  if [[ $pair != "T" ]]; then
      export nproc_exe1=$(( $NCPUOC * $NENS ))
      export nproc_exe2=$(( $NCPUATM * $NENS ))
  fi

  echo 'NCPU for Ocean: ' $nproc_exe1
  echo 'NCPU for Atmosphere: ' $nproc_exe2

  echo current script starts at year $SCRIPT_BEGIN_YEAR, month: $SCRIPT_BEGIN_MONTH, day: $SCRIPT_BEGIN_DAY
  echo current script ends   at year $SCRIPT_END_YEAR,   month: $SCRIPT_END_MONTH,   day: $SCRIPT_END_DAY
  export whichloop=${SCRIPT_END_YEAR}${SCRIPT_END_MONTH}${SCRIPT_END_DAY}
 
    # Actually below two lines is done in cpl_initial.inc, however, we need to do this manually.
  cd $RUN_DIR
  
  export LRESUME=T
  export OASIS_seqmode=concurrent 
  #testing----------------------------20190605
  #source ${BASE_DIR}run_scripts/include_oasis3mct/cpl_initial.inc.ucr
  ln -sf $FESOM_BASE_DIR/fesom.x fesom.x
  ln -sf $ECHAM6_BASE_DIR/bin/echam6 echam.x
  cp ${BASE_DIR}/run_scripts/weightgen/rmp* .

  cp $NAMDIR_FEOM/namelist.* .
  cp $NAMDIR_FEOM/cf_name_table.txt .

  if [ $OCERES == "CORE2" ] ; then
  cp $NAMDIR_FEOM/namelist.config.core namelist.config
  fi



#-----------------------------------------------------------------------------

fesom_startyear=$((BEGIN_YEAR))
export oce_stp_pday=$(( 86400 / $OCESTP ))

sed -e "s@<stp_per_day>@${oce_stp_pday}@g" \
    -e "s@<meshpath>@\'${FESOM_meshpath}\'@g" \
    -e "s@<opbndpath>@\'${FESOM_opbndpath}\'@g" \
    -e "s@<climatedatapath>@\'${FESOM_climatedatapath}\'@g" \
    -e "s@<forcingdatapath>@\'${FESOM_forcingdatapath}\'@g" \
    -e "s@<tideforcingpath>@\'${FESOM_tideforcingpath}\'@g" \
    -e "s@<resultpath>@\'${FESOM_resultpath}\'@g" \
    -e "s@<restart>@${fesom_restart}@g" \
    -e "s@<restart_unit>@${fesom_restart_unit}@g" \
    -e "s@<startyear>@${fesom_startyear}@g" \
    -e "s@<restartflag>@\'${FESOM_restartflag}\'@g" \
    namelist.config > toto
    mv toto namelist.config

  sed -i `grep -n ResultPath= namelist.config|cut -d ':' -f 1`"c ResultPath='${RUN_DIR}'" namelist.config

ln -s  ${INIECH}/rrtmg_lw.nc rrtmg_lw.nc
ln -s  ${INIECH}/rrtmg_sw.nc rrtmg_sw.nc
ln -s  ${INIECH}/ECHAM6_CldOptProps.nc ECHAM6_CldOptProps.nc
ln -s ${INIJSBFES}/jsbach_${RES}${OCERES}_11tiles_5layers_1976.nc jsbach.nc
ln -s ${ECHAM6_BASE_DIR}/util/lctlib_nlct21.def lctlib.def
ln -s  ${INIJSB}/HD/hdpara.nc hdpara.nc
ln -s  ${INIJSB}/HD/hdstart.nc hdstart.nc
ln -s  ${INIJSB}/HD/rmp_${RES}_to_hd.nc rmp_hd.nc
ln -s  ${INIECH}/${RES}/${RES}${LEVELS}_jan_spec.nc unit.23
ln -s  ${INIECHFES}/${RES}${OCERES}_jan_surf.nc unit.24
ln -s  ${INIECHFES}/${RES}${OCERES}_VLTCLIM.nc unit.90
ln -s  ${INIECHFES}/${RES}${OCERES}_VGRATCLIM.nc unit.91
ln -s  ${INIECH}/${RES}/${RES}_TSLCLIM.nc unit.92

  source ${BASE_DIR}/run_scripts/include_oasis3mct/cpl_pdaf.inc.ucr
  source ${BASE_DIR}/run_scripts/include_oasis3mct/cpl_restart.inc.ucr



exe1=fesom.x
exe2=echam.x
  
  if [[ $OASIS_restart = "T" ]]; then
      cp ${ECHAM_spinup_path}a2o_flux.201611 ${RUN_DIR}a2o_flux
      cp ${ECHAM_spinup_path}o2a_flux.201611 ${RUN_DIR}o2a_flux
  fi

  cp ${FESOM_spinup_path}fesom.2015.oce.nc  ${RUN_DIR}fesom.2015.oce.nc
  cp ${FESOM_spinup_path}fesom.2015.ice.nc  ${RUN_DIR}fesom.2015.ice.nc
  cp ${ECHAM_spinup_path}fesom.clock  .
  cp ${ECHAM_spinup_path}*restart*  .       
  cp ${ECHAM_spinup_path}*hdrestart*  .
done

cd $ENS_DIR

echo run started at realtime: 
date
echo `pwd`

# create MPMD configuration file
if [ -e mpmd.conf ];then
  rm mpmd.conf
fi
touch mpmd.conf

for((i=1;i<=$NENS;i++))
do
  ENSstr=`printf %02d $i`

  echo "-np $nproc_exe1 -wdir ${ENS_DIR}/${ENSstr} ${ENS_DIR}/${ENSstr}/fesom.x " >> mpmd.conf
  echo "-np $nproc_exe2 -wdir ${ENS_DIR}/${ENSstr} ${ENS_DIR}/${ENSstr}/echam.x " >> mpmd.conf
done
cat mpmd.conf

# Wait for 30 seconds to ensure that files are fully synced
sleep 30

# Now run
mpirun --configfile ${ENS_DIR}/mpmd.conf &> ${OCERES}-${RES}.out.${SLURM_JOB_ID}

export EXITSTATUS=$?
echo run ended at realtime: 
date
echo with exitstatus ${EXITSTATUS}


for((i=1;i<=$NENS;i++))
do
  ENSstr=`printf %02d $i`
  cd ${ENS_DIR}/${ENSstr}/
  source ${BASE_DIR}/run_scripts/include_oasis3mct/cpl_isbefore.inc
done

export FIRSTRUN="FALSE"

if test "x${is_before}" = "xFALSE"; then
exit
fi

