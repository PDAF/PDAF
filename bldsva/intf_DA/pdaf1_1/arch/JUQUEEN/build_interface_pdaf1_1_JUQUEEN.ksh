#! /bin/ksh
#

always_da(){
route "${cyellow}>> always_da${cnormal}"
route "${cyellow}<< always_da${cnormal}"
}

substitutions_da(){
route "${cyellow}>> substitutions_da${cnormal}"

  comment "   mkdir  $dadir/interface"
    mkdir -p $dadir/interface  >> $log_file 2>> $err_file
  check

  comment "   cp pdaf interface model to $dadir/interface"
    patch $rootdir/bldsva/intf_DA/pdaf1_1/model $dadir/interface 
  check

  comment "   cp pdaf interface framework to $dadir/interface"
    patch $rootdir/bldsva/intf_DA/pdaf1_1/framework $dadir/interface 
  check

route "${cyellow}<< substitutions_da${cnormal}"
}

configure_da(){
route "${cyellow}>> configure_da${cnormal}"
  export PDAF_DIR=$dadir
  export PDAF_ARCH=ibm_xlf_mpi

#PDAF part
  file=$dadir/make.arch/ibm_xlf_mpi.h
  
  comment "   cp pdaf config to $dadir"
    cp $rootdir/bldsva/intf_DA/pdaf1_1/arch/$platform/config/ibm_xlf_mpi.h $file >> $log_file 2>> $err_file
  check

  comment "   sed MPI dir to $file" 
    sed -i "s@__mpidir__@${mpiPath}@" $file >> $log_file 2>> $err_file  
  check

  comment "   sed LIBS to $file"
    sed -i "s@__LIBS__@ -L${lapackPath}/lib  -L/bgsys/local/lib -llapack -lesslbg @" $file >> $log_file 2>> $err_file
  check

  comment "   sed optimizations to $file"
    sed -i "s@__OPT__@${optComp}@" $file >> $log_file 2>> $err_file
  check

  comment "   cd to $dadir/src"
    cd $dadir/src >> $log_file 2>> $err_file
  check
  comment "   make clean pdaf"
    make clean >> $log_file 2>> $err_file
  check



#PDAF interface part
  file1=$dadir/interface/model/Makefile
  file2=$dadir/interface/framework/Makefile
  comment "   cp pdaf interface Makefiles to $dadir"
    cp $rootdir/bldsva/intf_DA/pdaf1_1/model/Makefile  $file1 >> $log_file 2>> $err_file
  check
    cp $rootdir/bldsva/intf_DA/pdaf1_1/framework/Makefile  $file2 >> $log_file 2>> $err_file
  check

  importFlags=" "
  importFlagsOAS=" "
  importFlagsPFL=" "
  importFlagsCLM=" "
  importFlagsCOS=" "
  importFlagsDA=" "
  cppdefs=" "
  fcppdefs=" "
  obj=' ' 
  libs=" "
  pf="-WF,"
 
  # Oasis include dirs
  importFlagsOAS+="-I$oasdir/$platform/build/lib/psmile.MPI1 "
  importFlagsOAS+="-I$oasdir/$platform/build/lib/scrip "

  # CLM include dirs
  importFlagsCLM+="-I$clmdir/build/ "

  # COSMO include dirs
  importFlagsCOS+="-I$cosdir/obj "

  # ParFlow include dirs
  importFlagsPFL+="-I$pfldir/pfsimulator/parflow_lib "
  importFlagsPFL+="-I$pfldir/pfsimulator/amps/oas3 "
  importFlagsPFL+="-I$pfldir/pfsimulator/amps/common "
  importFlagsPFL+="-I$pfldir/pfsimulator/include "
  
  # DA include dirs
  importFlagsDA+="-I$dadir/interface/model/common "
  if [[ $withPFL == "true" ]] ; then
    importFlagsDA+="-I$dadir/interface/model/${mList[3]} "
  fi

  if [[ $withOAS == "false" && $withPFL == "true" ]] ; then
     importFlags+=$importFlagsPFL
     importFlags+=$importFlagsDA
     cppdefs+=" -DPARFLOW_STAND_ALONE "
     fcppdefs+=" ${pf}-DPARFLOW_STAND_ALONE "
     libs+=" -L$hyprePath/lib -L$siloPath/lib -lparflow -lamps -lamps_common -lamps -lamps_common -lkinsol -lgfortran -lHYPRE -lsilo "
     obj+=' $(OBJPF) '
  fi

  if [[ $withOAS == "false" && $withCLM == "true" ]] ; then
     importFlags+=$importFlagsCLM
     importFlags+=$importFlagsDA
     fcppdefs+=" ${pf}-DCLMSA "
     libs+=" -lclm "
     obj+=' $(OBJCLM) print_update_clm.o'
  fi

  if [[ $withCLM == "true" && $withCOS == "true" && $withPFL == "false" ]] ; then
     importFlags+=$importFlagsCLM
     importFlags+=$importFlagsOAS
     importFlags+=$importFlagsCOS
     importFlags+=$importFlagsDA
     fcppdefs+=" ${pf}-Duse_comm_da ${pf}-DCOUP_OAS_COS ${pf}-DGRIBDWD ${pf}-DNETCDF ${pf}-DHYMACS ${pf}-DMAXPATCH_PFT=1 "
     if [[ $cplscheme == "true" ]] ; then ; fcppdefs+=" ${pf}-DCPL_SCHEME_F " ; fi
     if [[ $readCLM == "true" ]] ; then ; fcppdefs+=" ${pf}-DREADCLM " ; fi
     if [[ ${mList[2]} == "cosmo5_1" ]] ; then
       libs+=" -lclm -lcosmo -lpsmile.MPI1 -lmct -lmpeu -lscrip -L$gribPath/lib/ -leccodes_f90 -leccodes"
     else
       libs+=" -lclm -lcosmo -lpsmile.MPI1 -lmct -lmpeu -lscrip $grib1Path/libgrib1.a "
     fi
     obj+=' $(OBJCLM) $(OBJCOSMO) '
  fi 

  if [[ $withCLM == "true" && $withCOS == "false" && $withPFL == "true" ]] ; then
     importFlags+=$importFlagsCLM
     importFlags+=$importFlagsOAS
     importFlags+=$importFlagsPFL
     importFlags+=$importFlagsDA
     cppdefs+=" -Duse_comm_da -DMAXPATCH_PFT=1 -DCOUP_OAS_PFL "
     fcppdefs+=" ${pf}-Duse_comm_da ${pf}-DCOUP_OAS_PFL ${pf}-DMAXPATCH_PFT=1 "
     if [[ $readCLM == "true" ]] ; then ; cppdefs+=" -DREADCLM " ; fcppdefs+=" ${pf}-DREADCLM " ; fi
     if [[ $freeDrain == "true" ]] ; then ; cppdefs+=" -DFREEDRAINAGE " ; fcppdefs+=" ${pf}-DFREEDRAINAGE " ; fi
     libs+=" -lclm -lpsmile.MPI1 -lmct -lmpeu -lscrip -L$hyprePath/lib -L$siloPath/lib -lparflow -lamps -lamps_common -lamps -lamps_common -lkinsol -lgfortran -lHYPRE -lsilo "
     obj+=' $(OBJCLM) $(OBJPF) '
  fi
  if [[ $withCLM == "true" && $withCOS == "true" && $withPFL == "true" ]] ; then
     importFlags+=$importFlagsCLM
     importFlags+=$importFlagsOAS
     importFlags+=$importFlagsPFL
     importFlags+=$importFlagsCOS
     importFlags+=$importFlagsDA
     cppdefs+=" -Duse_comm_da -DCOUP_OAS_COS -DGRIBDWD -DNETCDF -DHYMACS -DMAXPATCH_PFT=1 -DCOUP_OAS_PFL "
     fcppdefs+=" ${pf}-Duse_comm_da ${pf}-DCOUP_OAS_COS ${pf}-DGRIBDWD ${pf}-DNETCDF ${pf}-DHYMACS ${pf}-DMAXPATCH_PFT=1 ${pf}-DCOUP_OAS_PFL "
     if [[ $cplscheme == "true" ]] ; then ; cppdefs+=" -DCPL_SCHEME_F " ; fcppdefs+=" ${pf}-DCPL_SCHEME_F " ; fi
     if [[ $readCLM == "true" ]] ; then ; cppdefs+=" -DREADCLM " ; fcppdefs+=" ${pf}-DREADCLM " ; fi
     if [[ $freeDrain == "true" ]] ; then ; cppdefs+=" -DFREEDRAINAGE " ; fcppdefs+=" ${pf}-DFREEDRAINAGE " ; fi
     if [[ ${mList[2]} == "cosmo5_1" ]] ; then
       libs+=" -lclm -lpsmile.MPI1 -lmct -lmpeu -lscrip -L$gribPath/lib/      -L$hyprePath -L$siloPath -lparflow -lamps -lamps_common -lamps -lamps_common -lkinsol -lgfortran -lHYPRE -lsilo -lcosmo -leccodes_f90 -leccodes"
     else
       libs+=" -lclm -lpsmile.MPI1 -lmct -lmpeu -lscrip $grib1Path/libgrib1.a -L$hyprePath -L$siloPath -lparflow -lamps -lamps_common -lamps -lamps_common -lkinsol -lgfortran -lHYPRE -lsilo -lcosmo $grib1Path/libgrib1.a"
     fi
     obj+=' $(OBJCLM) $(OBJCOSMO) $(OBJPF) '
  fi

  libs+=" -L$mpiPath/lib -lmpich -L$ncdfPath/lib/ -lnetcdf "

  comment "   sed bindir to Makefiles"
    sed -i "s,__bindir__,$bindir," $file1 $file2 >> $log_file 2>> $err_file
  check
  comment "   sed comp flags to Makefiles"
    sed -i "s,__fflags__, -qfree=f90 -qextname=flush -qfullpath -I$dadir/interface/model -I$ncdfPath/include $importFlags," $file1 $file2 >> $log_file 2>> $err_file
  check
    sed -i "s,__ccflags__, -I$dadir/interface/model -I$ncdfPath/include $importFlags," $file1 $file2 >> $log_file 2>> $err_file
  check
  comment "   sed preproc flags to Makefiles"
    sed -i "s@__cpp_defs__@$cppdefs@" $file1 $file2 >> $log_file 2>> $err_file
  check
    sed -i "s@__fcpp_defs__@$fcppdefs@" $file1 $file2 >> $log_file 2>> $err_file
  check
  comment "   sed libs to Makefiles"
    sed -i "s,__libs__,$libs," $file2 >> $log_file 2>> $err_file
  check
  comment "   sed obj to Makefiles"
    sed -i "s,__obj__,$obj," $file1 >> $log_file 2>> $err_file
  check
  comment "   sed -D prefix to Makefiles"
    sed -i "s@__pf__@$pf@" $file1 $file2 >> $log_file 2>> $err_file
  check
  comment "   sed clm directory to Makefiles"
    sed -i "s,__clmdir__,${mList[1]}," $file1 $file2 >> $log_file 2>> $err_file
  check
  comment "   sed cosmo directory to Makefiles"
    sed -i "s,__cosdir__,${mList[2]}," $file1 $file2 >> $log_file 2>> $err_file
  check
  comment "   sed parflow directory to Makefiles"
    sed -i "s,__pfldir__,${mList[3]}," $file1 $file2 >> $log_file 2>> $err_file
  check


  comment "   cd to $dadir/interface/model"
    cd $dadir/interface/model >> $log_file 2>> $err_file
  check
  comment "   make clean model"
    make clean >> $log_file 2>> $err_file
  check
  comment "   cd to $dadir/src/interface/framework"
    cd $dadir/interface/framework >> $log_file 2>> $err_file
  check
  comment "   make clean framework"
    make clean >> $log_file 2>> $err_file
  check


route "${cyellow}<< configure_da${cnormal}"
}

make_da(){
route "${cyellow}>> make_da${cnormal}"
  export PDAF_DIR=$dadir
  export PDAF_ARCH=ibm_xlf_mpi

  comment "   cd to $dadir/src"
    cd $dadir/src >> $log_file 2>> $err_file
  check
  comment "   make pdaf"
    make >> $log_file 2>> $err_file
  check

  comment "   cd to $dadir/interface/model"
    cd $dadir/interface/model >> $log_file 2>> $err_file
  check
  comment "   make pdaf model"
    make >> $log_file 2>> $err_file
  check

  comment "   cd to $dadir/interface/framework"
    cd $dadir/interface/framework >> $log_file 2>> $err_file
  check
  comment "   make pdaf framework"
    make >> $log_file 2>> $err_file
  check


route "${cyellow}<< make_da${cnormal}"
}


setup_da(){
route "${cyellow}>> setup_da${cnormal}"
  c_setup_pdaf
route "${cyellow}<< setup_da${cnormal}"
}

