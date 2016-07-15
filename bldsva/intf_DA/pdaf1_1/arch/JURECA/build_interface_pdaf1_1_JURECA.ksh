#! /bin/ksh
#

always_da(){
route "${cblue}>> always_da${cnormal}"
route "${cblue}<< always_da${cnormal}"
}

substitutions_da(){
route "${cblue}>> substitutions_da${cnormal}"
route "${cblue}<< substitutions_da${cnormal}"
}

configure_da(){
route "${cblue}>> configure_da${cnormal}"
  export PDAF_DIR=$dadir
  export PDAF_ARCH=linux_ifort_jureca

  comment "   cd to $dadir/src"
    cd $dadir/src >> $log_file 2>> $err_file
  check
  comment "   make clean pdaf"
    make clean >> $log_file 2>> $err_file
  check

  file=$dadir/make.arch/linux_ifort_jureca.h

  comment "   cp pdaf config to $dadir"
    cp $rootdir/bldsva/intf_DA/pdaf1_1/arch/$platform/config/linux_ifort_jureca.h $file >> $log_file 2>> $err_file
  check

  comment "   sed MPI dir to $file" 
    sed -i "s@__MPI_INC__@-I${mpiPath}/include@" $file >> $log_file 2>> $err_file  
  check

  comment "   sed LIBS to $file"
    sed -i "s@__LIBS__@ $lapackPath/mkl/lib/intel64/libmkl_intel_lp64.a $lapackPath/mkl/lib/intel64/libmkl_intel_thread.a $lapackPath/mkl/lib/intel64/libmkl_core.a -L${mpiPath}/lib64@" $file >> $log_file 2>> $err_file
  check

  comment "   sed optimizations to $file"
    sed -i "s@__OPT__@${optComp}@" $file >> $log_file 2>> $err_file
  check

route "${cblue}<< configure_da${cnormal}"
}

make_da(){
route "${cblue}>> make_da${cnormal}"
  export PDAF_DIR=$dadir
  export PDAF_ARCH=linux_ifort_jureca

  comment "   cd to $dadir/src"
    cd $dadir/src >> $log_file 2>> $err_file
  check
  comment "   make pdaf"
    make >> $log_file 2>> $err_file
  check
route "${cblue}<< make_da${cnormal}"
}


setup_da(){
route "${cblue}>> setup_da${cnormal}"
route "${cblue}<< setup_da${cnormal}"
}

