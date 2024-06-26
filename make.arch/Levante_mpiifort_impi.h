######################################################
# Include file with machine-specific definitions     #
# for building PDAF.                                 #
#                                                    #
# Variant for Levante@DKRZ                           #
# (Intel Fortran Compiler and Intel MPI)             #
######################################################

# Compiler, Linker, and Archiver
FC = mpiifort 
#FC = /sw/rhel6-x64/intel/impi/2017.3.196/compilers_and_libraries/linux/mpi/bin64/mpiifort 
LD = $(FC)
AR = ar
RANLIB = ranlib 

# C preprocessor
# (only required, if preprocessing is not performed via the compiler)
CPP = /usr/bin/cpp

# Definitions for CPP
# Define USE_PDAF to include PDAF
# (if the compiler does not support get_command_argument()
# from Fortran 2003 you should define F77 here.)
CPP_DEFS = -DUSE_PDAF

# Optimization specs for compiler
#   (You should explicitly define double precision for floating point
#   variables in the compilation)  
OPT= -qmkl -O3 -xHost -r8 -qopenmp

# Optimization specifications for Linker
OPT_LNK = $(OPT)

# Linking libraries (BLAS, LAPACK, if required: MPI)
# from email by Mathis and https://docs.dkrz.de/doc/levante/code-development/compiling-and-linking.html for netcdf-fortran
LINK_LIBS =-L/sw/spack-levante/intel-oneapi-mkl-2022.0.1-ttdktf/mkl/2022.0.1/lib/intel64 \
           -Wl,-rpath,/sw/spack-levante/intel-oneapi-mkl-2022.0.1-ttdktf/mkl/2022.0.1/lib/intel64 \
           -Wl,-rpath,/sw/spack-levante/netcdf-fortran-4.5.3-l2ulgp/lib


# Specifications for the archiver
AR_SPEC = 

# Specifications for ranlib
RAN_SPEC =

# Specification for directory holding modules (-module for Intel, -J for GNU)
MODULEOPT = -module

# Include path for MPI header file
MPI_INC =            # dummy MPI for PDAF versions < 2

# Object for nullMPI - if compiled without MPI library
OBJ_MPI =            # needed when running without MPI for PDAF versions < 2

# NetCDF 
# June 2022 on levante
NC_LIB   = -L/sw/spack-levante/netcdf-fortran-4.5.3-l2ulgp/lib -lnetcdff  #  gcc@11.2.0
NC_INC   = -I/sw/spack-levante/netcdf-fortran-4.5.3-l2ulgp/include


