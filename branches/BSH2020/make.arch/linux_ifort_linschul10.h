######################################################
# Include file with machine-specific definitions     #
# for building PDAF.                                 #
#                                                    #
# Variant for Linux with Intel Fortran Compiler      #
# without MPI at AWI                                 #
#                                                    #
# In the case of compilation without MPI, a dummy    #
# implementation of MPI, like provided in the        #
# directory nullmpi/ has to be linked when building  #
# an executable.                                     #
######################################################
# $Id: linux_ifort.h 1395 2013-05-03 13:44:37Z lnerger $

# Compiler, Linker, and Archiver
FC = mpiifort
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
OPT= -O2 -xHost -r8 -fp-model precise -mkl -traceback 

# Optimization specifications for Linker
OPT_LNK = $(OPT)

# Linking libraries (BLAS, LAPACK, if required: MPI)
LINK_LIBS = -L/software/Intel/compilers_and_libraries_2018.3.222/linux/mkl/lib/intel64/libmkl_intel_thread.a -L/software/intel/compilers_and_libraries_2018.3.222/linux/mkl/lib/intel64/libmkl_core.a 



# Specifications for the archiver
AR_SPEC = 

# Specifications for ranlib
RAN_SPEC =

# Include path for MPI header file
MPI_INC =  #-I/home/share/mpich3_ifort/include

# Object for nullMPI - if compiled without MPI library
OBJ_MPI = #nullmpi.o

# NetCDF (only required for Lorenz96)
NC_LIB   =-L/usr/local/netcdf-fortran-4.4.4-ifort/lib -lnetcdf -lnetcdff
NC_INC   =-I/usr/local/netcdf-fortran-4.4.4-ifort/include
