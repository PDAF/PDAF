######################################################
# Include file with machine-specific definitions     #
# for building PDAF.                                 #
#                                                    #
# Variant for HLRN-III                               #
# (Cray XC30 with CCE, LibSci, and mpich)            #
#                                                    #
# In the case of compilation without MPI, a dummy    #
# implementation of MPI, like provided in the        #
# directory nullmpi/ has to be linked when building  #
# an executable.                                     #
######################################################
# $Id: cray_cce_mpi.h 1645 2016-08-30 11:52:32Z lnerger $


# Compiler, Linker, and Archiver
FC = ftn
LD = $(FC)
AR = ar
RANLIB = ranlib 

# C preprocessor
# (only required, if preprocessing is not performed via the compiler)
CPP = /usr/bin/cpp

# Definitions for CPP
# Define USE_PDAF to include PDAF
# Define BLOCKING_MPI_EXCHANGE to use blocking MPI commands to exchange data between model and PDAF
# (if the compiler does not support get_command_argument()
# from Fortran 2003 you should define F77 here.)
CPP_DEFS = -DUSE_PDAF

# Optimization specs for compiler
#   (You should explicitly define double precision for floating point
#   variables in the compilation)  
OPT= -sreal64 -hnoomp

# Optimization specifications for Linker
OPT_LNK = $(OPT)

# Linking libraries (BLAS, LAPACK, if required: MPI)
LINK_LIBS =

# Specifications for the archiver
AR_SPEC = 

# Specifications for ranlib
RAN_SPEC =

# Include path for MPI header file
MPI_INC =  

# Object for nullMPI - if compiled without MPI library
OBJ_MPI = 

# NetCDF (only required for Lorenz96)
NC_LIB   = -L${NETCDF_DIR}/lib
NC_INC   = -I${NETCDF_DIR}/include
