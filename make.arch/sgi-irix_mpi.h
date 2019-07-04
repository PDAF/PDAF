######################################################
# Include file with machine-specific definitions     #
# for building PDAF.                                 #
#                                                    #
# Variant for SGI Origin 2000 at AWI using MPI.      #
#                                                    #
# In the case of compilation without MPI, a dummy    #
# implementation of MPI, like provided in the        #
# directory nullmpi/ has to be linked when building  #
# an executable.                                     #
######################################################
# $Id: sgi-irix_mpi.h 1565 2015-02-28 17:04:41Z lnerger $


# Compiler, Linker, and Archiver
FC = f90
LD = f90
AR = ar
RANLIB = ranlib 

# C preprocessor
# (only required, if preprocessing is not performed via the compiler)
CPP = /lib/cpp

# Definitions for CPP
# Define USE_PDAF to include PDAF
# Define BLOCKING_MPI_EXCHANGE to use blocking MPI commands to exchange data between model and PDAF
# (if the compiler does not support get_command_argument()
# from Fortran 2003 you should define F77 here.)
CPP_DEFS = -DF77

# Optimization specs for compiler
#   (You should explicitly define double precision for floating point
#   variables in the compilation)  
OPT= -r8 -O3 

# Optimization specifications for Linker
OPT_LNK = $(OPT)

# Linking libraries (BLAS, LAPACK, if required: MPI)
LINK_LIBS = -lblas -lscs -lmpi

# Specifications for the archiver
AR_SPEC = 

# Specifications for ranlib
RAN_SPEC =

# Include path for MPI header file
MPI_INC =

# Object for nullMPI - if compiled without MPI library
OBJ_MPI = 

# NetCDF (only required for Lorenz96)
NC_LIB   = 
NC_INC   = 
