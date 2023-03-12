######################################################
# Include file with machine-specific definitions     #
# for building PDAF.                                 #
#                                                    #
# Variant for Sun with Sun Studio 12 Compiler at AWI #
# without MPI.                                       #
#                                                    #
# In the case of compilation without MPI, a dummy    #
# implementation of MPI, like provided in the        #
# directory nullmpi/ has to be linked when building  #
# an executable.                                     #
######################################################
# $Id: sun_studio.h 850 2010-02-02 09:53:06Z lnerger $


# Compiler, Linker, and Archiver
FC = f90
LD = f90
AR = ar
RANLIB = ranlib 

# C preprocessor
# (only required, if preprocessing is not performed via the compiler)
CPP = /usr/ccs/lib/cpp

# Definitions for CPP
# Define USE_PDAF to include PDAF
# (if the compiler does not support get_command_argument()
# from Fortran 2003 you should define F77 here.)
CPP_DEFS = 

# Optimization specs for compiler
#   (You should explicitly define double precision for floating point
#   variables in the compilation)  
OPT= -fast -xtypemap=real:64 -m64

# Optimization specifications for Linker
OPT_LNK = $(OPT)

# Linking libraries (BLAS, LAPACK, if required: MPI)
LINK_LIBS = -xlic_lib=sunperf

# Specifications for the archiver
AR_SPEC = 

# Specifications for ranlib
RAN_SPEC =

# Include path for MPI header file
MPI_INC =  -Idummympi

# Object for nullMPI - if compiled without MPI library
OBJ_MPI = nullmpi.o

# NetCDF (only required for Lorenz96)
NC_LIB   = 
NC_INC   = 
