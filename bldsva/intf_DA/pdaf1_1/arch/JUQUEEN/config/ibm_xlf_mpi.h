######################################################
# Include file with machine-specific definitions     #
# for building PDAF.                                 #
#                                                    #
# Variant for IBM BladeCenter (Power6) at AWI        #
# with MPI parallelization!                          #
#                                                    #
# In the case of compilation without MPI, a dummy    #
# implementation of MPI, like provided in the        #
# directory nullmpi/ has to be linked when building  #
# an executable.                                     #
######################################################
# $Id: ibm_xlf_mpi.h 970 2010-04-16 13:40:54Z lnerger $


# Compiler, Linker, and Archiver
AR = ar
RANLIB = ranlib 
MPIDIR  = __mpidir__
CC      = __comCC__
FC      = __comFC__
LD      = $(FC)


# C preprocessor
# (only required, if preprocessing is not performed via the compiler)
CPP = /usr/ccs/lib/cpp

# Definitions for CPP
# Define USE_PDAF to include PDAF
# (if the compiler does not support get_command_argument()
# from Fortran 2003 you should define F77 here.)
CPP_DEFS = -WF,-DUSE_PDAF

# Optimization specs for compiler
#   (You should explicitly define double precision for floating point
#   variables in the compilation)  
OPT= -q64 __OPT__ -qflag=i:x -qsuffix=f=f90 -qrealsize=8

# Optimization specifications for Linker
OPT_LNK = $(OPT)

# Linking libraries (BLAS, LAPACK, if required: MPI)
#LINK_LIBS = -lessl
#LINK_LIBS  = -lesslbg
LINK_LIBS = -Wl,-allow-multiple-definition __LIBS__

# Specifications for the archiver
#AR_SPEC = -X64

# Specifications for ranlib
RAN_SPEC =

# Include path for MPI header file
MPI_INC = __MPI_INC__

# Object for nullMPI - if compiled without MPI library
OBJ_MPI =

# NetCDF (only required for Lorenz96)
NC_LIB   = 
NC_INC   = 

CC += $(OPT)
FC += $(OPT)
LD += $(OPT)
