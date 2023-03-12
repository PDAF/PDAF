######################################################
# Include file with machine-specific definitions     #
# for building PDAF.                                 #
#                                                    #
# Variant for Linux with PGI Compiler and OpenMPI    #
# at DKRZ (tornado.dkrz.de)                          #
#                                                    #
# In the case of compilation without MPI, a dummy    #
# implementation of MPI, like provided in the        #
# directory nullmpi/ has to be linked when building  #
# an executable.                                     #
######################################################
# $Id: linux_pgi_openmpi.h 1580 2015-06-07 08:15:58Z lnerger $


# Compiler, Linker, and Archiver
FC = mpif90
LD = $(FC)
AR = ar
RANLIB = ranlib 

# C preprocessor
# (only required, if preprocessing is not performed via the compiler)
CPP = /usr/bin/cpp

# Definitions for CPP
# Define BLOCKING_MPI_EXCHANGE to use blocking MPI commands to exchange data between model and PDAF
# (if the compiler does not support get_command_argument()
# from Fortran 2003 you should define F77 here.)
CPP_DEFS = -DUSE_PDAF

# Optimization specs for compiler
#   (You should explicitly define double precision for floating point
#   variables in the compilation)  
OPT= -O3 -fast -Mr8 -Mr8intrinsics -Mbyteswapio -Mpreprocess

# Optimization specifications for Linker
OPT_LNK = $(OPT)

# Linking libraries (BLAS, LAPACK, if required: MPI)
LINK_LIBS = -L/sw/sles9-x64/acml-3.6.0/pgi64/lib -lacml -lacml_mv -Bstatic

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
