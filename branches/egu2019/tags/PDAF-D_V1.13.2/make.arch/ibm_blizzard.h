######################################################
# Include file with machine-specific definitions     #
# for building PDAF with the dummy model example.    #
#                                                    #
# Variant for blizzard at DKRZ 			     #
#                                                    #
# Thanks to Jan Saynisch, GFZ Potsdam, for this.     #
#                                                    #
######################################################
#  $Id: ibm_blizzard.h 854 2010-02-02 13:56:53Z lnerger $

# Compiler, Linker, and Archiver
FC = mpxlf
LD = mpxlf
AR = ar
RANLIB = ranlib 

# C preprocessor
# (only required, if preprocessing is not performed via the compiler)
CPP = /usr/bin/cpp

# Definitions for CPP
# Define USE_PDAF to activate PDAF
# (if the compiler does not support get_command_argument()
# from Fortran 2003 you should define F77 here.)
CPP_DEFS =  -WF,-DUSE_PDAF


# Optimization specs for compiler
OPT= -q64 -O3 -qrealsize=8 -qextname -qlibessl

# Optimization specifications for Linker
OPT_LNK = $(OPT)

# Linking libraries (BLAS, LAPACK, if required: MPI)
LINK_LIBS = -L/sw/aix53/lapack-3.2.0/lib -llapack -lessl  

# Specifications for the archiver
AR_SPEC =  -X64

# Include path for MPI header file
MPI_INC = 

# Object for nullMPI - if compiled without MPI library
OBJ_MPI= 

# NetCDF (only required for Lorenz96)
NC_LIB   = 
NC_INC   = 

