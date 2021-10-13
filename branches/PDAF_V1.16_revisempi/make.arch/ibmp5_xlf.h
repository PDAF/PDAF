######################################################
# Include file with machine-specific definitions     #
# for building PDAF.                                 #
#                                                    #
# Variant for IBM p575 (tsuna) at AWI                #
# without MPI.                                       #
#                                                    #
# In the case of compilation without MPI, a dummy    #
# implementation of MPI, like provided in the        #
# directory nullmpi/ has to be linked when building  #
# an executable.                                     #
######################################################
# $Id: ibmp5_xlf.h 1187 2011-09-16 08:14:28Z lnerger $


# Compiler, Linker, and Archiver
FC = xlf90_r
LD = xlf90_r
AR = ar
RANLIB = ranlib 

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
OPT= -q64 -O3 -qtune=auto -qarch=pwr5 -qsuffix=f=f90 -qrealsize=8

# Optimization specifications for Linker
OPT_LNK = $(OPT)

# Linking libraries (BLAS, LAPACK, if required: MPI)
LINK_LIBS = -lessl -L/edvir1/soft/LAPACK -llapack_pwr4_64

# Specifications for the archiver
AR_SPEC = -X64

# Specifications for ranlib
RAN_SPEC =

# Include path for MPI header file
MPI_INC =  -Idummympi

# Object for nullMPI - if compiled without MPI library
OBJ_MPI = nullmpi.o

# NetCDF (only required for Lorenz96)
NC_LIB   = -L/edvir1/soft/netcdf/3.6.2/64/lib -lnetcdf
NC_INC   = -I/edvir1/soft/netcdf/3.6.2/64/include
