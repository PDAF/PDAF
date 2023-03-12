######################################################
# Include file with machine-specific definitions     #
# for building PDAF.                                 #
#                                                    #
# Variant for IBM BladeCenter (Power6) at AWI        #
# without MPI.                                       #
#                                                    #
# In the case of compilation without MPI, a dummy    #
# implementation of MPI, like provided in the        #
# directory nullmpi/ has to be linked when building  #
# an executable.                                     #
######################################################
# $Id: ibm_xlf.h 850 2010-02-02 09:53:06Z lnerger $


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
OPT= -q64 -O3 -qtune=auto -qsuffix=f=f90 -qrealsize=8

# Optimization specifications for Linker
OPT_LNK = $(OPT)

# Linking libraries (BLAS, LAPACK, if required: MPI)
LINK_LIBS = -lessl -L/iblade/soft/lapack/3.2.1 -llapack_pwr6

# Specifications for the archiver
AR_SPEC = -X64

# Specifications for ranlib
RAN_SPEC =

# Include path for MPI header file
MPI_INC =  -Idummympi

# Object for nullMPI - if compiled without MPI library
OBJ_MPI = nullmpi.o

# NetCDF (only required for Lorenz96)
NC_LIB   = -L/iblade/soft/netcdf/3.6.3/lib -lnetcdf
NC_INC   = -I/iblade/soft/netcdf/3.6.3/include
