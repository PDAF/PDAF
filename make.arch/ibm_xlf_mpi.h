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
# $Id: ibm_xlf_mpi.h 1565 2015-02-28 17:04:41Z lnerger $


# Compiler, Linker, and Archiver
FC = mpxlf90_r
LD = mpxlf90_r
AR = ar
RANLIB = ranlib 

# C preprocessor
# (only required, if preprocessing is not performed via the compiler)
CPP = /usr/ccs/lib/cpp

# Definitions for CPP
# Define USE_PDAF to include PDAF
# Define BLOCKING_MPI_EXCHANGE to use blocking MPI commands to exchange data between model and PDAF
# (if the compiler does not support get_command_argument()
# from Fortran 2003 you should define F77 here.)
CPP_DEFS = -WF,-DUSE_PDAF

# Optimization specs for compiler
#   (You should explicitly define double precision for floating point
#   variables in the compilation)  
OPT= -q64 -O3 -qtune=auto -qarch=pwr4 -qsuffix=f=f90 -qrealsize=8

# Optimization specifications for Linker
OPT_LNK = $(OPT)

# Linking libraries (BLAS, LAPACK, if required: MPI)
LINK_LIBS = -lessl -L/iblade/soft/lapack/3.2.1 -llapack_pwr6

# Specifications for the archiver
AR_SPEC = -X64

# Specifications for ranlib
RAN_SPEC =

# Specification for directory holding modules (-module for Intel, -J for GNU)
MODULEOPT = -I_DONT_KNOW

# Include path for MPI header file
MPI_INC =  -I/usr/lpp/ppe.poe/include

# Object for nullMPI - if compiled without MPI library
OBJ_MPI =

# NetCDF (only required for Lorenz96)
NC_LIB   = 
NC_INC   = 
