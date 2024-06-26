######################################################
# Include file with machine-specific definitions     #
# for building PDAF.                                 #
#                                                    #
# Variant for NEC SX8-R with MPI at AWI.             #
#                                                    #
# In the case of compilation without MPI, a dummy    #
# implementation of MPI, like provided in the        #
# directory nullmpi/ has to be linked when building  #
# an executable.                                     #
######################################################
# $Id: nec-sx8_mpi.h 1187 2011-09-16 08:14:28Z lnerger $


# Compiler, Linker, and Archiver
FC = sxf90
LD = sxf90
AR = sxar
RANLIB = ls #ranlib 

# C preprocessor
# (only required, if preprocessing is not performed via the compiler)
CPP = /usr/bin/cpp

# Definitions for CPP
# Define USE_PDAF to include PDAF
# (if the compiler does not support get_command_argument()
# from Fortran 2003 you should define F77 here.)
CPP_DEFS =  -DUSE_PDAF 

# Optimization specs for compiler
#   (You should explicitly define double precision for floating point
#   variables in the compilation)  
OPT= -sxace -Wf"-A dbl" -ftrace 

# Optimization specifications for Linker
OPT_LNK = $(OPT)

# Linking libraries (BLAS, LAPACK, if required: MPI)
LINK_LIBS = -L/SX/opt/mathkeisan/MK4_0_1/lib0 -llapack -lblas

# Specifications for the archiver
AR_SPEC = 

# Specifications for ranlib
RAN_SPEC =

# Specification for directory holding modules (-module for Intel, -J for GNU)
MODULEOPT = -module

# Include path for MPI header file
MPI_INC =  -Idummympi

# Object for nullMPI - if compiled without MPI library
OBJ_MPI = nullmpi.o

# NetCDF (only required for Lorenz96)
NC_LIB   = -L/home/repo/nodes/netcdf/3.6.1_32_dw/SX/lib -lnetcdf
NC_INC   = -I/home/repo/nodes/netcdf/3.6.1_32_dw/SX/include
