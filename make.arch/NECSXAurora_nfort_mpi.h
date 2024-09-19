######################################################
# Include file with machine-specific definitions     #
# for building PDAF.                                 #
#                                                    #
# Variant for NEC Aurora with nfort and MPI          #
######################################################

# Compiler, Linker, and Archiver
FC = mpinfort
LD = $(FC)
AR = ar
RANLIB = ranlib

# C preprocessor
# (only required, if preprocessing is not performed via the compiler)
CPP = /usr/bin/cpp

# Definitions for CPP
# Define USE_PDAF to include PDAF
# Define BLOCKING_MPI_EXCHANGE to use blocking MPI commands to exchange data between model and PDAF
# Define PDAF_NO_UPDATE to deactivate the analysis step of the filter
# (if the compiler does not support get_command_argument()
# from Fortran 2003 you should define F77 here.)
CPP_DEFS = -DUSE_PDAF

# Optimization specs for compiler
# To use OpenMP parallelization in PDAF, specify it here (-fopenmp (gfortran) or -openmp (ifort))
# For gfortran 10 and later, you probably need to set -fallow-argument-mismatch
#   (You should explicitly define double precision for floating point
#   variables in the compilation)
OPT = -O3 -fdefault-real=8 -shared_mpi 

# Optimization specifications for Linker
OPT_LNK = $(OPT)

# Linking libraries (BLAS, LAPACK, if required: MPI)
LINK_LIBS = /opt/nec/ve/nlc/2.3.0/lib/liblapack.a /opt/nec/ve/nlc/2.3.0/lib/libblas_sequential.a -static-nec

# Specifications for the archiver
AR_SPEC =

# Specifications for ranlib
RAN_SPEC =

# Include path for MPI header file
MPI_INC =

# Object for nullMPI - if compiled without MPI library
OBJ_MPI =

MODULEOPT = -module

# NetCDF (only required for Lorenz96)
NC_LIB   = -L$(NETCDF_DIR)/lib -lnetcdff -lnetcdf -L/hpc/sw/hdf5/1.10.5/VE/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lsz -lz
NC_INC   = -I$(NETCDF_DIR)/include -I$(HDF5ROOT)/include
