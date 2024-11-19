######################################################
# Include file with machine-specific definitions     #
# for building PDAF.                                 #
#                                                    #
# Variant for CMake-based TSMP2-PDAF                 #
#                                                    #
# See: https://github.com/HPSCTerrSys/TSMP2          #
#                                                    #
# Environment variables are set in one of            #
# - environment file                                 #
#   - `env/[machine].[year]_[compiler].sh`           #
# - PDAF-related CMake scripts                       #
#   - `cmake/BuildPDAF.cmake`                        #
#   - `cmake/BuildPDAFMODEL.cmake`                   #
#   - `cmake/BuildPDAFFRAMEWORK.cmake`               #
######################################################

# Compiler, Linker, and Archiver
# FC = # Using environment default
LD = $(FC)
# CC = # Using environment default
AR = ar
RANLIB = ranlib 

# C preprocessor
# (only required, if preprocessing is not performed via the compiler)
CPP = 

# Definitions for CPP
# Define USE_PDAF to include PDAF
# Define BLOCKING_MPI_EXCHANGE to use blocking MPI commands to exchange data between model and PDAF
# (if the compiler does not support get_command_argument()
# from Fortran 2003 you should define F77 here.)
CPP_DEFS = ${TSMPPDAFCPP_DEFS}

# Optimization specs for compiler
#   (You should explicitly define double precision for floating point
#   variables in the compilation)  
OPT= ${TSMPPDAFFOPT} ${TSMPPDAFDOUBLEPRECISION}

# Optimization specifications for Linker
OPT_LNK = $(OPT)

# Linking libraries (BLAS, LAPACK, if required: MPI)
LINK_LIBS = ${TSMPPDAFLINK_LIBS}


# Specifications for the archiver
AR_SPEC = 

# Specifications for ranlib
RAN_SPEC =

# Specification for directory holding modules (-module for Intel, -J for GNU)
MODULEOPT = -module

# Include path for MPI header file
MPI_INC = ${TSMPPDAFMPI_INC}

# Object for nullMPI - if compiled without MPI library
#OBJ_MPI = nullmpi.o

# NetCDF (only required for Lorenz96)
NC_LIB   = 
NC_INC   = 
