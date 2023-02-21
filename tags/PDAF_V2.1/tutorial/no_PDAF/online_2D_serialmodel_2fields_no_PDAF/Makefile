# $Id: Makefile 746 2009-08-04 12:16:28Z lnerger $

#######################################################
# Generic Makefile for to build PDAF with dummy model #
# To choose the architecture set $PDAF_ARCH           #
#######################################################

######################################################


# User specifications
# 1. Set BASEDIR, the directory where the PDAF package resides
# 2. Set PDAF_ARCH to include compile definitions
#    (See directory BASEDIR/make.arch for files. PDAF_ARCH is filename without .h)

# Root directory of PDAF package
BASEDIR = ../..

# Include machine-specific definitions
# For available include files see directory make.arch
# To choose a file, set PDAF_ARCH either here or by an
# environment variable.
include $(BASEDIR)/make.arch/$(PDAF_ARCH).h


# End of user specifications
######################################################

.SUFFIXES: .F90 .o 

# Modules used for the model part
MODULES =  	mod_model.o

# Model routines
OBJ_MODEL =	main.o\
		initialize.o \
		integrate.o

######################################################

model : $(MODULES) $(OBJ_MODEL)
	$(LD)  $(OPT_LNK)  -o $@  $(MODULES) $(OBJ_MODEL)


######################################################

.F90.o :
	$(FC) $(OPT) $(CPP_DEFS) -I$(BASEDIR)/include -c $*.F90

# For older compilers one might need to separate the
# preprocessing from the compilation as defined below:
#.F90.o :
#	$(CPP) -P -C $(MPI_INC) $(CPP_DEFS) $*.F90 $*.f90
#	$(FC) $(OPT) $(MPI_INC) -c $*.F90
#	@rm -f  $*.f90

######################################################
# Cleans

cleanall : clean cleandata

clean :
	rm -f *.o *.mod *.a model

cleandata:
	rm -f true*txt

cleandataq:
	@rm -f true*txt

######################################################
# List arch files

listarch:
	@echo Available architecture-specific input files for PDAF_ARCH
	@echo ---------------------------------------------------------
	@ls -1 ../../make.arch | cut -d"." -f1

