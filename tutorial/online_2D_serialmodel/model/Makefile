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
BASEDIR = ../../..

# Include machine-specific definitions
# For available include files see directory make.arch
# To choose a file, set PDAF_ARCH either here or by an
# environment variable.
include $(BASEDIR)/make.arch/$(PDAF_ARCH).h


# End of user specifications
######################################################

.SUFFIXES: .F90 .o

.F90.o :
	$(FC) $(OPT) $(CPP_DEFS)  -c $*.F90

# Modules used for the model part
MODULES =  	mod_model.o

# Moduls used for PDAF
#MOD_ASSIM =    mod_parallel_pdaf.o \
#			parser_mpi.o

# Model routines
OBJ_MODEL =	main.o\
		initialize.o \
		integrate.o

# Interface to PDAF - model sided
#OBJ_PDAF_INT =  init_parallel_pdaf.o

# User-supplied routines for PDAF
#OBJ_PDAF_USER =

model : $(MODULES) $(OBJ_MODEL) #$(MOD_ASSIM) $(OBJ_PDAF_INT) $(OBJ_PDAF_USER)
	$(LD)  $(OPT_LNK)  -o $@  $(MODULES) $(OBJ_MODEL) #$(MOD_ASSIM) $(OBJ_PDAF_INT) $(OBJ_PDAF_USER)

clean :
	rm -f *.o *.mod *.a model

cleandata :
	rm -f true*.txt
