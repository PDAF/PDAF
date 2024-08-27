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


# Name of executable
EXE = PDAF_online

# End of user specifications
######################################################

.SUFFIXES: .F90 .o 

# Modules used for the model part
MODULES = 	mod_parallel_pdaf.o \
		mod_model.o \
		parser_mpi.o

# Model routines
OBJ_MODEL =	main.o

# Moduls used for PDAF
MOD_ASSIM = 	mod_assimilation.o

# Routines of observation handling (PDAF-OMI)
OBJ_USER_PDAFOMI = obs_OBSTYPE_pdafomi_TEMPLATE.o \
		callback_obs_pdafomi.o

# Interface to PDAF - model sided
OBJ_PDAF_INT =  init_parallel_pdaf.o \
		finalize_pdaf.o \
		init_pdaf.o \
		init_pdaf_parse.o \
		init_pdaf_info.o \
		assimilate_pdaf.o 

# Generic user-supplied routines
OBJ_USER_GEN  = init_ens_pdaf.o \
		prepoststep_pdaf.o \
		next_observation_pdaf.o \
		distribute_state_pdaf.o \
		collect_state_pdaf.o

# User-supplied routines for localized analysis
OBJ_USER_LOCAL = init_n_domains_pdaf.o \
		init_dim_l_pdaf.o

# Full list of user-supplied routines for online modes
OBJ_PDAF_USER = $(OBJ_USER_PDAFOMI) $(OBJ_PDAF_INT) $(OBJ_USER_GEN) $(OBJ_USER_LOCAL)


######################################################

$(EXE) : libpdaf-d.a \
	$(MODULES) $(MOD_ASSIM) $(OBJ_MODEL) $(OBJ_MODEL_PDAF) $(OBJ_PDAF_USER)
	$(LD)  $(OPT_LNK)  -o $@ \
	$(MODULES) $(MOD_ASSIM) $(OBJ_MODEL) $(OBJ_MODEL_PDAF) $(OBJ_PDAF_USER) \
	-L$(BASEDIR)/lib -lpdaf-d $(LINK_LIBS)
#	$(MODULES) $(MOD_ASSIM) $(OBJ_MODEL) $(OBJ_PDAF_USER) $(OBJ_PDAF_INT) $(OBJ_MPI)
#	$(LD)  $(OPT_LNK)  -o $@  \
#	$(MODULES) $(MOD_ASSIM) $(OBJ_MODEL) $(OBJ_PDAF_USER) $(OBJ_PDAF_INT) $(OBJ_MPI) \
#	-L$(BASEDIR)/lib -lpdaf-d \
#	$(LINK_LIBS)
	@echo "++++++ Done ++++++"

######################################################

libpdaf-d.a: 
	@echo "++++++ Generate Filter library ++++++"
	@cd $(BASEDIR); make;

######################################################

.F90.o :
	$(FC) $(OPT) $(CPP_DEFS) -I$(BASEDIR)/include $(MPI_INC) -c $*.F90

# For older compilers one might need to separate the
# preprocessing from the compilation as defined below:
#.F90.o :
#	$(CPP) -P -C $(MPI_INC) $(CPP_DEFS) $*.F90 $*.f90
#	$(FC) $(OPT) $(MPI_INC) -c $*.F90
#	@rm -f  $*.f90

cleanall : cleanpdaf clean

clean :
	rm -f *.o *.mod *.a $(EXE)

cleanpdaf:
	@echo "+++ Clean up PDAF directory"
	cd $(BASEDIR); make clean
