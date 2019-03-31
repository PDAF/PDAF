# $Id: Makefile 714 2009-05-27 12:00:25Z lnerger $

##########################################################
# Generic Makefile for to build PDAF with Lorenz63 model #
# To choose the architecture set $PDAF_ARCH              #
##########################################################

.SUFFIXES: .F90 .o 

# include object lists
#include make.objects

# Modules used for the model part
MODULES =  	mod_model.o

# Module required for assimilation
MOD_ASSIM = 	mod_assimilation.o

# Model routines
OBJ_MODEL =  	initialize.o \
		integration.o

# Interface to PDAF - model sided
OBJ_PDAF_INT =  init_parallel_pdaf_D.o \
		init_pdaf_dummy_D.o \
		assimilation_pdaf_dummy_D.o

# Generic user-supplied routines
OBJ_USER_GEN  = next_observation_dummy_D.o \
		distribute_state_dummy_D.o \
		collect_state_dummy_D.o \
		init_dim_obs_dummy_D.o \
		init_obs_dummy_D.o \
		obs_op_dummy_D.o \
		distribute_stateinc_dummy_D.o

# User-supplied routines for SEEK
OBJ_USER_SEEK = init_seek_dummy_D.o \
		prepoststep_seek_dummy_D.o

# User-supplied routines for SEIK
OBJ_USER_SEIK = init_seik_dummy_D.o \
		prepoststep_seik_dummy_D.o

# User-supplied routines for EnKF
OBJ_USER_ENKF = init_enkf_dummy_D.o \
		add_obs_error_dummy_D.o \
		init_obscovar_dummy_D.o \
		prepoststep_enkf_dummy_D.o \

# User-supplied routine for SEEK and SEIK
OBJ_USER_SEEKSEIK = prodrinva_dummy_D.o \
		init_obsvar_dummy_D.o

# User-supplied routines for localized analysis (LSEIK)
# Variant for assimilation of observations within cut-off distance
OBJ_USER_LOCAL = init_n_domains_dummy_D.o \
		init_dim_local_dummy_D.o \
		global2local_state_dummy_D.o \
		local2global_state_dummy_D.o \
		prodrinva_local_dummy_D.o \
		init_obs_local_domain_dummy_D.o \
		init_dim_obs_local_domain_dummy_D.o \
		global2local_obs_domain_dummy_D.o \
		init_obsvar_local_dummy_D.o \
		init_dim_obs_full_dummy_D.o \
		init_obs_full_dummy_D.o \
		obs_op_full_dummy_D.o
# Variant for assimilation of global observation vector
#OBJ_USER_LOCAL = init_n_domains_dummy_D.o \
#		init_dim_local_dummy_D.o \
#		global2local_state_dummy_D.o \
#		local2global_state_dummy_D.o \
#		prodrinva_local_dummy_D.o \
#		init_obs_local_global_dummy_D.o \
#		init_dim_obs_local_global_dummy_D.o \
#		global2local_obs_global_dummy_D.o \
#		init_obsvar_local_dummy_D.o \
#		init_obs_full_dummy_D.o \
#		init_dim_obs_full_dummy_D.o \
#		obs_op_full_dummy_D.o

# Full list of user-supplied routines for online modes
OBJ_PDAF_USER = $(OBJ_USER_GEN) $(OBJ_USER_SEEK) $(OBJ_USER_SEIK) $(OBJ_USER_SEEKSEIK) $(OBJ_USER_ENKF) $(OBJ_USER_LOCAL)

######################################################

# Include machine-specific definitions
# For available include files see directory make.arch
# To choose a file, set PDAF_ARCH either here or by an
# environment variable.
include ../../../make.arch/$(PDAF_ARCH).h

######################################################

lorenz63.a : $(MODULES) $(OBJ_MODEL)
	@echo "++++++ Generate archive for model lorenz63 ++++++"
	$(AR) -r $(AR_SPEC) $@ $(MODULES) $(OBJ_MODEL)
	@echo "++++++ Done ++++++"

lorenz63_pdaf.a : $(MODULES) $(MOD_ASSIM) $(OBJ_MODEL) $(OBJ_PDAF_INT) $(OBJ_PDAF_USER)
	@echo "++++++ Generate archive for model lorenz63 ++++++"
	$(AR) -r $(AR_SPEC) $@ $(MODULES) $(MOD_ASSIM) $(OBJ_MODEL) $(OBJ_PDAF_INT) $(OBJ_PDAF_USER)
	@echo "++++++ Done ++++++"

.F90.o :
	$(FC) $(OPT) $(CPP_DEFS) -I../main -c $*.F90

# For older compilers one might need to separate the
# preprocessing from the compilation as defined below:
#.F90.o :
#	$(CPP) -P -C $(MPI_INC) $(CPP_DEFS) $*.F90 $*.f90
#	$(FC) $(OPT) $(MPI_INC) -c $*.F90
#	@rm -f  $*.f90

clean :
	rm -f *.o *.mod *.a


$(OBJ_MODEL) : $(MODULES)
