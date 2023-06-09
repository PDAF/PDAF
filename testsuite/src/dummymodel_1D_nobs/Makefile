# $Id: Makefile 746 2009-08-04 12:16:28Z lnerger $

#######################################################
# Generic Makefile for to build PDAF with dummy model #
# To choose the architecture set $PDAF_ARCH           #
#######################################################

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
OBJ_PDAF_INT =  init_parallel_pdaf.o \
		init_pdaf.o \
		init_pdaf_parse.o \
		init_pdaf_info.o \
		assimilation_pdaf.o

# Generic user-supplied routines
OBJ_USER_GEN  = next_observation_pdaf.o \
		distribute_state_pdaf.o \
		collect_state_pdaf.o \
		init_dim_obs_pdaf.o \
		init_obs_pdaf.o \
		obs_op_pdaf.o \
		distribute_stateinc_pdaf.o

# User-supplied routines for SEEK
OBJ_USER_SEEK = init_seek_pdaf.o \
		prepoststep_seek_pdaf.o

# User-supplied routines for SEIK/ETKF
OBJ_USER_ENS  = init_seik_pdaf_large.o \
		prepoststep_ens_pdaf_large.o

# User-supplied routines for EnKF
OBJ_USER_ENKF = init_enkf_pdaf.o \
		add_obs_error_pdaf.o \
		init_obscovar_pdaf.o

# User-supplied routine for SEEK/SEIK/ETKF
OBJ_USER_SEEKSEIK = prodrinva_pdaf.o \
		init_obsvar_pdaf.o

# User-supplied routines for localized analysis (LSEIK)
OBJ_USER_LOCAL = init_n_domains_pdaf.o \
		init_dim_l_pdaf.o \
		g2l_state_pdaf.o \
		l2g_state_pdaf.o \
		prodrinva_l_pdaf.o \
		init_obs_l_pdaf.o \
		init_dim_obs_l_pdaf.o \
		g2l_obs_pdaf.o \
		init_obsvar_l_pdaf.o \
		init_dim_obs_f_pdaf.o \
		init_obs_f_pdaf.o \
		obs_op_f_pdaf.o

# Full list of user-supplied routines for online modes
OBJ_PDAF_USER = $(OBJ_USER_GEN) $(OBJ_USER_SEEK) $(OBJ_USER_ENS) \
	$(OBJ_USER_SEEKSEIK) $(OBJ_USER_ENKF) $(OBJ_USER_LOCAL)

######################################################

# Include machine-specific definitions
# For available include files see directory make.arch
# To choose a file, set PDAF_ARCH either here or by an
# environment variable.
include ../../../make.arch/$(PDAF_ARCH).h

######################################################

dummymodel_1D.a : $(MODULES) $(OBJ_MODEL)
	@echo "++++++ Generate archive for model dummymodel_1D ++++++"
	$(AR) -r $(AR_SPEC) $@ $(MODULES) $(OBJ_MODEL)
	@echo "++++++ Done ++++++"

dummymodel_1D_pdaf_nobs.a : $(MODULES) $(MOD_ASSIM) $(OBJ_MODEL) $(OBJ_PDAF_INT) $(OBJ_PDAF_USER)
	@echo "++++++ Generate archive for model dummymodel_1D ++++++"
	$(AR) -r $(AR_SPEC) $@ $(MODULES) $(MOD_ASSIM) $(OBJ_MODEL) $(OBJ_PDAF_INT) $(OBJ_PDAF_USER)
	@echo "++++++ Done ++++++"

.F90.o :
	$(FC) $(OPT) $(MPI_INC) $(CPP_DEFS) -I../main -c $*.F90

# For older compilers one might need to separate the
# preprocessing from the compilation as defined below:
#.F90.o :
#	$(CPP) -P -C $(MPI_INC) $(CPP_DEFS) $*.F90 $*.f90
#	$(FC) $(OPT) $(MPI_INC) -c $*.F90
#	@rm -f  $*.f90

clean :
	rm -f *.o *.mod *.a
