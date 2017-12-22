# $Id: Makefile 714 2009-05-27 12:00:25Z lnerger $

##########################################################
# Generic Makefile for to build PDAF with Lorenz96 model #
# To choose the architecture set $PDAF_ARCH              #
##########################################################

.SUFFIXES: .F90 .o 

# include object lists
#include make.objects

# Modules used for the model part
MODULES =  	mod_model.o \
		output_netcdf.o

# Module required for assimilation
MOD_ASSIM = 	mod_assimilation.o \
		output_netcdf_asml.o

# Model routines
OBJ_MODEL =  	initialize.o \
		integration.o

# Interface to PDAF - model sided
OBJ_PDAF_INT =  init_parallel_pdaf_D.o \
		init_pdaf.o \
		init_pdaf_parse.o \
		assimilation_pdaf.o

# Generic user-supplied routines
OBJ_USER_GEN  = next_observation.o \
		distribute_state.o \
		collect_state.o \
		init_dim_obs.o \
		init_obs.o \
		obs_op.o \
		distribute_stateinc.o \
		compute_truermse.o \
		compute_rms_smoother.o \
		init_ens_eof.o \
		seik_omega.o \
		init_ens_rnd.o \
		init_obs_mask.o \
		add_model_noise.o

# User-supplied routines for SEEK
#OBJ_USER_SEEK = init_seek_dummy_D.o \
#		prepoststep_seek_dummy_D.o

# User-supplied routines for SEIK
OBJ_USER_SEIK = init_seik.o \
		prepoststep_seik.o

# User-supplied routines for EnKF
OBJ_USER_ENKF = add_obs_error.o \
		init_obscovar.o

# User-supplied routine for SEEK and SEIK
OBJ_USER_SEEKSEIK = prodrinva.o \
		init_obsvar.o

# User-supplied routines for localized analysis (LSEIK)
# Variant for assimilation of observations within cut-off distance
OBJ_USER_LOCAL = init_n_domains.o \
		init_dim_local.o \
		global2local_state.o \
		local2global_state.o \
		prodrinva_local.o \
		init_obs_local.o \
		init_dim_obs_local.o \
		global2local_obs.o \
		init_obsvar_local.o \
		init_dim_obs_full.o \
		init_obs_full.o \
		obs_op_full.o
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

# User-supplied routines for NETF
OBJ_USER_NETF = likelihood.o

# User-supplied routines for localized NETF
OBJ_USER_LNETF = likelihood_local.o

# Full list of user-supplied routines for online modes
OBJ_PDAF_USER = $(OBJ_USER_GEN) $(OBJ_USER_SEEK) $(OBJ_USER_SEIK) $(OBJ_USER_SEEKSEIK) \
	$(OBJ_USER_ENKF) $(OBJ_USER_LOCAL) $(OBJ_USER_LNETF) $(OBJ_USER_NETF)

######################################################

# Include machine-specific definitions
# For available include files see directory make.arch
# To choose a file, set PDAF_ARCH either here or by an
# environment variable.
include ../../../make.arch/$(PDAF_ARCH).h

######################################################

lorenz96.a : $(MODULES) $(OBJ_MODEL)
	@echo "++++++ Generate archive for model Lorenz96 ++++++"
	$(AR) -r $(AR_SPEC) $@ $(MODULES) $(OBJ_MODEL)
	@echo "++++++ Done ++++++"

lorenz96_pdaf.a : $(MODULES) $(MOD_ASSIM) $(OBJ_MODEL) $(OBJ_PDAF_INT) $(OBJ_PDAF_USER)
	@echo "++++++ Generate archive for model Lorenz96 ++++++"
	$(AR) -r $(AR_SPEC) $@ $(MODULES) $(MOD_ASSIM) $(OBJ_MODEL) $(OBJ_PDAF_INT) $(OBJ_PDAF_USER)
	@echo "++++++ Done ++++++"

.F90.o :
	$(FC) $(OPT) $(CPP_DEFS) -I../main -I../../../src $(NC_INC) -c $*.F90 

# For older compilers one might need to separate the
# preprocessing from the compilation as defined below:
#.F90.o :
#	$(CPP) -P -C $(MPI_INC) $(CPP_DEFS) $*.F90 $*.f90
#	$(FC) $(OPT) $(MPI_INC) -c $*.F90
#	@rm -f  $*.f90

clean :
	rm -f *.o *.mod *.a


$(OBJ_MODEL) : $(MODULES)
