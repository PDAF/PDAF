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

# Modules used for the model part
MODULES =  	mod_model.o

# Model routines
OBJ_MODEL =	main.o\
		initialize.o \
		integrate.o

# Moduls used for PDAF
MOD_ASSIM = 	mod_parallel_pdaf.o \
		mod_assimilation.o \
		parser_mpi.o

# Model routines used with PDAF
OBJ_MODEL_PDAF =main_pdaf.o\
		initialize.o \
		integrate_pdaf.o

# Interface to PDAF - model sided
OBJ_PDAF_INT =  init_parallel_pdaf.o \
		finalize_pdaf.o \
		init_pdaf.o \
		init_pdaf_parse.o \
		init_pdaf_info.o \
		assimilate_pdaf.o

# Generic user-supplied routines
OBJ_USER_GEN  = init_ens_pdaf.o \
		prepoststep_ens_pdaf.o \
		next_observation_pdaf.o \
		distribute_state_pdaf.o \
		collect_state_pdaf.o \
		init_obs_pdaf.o \
		init_dim_obs_pdaf.o \
		obs_op_pdaf.o

# User-supplied routines for ESTKF/SEIK/ETKF
OBJ_USER_SQRTF = prodrinva_pdaf.o \
		init_obsvar_pdaf.o

# User-supplied routines for EnKF
OBJ_USER_ENKF = add_obs_error_pdaf.o \
		init_obscovar_pdaf.o 

# User-supplied routines for localized analysis (LESTKF/LSEIK/LETKF)
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
OBJ_PDAF_USER = $(OBJ_USER_GEN) $(OBJ_USER_SQRTF) $(OBJ_USER_ENKF) $(OBJ_USER_LOCAL)

######################################################

info:
	@echo "Makefile to build PDAF tutorial online implementation";
	@echo "Example: 2D serial model (without parallelization)";
	@echo "------------------------------------------------------------------";
	@echo "Use as   make TARGET   where TARGET can be:";
	@echo "  model         - just the forward model (use without -DUSE_PDAF!)";
	@echo "  model_pdaf    - assimilation model     (use with -DUSE_PDAF!)";
	@echo "  (if compiled without -DUSE_PDAF, both programs behave equally)"
	@echo "";
	@echo "Depending on whether PDAF should be active, one has to define";
	@echo "USE_PDAF in the include file in make.arch/ file for preprocessing.";
	@echo "List available settings for PDAF_ARCH with 'make listarch'";
	@echo "Use machine-specific definitions for PDAF_ARCH: " $(PDAF_ARCH);
	@echo "------------------------------------------------------------------";

model : $(MODULES) $(OBJ_MODEL)
	$(LD)  $(OPT_LNK)  -o $@  $(MODULES) $(OBJ_MODEL)

model_pdaf : libpdaf-d.a \
	$(MODULES) $(MOD_ASSIM) $(OBJ_MODEL_PDAF) $(OBJ_PDAF_INT) $(OBJ_PDAF_USER)
	$(LD)  $(OPT_LNK)  -o $@ \
	$(MODULES) $(MOD_ASSIM) $(OBJ_MODEL_PDAF) $(OBJ_PDAF_INT) $(OBJ_PDAF_USER) \
	-L$(BASEDIR)/lib -lpdaf-d $(LINK_LIBS)

######################################################

libpdaf-d.a: 
	@echo "++++++ Generate Filter library ++++++"
	@cd $(BASEDIR); make;


######################################################

.F90.o :
	$(FC) $(OPT) $(CPP_DEFS)  -c $*.F90

# For older compilers one might need to separate the
# preprocessing from the compilation as defined below:
#.F90.o :
#	$(CPP) -P -C $(MPI_INC) $(CPP_DEFS) $*.F90 $*.f90
#	$(FC) $(OPT) $(MPI_INC) -c $*.F90
#	@rm -f  $*.f90

######################################################
# Cleans

cleanall : cleanpdaf clean cleandata

clean :
	rm -f *.o *.mod *.a model model_pdaf

cleanpdaf:
	@echo "+++ Clean up PDAF directory"
	cd $(BASEDIR); make clean

cleandata:
	rm -f ens*ana.txt ens*for.txt state_*ana.txt state_*for.txt true*.txt

cleandataq:
	@rm -f ens*ana.txt ens*for.txt state_*ana.txt state_*for.txt true*.txt

######################################################
# List arch files

listarch:
	@echo Available architecture-specific input files for PDAF_ARCH
	@echo ---------------------------------------------------------
	@ls -1 ../../make.arch | cut -d"." -f1

