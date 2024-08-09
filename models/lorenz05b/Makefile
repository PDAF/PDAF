# $Id: Makefile 714 2009-05-27 12:00:25Z lnerger $

##########################################################
# Generic Makefile for to build PDAF with Lorenz96 model #
# To choose the architecture set $PDAF_ARCH              #
##########################################################

# Root directory of PDAF package
BASEDIR = ../..

# Include machine-specific definitions
# For available include files see directory make.arch
# To choose a file, set PDAF_ARCH either here or by an
# environment variable.
include $(BASEDIR)/make.arch/$(PDAF_ARCH).h


######################################################

.SUFFIXES: .F90 .o 

# Modules used for the model part
MODULES =  	mod_parallel.o \
		mod_modeltime.o \
		timer.o \
		mod_memcount.o \
		mod_model.o \
		output_netcdf.o \
		parser_mpi.o 

# Model routines
OBJ_MODEL =  	main.o \
		initialize.o \
		integration.o

# Module required for assimilation
MOD_ASSIM = 	mod_assimilation.o \
		output_netcdf_asml.o

# Routines of observation handling (PDAF-OMI)
OBJ_USER_PDAFOMI = obs_gp_pdafomi.o \
		callback_obs_pdafomi.o

# Interface to PDAF - model sided
OBJ_PDAF_INT =  init_parallel_pdaf.o \
		init_pdaf.o \
		init_pdaf_parse.o \
		assimilation_pdaf.o

# Generic user-supplied routines
OBJ_USER_GEN  = init_ens_pdaf.o \
		init_ens_rnd.o \
		init_ens_eof.o \
		prepoststep_pdaf.o \
		next_observation_pdaf.o \
		distribute_state_pdaf.o \
		collect_state_pdaf.o \
		distribute_stateinc_pdaf.o \
		compute_truermse.o \
		compute_rms_smoother.o \
		init_obs_mask.o \
		add_model_noise.o

# User-supplied routines for localized analysis (LESKTF/LSEIK/LETKF/LNETF)
# Variant for assimilation of observations within cut-off distance
OBJ_USER_LOCAL = init_n_domains_pdaf.o \
		init_dim_l_pdaf.o \
		g2l_state_pdaf.o \
		l2g_state_pdaf.o

# User-suppled routines for generation of observations
OBJ_USER_OBSGEN = get_obs_f_pdaf.o \
		readwrite_obs.o

# Full list of user-supplied routines for online modes
OBJ_PDAF_USER = $(OBJ_USER_PDAFOMI) $(OBJ_USER_GEN) \
	$(OBJ_USER_LOCAL) $(OBJ_USER_OBSGEN)

######################################################

info:
	@echo "Makefile to build PDAF examples";
	@echo  "------------------------------------------------------------------";
	@echo  "Use as   make TARGET   where TARGET can be:";
	@echo  "  ----- BUILD -----------------------------------------------------";
	@echo  "  info               - This text";
	@echo  "  lorenz_05b         - Lorenz 2005 model II - forward (use without -DUSE_PDAF)";
	@echo  "  pdaf_lorenz_05b    - Lorenz 2005 model II - PDAF online  (use with -DUSE_PDAF)";
	@echo  "  ----- CLEAN -----------------------------------------------------";
	@echo  "  clean              - Clean up objects in this directory";
	@echo  "  cleandata          - Remove all .nc output files";
	@echo  "  cleanpdaf          - Clean up directory of PDAF library";
	@echo  "  cleanall           - Combined cleanup (cleans all)";
	@echo  "";
	@echo  "Depending on whether PDAF should be active, one has to define";
	@echo  "USE_PDAF in the include file in make.arch/ file for preprocessing.";
	@echo  "List available settings for PDAF_ARCH with 'make listarch'";
	@echo  "Use machine-specific definitions for PDAF_ARCH: " $(PDAF_ARCH);
	@echo  "------------------------------------------------------------------";


lorenz_05b : $(MODULES) $(OBJ_MPI) $(OBJ_MODEL)
	@echo "++++++ Build Lorenz05ii model without assimilation ++++++"
	$(LD)  $(OPT_LNK) -o $@  $(MODULES) $(OBJ_MPI) $(OBJ_MODEL) $(NC_LIB) $(LINK_LIBS)
	@echo "++++++ Done ++++++"

pdaf_lorenz_05b : libpdaf-d.a \
	$(MODULES) $(MOD_ASSIM) $(OBJ_MPI) $(OBJ_MODEL) $(OBJ_PDAF_USER) $(OBJ_PDAF_INT) 
	@echo "++++++  Build Lorenz96 model wit PDAF assimilation ++++++"
	$(LD) $(OPT_LNK) -o $@  \
	$(MODULES) $(MOD_ASSIM) $(OBJ_MPI) $(OBJ_MODEL) $(OBJ_PDAF_USER) $(OBJ_PDAF_INT) \
	-L$(BASEDIR)/lib -lpdaf-d $(NC_LIB) $(LINK_LIBS)
	@echo "++++++ Done ++++++"

######################################################

libpdaf-d.a: 
	@echo "++++++ Generate Filter library ++++++"
	@cd $(BASEDIR); make;


######################################################

.F90.o :
	$(FC) $(OPT) $(CPP_DEFS) -I$(BASEDIR)/src $(MPI_INC) $(NC_INC) -c $*.F90 

# For older compilers one might need to separate the
# preprocessing from the compilation as defined below:
#.F90.o :
#	$(CPP) -P -C $(MPI_INC) $(CPP_DEFS) $*.F90 $*.f90
#	$(FC) $(OPT) $(MPI_INC) -c $*.F90
#	@rm -f  $*.f90

######################################################
# Cleans

clean :
	rm -f *.o *.mod *.a lorenz_05b pdaf_lorenz_05b

cleanall : cleanpdaf clean cleandata

cleandata:
	rm -f *.nc

cleanpdaf:
	@echo "+++ Clean up PDAF directory"
	cd $(BASEDIR); make clean


######################################################
# List arch files

listarch:
	@echo Available architecture-specific input files for PDAF_ARCH
	@echo ---------------------------------------------------------
	@ls -1 ../../make.arch | cut -d"." -f1



$(OBJ_MODEL) : $(MODULES)
