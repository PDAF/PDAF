# $Id: Makefile 746 2009-08-04 12:16:28Z lnerger $

#######################################################
# Generic Makefile for to build PDAF with dummy model #
# To choose the architecture set $PDAF_ARCH           #
#######################################################

#######################################################
# user specifications
######################################################

#PROG=tsmp_pdaf
PROG=$(TSMPPDAFBINDIR)/tsmp-pdaf
# User specifications
# 1. Set BASEDIR, the directory where the PDAF package resides
# 2. Set PDAF_ARCH to include compile definitions
#    (See directory BASEDIR/make.arch for files. PDAF_ARCH is filename without .h)

# Root directory of PDAF package
BASEDIR = $(PDAF_DIR)
# Include machine-specific definitions
# For available include files see directory make.arch
# To choose a file, set PDAF_ARCH either here or by an
# environment variable.
include $(BASEDIR)/make.arch/$(PDAF_ARCH).h

FFLAGS = -cpp $(TSMPPDAFIMPORTFLAGS)
MODELDIR = $(TSMPPDAFLIBDIR)
LIBMODEL = $(MODELDIR)/libmodel.a

# Optimization specs for C/C++ compiler
#   C/C++ compiler do not get the floating point definition
COPT= ${TSMPPDAFCOPT}

CPP_DEFS += ${TSMPPDAFCPPDEFS}


# End of user specifications

######################################################
# case dependent settings for TSMP-PDAF
######################################################
## template
#ifeq ($(COUPLE_TYPE), parflow_stand_alone)
## parflow stand alone
#else
#ifeq ($(COUPLE_TYPE), clm_stand_alone)
## clm stand alone
#else
## terrsysmp
#endif
#endif

# optimisation flags


# case dependent libraries

LIBS += -L$(MODELDIR) $(TSMPPDAFLIBS)



##################################################
# Rules et al.
######################################################

.SUFFIXES: .c .F90 .o 

# Modules used for the model part
MODULES =  	mod_tsmp.o

# Moduls used for PDAF
MOD_ASSIM = 	mod_parallel_pdaf.o \
		mod_assimilation.o \
		parser_mpi.o \
		mod_read_obs.o

# Model routines used with PDAF
OBJ_MODEL_PDAF =pdaf_terrsysmp.o\
		integrate_pdaf.o

# Interface to PDAF - model sided
OBJ_PDAF_INT =  init_parallel_pdaf.o \
		finalize_pdaf.o \
		init_pdaf.o \
		init_pdaf_parse.o \
		init_pdaf_info.o \
		assimilate_pdaf.o

# Generic user-supplied routines
OBJ_USER_GEN  = init_ens.o \
		prepoststep_ens_pdaf.o \
		next_observation_pdaf.o \
		distribute_state_pdaf.o \
		collect_state_pdaf.o \
		init_obs_pdaf.o \
		init_dim_obs_pdaf.o \
		obs_op_pdaf.o

# User-supplied routines for ESTKF/SEIK/ETKF
OBJ_USER_SEIK = prodrinva_pdaf.o \
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
		obs_op_f_pdaf.o \
		localize_covar_pdaf.o

# Full list of user-supplied routines for online modes
OBJ_PDAF_USER = $(OBJ_USER_GEN) $(OBJ_USER_SEIK) $(OBJ_USER_ENKF) $(OBJ_USER_LOCAL)

######################################################

all: $(PROG)

info:
	@echo "Makefile to build PDAF tutorial online implementation";
	@echo "Example: 2D model with parallelization";
	@echo "------------------------------------------------------------------";
	@echo "Use as   make TARGET   where TARGET can be:";
	@echo "  model         - just the forward model";
	@echo "  model_pdaf    - assimilation model (use with -DUSE_PDAF!)";
	@echo "  (if compiled without -DUSE_PDAF, both programs behave equally)"
	@echo "";
	@echo "Depending on whether PDAF should be active, one has to define";
	@echo "USE_PDAF in the include file in make.arch/ file for preprocessing.";
	@echo "Use machine-specific definitions for PDAF_ARCH: " $(PDAF_ARCH);
	@echo "------------------------------------------------------------------";


$(PROG) : $(LIBMODEL) libpdaf-d.a \
	$(MODULES) $(MOD_ASSIM) $(OBJ_MODEL_PDAF) $(OBJ_PDAF_INT) $(OBJ_PDAF_USER) 
	$(PREP_C) $(LD)  $(OPT_LNK)  -o $@ \
	$(MODULES) $(MOD_ASSIM) $(OBJ_MODEL_PDAF) $(OBJ_PDAF_INT) $(OBJ_PDAF_USER) \
	-L$(BASEDIR)/lib -lpdaf-d $(LINK_LIBS) $(LIBMODEL) $(LIBS)

######################################################

libpdaf-d.a: 
	@echo "++++++ Generate Filter library +++++++"

#	@cd $(BASEDIR)/src; make;

######################################################
%.o : %.c
	$(PREP_C) $(CC) $(COPT) $(CPP_DEFS) $(CFLAGS) -c $< -o $@

%.o %.mod: %.F90
	$(PREP_C) $(FC) $(OPT) $(CPP_DEFS) $(FFLAGS) -c $*.F90

.PHONY: clean depend
cleanall : cleanpdaf clean

clean :
	rm -f *.o *.mod *.a $(PROG)

cleanpdaf:
	@echo "+++ Clean up PDAF directory"
	cd $(BASEDIR)/src; make clean

depend .depend:
	mpif90 -M *.F90 > .depend 
