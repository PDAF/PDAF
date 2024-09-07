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


# Name of executable
EXE = PDAF_offline

# End of user specifications
######################################################

.SUFFIXES: .F90 .o 

# Modules used for the model part
MODULES =  	mod_parallel.o \
		timer.o \
		mod_memcount.o \
		parser_mpi.o

# Module required for assimilation
MOD_ASSIM = 	mod_assimilation.o

# Model routines
OBJ_MODEL =	main_offline.o\
		initialize.o

# Interface to PDAF - model sided
OBJ_PDAF_INT =  init_parallel_pdaf.o \
		init_pdaf_offline.o \
		init_pdaf_parse.o \
		init_pdaf_info.o \
		assimilation_pdaf_offline.o \
		finalize_pdaf.o

# Generic user-supplied routines
OBJ_USER_GEN  = init_ens_offline.o \
		collect_state_pdaf_offline.o \
		init_obs_pdaf.o \
		init_dim_obs_pdaf.o \
		obs_op_pdaf.o \
		prepoststep_ens_offline.o

# User-supplied routines for SEIK
OBJ_USER_SEEKSEIK = prodrinva_pdaf.o \
		init_obsvar_pdaf.o

# User-supplied routines for EnKF
OBJ_USER_ENKF = add_obs_error_pdaf.o \
		init_obscovar_pdaf.o 

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

# User-supplied routines for LEnKF
OBJ_USER_LENKF = localize_covar_pdaf.o

# User-supplied routines for NETF
OBJ_USER_NETF = likelihood_pdaf.o

# User-supplied routines for localized NETF
OBJ_USER_LNETF = likelihood_l_pdaf.o

# User-supplied routines for SEIK
OBJ_USER_LKNETF = prodrinva_hyb_l_pdaf.o \
		likelihood_hyb_l_pdaf.o

# Full list of user-supplied routines for online modes
OBJ_PDAF_USER = $(OBJ_USER_GEN) $(OBJ_USER_SEEK) $(OBJ_USER_ENS) \
	$(OBJ_USER_SEEKSEIK) $(OBJ_USER_ENKF) $(OBJ_USER_LOCAL) \
	$(OBJ_USER_LENKF) $(OBJ_USER_LNETF) $(OBJ_USER_NETF) \
	$(OBJ_USER_LKNETF)

######################################################

$(EXE) : libpdaf-d.a \
	$(MODULES) $(MOD_ASSIM) $(OBJ_MODEL) $(OBJ_PDAF_INT) $(OBJ_PDAF_USER) $(OBJ_MPI)
	$(LD)  $(OPT_LNK)  -o $@  \
	$(MODULES) $(MOD_ASSIM) $(OBJ_MODEL) $(OBJ_PDAF_INT) $(OBJ_PDAF_USER) $(OBJ_MPI) \
	-L$(BASEDIR)/lib -lpdaf-d \
	$(LINK_LIBS)
	@echo "++++++ Done ++++++"

######################################################

libpdaf-d.a: 
	@echo "++++++ Generate Filter library ++++++"
	@cd $(BASEDIR); make;

######################################################

.F90.o :
	$(FC) $(OPT) $(CPP_DEFS) $(MPI_INC) -I$(BASEDIR)/include -c $*.F90

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
