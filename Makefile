# $Id: Makefile 1856 2017-12-06 08:36:03Z lnerger $

#######################################################
# Generic Makefile for to build PDAF library          #
# To choose the architecture set $PDAF_ARCH           #
#                                                     #
# Armin Corbin - University of Bonn                   #
#######################################################

#######################################################
# definitions for highlighting outputs
bold := $(shell tput bold)
sgr0 := $(shell tput sgr0)

#######################################################

# Include machine-specific definitions
# For available include files see directory make.arch
# To choose a file, set PDAF_ARCH either here or by an
# environment variable.
include ./make.arch/$(PDAF_ARCH).h
$(info $(bold)use machine-specific definitions in $(PDAF_ARCH).h$(sgr0))

#######################################################

OBJDIR:=build
INCDIR:=include
LIBDIR:=lib

SRCDIR:=src
EXTDIR:=external

#######################################################
# List of sources for PDAF library
#######################################################

# Modules used in PDAF
SRC_MOD_PDAF =  PDAF_timer.F90 \
		PDAF_memcount.F90 \
		PDAF_mod_filtermpi.F90 \
		PDAF_mod_filter.F90 \
		PDAFobs.F90

# Module file with interface definitions
SRC_MOD_INTERFACE = PDAF_da.F90 \
		PDAF.F90 \
		PDAF_forecast.F90 \
		PDAF_iau.F90 \
		PDAF_assimilate_ens.F90 \
		PDAF_put_state_ens.F90 \
		PDAF_assimilate_ens_nondiagR.F90 \
		PDAF_put_state_ens_nondiagR.F90

# Generic routines in PDAF
SRC_PDAF_GEN = 	PDAF_init_mod.F90 \
		PDAF_init.F90 \
		PDAF_utils.F90 \
		PDAF_get.F90 \
		PDAF_set.F90 \
		PDAF_info.F90 \
		PDAF_utils_filters.F90 \
		PDAF_analysis_utils.F90 \
		PDAF_communicate_ens.F90 \
		PDAF_get_state_mod.F90 \
		PDAF_get_state.F90 \
		PDAF_smoother.F90 \
		PDAF_put_state_prepost.F90 \
		PDAF_assimilate_prepost.F90 \
		PDAF_prepost.F90 \
		PDAF_sample.F90 \
		PDAF_diag.F90 \
		PDAF_comm_obs.F90 \
		PDAFlocal.F90 \
		PDAFlocal_callback.F90 \
		PDAFlocal_assimilate_ens.F90 \
		PDAFlocal_put_state_ens.F90 \
		PDAFlocal_assimilate_3dvars.F90 \
		PDAFlocal_put_state_3dvars.F90

# Specific PDAF-routines for SEIK
SRC_SEIK =	PDAF_seik.F90 \
		PDAF_seik_analysis.F90 \
		PDAF_seik_analysis_newT.F90 \
		PDAF_seik_analysis_trans.F90 \
		PDAF_seik_update.F90 \
		PDAF_put_state_seik.F90 \
		PDAF_assimilate_seik.F90

# Specific PDAF-routines for local SEIK
SRC_LSEIK =     PDAF_lseik.F90 \
		PDAF_lseik_analysis.F90 \
		PDAF_lseik_analysis_trans.F90 \
		PDAF_lseik_update.F90 \
		PDAF_put_state_lseik.F90 \
		PDAF_assimilate_lseik.F90

# Specific PDAF-routines for SEEK
SRC_SEEK =      PDAF_seek.F90 \
		PDAF_seek_analysis.F90 \
		PDAF_seek_update.F90 \
		PDAF_put_state_seek.F90 \
		PDAF_assimilate_seek.F90

# Specific PDAF-routines for EnKF
SRC_ENKF =	PDAF_enkf.F90 \
		PDAF_enkf_analysis_rlm.F90 \
		PDAF_enkf_analysis_rsm.F90 \
		PDAF_enkf_update.F90 \
		PDAF_put_state_enkf.F90 \
		PDAF_assimilate_enkf.F90

# Specific PDAF-routines for ETKF
SRC_ETKF =	PDAF_etkf.F90 \
		PDAF_etkf_analysis.F90 \
		PDAF_etkf_analysis_T.F90 \
		PDAF_etkf_analysis_fixed.F90 \
		PDAF_etkf_update.F90 \
		PDAF_put_state_etkf.F90 \
		PDAF_assimilate_etkf.F90

# Specific PDAF-routines for LETKF
SRC_LETKF =     PDAF_letkf.F90 \
		PDAF_letkf_analysis.F90 \
		PDAF_letkf_analysis_T.F90 \
		PDAF_letkf_analysis_fixed.F90 \
		PDAF_letkf_update.F90 \
		PDAF_put_state_letkf.F90 \
		PDAF_assimilate_letkf.F90

# Specific PDAF-routines for ESTKF
SRC_ESTKF =	PDAF_estkf.F90 \
		PDAF_estkf_analysis.F90 \
		PDAF_estkf_analysis_fixed.F90 \
		PDAF_estkf_update.F90 \
		PDAF_put_state_estkf.F90 \
		PDAF_assimilate_estkf.F90

# Specific PDAF-routines for LESTKF
SRC_LESTKF =	PDAF_lestkf.F90 \
		PDAF_lestkf_analysis.F90 \
		PDAF_lestkf_analysis_fixed.F90 \
		PDAF_lestkf_update.F90 \
		PDAF_put_state_lestkf.F90 \
		PDAF_assimilate_lestkf.F90

# Specific PDAF-routines for LEnKF
SRC_LENKF =	PDAF_lenkf.F90 \
		PDAF_put_state_lenkf.F90 \
		PDAF_assimilate_lenkf.F90 \
		PDAF_lenkf_update.F90 \
		PDAF_lenkf_analysis_rsm.F90

# Specific PDAF-routines for NETF
SRC_NETF =	PDAF_netf.F90 \
		PDAF_netf_analysis.F90 \
		PDAF_netf_update.F90 \
		PDAF_put_state_netf.F90 \
		PDAF_assimilate_netf.F90

# Specific PDAF-routines for LNETF
SRC_LNETF =	PDAF_lnetf.F90 \
		PDAF_lnetf_analysis.F90 \
		PDAF_lnetf_update.F90 \
		PDAF_put_state_lnetf.F90 \
		PDAF_assimilate_lnetf.F90

# Specific PDAF-routines for PF
SRC_PF =	PDAF_pf.F90 \
		PDAF_pf_analysis.F90 \
		PDAF_pf_update.F90 \
		PDAF_put_state_pf.F90 \
		PDAF_assimilate_pf.F90

# Specific PDAF-routines for LKNETF
SRC_LKNETF =	PDAF_lknetf.F90 \
		PDAF_lknetf_analysis_step.F90 \
		PDAF_lknetf_update_step.F90 \
		PDAF_lknetf_analysis_sync.F90 \
		PDAF_lknetf_update_sync.F90 \
		PDAF_lknetf_reset_gamma.F90 \
		PDAF_put_state_lknetf.F90 \
		PDAF_assimilate_lknetf.F90

# Specific PDAF-routines for ENSRF with serial observation processing
SRC_ENSRF =	PDAF_ensrf.F90 \
		PDAF_ensrf_analysis.F90 \
		PDAF_ensrf_update.F90 \
		PDAF_put_state_ensrf.F90 \
		PDAF_assimilate_ensrf.F90

# Specific PDAF-routines for generating observations
SRC_GENOBS =	PDAF_genobs.F90 \
		PDAF_generate_obs_update.F90 \
		PDAF_put_state_generate_obs.F90 \
		PDAF_generate_obs.F90

# Specific PDAF-routines for 3DVAR
SRC_3DVAR =	PDAF_3dvar.F90 \
		PDAF_3dvar_optim.F90 \
		PDAF_3dvar_analysis_cvt.F90 \
		PDAF_3dvar_update.F90 \
		PDAF_en3dvar_optim.F90 \
		PDAF_en3dvar_analysis_cvt.F90 \
		PDAF_en3dvar_update.F90 \
		PDAF_hyb3dvar_optim.F90 \
		PDAF_hyb3dvar_analysis_cvt.F90\
		PDAF_hyb3dvar_update.F90 \
		PDAF_assimilate_3dvars.F90 \
		PDAF_assimilate_3dvars_nondiagR.F90 \
		PDAF_put_state_3dvars.F90 \
		PDAF_put_state_3dvars_nondiagR.F90 \
		PDAF_put_state_3dvar.F90 \
		PDAF_assimilate_3dvar.F90 \
		PDAF_put_state_en3dvar_lestkf.F90 \
		PDAF_put_state_en3dvar_estkf.F90 \
		PDAF_assimilate_en3dvar_lestkf.F90 \
		PDAF_assimilate_en3dvar_estkf.F90 \
		PDAF_put_state_hyb3dvar_estkf.F90 \
		PDAF_put_state_hyb3dvar_lestkf.F90 \
		PDAF_assimilate_hyb3dvar_lestkf.F90 \
		PDAF_assimilate_hyb3dvar_estkf.F90

# Routines for PDAF-OMI
SRC_PDAFOMI =	PDAFomi_obs_f.F90 \
		PDAFomi_obs_l.F90 \
		PDAFomi_dim_obs_l.F90 \
		PDAFomi_obs_op.F90 \
		PDAFomi.F90 \
		PDAFomi_callback.F90 \
		PDAFomi_assimilate_ens.F90 \
		PDAFomi_assimilate_ens_nondiagR.F90 \
		PDAFomi_put_state_ens.F90 \
		PDAFomi_put_state_ens_nondiagR.F90 \
		PDAFomi_assimilate_3dvars.F90 \
		PDAFomi_put_state_3dvars.F90 \
		PDAFlocalomi_assimilate_ens.F90 \
		PDAFlocalomi_put_state_ens.F90 \
		PDAFlocalomi_assimilate_3dvars.F90 \
		PDAFlocalomi_put_state_3dvars.F90


# collect all PDAF sources
SRC_PDAF =  $(SRC_PDAFOMI) $(SRC_PDAF_GEN) \
	    $(SRC_SEIK) $(SRC_LSEIK) \
	    $(SRC_ENKF) $(SRC_ETKF) $(SRC_LETKF) \
	    $(SRC_ESTKF) $(SRC_LESTKF) $(SRC_LENKF) $(SRC_NETF) $(SRC_LNETF) \
	    $(SRC_LKNETF) $(SRC_PF) $(SRC_ENSRF) $(SRC_GENOBS) $(SRC_3DVAR) \
	    $(SRC_MOD_PDAF) $(SRC_MOD_INTERFACE)

# external sources
SRC_SANGOMA = $(EXTDIR)/SANGOMA/SANGOMA_quicksort.F90

#######################################################
# object files
#######################################################
OBJ_PDAF := $(SRC_PDAF:%.F90=$(OBJDIR)/%.o)

OBJ_PDAF_VAR = $(SRC_3DVAR:%.F90=$(OBJDIR)/%.o)

OBJ_SANGOMA = $(OBJDIR)/SANGOMA_quicksort.o

# External optimizer libraries implicitly build by make
OBJ_OPTIM = $(EXTDIR)/CG+_mpi/cgfam.o $(EXTDIR)/CG+_mpi/cgsearch.o \
	$(EXTDIR)/CG+/cgfam.o $(EXTDIR)/CG+/cgsearch.o \
	$(EXTDIR)/LBFGS/lbfgsb.o $(EXTDIR)/LBFGS/linpack.o \
	$(EXTDIR)/LBFGS/timer.o

#######################################################
# compiler instructions

COMPILE.f90 = $(FC) $(OPT) $(MPI_INC) $(CPP_DEFS) -c -o $@ $(MODULEOPT) $(INCDIR)

#######################################################
.PHONY: all
all: directories libpdaf libpdafvar

.PHONY: libpdaf
libpdaf: $(LIBDIR)/libpdaf-d.a

.PHONY: libpdafvar
libpdafvar: $(LIBDIR)/libpdaf-var.a
#######################################################

$(LIBDIR)/libpdaf-d.a: $(OBJ_PDAF) $(OBJ_SANGOMA)
	$(info $(bold)Generate Filter library$(sgr0))
	$(AR) rs $(AR_SPEC) $(LIBDIR)/libpdaf-d.a $(OBJ_PDAF) $(OBJ_SANGOMA)

$(LIBDIR)/libpdaf-var.a:  $(OBJ_PDAF) $(OBJ_PDAF_VAR) $(OBJ_OPTIM)
	$(info $(bold)Generate Var Filter library$(sgr0))
	$(AR) rs $(AR_SPEC) $(LIBDIR)/libpdaf-var.a $(OBJ_PDAF) $(OBJ_PDAF_VAR) $(OBJ_OPTIM)

# use pattern rule to create rules for all object files
$(OBJDIR)/%.o: $(SRCDIR)/%.F90
	$(info $(bold)compile $<$(sgr0))
	$(COMPILE.f90) $<

# explicite rule for sangoma
$(OBJ_SANGOMA) : $(SRC_SANGOMA)
	$(info $(bold)compile external dependency $<$(sgr0))
	$(COMPILE.f90) $<

#######################################################
MISSINGDIRS:= $(INCDIR) $(OBJDIR) $(LIBDIR)
.PHONY: directories
directories: $(MISSINGDIRS)

$(MISSINGDIRS):
	mkdir -p $@

#######################################################
.PHONY: clean
clean :
	$(info $(bold)clean PDAF$(sgr0))
	rm -rf $(OBJDIR)/*.o
	rm -rf $(INCDIR)/*.mod
	rm -rf $(LIBDIR)/libpdaf-d.a
	rm -rf $(LIBDIR)/libpdaf-var.a
	rm -rf $(EXTDIR)/CG+/*.o
	rm -rf $(EXTDIR)/CG+_mpi/*.o
	rm -rf $(EXTDIR)/LBFGS/*.o

#######################################################
.PHONY: listarch
listarch:
	@echo Available architecture-specific input files for PDAF_ARCH
	@echo ---------------------------------------------------------
	@ls -1 ../make.arch | cut -d"." -f1

#######################################################
# include file containing dependencies for parallel
# execution of make
# created with ./external/mkdepends/mkdepends ./src ./external/* '$(OBJDIR)'
include Depends
