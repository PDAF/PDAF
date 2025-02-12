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
SRC_MOD_INTERFACE = PDAF_interfaces_module.F90 \
		PDAFlocal_interfaces.F90 

# Generic routines in PDAF
SRC_PDAF_GEN = 	PDAF_utils_filters.F90 \
		PDAF_init.F90 \
		PDAF_init_si.F90 \
		PDAF_alloc.F90 \
		PDAF_print_version.F90 \
		PDAF_print_info.F90 \
		PDAF_analysis_utils.F90 \
		PDAF_set_iparam.F90 \
		PDAF_set_rparam.F90 \
		PDAF_communicate_ens.F90 \
		PDAF_set_comm_pdaf.F90 \
		PDAF_get_state.F90 \
		PDAF_get_state_si.F90 \
		PDAF_incremental.F90 \
		PDAF_incremental_si.F90 \
		PDAF_add_increment.F90 \
		PDAF_generate_rndmat.F90 \
		PDAF_local_weights.F90 \
		PDAF_local_weight.F90 \
		PDAF_force_analysis.F90 \
		PDAF_set_seedset.F90 \
		PDAF_set_memberid.F90 \
		PDAF_get_memberid.F90 \
		PDAF_get_obsmemberid.F90 \
		PDAF_smoother.F90 \
		PDAF_set_smootherens.F90 \
		PDAF_get_smootherens.F90 \
		PDAF_set_ens_pointer.F90 \
		PDAF_put_state_prepost.F90 \
		PDAF_put_state_prepost_si.F90 \
		PDAF_assimilate_prepost.F90 \
		PDAF_assimilate_prepost_si.F90 \
		PDAF_prepost.F90 \
		PDAF_prepost_si.F90 \
		PDAF_sampleens.F90 \
		PDAF_mvnormalize.F90 \
		PDAF_eofcovar.F90 \
		PDAF_diag_histogram.F90 \
		PDAF_diag_ensstats.F90 \
		PDAF_diag_effsample.F90 \
		PDAF_diag_crps.F90 \
		PDAF_gather_dim_obs_f.F90 \
		PDAF_gather_obs_f.F90 \
		PDAF_gather_obs_f2.F90 \
		PDAF_gather_obs_f_flex.F90 \
		PDAF_gather_obs_f2_flex.F90 \
		PDAF_allreduce.F90 \
		PDAF_deallocate.F90 \
		PDAF_get_assim_flag.F90 \
		PDAF_get_localfilter.F90 \
		PDAF_get_globalobs.F90 \
		PDAFomi_put_state_global.F90 \
		PDAFomi_put_state_global_si.F90 \
		PDAFomi_put_state_global_nondiagR.F90 \
		PDAFomi_put_state_global_nondiagR_si.F90 \
		PDAFomi_put_state_nonlin_nondiagR.F90 \
		PDAFomi_put_state_nonlin_nondiagR_si.F90 \
		PDAFomi_put_state_local.F90 \
		PDAFomi_put_state_local_si.F90 \
		PDAFomi_put_state_local_nondiagR.F90 \
		PDAFomi_put_state_local_nondiagR_si.F90 \
		PDAFomi_assimilate_global.F90 \
		PDAFomi_assimilate_global_si.F90 \
		PDAFomi_assimilate_global_nondiagR.F90 \
		PDAFomi_assimilate_global_nondiagR_si.F90 \
		PDAFomi_assimilate_nonlin_nondiagR.F90 \
		PDAFomi_assimilate_nonlin_nondiagR_si.F90 \
		PDAFomi_assimilate_local.F90 \
		PDAFomi_assimilate_local_si.F90 \
		PDAFomi_assimilate_local_nondiagR.F90 \
		PDAFomi_assimilate_local_nondiagR_si.F90 \
		PDAF_reset_forget.F90 \
		PDAF_reset_dim_p.F90 \
		PDAF_reset_dim_ens.F90 \
		PDAF_get_ensstats.F90 \
		PDAF_set_debug_flag.F90 \
		PDAF_set_offline_mode.F90 \
		PDAFlocal.F90 \
		PDAFlocal_set_indices.F90 \
		PDAFlocal_set_increment_weights.F90 \
		PDAFlocal_clear_increment_weights.F90 \
		PDAFlocal_g2l_cb.F90 \
		PDAFlocal_l2g_cb.F90 \
		PDAFlocalomi_assimilate.F90 \
		PDAFlocalomi_assimilate_nondiagR.F90 \
		PDAFlocalomi_assimilate_nondiagR_si.F90 \
		PDAFlocalomi_assimilate_si.F90 \
		PDAFlocalomi_put_state.F90 \
		PDAFlocalomi_put_state_nondiagR.F90 \
		PDAFlocalomi_put_state_nondiagR_si.F90 \
		PDAFlocalomi_put_state_si.F90 \
		PDAF_correlation_function.F90

# Specific PDAF-routines for SEIK
SRC_SEIK =	PDAF_seik.F90 \
		PDAF_seik_analysis.F90 \
		PDAF_seik_analysis_newT.F90 \
		PDAF_seik_analysis_trans.F90 \
		PDAF_seik_update.F90 \
		PDAF_put_state_seik.F90 \
		PDAF_put_state_seik_si.F90 \
		PDAF_assimilate_seik.F90 \
		PDAF_assimilate_seik_si.F90

# Specific PDAF-routines for local SEIK
SRC_LSEIK =     PDAF_lseik.F90 \
		PDAF_lseik_analysis.F90 \
		PDAF_lseik_analysis_trans.F90 \
		PDAF_lseik_update.F90 \
		PDAF_put_state_lseik.F90 \
		PDAF_put_state_lseik_si.F90 \
		PDAF_assimilate_lseik.F90 \
		PDAF_assimilate_lseik_si.F90 \
		PDAFlocal_put_state_lseik.F90 \
		PDAFlocal_put_state_lseik_si.F90 \
		PDAFlocal_assimilate_lseik.F90 \
		PDAFlocal_assimilate_lseik_si.F90

# Specific PDAF-routines for SEEK
SRC_SEEK =      PDAF_seek.F90 \
		PDAF_seek_analysis.F90 \
		PDAF_seek_update.F90 \
		PDAF_put_state_seek.F90 \
		PDAF_put_state_seek_si.F90 \
		PDAF_assimilate_seek.F90 \
		PDAF_assimilate_seek_si.F90

# Specific PDAF-routines for EnKF
SRC_ENKF =	PDAF_enkf.F90 \
		PDAF_enkf_analysis_rlm.F90 \
		PDAF_enkf_analysis_rsm.F90 \
		PDAF_enkf_update.F90 \
		PDAF_put_state_enkf.F90 \
		PDAF_put_state_enkf_si.F90 \
		PDAF_assimilate_enkf.F90 \
		PDAF_assimilate_enkf_si.F90 \
		PDAFomi_put_state_enkf_nondiagR.F90 \
		PDAFomi_put_state_enkf_nondiagR_si.F90 \
		PDAFomi_assimilate_enkf_nondiagR.F90 \
		PDAFomi_assimilate_enkf_nondiagR_si.F90

# Specific PDAF-routines for ETKF
SRC_ETKF =	PDAF_etkf.F90 \
		PDAF_etkf_analysis.F90 \
		PDAF_etkf_analysis_T.F90 \
		PDAF_etkf_analysis_fixed.F90 \
		PDAF_etkf_update.F90 \
		PDAF_put_state_etkf.F90 \
		PDAF_put_state_etkf_si.F90 \
		PDAF_assimilate_etkf.F90 \
		PDAF_assimilate_etkf_si.F90

# Specific PDAF-routines for LETKF
SRC_LETKF =     PDAF_letkf.F90 \
		PDAF_letkf_analysis.F90 \
		PDAF_letkf_analysis_T.F90 \
		PDAF_letkf_analysis_fixed.F90 \
		PDAF_letkf_update.F90 \
		PDAF_put_state_letkf.F90 \
		PDAF_put_state_letkf_si.F90 \
		PDAF_assimilate_letkf.F90 \
		PDAF_assimilate_letkf_si.F90 \
		PDAFlocal_put_state_letkf.F90 \
		PDAFlocal_put_state_letkf_si.F90 \
		PDAFlocal_assimilate_letkf.F90 \
		PDAFlocal_assimilate_letkf_si.F90

# Specific PDAF-routines for ESTKF
SRC_ESTKF =	PDAF_estkf.F90 \
		PDAF_estkf_analysis.F90 \
		PDAF_estkf_analysis_fixed.F90 \
		PDAF_estkf_update.F90 \
		PDAF_put_state_estkf.F90 \
		PDAF_put_state_estkf_si.F90 \
		PDAF_assimilate_estkf.F90 \
		PDAF_assimilate_estkf_si.F90

# Specific PDAF-routines for LESTKF
SRC_LESTKF =	PDAF_lestkf.F90 \
		PDAF_lestkf_analysis.F90 \
		PDAF_lestkf_analysis_fixed.F90 \
		PDAF_lestkf_update.F90 \
		PDAF_put_state_lestkf.F90 \
		PDAF_put_state_lestkf_si.F90 \
		PDAF_assimilate_lestkf.F90 \
		PDAF_assimilate_lestkf_si.F90 \
		PDAFlocal_put_state_lestkf.F90 \
		PDAFlocal_put_state_lestkf_si.F90 \
		PDAFlocal_assimilate_lestkf.F90 \
		PDAFlocal_assimilate_lestkf_si.F90

# Specific PDAF-routines for LEnKF
SRC_LENKF =	PDAF_lenkf.F90 \
		PDAF_put_state_lenkf.F90 \
		PDAF_put_state_lenkf_si.F90 \
		PDAFomi_put_state_lenkf.F90 \
		PDAFomi_put_state_lenkf_si.F90 \
		PDAFomi_put_state_lenkf_nondiagR.F90 \
		PDAFomi_put_state_lenkf_nondiagR_si.F90 \
		PDAF_assimilate_lenkf.F90 \
		PDAF_assimilate_lenkf_si.F90 \
		PDAFomi_assimilate_lenkf.F90 \
		PDAFomi_assimilate_lenkf_si.F90 \
		PDAFomi_assimilate_lenkf_nondiagR.F90 \
		PDAFomi_assimilate_lenkf_nondiagR_si.F90 \
		PDAF_lenkf_update.F90 \
		PDAF_lenkf_analysis_rsm.F90
# Additional objects used by LEnKF but already specified for EnKF
#		PDAF_enkf_gather_resid.F90
#		PDAF_enkf_obs_ensemble.F90
#		PDAF_enkf_omega.F90
#		PDAF_enkf_Tleft.F90

# Specific PDAF-routines for NETF
SRC_NETF =	PDAF_netf.F90 \
		PDAF_netf_analysis.F90 \
		PDAF_netf_update.F90 \
		PDAF_put_state_netf.F90 \
		PDAF_put_state_netf_si.F90 \
		PDAF_assimilate_netf.F90 \
		PDAF_assimilate_netf_si.F90 \

# Specific PDAF-routines for LNETF
SRC_LNETF =	PDAF_lnetf.F90 \
		PDAF_lnetf_analysis.F90 \
		PDAF_lnetf_update.F90 \
		PDAF_put_state_lnetf.F90 \
		PDAF_put_state_lnetf_si.F90 \
		PDAF_assimilate_lnetf.F90 \
		PDAF_assimilate_lnetf_si.F90 \
		PDAFomi_put_state_lnetf_nondiagR.F90 \
		PDAFomi_put_state_lnetf_nondiagR_si.F90 \
		PDAFomi_assimilate_lnetf_nondiagR.F90 \
		PDAFomi_assimilate_lnetf_nondiagR_si.F90 \
		PDAFlocal_put_state_lnetf.F90 \
		PDAFlocal_put_state_lnetf_si.F90 \
		PDAFlocal_assimilate_lnetf.F90 \
		PDAFlocal_assimilate_lnetf_si.F90 \
		PDAFlocalomi_assimilate_lnetf_nondiagR.F90 \
		PDAFlocalomi_assimilate_lnetf_nondiagR_si.F90 \
		PDAFlocalomi_put_state_lnetf_nondiagR.F90 \
		PDAFlocalomi_put_state_lnetf_nondiagR_si.F90

# Specific PDAF-routines for PF
SRC_PF =	PDAF_pf.F90 \
		PDAF_pf_analysis.F90 \
		PDAF_pf_update.F90 \
		PDAF_put_state_pf.F90 \
		PDAF_put_state_pf_si.F90 \
		PDAF_assimilate_pf.F90 \
		PDAF_assimilate_pf_si.F90

# Specific PDAF-routines for LKNETF
SRC_LKNETF =	PDAF_lknetf.F90 \
		PDAF_lknetf_analysis_step.F90 \
		PDAF_lknetf_update_step.F90 \
		PDAF_lknetf_analysis_sync.F90 \
		PDAF_lknetf_update_sync.F90 \
		PDAF_lknetf_reset_gamma.F90 \
		PDAF_put_state_lknetf.F90 \
		PDAF_put_state_lknetf_si.F90 \
		PDAF_assimilate_lknetf.F90 \
		PDAF_assimilate_lknetf_si.F90 \
		PDAFomi_put_state_lknetf_nondiagR.F90 \
		PDAFomi_put_state_lknetf_nondiagR_si.F90 \
		PDAFomi_assimilate_lknetf_nondiagR.F90 \
		PDAFomi_assimilate_lknetf_nondiagR_si.F90 \
		PDAFlocal_put_state_lknetf.F90 \
		PDAFlocal_put_state_lknetf_si.F90 \
		PDAFlocal_assimilate_lknetf.F90 \
		PDAFlocal_assimilate_lknetf_si.F90 \
		PDAFlocalomi_assimilate_lknetf_nondiagR.F90 \
		PDAFlocalomi_assimilate_lknetf_nondiagR_si.F90 \
		PDAFlocalomi_put_state_lknetf_nondiagR.F90 \
		PDAFlocalomi_put_state_lknetf_nondiagR_si.F90

# Specific PDAF-routines for generating observations
SRC_GENOBS =	PDAF_genobs.F90 \
		PDAF_generate_obs_update.F90 \
		PDAF_put_state_generate_obs.F90 \
		PDAF_put_state_generate_obs_si.F90 \
		PDAFomi_put_state_generate_obs.F90 \
		PDAF_generate_obs.F90 \
		PDAF_generate_obs_si.F90 \
		PDAFomi_generate_obs.F90

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
		PDAF_put_state_3dvar.F90 \
		PDAF_assimilate_3dvar.F90 \
		PDAF_put_state_en3dvar_lestkf.F90 \
		PDAF_put_state_en3dvar_estkf.F90 \
		PDAF_assimilate_en3dvar_lestkf.F90 \
		PDAF_assimilate_en3dvar_estkf.F90 \
		PDAF_put_state_hyb3dvar_estkf.F90 \
		PDAF_put_state_hyb3dvar_lestkf.F90 \
		PDAF_assimilate_hyb3dvar_lestkf.F90 \
		PDAF_assimilate_hyb3dvar_estkf.F90 \
		PDAFlocal_put_state_en3dvar_lestkf.F90 \
		PDAFlocal_put_state_hyb3dvar_lestkf.F90 \
		PDAFlocal_assimilate_en3dvar_lestkf.F90 \
		PDAFlocal_assimilate_hyb3dvar_lestkf.F90 \
		PDAFomi_assimilate_3dvar.F90 \
		PDAFomi_assimilate_en3dvar_estkf.F90 \
		PDAFomi_assimilate_en3dvar_lestkf.F90 \
		PDAFomi_assimilate_hyb3dvar_estkf.F90 \
		PDAFomi_assimilate_hyb3dvar_lestkf.F90 \
		PDAFomi_assimilate_3dvar_nondiagR.F90 \
		PDAFomi_assimilate_en3dvar_estkf_nondiagR.F90 \
		PDAFomi_assimilate_en3dvar_lestkf_nondiagR.F90 \
		PDAFomi_assimilate_hyb3dvar_estkf_nondiagR.F90 \
		PDAFomi_assimilate_hyb3dvar_lestkf_nondiagR.F90 \
		PDAFlocalomi_assimilate_en3dvar_lestkf.F90 \
		PDAFlocalomi_assimilate_en3dvar_lestkf_nondiagR.F90 \
		PDAFlocalomi_assimilate_hyb3dvar_lestkf.F90 \
		PDAFlocalomi_assimilate_hyb3dvar_lestkf_nondiagR.F90 \
		PDAFomi_put_state_3dvar.F90 \
		PDAFomi_put_state_en3dvar_estkf.F90 \
		PDAFomi_put_state_en3dvar_lestkf.F90 \
		PDAFomi_put_state_hyb3dvar_estkf.F90 \
		PDAFomi_put_state_hyb3dvar_lestkf.F90 \
		PDAFomi_put_state_3dvar_nondiagR.F90 \
		PDAFomi_put_state_en3dvar_estkf_nondiagR.F90 \
		PDAFomi_put_state_en3dvar_lestkf_nondiagR.F90 \
		PDAFomi_put_state_hyb3dvar_estkf_nondiagR.F90 \
		PDAFomi_put_state_hyb3dvar_lestkf_nondiagR.F90 \
		PDAFlocalomi_put_state_en3dvar_lestkf.F90 \
		PDAFlocalomi_put_state_en3dvar_lestkf_nondiagR.F90 \
		PDAFlocalomi_put_state_hyb3dvar_lestkf.F90 \
		PDAFlocalomi_put_state_hyb3dvar_lestkf_nondiagR.F90

# Routines for PDAF-OMI
SRC_PDAFOMI =	PDAFomi_obs_f.F90 \
		PDAFomi_obs_l.F90 \
		PDAFomi_dim_obs_l.F90 \
		PDAFomi_obs_op.F90 \
		PDAFomi.F90 \
		PDAFomi_callback.F90

# collect all PDAF sources
SRC_PDAF =  $(SRC_PDAFOMI) $(SRC_PDAF_GEN) \
	    $(SRC_SEIK) $(SRC_LSEIK) $(SRC_SEEK) \
	    $(SRC_ENKF) $(SRC_ETKF) $(SRC_LETKF) \
	    $(SRC_ESTKF) $(SRC_LESTKF) $(SRC_LENKF) $(SRC_NETF) $(SRC_LNETF) \
	    $(SRC_LKNETF) $(SRC_PF) $(SRC_GENOBS) $(SRC_3DVAR) \
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
