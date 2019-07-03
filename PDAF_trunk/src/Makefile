# $Id: Makefile 1856 2017-12-06 08:36:03Z lnerger $

#######################################################
# Generic Makefile for to build PDAF library          #
# To choose the architecture set $PDAF_ARCH           #
#######################################################

.SUFFIXES: .F90 .o 

######################################################

# Include machine-specific definitions
# For available include files see directory make.arch
# To choose a file, set PDAF_ARCH either here or by an
# environment variable.
include ../make.arch/$(PDAF_ARCH).h


######################################################
# Define objects for PDAF library
######################################################

# Modules used in PDAF
MOD_PDAF =    	PDAF_timer.o \
		PDAF_memcount.o \
		PDAF-D_mod_filtermpi.o \
		PDAF-D_mod_filter.o

# Module file with interface definitions
MOD_INTERFACE = PDAF_interfaces_module.o

# Generic routines in PDAF
OBJ_PDAF_GEN = 	PDAF-D_init.o \
		PDAF-D_init_si.o \
		PDAF-D_init_filters.o \
		PDAF-D_alloc_filters.o \
		PDAF-D_print_info.o \
		PDAF-D_print_version.o \
		PDAF_options_filters.o \
		PDAF-D_get_state.o \
		PDAF-D_get_state_si.o \
		PDAF-D_incremental.o \
		PDAF-D_incremental_si.o \
		PDAF-D_set_forget.o \
		PDAF-D_set_forget_local.o \
		PDAF-D_add_increment.o \
		PDAF_generate_rndmat.o \
		PDAF_local_weights.o \
		PDAF_local_weight.o \
		PDAF-D_gather_ens.o \
		PDAF-D_get_memberid.o \
		PDAF-D_get_obsmemberid.o \
		PDAF-D_smoother_shift.o \
		PDAF-D_smoother.o \
		PDAF-D_smoother_local.o \
		PDAF-D_set_smootherens.o \
		PDAF-D_get_smootherens.o \
		PDAF_set_ens_pointer.o \
		PDAF-D_put_state_prepost.o \
		PDAF-D_put_state_prepost_si.o \
		PDAF-D_prepost.o \
		PDAF-D_prepost_si.o \
		PDAF_inflate_ens.o \
		PDAF_sampleens.o \
		PDAF_mvnormalize.o \
		PDAF_eofcovar.o \
		PDAF_diag_histogram.o \
		PDAF_diag_ensstats.o \
		PDAF_diag_effsample.o \
		PDAF-D_gather_dim_obs_f.o \
		PDAF-D_gather_obs_f.o \
		PDAF-D_gather_obs_f2.o \
		PDAF-D_gather_obs_f_flex.o \
		PDAF-D_gather_obs_f2_flex.o \
		PDAF-D_allreduce.o \
		PDAF-D_deallocate.o \
		PDAF-D_get_assim_flag.o

# Specific PDAF-routines for SEIK
OBJ_SEIK =	PDAF-D_seik_init.o \
		PDAF-D_seik_alloc.o \
		PDAF_seik_options.o \
		PDAF-D_seik_memtime.o \
		PDAF-D_put_state_seik.o \
		PDAF-D_put_state_seik_si.o \
		PDAF-D_assimilate_seik.o \
		PDAF-D_assimilate_seik_si.o \
		PDAF-D_seik_update.o \
		PDAF-D_seik_analysis.o \
		PDAF-D_seik_resample.o \
		PDAF-D_seik_analysis_newT.o \
		PDAF-D_seik_resample_newT.o \
		PDAF-D_seik_analysis_trans.o \
		PDAF-D_seik_matrixT.o \
		PDAF-D_seik_uinv.o \
		PDAF_seik_omega.o \
		PDAF_seik_TtimesA.o

# Specific PDAF-routines for local SEIK
OBJ_LSEIK =     PDAF-D_lseik_init.o \
		PDAF-D_lseik_alloc.o \
		PDAF_lseik_options.o \
		PDAF-D_lseik_memtime.o \
		PDAF-D_put_state_lseik.o \
		PDAF-D_put_state_lseik_si.o \
		PDAF-D_assimilate_lseik.o \
		PDAF-D_assimilate_lseik_si.o \
		PDAF-D_lseik_update.o \
		PDAF-D_lseik_analysis.o \
		PDAF-D_lseik_resample.o \
		PDAF-D_lseik_analysis_trans.o

# Specific PDAF-routines for SEEK
OBJ_SEEK =      PDAF-D_seek_init.o \
		PDAF-D_seek_alloc.o \
		PDAF_seek_options.o \
		PDAF-D_seek_memtime.o \
		PDAF-D_put_state_seek.o \
		PDAF-D_put_state_seek_si.o \
		PDAF-D_assimilate_seek.o \
		PDAF-D_assimilate_seek_si.o \
		PDAF-D_seek_update.o \
		PDAF-D_seek_analysis.o \
		PDAF-D_seek_rediag.o

# Specific PDAF-routines for EnKF
OBJ_ENKF =	PDAF-D_enkf_init.o \
		PDAF-D_enkf_alloc.o \
		PDAF_enkf_options.o \
		PDAF-D_enkf_memtime.o \
		PDAF-D_put_state_enkf.o \
		PDAF-D_put_state_enkf_si.o \
		PDAF-D_assimilate_enkf.o \
		PDAF-D_assimilate_enkf_si.o \
		PDAF-D_enkf_update.o \
		PDAF-D_enkf_obs_ensemble.o \
		PDAF-D_enkf_gather_resid.o \
		PDAF-D_enkf_analysis_rlm.o \
		PDAF-D_enkf_analysis_rsm.o \
		PDAF_enkf_omega.o \
		PDAF_enkf_Tleft.o \
		PDAF-D_smoother_enkf.o

# Specific PDAF-routines for ETKF
OBJ_ETKF =	PDAF-D_etkf_init.o \
		PDAF-D_etkf_alloc.o \
		PDAF_etkf_options.o \
		PDAF-D_etkf_memtime.o \
		PDAF-D_put_state_etkf.o \
		PDAF-D_put_state_etkf_si.o \
		PDAF-D_assimilate_etkf.o \
		PDAF-D_assimilate_etkf_si.o \
		PDAF-D_etkf_update.o \
		PDAF-D_etkf_analysis.o \
		PDAF-D_etkf_analysis_T.o \
		PDAF-D_etkf_analysis_fixed.o \
		PDAF-D_etkf_Tright.o \
		PDAF_etkf_Tleft.o

# Specific PDAF-routines for LETKF
OBJ_LETKF =     PDAF-D_letkf_init.o \
		PDAF-D_letkf_alloc.o \
		PDAF_letkf_options.o \
		PDAF-D_letkf_memtime.o \
		PDAF-D_put_state_letkf.o \
		PDAF-D_put_state_letkf_si.o \
		PDAF-D_assimilate_letkf.o \
		PDAF-D_assimilate_letkf_si.o \
		PDAF-D_letkf_update.o \
		PDAF-D_letkf_analysis.o \
		PDAF-D_letkf_analysis_T.o \
		PDAF-D_letkf_analysis_fixed.o

# Specific PDAF-routines for ESTKF
OBJ_ESTKF =	PDAF-D_estkf_init.o \
		PDAF-D_estkf_alloc.o \
		PDAF_estkf_options.o \
		PDAF-D_estkf_memtime.o \
		PDAF-D_put_state_estkf.o \
		PDAF-D_put_state_estkf_si.o \
		PDAF-D_assimilate_estkf.o \
		PDAF-D_assimilate_estkf_si.o \
		PDAF-D_estkf_update.o \
		PDAF-D_estkf_analysis.o \
		PDAF-D_estkf_analysis_fixed.o \
		PDAF-D_estkf_AOmega.o \
		PDAF-D_estkf_OmegaA.o

# Specific PDAF-routines for LESTKF
OBJ_LESTKF =	PDAF-D_lestkf_init.o \
		PDAF-D_lestkf_alloc.o \
		PDAF_lestkf_options.o \
		PDAF-D_lestkf_memtime.o \
		PDAF-D_put_state_lestkf.o \
		PDAF-D_put_state_lestkf_si.o \
		PDAF-D_assimilate_lestkf.o \
		PDAF-D_assimilate_lestkf_si.o \
		PDAF-D_lestkf_update.o \
		PDAF-D_lestkf_analysis.o \
		PDAF-D_lestkf_analysis_fixed.o

# Specific PDAF-routines for EnKF
OBJ_LENKF =	PDAF-D_lenkf_init.o \
		PDAF-D_lenkf_alloc.o \
		PDAF_lenkf_options.o \
		PDAF-D_lenkf_memtime.o \
		PDAF-D_put_state_lenkf.o \
		PDAF-D_put_state_lenkf_si.o \
		PDAF-D_assimilate_lenkf.o \
		PDAF-D_assimilate_lenkf_si.o \
		PDAF-D_lenkf_update.o \
		PDAF-D_enkf_obs_ensemble.o \
		PDAF-D_enkf_gather_resid.o \
		PDAF-D_lenkf_analysis_rsm.o \
		PDAF_enkf_omega.o \
		PDAF_enkf_Tleft.o

# Specific PDAF-routines for NETF
OBJ_NETF =	PDAF-D_netf_init.o \
		PDAF-D_netf_alloc.o \
		PDAF_netf_options.o \
		PDAF-D_netf_memtime.o \
		PDAF-D_put_state_netf.o \
		PDAF-D_put_state_netf_si.o \
		PDAF-D_assimilate_netf.o \
		PDAF-D_assimilate_netf_si.o \
		PDAF-D_netf_update.o \
		PDAF-D_netf_smootherT.o \
		PDAF-D_netf_analysis.o \
		PDAF-D_smoother_netf.o

# Specific PDAF-routines for LNETF
OBJ_LNETF =	PDAF-D_lnetf_init.o \
		PDAF-D_lnetf_alloc.o \
		PDAF_lnetf_options.o \
		PDAF-D_lnetf_memtime.o \
		PDAF-D_put_state_lnetf.o \
		PDAF-D_put_state_lnetf_si.o \
		PDAF-D_assimilate_lnetf.o \
		PDAF-D_assimilate_lnetf_si.o \
		PDAF-D_lnetf_update.o \
		PDAF-D_lnetf_analysis.o \
		PDAF-D_lnetf_smootherT.o \
		PDAF-D_smoother_lnetf.o

# Specific PDAF-routines for 
OBJ_PF =	PDAF-D_pf_init.o \
		PDAF-D_pf_alloc.o \
		PDAF_pf_options.o \
		PDAF-D_pf_memtime.o \
		PDAF-D_put_state_pf.o \
		PDAF-D_put_state_pf_si.o \
		PDAF-D_assimilate_pf.o \
		PDAF-D_assimilate_pf_si.o \
		PDAF-D_pf_update.o \
		PDAF-D_pf_analysis.o \
		PDAF-D_pf_resampling.o \
		PDAF-D_pf_add_noise.o

# Specific PDAF-routines for generating observations
OBJ_OBSGEN =	PDAF-D_genobs_init.o \
		PDAF-D_genobs_alloc.o \
		PDAF_genobs_options.o \
		PDAF-D_put_state_generate_obs.o \
		PDAF-D_put_state_generate_obs_si.o \
		PDAF-D_generate_obs.o \
		PDAF-D_generate_obs_si.o \
		PDAF-D_gen_obs.o

OBJ_PDAF =   $(OBJ_PDAF_GEN) $(OBJ_SEIK) $(OBJ_LSEIK) $(OBJ_SEEK) \
		$(OBJ_ENKF) $(OBJ_ETKF) $(OBJ_LETKF) \
		$(OBJ_ESTKF) $(OBJ_LESTKF) $(OBJ_LENKF) $(OBJ_NETF) $(OBJ_LNETF) \
		$(OBJ_PF) $(OBJ_OBSGEN)

######################################################

../lib/libpdaf-d.a: 	$(MOD_PDAF) $(MOD_INTERFACE) $(OBJ_PDAF)
	@echo "++++++ Generate Filter library ++++++"
	$(AR) -r $(AR_SPEC) $@ \
	$(MOD_PDAF) $(MOD_INTERFACE) $(OBJ_PDAF)
	$(RANLIB) ../lib/libpdaf-d.a

.F90.o :
	$(FC) $(OPT) $(MPI_INC) $(CPP_DEFS) -c $*.F90

# For older compilers one might need to separate the
# preprocessing from the compilation as defined below:
#.F90.o :
#	$(CPP) -P -C $(MPI_INC) $(CPP_DEFS) $*.F90 $*.f90
#	$(FC) $(OPT) $(MPI_INC) -c $*.F90
#	@rm -f  $*.f90

######################################################
# Cleans

clean :
	rm -f *.o *.mod ../lib/libpdaf-d.a

######################################################
# List arch files

listarch:
	@echo Available architecture-specific input files for PDAF_ARCH
	@echo ---------------------------------------------------------
	@ls -1 ../make.arch | cut -d"." -f1
