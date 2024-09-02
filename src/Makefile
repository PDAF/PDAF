# $Id: Makefile 1856 2017-12-06 08:36:03Z lnerger $

#######################################################
# Generic Makefile for to build PDAF library          #
# To choose the architecture set $PDAF_ARCH           #
#######################################################

.SUFFIXES: .F90 .f .o

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
MOD_PDAF =    	PDAF_timer_mpi.o \
		PDAF_memcount.o \
		PDAF_mod_filtermpi.o \
		PDAF_mod_filter.o

# Module file with interface definitions
MOD_INTERFACE = PDAFlocal_interfaces.o \
		PDAF_interfaces_module.o

# Generic routines in PDAF
OBJ_PDAF_GEN = 	PDAF_analysis_utils.o \
		PDAF_init.o \
		PDAF_init_si.o \
		PDAF_init_filters.o \
		PDAF_alloc_filters.o \
		PDAF_print_info.o \
		PDAF_print_version.o \
		PDAF_communicate_ens.o \
		PDAF_set_comm_pdaf.o \
		PDAF_options_filters.o \
		PDAF_get_state.o \
		PDAF_get_state_si.o \
		PDAF_incremental.o \
		PDAF_incremental_si.o \
		PDAF_set_forget.o \
		PDAF_set_forget_local.o \
		PDAF_add_increment.o \
		PDAF_generate_rndmat.o \
		PDAF_local_weights.o \
		PDAF_local_weight.o \
		PDAF_force_analysis.o \
		PDAF_set_memberid.o \
		PDAF_get_memberid.o \
		PDAF_get_obsmemberid.o \
		PDAF_smoother_shift.o \
		PDAF_smoother.o \
		PDAF_smoother_local.o \
		PDAF_set_smootherens.o \
		PDAF_get_smootherens.o \
		PDAF_set_ens_pointer.o \
		PDAF_put_state_prepost.o \
		PDAF_put_state_prepost_si.o \
		PDAF_assimilate_prepost.o \
		PDAF_assimilate_prepost_si.o \
		PDAF_prepost.o \
		PDAF_prepost_si.o \
		PDAF_inflate_ens.o \
		PDAF_sampleens.o \
		PDAF_mvnormalize.o \
		PDAF_eofcovar.o \
		PDAF_diag_histogram.o \
		PDAF_diag_ensstats.o \
		PDAF_diag_effsample.o \
		PDAF_diag_crps.o \
		PDAF_gather_dim_obs_f.o \
		PDAF_gather_obs_f.o \
		PDAF_gather_obs_f2.o \
		PDAF_gather_obs_f_flex.o \
		PDAF_gather_obs_f2_flex.o \
		PDAF_allreduce.o \
		PDAF_deallocate.o \
		PDAF_get_assim_flag.o \
		PDAF_get_localfilter.o \
		PDAF_get_globalobs.o \
		PDAF_inflate_weights.o \
		PDAFomi_put_state_global.o \
		PDAFomi_put_state_global_si.o \
		PDAFomi_put_state_global_nondiagR.o \
		PDAFomi_put_state_global_nondiagR_si.o \
		PDAFomi_put_state_nonlin_nondiagR.o \
		PDAFomi_put_state_nonlin_nondiagR_si.o \
		PDAFomi_put_state_local.o \
		PDAFomi_put_state_local_si.o \
		PDAFomi_put_state_local_nondiagR.o \
		PDAFomi_put_state_local_nondiagR_si.o \
		PDAFomi_assimilate_global.o \
		PDAFomi_assimilate_global_si.o \
		PDAFomi_assimilate_global_nondiagR.o \
		PDAFomi_assimilate_global_nondiagR_si.o \
		PDAFomi_assimilate_nonlin_nondiagR.o \
		PDAFomi_assimilate_nonlin_nondiagR_si.o \
		PDAFomi_assimilate_local.o \
		PDAFomi_assimilate_local_si.o \
		PDAFomi_assimilate_local_nondiagR.o \
		PDAFomi_assimilate_local_nondiagR_si.o \
		PDAF_reset_forget.o \
		PDAF_get_ensstats.o \
		PDAF_set_debug_flag.o \
		PDAF_set_offline_mode.o \
		PDAFlocal.o \
		PDAFlocal_set_indices.o \
		PDAFlocal_set_increment_weights.o \
		PDAFlocal_clear_increment_weights.o \
		PDAFlocal_g2l_cb.o \
		PDAFlocal_l2g_cb.o \
		PDAFlocalomi_assimilate.o \
		PDAFlocalomi_assimilate_nondiagR.o \
		PDAFlocalomi_assimilate_nondiagR_si.o \
		PDAFlocalomi_assimilate_si.o \
		PDAFlocalomi_put_state.o \
		PDAFlocalomi_put_state_nondiagR.o \
		PDAFlocalomi_put_state_nondiagR_si.o \
		PDAFlocalomi_put_state_si.o \
		PDAF_correlation_function.o


# Specific PDAF-routines for SEIK
OBJ_SEIK =	PDAF_seik_init.o \
		PDAF_seik_alloc.o \
		PDAF_seik_options.o \
		PDAF_seik_memtime.o \
		PDAF_put_state_seik.o \
		PDAF_put_state_seik_si.o \
		PDAF_assimilate_seik.o \
		PDAF_assimilate_seik_si.o \
		PDAF_seik_update.o \
		PDAF_seik_analysis.o \
		PDAF_seik_resample.o \
		PDAF_seik_analysis_newT.o \
		PDAF_seik_resample_newT.o \
		PDAF_seik_analysis_trans.o \
		PDAF_seik_matrixT.o \
		PDAF_seik_uinv.o \
		PDAF_seik_omega.o \
		PDAF_seik_TtimesA.o

# Specific PDAF-routines for local SEIK
OBJ_LSEIK =     PDAF_lseik_init.o \
		PDAF_lseik_alloc.o \
		PDAF_lseik_options.o \
		PDAF_lseik_memtime.o \
		PDAF_put_state_lseik.o \
		PDAF_put_state_lseik_si.o \
		PDAF_assimilate_lseik.o \
		PDAF_assimilate_lseik_si.o \
		PDAF_lseik_update.o \
		PDAF_lseik_analysis.o \
		PDAF_lseik_resample.o \
		PDAF_lseik_analysis_trans.o \
		PDAFlocal_put_state_lseik.o \
		PDAFlocal_put_state_lseik_si.o \
		PDAFlocal_assimilate_lseik.o \
		PDAFlocal_assimilate_lseik_si.o

# Specific PDAF-routines for SEEK
OBJ_SEEK =      PDAF_seek_init.o \
		PDAF_seek_alloc.o \
		PDAF_seek_options.o \
		PDAF_seek_memtime.o \
		PDAF_put_state_seek.o \
		PDAF_put_state_seek_si.o \
		PDAF_assimilate_seek.o \
		PDAF_assimilate_seek_si.o \
		PDAF_seek_update.o \
		PDAF_seek_analysis.o \
		PDAF_seek_rediag.o

# Specific PDAF-routines for EnKF
OBJ_ENKF =	PDAF_enkf_init.o \
		PDAF_enkf_alloc.o \
		PDAF_enkf_options.o \
		PDAF_enkf_memtime.o \
		PDAF_put_state_enkf.o \
		PDAF_put_state_enkf_si.o \
		PDAF_assimilate_enkf.o \
		PDAF_assimilate_enkf_si.o \
		PDAF_enkf_update.o \
		PDAF_enkf_obs_ensemble.o \
		PDAF_enkf_gather_resid.o \
		PDAF_enkf_analysis_rlm.o \
		PDAF_enkf_analysis_rsm.o \
		PDAF_enkf_omega.o \
		PDAF_enkf_Tleft.o \
		PDAF_smoother_enkf.o \
		PDAFomi_put_state_enkf_nondiagR.o \
		PDAFomi_put_state_enkf_nondiagR_si.o \
		PDAFomi_assimilate_enkf_nondiagR.o \
		PDAFomi_assimilate_enkf_nondiagR_si.o

# Specific PDAF-routines for ETKF
OBJ_ETKF =	PDAF_etkf_init.o \
		PDAF_etkf_alloc.o \
		PDAF_etkf_options.o \
		PDAF_etkf_memtime.o \
		PDAF_put_state_etkf.o \
		PDAF_put_state_etkf_si.o \
		PDAF_assimilate_etkf.o \
		PDAF_assimilate_etkf_si.o \
		PDAF_etkf_update.o \
		PDAF_etkf_analysis.o \
		PDAF_etkf_analysis_T.o \
		PDAF_etkf_analysis_fixed.o \
		PDAF_etkf_Tright.o \
		PDAF_etkf_Tleft.o

# Specific PDAF-routines for LETKF
OBJ_LETKF =     PDAF_letkf_init.o \
		PDAF_letkf_alloc.o \
		PDAF_letkf_options.o \
		PDAF_letkf_memtime.o \
		PDAF_put_state_letkf.o \
		PDAF_put_state_letkf_si.o \
		PDAF_assimilate_letkf.o \
		PDAF_assimilate_letkf_si.o \
		PDAF_letkf_update.o \
		PDAF_letkf_analysis.o \
		PDAF_letkf_analysis_T.o \
		PDAF_letkf_analysis_fixed.o \
		PDAFlocal_put_state_letkf.o \
		PDAFlocal_put_state_letkf_si.o \
		PDAFlocal_assimilate_letkf.o \
		PDAFlocal_assimilate_letkf_si.o

# Specific PDAF-routines for ESTKF
OBJ_ESTKF =	PDAF_estkf_init.o \
		PDAF_estkf_alloc.o \
		PDAF_estkf_options.o \
		PDAF_estkf_memtime.o \
		PDAF_put_state_estkf.o \
		PDAF_put_state_estkf_si.o \
		PDAF_assimilate_estkf.o \
		PDAF_assimilate_estkf_si.o \
		PDAF_estkf_update.o \
		PDAF_estkf_analysis.o \
		PDAF_estkf_analysis_fixed.o \
		PDAF_estkf_AOmega.o \
		PDAF_estkf_OmegaA.o

# Specific PDAF-routines for LESTKF
OBJ_LESTKF =	PDAF_lestkf_init.o \
		PDAF_lestkf_alloc.o \
		PDAF_lestkf_options.o \
		PDAF_lestkf_memtime.o \
		PDAF_put_state_lestkf.o \
		PDAF_put_state_lestkf_si.o \
		PDAF_assimilate_lestkf.o \
		PDAF_assimilate_lestkf_si.o \
		PDAF_lestkf_update.o \
		PDAF_lestkf_analysis.o \
		PDAF_lestkf_analysis_fixed.o \
		PDAFlocal_put_state_lestkf.o \
		PDAFlocal_put_state_lestkf_si.o \
		PDAFlocal_assimilate_lestkf.o \
		PDAFlocal_assimilate_lestkf_si.o

# Specific PDAF-routines for LEnKF
OBJ_LENKF =	PDAF_lenkf_init.o \
		PDAF_lenkf_alloc.o \
		PDAF_lenkf_options.o \
		PDAF_lenkf_memtime.o \
		PDAF_put_state_lenkf.o \
		PDAF_put_state_lenkf_si.o \
		PDAFomi_put_state_lenkf.o \
		PDAFomi_put_state_lenkf_si.o \
		PDAFomi_put_state_lenkf_nondiagR.o \
		PDAFomi_put_state_lenkf_nondiagR_si.o \
		PDAF_assimilate_lenkf.o \
		PDAF_assimilate_lenkf_si.o \
		PDAFomi_assimilate_lenkf.o \
		PDAFomi_assimilate_lenkf_si.o \
		PDAFomi_assimilate_lenkf_nondiagR.o \
		PDAFomi_assimilate_lenkf_nondiagR_si.o \
		PDAF_lenkf_update.o \
		PDAF_lenkf_analysis_rsm.o
# Additional objects used by LEnKF but already specified for EnKF
#		PDAF_enkf_gather_resid.o
#		PDAF_enkf_obs_ensemble.o
#		PDAF_enkf_omega.o
#		PDAF_enkf_Tleft.o

# Specific PDAF-routines for NETF
OBJ_NETF =	PDAF_netf_init.o \
		PDAF_netf_alloc.o \
		PDAF_netf_options.o \
		PDAF_netf_memtime.o \
		PDAF_put_state_netf.o \
		PDAF_put_state_netf_si.o \
		PDAF_assimilate_netf.o \
		PDAF_assimilate_netf_si.o \
		PDAF_netf_update.o \
		PDAF_netf_smootherT.o \
		PDAF_netf_analysis.o \
		PDAF_smoother_netf.o

# Specific PDAF-routines for LNETF
OBJ_LNETF =	PDAF_lnetf_init.o \
		PDAF_lnetf_alloc.o \
		PDAF_lnetf_options.o \
		PDAF_lnetf_memtime.o \
		PDAF_put_state_lnetf.o \
		PDAF_put_state_lnetf_si.o \
		PDAF_assimilate_lnetf.o \
		PDAF_assimilate_lnetf_si.o \
		PDAF_lnetf_update.o \
		PDAF_lnetf_analysis.o \
		PDAF_lnetf_smootherT.o \
		PDAF_smoother_lnetf.o \
		PDAFomi_put_state_lnetf_nondiagR.o \
		PDAFomi_put_state_lnetf_nondiagR_si.o \
		PDAFomi_assimilate_lnetf_nondiagR.o \
		PDAFomi_assimilate_lnetf_nondiagR_si.o \
		PDAFlocal_put_state_lnetf.o \
		PDAFlocal_put_state_lnetf_si.o \
		PDAFlocal_assimilate_lnetf.o \
		PDAFlocal_assimilate_lnetf_si.o \
		PDAFlocalomi_assimilate_lnetf_nondiagR.o \
		PDAFlocalomi_assimilate_lnetf_nondiagR_si.o \
		PDAFlocalomi_put_state_lnetf_nondiagR.o \
		PDAFlocalomi_put_state_lnetf_nondiagR_si.o

# Specific PDAF-routines for PF
OBJ_PF =	PDAF_pf_init.o \
		PDAF_pf_alloc.o \
		PDAF_pf_options.o \
		PDAF_pf_memtime.o \
		PDAF_put_state_pf.o \
		PDAF_put_state_pf_si.o \
		PDAF_assimilate_pf.o \
		PDAF_assimilate_pf_si.o \
		PDAF_pf_update.o \
		PDAF_pf_analysis.o \
		PDAF_pf_resampling.o \
		PDAF_pf_add_noise.o

# Specific PDAF-routines for LKNETF
OBJ_LKNETF =	PDAF_lknetf_init.o \
		PDAF_lknetf_alloc.o \
		PDAF_lknetf_options.o \
		PDAF_lknetf_memtime.o \
		PDAF_put_state_lknetf.o \
		PDAF_put_state_lknetf_si.o \
		PDAF_assimilate_lknetf.o \
		PDAF_assimilate_lknetf_si.o \
		PDAF_lknetf_update.o \
		PDAF_lknetf_analysis_T.o \
		PDAF_lknetf_step_update.o \
		PDAF_lknetf_ana_lnetf.o \
		PDAF_lknetf_ana_letkfT.o \
		PDAF_lknetf_compute_gamma.o \
		PDAF_lknetf_set_gamma.o \
		PDAF_lknetf_alpha_neff.o \
		PDAF_lknetf_reset_gamma.o \
		PDAFomi_put_state_lknetf_nondiagR.o \
		PDAFomi_put_state_lknetf_nondiagR_si.o \
		PDAFomi_assimilate_lknetf_nondiagR.o \
		PDAFomi_assimilate_lknetf_nondiagR_si.o \
		PDAFlocal_put_state_lknetf.o \
		PDAFlocal_put_state_lknetf_si.o \
		PDAFlocal_assimilate_lknetf.o \
		PDAFlocal_assimilate_lknetf_si.o \
		PDAFlocalomi_assimilate_lknetf_nondiagR.o \
		PDAFlocalomi_assimilate_lknetf_nondiagR_si.o \
		PDAFlocalomi_put_state_lknetf_nondiagR.o \
		PDAFlocalomi_put_state_lknetf_nondiagR_si.o

# Specific PDAF-routines for generating observations
OBJ_OBSGEN =	PDAF_genobs_init.o \
		PDAF_genobs_alloc.o \
		PDAF_genobs_options.o \
		PDAF_put_state_generate_obs.o \
		PDAF_put_state_generate_obs_si.o \
		PDAFomi_put_state_generate_obs.o \
		PDAF_generate_obs.o \
		PDAF_generate_obs_si.o \
		PDAFomi_generate_obs.o \
		PDAF_gen_obs.o

# Specific PDAF-routines for 3DVAR initialization part
OBJ_3DVAR_INI =	PDAF_3dvar_init.o \
		PDAF_3dvar_alloc.o \
		PDAF_3dvar_options.o \
		PDAF_3dvar_memtime.o

# Specific PDAF-routines for 3DVAR
OBJ_3DVAR =	PDAF_put_state_3dvar.o \
		PDAF_assimilate_3dvar.o \
		PDAF_3dvar_update.o \
		PDAF_3dvar_analysis_cvt.o \
		PDAF_3dvar_optim_lbfgs.o \
		PDAF_3dvar_optim_cgplus.o \
		PDAF_3dvar_costf_cvt.o \
		PDAF_3dvar_costf_cg_cvt.o \
		PDAF_3dvar_optim_cg.o \
		PDAF_put_state_en3dvar_lestkf.o \
		PDAF_assimilate_en3dvar_lestkf.o \
		PDAF_en3dvar_update_lestkf.o \
		PDAF_put_state_en3dvar_estkf.o \
		PDAF_assimilate_en3dvar_estkf.o \
		PDAF_en3dvar_update_estkf.o \
		PDAF_en3dvar_analysis_cvt.o \
		PDAF_en3dvar_optim_lbfgs.o \
		PDAF_en3dvar_optim_cgplus.o \
		PDAF_en3dvar_optim_cg.o \
		PDAF_en3dvar_costf_cvt.o \
		PDAF_en3dvar_costf_cg_cvt.o \
		PDAF_put_state_hyb3dvar_lestkf.o \
		PDAF_assimilate_hyb3dvar_lestkf.o \
		PDAF_hyb3dvar_update_lestkf.o \
		PDAF_hyb3dvar_analysis_cvt.o\
		PDAF_put_state_hyb3dvar_estkf.o \
		PDAF_assimilate_hyb3dvar_estkf.o \
		PDAF_hyb3dvar_update_estkf.o \
		PDAF_hyb3dvar_optim_lbfgs.o \
		PDAF_hyb3dvar_optim_cgplus.o \
		PDAF_hyb3dvar_optim_cg.o \
		PDAF_hyb3dvar_costf_cvt.o \
		PDAF_hyb3dvar_costf_cg_cvt.o \
		PDAFlocal_put_state_en3dvar_lestkf.o \
		PDAFlocal_put_state_hyb3dvar_lestkf.o \
		PDAFlocal_assimilate_en3dvar_lestkf.o \
		PDAFlocal_assimilate_hyb3dvar_lestkf.o \
		PDAFomi_assimilate_3dvar.o \
		PDAFomi_assimilate_en3dvar_estkf.o \
		PDAFomi_assimilate_en3dvar_lestkf.o \
		PDAFomi_assimilate_hyb3dvar_estkf.o \
		PDAFomi_assimilate_hyb3dvar_lestkf.o \
		PDAFomi_assimilate_3dvar_nondiagR.o \
		PDAFomi_assimilate_en3dvar_estkf_nondiagR.o \
		PDAFomi_assimilate_en3dvar_lestkf_nondiagR.o \
		PDAFomi_assimilate_hyb3dvar_estkf_nondiagR.o \
		PDAFomi_assimilate_hyb3dvar_lestkf_nondiagR.o \
		PDAFlocalomi_assimilate_en3dvar_lestkf.o \
		PDAFlocalomi_assimilate_en3dvar_lestkf_nondiagR.o \
		PDAFlocalomi_assimilate_hyb3dvar_lestkf.o \
		PDAFlocalomi_assimilate_hyb3dvar_lestkf_nondiagR.o \
		PDAFomi_put_state_3dvar.o \
		PDAFomi_put_state_en3dvar_estkf.o \
		PDAFomi_put_state_en3dvar_lestkf.o \
		PDAFomi_put_state_hyb3dvar_estkf.o \
		PDAFomi_put_state_hyb3dvar_lestkf.o \
		PDAFomi_put_state_3dvar_nondiagR.o \
		PDAFomi_put_state_en3dvar_estkf_nondiagR.o \
		PDAFomi_put_state_en3dvar_lestkf_nondiagR.o \
		PDAFomi_put_state_hyb3dvar_estkf_nondiagR.o \
		PDAFomi_put_state_hyb3dvar_lestkf_nondiagR.o \
		PDAFlocalomi_put_state_en3dvar_lestkf.o \
		PDAFlocalomi_put_state_en3dvar_lestkf_nondiagR.o \
		PDAFlocalomi_put_state_hyb3dvar_lestkf.o \
		PDAFlocalomi_put_state_hyb3dvar_lestkf_nondiagR.o
# Additional file for 3DVar already specified in OBJ_3DVAR_ini
#		PDAF_3dvar_memtime.o

# Routines for PDAF-OMI
OBJ_PDAFOMI =	PDAFomi_obs_f.o \
		PDAFomi_obs_l.o \
		PDAFomi_dim_obs_l.o \
		PDAFomi_obs_op.o \
		PDAFomi.o \
		PDAFomi_callback.o

OBJ_PDAF =   $(OBJ_PDAFOMI) $(OBJ_PDAF_GEN) $(OBJ_SEIK) $(OBJ_LSEIK) $(OBJ_SEEK) \
		$(OBJ_ENKF) $(OBJ_ETKF) $(OBJ_LETKF) \
		$(OBJ_ESTKF) $(OBJ_LESTKF) $(OBJ_LENKF) $(OBJ_NETF) $(OBJ_LNETF) \
		$(OBJ_LKNETF) $(OBJ_PF) $(OBJ_OBSGEN) $(OBJ_3DVAR_INI)

OBJ_PDAF_VAR = $(OBJ_PDAF) $(OBJ_3DVAR)

# External Routines for SANGOMA tools
OBJ_SANGOMA = ../external/SANGOMA/SANGOMA_quicksort.o

# External optimizer libraries
OBJ_OPTIM = ../external/CG+_mpi/cgfam.o ../external/CG+_mpi/cgsearch.o \
	../external/CG+/cgfam.o ../external/CG+/cgsearch.o \
	../external/LBFGS/lbfgsb.o ../external/LBFGS/linpack.o \
	../external/LBFGS/timer.o

######################################################

../lib/libpdaf-d.a: 	$(MOD_PDAF) $(OBJ_SANGOMA) $(MOD_INTERFACE) $(OBJ_PDAF)
	@echo "++++++ Generate Filter library ++++++"
	$(AR) -r $(AR_SPEC) $@ \
	$(MOD_PDAF) $(MOD_INTERFACE) $(OBJ_PDAF) $(OBJ_SANGOMA)
	$(RANLIB) ../lib/libpdaf-d.a
	@cp *.mod ../include

../lib/libpdaf-var.a: $(MOD_PDAF)  $(OBJ_SANGOMA) $(MOD_INTERFACE) $(OBJ_PDAF_VAR) $(OBJ_OPTIM)
	@echo "++++++ Generate Filter library ++++++"
	$(AR) -r $(AR_SPEC) $@ \
	$(MOD_PDAF) $(MOD_INTERFACE) $(OBJ_PDAF) $(OBJ_PDAF_VAR) $(OBJ_OPTIM) $(OBJ_SANGOMA)
	$(RANLIB) ../lib/libpdaf-var.a
	@cp *.mod ../include

# Short alias for ../lib/libpdaf-var.a
pdaf-var: ../lib/libpdaf-var.a


.F90.o :
	$(FC) $(OPT) $(MPI_INC) $(CPP_DEFS) -o $*.o -c $*.F90

.f.o :
	$(FC) -O3 $(MPI_INC) -o $*.o -c $*.f

# For older compilers one might need to separate the
# preprocessing from the compilation as defined below:
#.F90.o :
#	$(CPP) -P -C $(MPI_INC) $(CPP_DEFS) $*.F90 $*.f90
#	$(FC) $(OPT) $(MPI_INC) -c $*.F90
#	@rm -f  $*.f90

######################################################
# Cleans

clean :
	rm -f *.o *.mod ../lib/libpdaf-d.a ../lib/libpdaf-var.a ../include/*.mod
	cd ../external/SANGOMA; rm -f *.o *.mod; cd -
	cd ../external/CG+; rm -f *.o; cd -
	cd ../external/CG+_mpi; rm -f *.o; cd -
	cd ../external/LBFGS; rm -f *.o; cd -

######################################################
# List arch files

listarch:
	@echo Available architecture-specific input files for PDAF_ARCH
	@echo ---------------------------------------------------------
	@ls -1 ../make.arch | cut -d"." -f1
