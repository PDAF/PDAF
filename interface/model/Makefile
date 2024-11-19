BASEDIR = $(PDAF_DIR)
include $(BASEDIR)/make.arch/$(PDAF_ARCH).h
## includes

FFLAGS = -cpp $(TSMPPDAFIMPORTFLAGS)
CCFLAGS = $(TSMPPDAFIMPORTFLAGS)
CPP_FLAGS =  ${TSMPPDAFCPPDEFS}
FCPP_FLAGS= ${TSMPPDAFCPPDEFS}

# Optimization specs for C/C++ compiler
#   C/C++ compiler do not get the floating point definition
COPT= ${TSMPPDAFCOPT}

## common object files
OBJ =  dictionary.o\
	   iniparser.o\
	   read_enkfpar.o\
	   wrapper_tsmp.o\

## clm object files
OBJCLM = enkf_clm_mod.o\
		 mod_clm_statistics.o\
		 print_update_clm.o\
		 enkf_clm.o\
## clm 5 object
OBJCLM5 = enkf_clm_mod_5.o\
		mod_clm_statistics_5.o\
		print_update_clm_5.o\
		enkf_clm_5.o\

## parflow object files
OBJPF = enkf_parflow.o\
		enkf_ensemblestatistics.o\
		problem_saturationtopressure.o\

## cosmo object files
OBJCOSMO = enkf_cosmo_mod.o\
		   enkf_cosmo.o\

## object files
ifdef $(TSMPPDAFOBJ)

# Object files supplied in environmentvariable (Bash)
OBJ = $(TSMPPDAFOBJ)
$(info ************  TSMPPDAFOBJ ************)

else

# Object files set from input sources (CMake)
ifdef CLM35_SRC
OBJ += $(OBJCLM)
$(info ************  CLM35_SRC ************)
endif
ifdef eCLM_SRC
OBJ += $(OBJCLM5)
$(info ************  eCLM_SRC ************)
endif
ifdef PARFLOW_SRC
$(info ************  PARFLOW_SRC ************)
OBJ += $(OBJPF)
endif
ifdef COSMO_SRC
$(info ************  COSMO_SRC ************)
OBJ += $(OBJCOSMO)
endif

endif

$(info $$OBJ is [${OBJ}])

## target + rules
all: libmodel.a

libmodel.a: $(OBJ)
	ar rcs $@ *.o 
	ranlib $@
	mv $@ $(TSMPPDAFLIBDIR)
	@echo "library compilation suceeded"

wrapper_tsmp.o:
	$(PREP_C) $(CC) $(COPT) $(CPP_FLAGS) $(CCFLAGS) $(CINCS) -c wrapper_tsmp.c

%.o: common/%.c
	$(PREP_C) $(CC) $(COPT) $(CPP_FLAGS) $(CCFLAGS) $(CINCS) -c $<  -o $@

%.o: $(TSMPPDAFPFLDIR)/%.c
	$(PREP_C) $(CC) $(COPT) $(CPP_FLAGS) $(CCFLAGS) $(CINCS) -c $<  -o $@

%.o:$(TSMPPDAFCLMDIR)/%.F90
	$(PREP_C) $(FC) $(OPT) $(FCPP_FLAGS) $(FFLAGS) $(FINCS) -c $< -o $@

%.o: $(TSMPPDAFCOSDIR)/%.F90
	$(PREP_C) $(FC) $(OPT) $(FCPP_FLAGS) $(FFLAGS) $(FINCS) -c $< -o $@

clean:
	rm -f *.o *.mod $(PROG) libmodel.a 
