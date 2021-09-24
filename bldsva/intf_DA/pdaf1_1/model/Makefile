BASEDIR = $(PDAF_DIR)
include $(BASEDIR)/make.arch/$(PDAF_ARCH).h
PF =__pf__
## includes

FFLAGS =  __fflags__
CCFLAGS = __ccflags__
CPP_FLAGS = -Duse_libMPI -Duse_netCDF -Duse_comm_MPI1 -DVERBOSE -DDEBUG -DTREAT_OVERLAY -DFORTRANUNDERSCORE -DOFFLINE -DSPMD -DLINUX __cpp_defs__
FCPP_FLAGS= $(PF)-Duse_libMPI $(PF)-Duse_netCDF $(PF)-Duse_comm_MPI1 $(PF)-DVERBOSE $(PF)-DDEBUG $(PF)-DTREAT_OVERLAY $(PF)-DFORTRANUNDERSCORE $(PF)-DOFFLINE $(PF)-DSPMD $(PF)-DLINUX __fcpp_defs__


## common object files
OBJ =  dictionary.o\
	   iniparser.o\
	   read_enkfpar.o\
	   wrapper_tsmp.o\
       #print_memusage.o\

## clm object files
OBJCLM = enkf_clm_mod.o\
		 mod_clm_statistics.o\
		 enkf_clm.o\

## parflow object files
OBJPF = enkf_parflow.o\
		enkf_ensemblestatistics.o\
		problem_saturationtopressure.o\

## cosmo object files
OBJCOSMO = enkf_cosmo_mod.o\
		   enkf_cosmo.o\



OBJ += __obj__

## target + rules
all: libmodel.a

libmodel.a: $(OBJ)
	ar rcs $@ *.o 
	ranlib $@
	mv $@ __bindir__/libs
	@echo "library compilation suceeded"

%.o: %.c
	$(PREP_C) $(CC) $(CCFLAGS) $(CPP_FLAGS) $(CINCS) -c $<  -o $@

%.o: __pfldir__/%.c
	$(PREP_C) $(CC) $(CCFLAGS) $(CPP_FLAGS) $(CINCS) -c $<  -o $@

%.o: __clmdir__/%.F90
	$(PREP_C) $(FC) $(FCPP_FLAGS) $(FFLAGS) $(FINCS) -c $< -o $@

%.o: __cosdir__/%.F90
	$(PREP_C) $(FC) $(FCPP_FLAGS) $(FFLAGS) $(FINCS) -c $< -o $@

clean:
	rm -f *.o *.mod $(PROG) libmodel.a 
