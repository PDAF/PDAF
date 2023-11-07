BASEDIR = $(PDAF_DIR)
include $(BASEDIR)/make.arch/$(PDAF_ARCH).h
## includes

FFLAGS =  -cpp -I$(PDAF_DIR)/interface/model -I$(EBROOTNETCDFMINFORTRAN)/include $(TSMPPDAFIMPORTFLAGS)
CCFLAGS = -I$(PDAF_DIR)/interface/model -I$(EBROOTNETCDFMINFORTRAN)/include $(TSMPPDAFIMPORTFLAGS)
CPP_FLAGS = -Duse_libMPI -Duse_netCDF -Duse_comm_MPI1 -DVERBOSE -DDEBUG -DTREAT_OVERLAY -DFORTRANUNDERSCORE -DOFFLINE -DSPMD -DLINUX ${TSMPPDAFCPPDEFS}
FCPP_FLAGS= -Duse_libMPI -Duse_netCDF -Duse_comm_MPI1 -DVERBOSE -DDEBUG -DTREAT_OVERLAY -DFORTRANUNDERSCORE -DOFFLINE -DSPMD -DLINUX ${TSMPPDAFCPPDEFS}


## object files
OBJ = $(TSMPPDAFOBJ)

## target + rules
all: libmodel.a

libmodel.a: $(OBJ)
	ar rcs $@ *.o 
	ranlib $@
	mv $@ $(TSMPPDAFLIBDIR)
	@echo "library compilation suceeded"

wrapper_tsmp.o:
	$(PREP_C) $(CC) $(OPT) $(CCFLAGS) $(CPP_FLAGS) $(CINCS) -c wrapper_tsmp.c

%.o: common/%.c
	$(PREP_C) $(CC) $(OPT) $(CCFLAGS) $(CPP_FLAGS) $(CINCS) -c $<  -o $@

%.o: $(TSMPPDAFPFLDIR)/%.c
	$(PREP_C) $(CC) $(OPT) $(CCFLAGS) $(CPP_FLAGS) $(CINCS) -c $<  -o $@

%.o:$(TSMPPDAFCLMDIR)/%.F90
	$(PREP_C) $(FC) $(OPT) $(FCPP_FLAGS) $(FFLAGS) $(FINCS) -c $< -o $@

%.o: $(TSMPPDAFCOSDIR)/%.F90
	$(PREP_C) $(FC) $(OPT) $(FCPP_FLAGS) $(FFLAGS) $(FINCS) -c $< -o $@

clean:
	rm -f *.o *.mod $(PROG) libmodel.a 
