# $Id: Makefile 714 2009-05-27 12:00:25Z lnerger $

#######################################################
# Generic Makefile for to build PDAF with dummy model #
# To choose the architecture set $PDAF_ARCH           #
#######################################################

.SUFFIXES: .F90 .o 

######################################################

# Include machine-specific definitions
# For available include files see directory make.arch
# To choose a file, set PDAF_ARCH either here or by an
# environment variable.
include ../PDAF-D_V1.12/make.arch/$(PDAF_ARCH).h

######################################################

OBJ_GENOBS = generate_obs_AWI-CM.o
OBJ_COVAR  = generate_covar_AWI-CM.o

### Name of Executables ##############################################

EXE_GENOBS = generate_obs_AWI-CM
EXE_COVAR  = generate_covar_AWI-CM

######################################################

info:
	@echo "Makefile to build tool programs for the AWI-CM model";
	@echo  "-------------------------------------------------------------------------------";
	@echo  "Use as   make TARGET   where TARGET can be:";
	@echo  "  info           - This text";
	@echo  "  generate_obs_AWI-CM     - Generate an observation file by distubing a state trajectory";
	@echo  "  generate_covar_AWI-CM   - Generate a file for the covariance matrix from a trajectory";
	@echo  "  all            - Generate the two executables listed above";
	@echo  "  clean          - Clean up";
	@echo  "";
	@echo  "Use machine-specific definitions for PDAF_ARCH: " $(PDAF_ARCH);
	@echo  "-------------------------------------------------------------------------------";

$(EXE_GENOBS) : $(OBJ_GENOBS)
	rm -f $@
	$(LD)  $(OPT)  -o $@ $(MODULES) $(OBJ_GENOBS) $(NC_LIB) $(LINK_LIBS)

$(EXE_COVAR) : libpdaf-d.a $(OBJ_COVAR)
	rm -f $@
	$(LD)  $(OPT)  -o $@ $(MODULES) $(OBJ_COVAR) $(NC_LIB) $(LINK_LIBS) -L../PDAF-D_V1.12/lib -lpdaf-d

all: $(EXE_GENOBS) $(EXE_COVAR)

.F90.o :
	$(FC) $(OPT) $(MPI_INC) $(CPP_DEFS) $(NC_INC) -c $*.F90 

# For older compilers one might need to separate the
# preprocessing from the compilation as defined below:
#.F90.o :
#	$(CPP) -P -C $(MPI_INC) $(CPP_DEFS) $*.F90 $*.f90
#	$(FC) $(OPT) $(MPI_INC) -c $*.F90
#	@rm -f  $*.f90

libpdaf-d.a: 
	@echo "++++++ Generate Filter library ++++++"
	@cd ../PDAF-D_V1.12/src; make;


clean :
	rm -f *.o *.mod *.a $(EXE_GENOBS) $(EXE_COVAR)
